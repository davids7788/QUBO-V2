import time
import numpy as np

from qiskit.utils import algorithm_globals
from qiskit.algorithms import NumPyMinimumEigensolver
from qiskit_optimization.algorithms import MinimumEigenOptimizer

from qubo.optimisation import make_impact_list, bit_flip_optimisation, make_connection_list, make_paired_list, \
                              impact_without_conflicts
from qubo.solver import Solver
from qubo.hamiltonian import Hamiltonian
from qubo.ansatz import Ansatz
from qubo.qubo_logging import QuboLogging

from utility.time_tracking import hms_string

algorithm_globals.massive = True


class QuboProcessing:
    def __init__(self,
                 triplet_list_file: str,
                 config: dict,
                 solver: Solver,
                 ansatz: Ansatz,
                 qubo_logging: QuboLogging,
                 save_folder: str,
                 verbose=1):
        """Processes the Qubo and provides solving functions like a solving process via a chosen optimisation strategy.
        :param triplet_list_file: .npy file with triplet objects
        :param config: dictionary with configuration parameters
        :param solver: solver object
        :param ansatz: ansatz circuit
        :param qubo_logging: object for handling information about the solving process of class QuboLogging
        :param save_folder: folder to which results are stored
        :param verbose: 0: only energy output, 1: showing sub-QUBO process
        """
        self.verbose = verbose
        self.triplets = np.load(triplet_list_file, allow_pickle=True)
        self.t_mapping = {key: value for key, value in zip([t.triplet_id for t in self.triplets],
                                                           np.arange(len(self.triplets)))}

        self.config = config
        self.solver = solver
        self.ansatz = ansatz
        if self.solver is not None:
            self.solver.configure_solver(self.ansatz)
        self.qubo_logging = qubo_logging
        if "impact list" in self.config["qubo"]["optimisation strategy"]:
            self.optimisation_strategy = make_impact_list
        if "connection list" in self.config["qubo"]["optimisation strategy"]:
            self.optimisation_strategy = make_connection_list
        if "paired list" in self.config["qubo"]["optimisation strategy"]:
            self.optimisation_strategy = make_paired_list
        if "impact without conflicts" in self.config["qubo"]["optimisation strategy"]:
            self.optimisation_strategy = impact_without_conflicts
        self.save_folder = save_folder

        # Log truth minimum energy state and energy
        self.log_truth_energy()

        # Setting initial binary vector, providing variable for current solution candidate and its energy
        self.solution_candidate = self.set_initial_solution()
        self.energy_candidate = self.hamiltonian_energy(self.solution_candidate)
        self.qubo_logging.set_value("computed solution vector", self.solution_candidate)
        self.qubo_logging.set_value("computed minimum energy", self.energy_candidate)
        self.qubo_logging.add_entry("solution vector", 0, self.solution_candidate)
        self.qubo_logging.add_entry("energy", 0, self.energy_candidate)

        self.iteration = 1   # global counter for keeping track of iterations / epochs

        # Create solver object
        self.exact_solver = MinimumEigenOptimizer(NumPyMinimumEigensolver())

        # Performing initial bit flip optimisation
        if self.config["bit flip optimisation"]["iterations"] is not None:
            start_bit_flip = time.process_time()
            for i in range(self.config["bit flip optimisation"]["iterations"]):
                new_solution_candidate, energy_change = bit_flip_optimisation(self.triplets,
                                                                              self.solution_candidate,
                                                                              self.t_mapping,
                                                                              self.optimisation_strategy(
                                                                                   self.triplets,
                                                                                   self.t_mapping,
                                                                                   self.solution_candidate),
                                                                              self.config["bit flip optimisation"]
                                                                              ["reverse"])
                self.energy_candidate += energy_change
                self.qubo_logging.add_entry("solution vector", self.iteration, self.solution_candidate)
                self.qubo_logging.add_entry("energy", self.iteration, self.energy_candidate)
                self.iteration += 1
                print(f"Energy after performing {i + 1} bit flip optimisation(s): "
                      f"{np.around(self.energy_candidate, 2)}")
            if self.config["qubo"]["search depth"] is None:
                self.qubo_logging.save_results(self.save_folder)
                self.write_setup_info_to_file()
            end_bit_flip = time.process_time()
            print(f"Bit flip optimisation needed {QuboProcessing.hms_string((end_bit_flip - start_bit_flip))}")

        # Here starts the sub qubo algorithm
        self.loop_count = 0
        self.pass_count = 0
        self.sub_qubo_counter = 0

    def qubo_processing(self):
        """Solves the QUBO with the set optimisation strategy.
        """
        # timer
        start_solving_process_time = time.process_time()

        # global optimisation iterations
        while self.pass_count < self.config["qubo"]["search depth"]:
            start_loop = time.process_time()
            start_quantum_part = time.process_time()
            triplet_ordering = [i for i in range(len(self.triplets))]

            # triplet list ordering determines also how many sub-QUBOs are built within one iteration,
            # so it is possible to define a function with repeating indices resulting in a longer triplet ordering list
            # e.g. [1, 5, 3, 2, 4, 0], but also [1, 5, 1, 5, 3, 2, 4, 1, 5, ...] would be possible

            if "impact" in self.config["qubo"]["optimisation strategy"]:
                triplet_ordering = self.optimisation_strategy(self.triplets,
                                                              self.t_mapping,
                                                              self.solution_candidate,
                                                              False)
            if "connection list" in self.config["qubo"]["optimisation strategy"]:
                triplet_ordering = self.optimisation_strategy(self.triplets,
                                                              self.t_mapping)
            if "paired list" in self.config["qubo"]["optimisation strategy"]:
                triplet_ordering = self.optimisation_strategy(self.triplets,
                                                              self.t_mapping)
            if "reverse" in self.config["qubo"]["optimisation strategy"]:
                triplet_ordering.reverse()

            # Calculate number of sub qubos:
            if len(triplet_ordering) % self.config["qubo"]["num qubits"] == 0:
                num_sub_qubos = int(len(triplet_ordering) / self.config["qubo"]["num qubits"])
            else:
                num_sub_qubos = int(len(triplet_ordering) / self.config["qubo"]["num qubits"]) + 1

            # Sub-qubos with full size
            for i in range(int(len(self.triplets) / self.config["qubo"]["num qubits"])):
                if self.verbose == 1:
                    print(f"Processing Sub-QUBO {i + 1} of {num_sub_qubos}", end="\r")
                result = self.solve_subqubos([self.triplets[i] for i in
                                             triplet_ordering[self.config["qubo"]["num qubits"] * i:
                                                              self.config["qubo"]["num qubits"] * (i + 1)]])
                if result is None:
                    continue
                for k, entry in enumerate(triplet_ordering[self.config["qubo"]["num qubits"] * i:
                                                           self.config["qubo"]["num qubits"] * (i + 1)]):
                    self.solution_candidate[entry] = result[k]

            # Check if triplets are left -> take the last <sub_qubo_size> triplets and perform a vqe
            if len(triplet_ordering) % self.config["qubo"]["num qubits"] != 0:
                if self.verbose == 1:
                    print(f"Processing Sub-QUBO {num_sub_qubos} of {num_sub_qubos}", end="\r")
                result = self.solve_subqubos([self.triplets[i] for i in
                                              triplet_ordering[- self.config["qubo"]["num qubits"]:]])
                if result is None:
                    pass
                else:
                    for k, entry in enumerate(triplet_ordering[- self.config["qubo"]["num qubits"]:]):
                        self.solution_candidate[entry] = result[k]

            # Time tracking quantum part
            end_quantum_part = time.process_time()
            quantum_time = end_quantum_part - start_quantum_part
            self.qubo_logging.add_entry("time tracking quantum",
                                        str(self.loop_count),
                                        quantum_time)
            # update energy
            self.energy_candidate = self.hamiltonian_energy(self.solution_candidate)

            self.qubo_logging.add_entry("solution vector",
                                        str(self.loop_count),
                                        self.solution_candidate)
            self.qubo_logging.add_entry("energy",
                                        str(self.loop_count),
                                        self.energy_candidate)

            # check if new solution candidate has a lower energy value than the one before
            if self.energy_candidate < self.qubo_logging.qubo_log["computed minimum energy"]:
                self.qubo_logging.set_value("computed minimum energy",
                                            self.energy_candidate)
                self.qubo_logging.set_value("computed solution vector",
                                            self.solution_candidate)

            # increasing pass count if a lower energy was not found in this iteration
                self.pass_count = 0
            else:
                self.pass_count += 1

            # increasing overall loop count
            self.loop_count += 1

            # time tracking iteration
            end_loop = time.process_time()
            loop_time = end_loop - start_loop
            self.qubo_logging.add_entry("time tracking qubo iteration",
                                        str(self.loop_count),
                                        loop_time)

            print(f"Energy after iteration {self.loop_count} at pass count {self.pass_count}: "
                  f"{np.around(self.energy_candidate, 2)}  |  "
                  f"Time elapsed: {QuboProcessing.hms_string(loop_time)}")

            self.qubo_logging.save_results(self.save_folder)

        end_solving_process_time = time.process_time()
        solving_process_time = end_solving_process_time - start_solving_process_time
        self.qubo_logging.add_entry("time tracking complete",
                                    "complete run",
                                    solving_process_time)
        self.qubo_logging.save_results(self.save_folder)
        self.write_setup_info_to_file()
        if self.config["qubo"]["search depth"] > 0:
            print(f"Qubo solving process needed {QuboProcessing.hms_string(solving_process_time)}")

    def minimal_energy_and_solution(self):
        """Calculates the minimum energy state and value.
        :return:
            minimum energy state, minimum energy value """
        minimum_energy_state = []
        for t in self.triplets:
            if t.from_same_particle():
                minimum_energy_state.append(1)
            else:
                minimum_energy_state.append(0)
        return minimum_energy_state, self.hamiltonian_energy(minimum_energy_state)

    def hamiltonian_energy(self, binary_vector, triplet_subset=None):
        """Calculates the energy according to a binary vector matching the triplet list
        :param binary_vector: binary solution candidate vector
        :param triplet_subset: subset of triplet, represented by indices
        :return:
            energy value
        """

        hamiltonian_energy = 0
        if triplet_subset is None:
            for i, b1 in enumerate(binary_vector):
                if b1 == 1:
                    hamiltonian_energy += self.triplets[i].quality
                for j in self.triplets[i].interactions.keys():
                    if self.t_mapping[j] < i:
                        continue
                    else:
                        pass
                        if binary_vector[i] == binary_vector[self.t_mapping[j]] == 1:
                            hamiltonian_energy += self.triplets[i].interactions[j]
            return hamiltonian_energy

        else:
            for i, b1 in enumerate(binary_vector):
                if b1 == 1:
                    hamiltonian_energy += self.triplets[triplet_subset[i]].quality

                for j in self.triplets[triplet_subset[i]].interactions.keys():
                    if j < triplet_subset[i]:
                        continue
                    if j not in triplet_subset:
                        continue
                    position = triplet_subset.index(j)
                    if binary_vector[i] == binary_vector[position] == 1:
                        hamiltonian_energy += self.triplets[triplet_subset[i]].interactions[j]
            return hamiltonian_energy

    def log_truth_energy(self):
        """Obtaining minimal energy solution and minimal energy printing information about ideal solution and energy.
        """
        min_energy_state, min_energy = self.minimal_energy_and_solution()
        self.qubo_logging.set_value("truth solution vector", min_energy_state)
        self.qubo_logging.set_value("truth minimum energy", min_energy)
        print(f"Minimal Energy of the problem: {np.around(min_energy, 3)}")

    def set_initial_solution(self):
        """Setting the initial vector according to the configuration file.
        """
        if self.config["qubo"]["initial binary vector"] == "zeros":
            solution_candidate = np.zeros(len(self.triplets))
            print("Starting with every triplet set to 0\n")
        elif self.config["qubo"]["initial binary vector"] == "ones":
            solution_candidate = np.ones(len(self.triplets))
            print("Starting with every triplet set to 1\n")
        else:
            solution_candidate = np.random.randint(2, size=len(self.triplets))
            print("Starting with a random vector\n")
        return solution_candidate

    def solve_quantum(self,
                      hamiltonian):
        """Solves the sub qubo. Returns the result in the correct order!  -> qiskit qubit ordering
        :param hamiltonian: hamiltonian in the form of a PauliSumOp list from qiskit
        :return
            result of the solving process"""
        result_quantum = dict(self.solver.quantum_algorithm.compute_minimum_eigenvalue(hamiltonian).eigenstate)

        def get_result_in_correct_order(res):
            """Returns the result with the highest probability in the correct ordering of the qubits.
            :param res: result obtained via quantum algorithm
            :return
                result with the highest probability in the correct ordering of the qubits"""
            highest_prob = 0
            key_best_vqe = ""
            for key, prob, in zip(res.keys(), res.values()):
                if prob > highest_prob:
                    highest_prob = prob
                    key_best_vqe = key
            res_out = [0] * (len(key_best_vqe))
            for i, bit in enumerate(key_best_vqe):
                res_out[- i - 1] = bit
            return res_out

        return get_result_in_correct_order(result_quantum)

    def solve_subqubos(self,
                       triplet_slice,
                       subsequent_list=None):
        """Function for solving a SubQUBO
        :param triplet_slice: slice of triplets used for the SubQUBO
        :param subsequent_list if only triplets in a certain area should be considered
        :return
            result of computation
        """
        result = None
        sub_qubo_time = None

        # start timer sub-qubo creation and solving

        hamiltonian = Hamiltonian(triplet_slice=triplet_slice,
                                  t_mapping=self.t_mapping,
                                  solution_candidate=self.solution_candidate,
                                  rescaling=self.config["qubo"]["hamiltonian rescaling"],
                                  only_specified_connections=subsequent_list)
        if self.loop_count == 0:  # Uses a lot of disk space, so just results of first iteration
            self.qubo_logging.add_entry("hamiltonian",
                                        self.sub_qubo_counter,
                                        [hamiltonian.linear_term(), hamiltonian.quadratic_term()])

        # Create hamiltonian and check if no linear and quadratic terms --> do nothing
        if hamiltonian.qubo_representation_quantum() is None:
            return

        if self.config["solver"]["algorithm"] == "VQE" or self.config["solver"]["algorithm"] == "QAOA":

            # Overwrite ansatz for vqe if HamiltonianDriven

            if self.config["ansatz"]["layout"] == "HamiltonianDriven":
                triplet_slice = triplet_slice
                self.ansatz.set_hamiltonian_driven(triplet_slice)
                self.solver.quantum_instance.ansatz = self.ansatz

            # Timestamps and solving the SubQUBO
            start_sub_qubo = time.process_time()
            result = self.solve_quantum(hamiltonian=hamiltonian.qubo_representation_quantum())
            end_sub_qubo = time.process_time()
            sub_qubo_time = end_sub_qubo - start_sub_qubo

            # Check if comparing to Eigensolver
            if self.config["qubo"]["compare to eigensolver"]:
                exact_result = self.exact_solver.solve(hamiltonian.qubo_representation()).x
                self.qubo_logging.add_entry("compare to analytical solution",
                                            self.sub_qubo_counter,
                                            list([float(r) for r in result] == exact_result).
                                            count(True) == self.config["qubo"]["num qubits"])

        elif self.config["solver"]["algorithm"] == "Numpy Eigensolver":
            start_sub_qubo = time.process_time()
            result = self.exact_solver.solve(hamiltonian.qubo_representation()).x
            end_sub_qubo = time.process_time()
            sub_qubo_time = end_sub_qubo - start_sub_qubo

        # update solution

        self.qubo_logging.add_entry("time tracking subQUBOs",
                                    self.sub_qubo_counter,
                                    sub_qubo_time)

        self.sub_qubo_counter += 1
        return result

    @staticmethod
    def hms_string(sec_elapsed):
        """Nicely formatted time string.
        """
        h = int(sec_elapsed / (60 * 60))
        m = int((sec_elapsed % (60 * 60)) / 60)
        s = sec_elapsed % 60
        return "{}:{:>02}:{:>05.2f}".format(h, m, s)

    def write_setup_info_to_file(self):
        """Writes configuration details and track efficiency and fake rate results into a .txt file.
        """
        with open(self.save_folder + "/qubo_info.txt", "w") as f:
            f.write("Qubo solved with the following configuration: \n")
            f.write("---\n")
            for outer_key in self.config.keys():
                for inner_key, value in self.config[outer_key].items():
                    f.write(f"\n{inner_key}: {value}")
                f.write("\n")
            f.write("\n\n")
            f.write("---\n")

    def get_kept_triplets(self):
        """Collects the kept triplets from the triplet list
        :return
            kept triplets selected by the optimisation process
        """
        kept_triplets = []
        for triplet, entry in zip(self.triplets, self.solution_candidate):
            if int(entry) == 1:
                kept_triplets.append(triplet)
        return kept_triplets
