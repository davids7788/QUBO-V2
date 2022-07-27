import copy
import time
import numpy as np

from qiskit.utils import algorithm_globals
from qiskit.algorithms import NumPyMinimumEigensolver
from qiskit.optimization.algorithms import MinimumEigenOptimizer

from hamiltonian import Hamiltonian
from sub_qubo import SubQubo

algorithm_globals.massive = True


class QuboProcessing:
    def __init__(self,
                 triplet_list_file: str,
                 config: dict,
                 solver: Solver,
                 ansatz: Ansatz,
                 qubo_logging: QuboLogging):
        """Processes the Qubo and provides solving functions like a solving process via an impact list.
        :param triplet_list_file: .npy file with triplet objects
        :param config: dictionary with configuration parameters
        :param solver: solver object
        :param ansatz: ansatz circuit
        :param qubo_logging: object for handling information about the solving process of class QuboLogging
        """
        self.triplets = np.load(triplet_list_file, allow_pickle=True)
        self.config = config
        self.solver = solver
        self.solver.configure_solver()
        self.ansatz = ansatz
        self.error_mitigation = None
        self.qubo_logging = qubo_logging

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

        # Performing initial bit flip optimization
        if self.config["bit flip optimization"]["initial"] > 0:
            for i in range(self.config["bit flip optimization"]["initial"]):
                new_solution_candidate = self.bit_flip_optimization(self.triplets,
                                                                    self.solution_candidate,
                                                                    self.make_impact_list())
                self.energy_candidate = self.hamiltonian_energy(new_solution_candidate)
                self.qubo_logging.add_entry("solution vector", self.iteration, self.solution_candidate)
                self.qubo_logging.add_entry("energy", self.iteration, self.energy_candidate)
                self.iteration += 1
                print(f"Energy after performing {i + 1} bit flip tabu search(es): "
                      f"{np.around(self.energy_candidate, 4)}")

        # Here starts the sub qubo algorithm
        self.loop_count = 0
        self.pass_count = 0
        self.sub_qubo_counter = 0

    def impact_list_solve(self):
        """Performs impact list algorithm"""
        self.qubo_process(optimization_strategy=self.make_impact_list)

    def qubo_process(self,
                     optimization_strategy):
        """Solves the QUBO with
        :param optimization_strategy: function which is used to determine triplet ordering
        """
        # timer
        start_solving_process_time = time.time()

        # global optimization iterations
        while self.pass_count < search_depth:
            start_loop = time.time()
            start_quantum_part = time.time()

            # triplet list ordering determines also how many subqubos are built within one iteration
            # so it is possible to define a function with repeating indices resulting in a longer triplet ordering list
            # e.g. [1, 5, 3, 2, 4, 0], but also [1, 5, 1, 5, 3, 2, 4, 1, 5, ...]
            triplet_ordering = optimization_strategy

            # Calculate number of sub qubos:
            if len(triplet_ordering) % self.config["solver"]["num qubits"] == 0:
                num_sub_qubos = int(len(triplet_ordering) / self.config["solver"]["num qubits"])
            else:
                num_sub_qubos = int(len(triplet_ordering) / self.config["solver"]["num qubits"]) + 1
            print(f"Number of Subqubos: {num_sub_qubos} for iteration {self.loop_count}")

            # Sub-qubos with full size
            for i in range(int(len(self.triplets) / sub_qubo_size)):
                solve_subqubos([self.triplets[i] for i in triplet_ordering[sub_qubo_size * i: sub_qubo_size * (i + 1)]],
                               i)

            # Check if triplets are left -> take the last <sub_qubo_size> triplets and perform a vqe
            if len(triplet_ordering) % self.config["solver"]["num qubits"] != 0:
                solve_subqubos([self.triplets[i] for i in triplet_ordering[- sub_qubo_size:]],
                               int(len(triplet_ordering) / self.config["solver"]["num qubits"]) + 1)

            # Time tracking quantum part
            end_quantum_part = time.time()
            quantum_time = end_quantum_part - start_quantum_part
            self.qubo_logging.add_qubo_log_entry("time tracking quantum",
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
            end_loop = time.time()
            loop_time = end_loop - start_loop
            self.qubo_logging.add_qubo_log_entry("time tracking qubo iteration",
                                                 str(self.loop_count),
                                                 loop_time)

            self.qubo_logging.save_results(self.save_folder)

        end_solving_process_time = time.time()
        solving_process_time = end_solving_process_time - start_solving_process_time
        self.qubo_logging.add_qubo_log_entry("time tracking complete",
                                             "complete run",
                                             solving_process_time)
        print(f"Solving process for all solver needed {hms_string(solving_process_time)}")

    def make_impact_list(self):
        """
        Creates an impact list based on how much influence on the energy a bit flip has
        :return: list of indices ordered from lowest to highest impact of triplets in triplet list
        """
        impact_list_values = []
        for triplet, t_i in zip(self.triplets, self.solution_candidate):
            energy_change = 0
            if t_i == 0:
                energy_change += triplet.quality
            else:
                energy_change -= triplet.quality

            for interaction in triplet.interactions.keys():
                if t_i == 0 and self.solution_candidate[interaction] == 0:
                    pass
                elif t_i == 0 and self.solution_candidate[interaction] == 1:
                    energy_change += triplet.interactions[interaction]
                elif t_i == 1 and self.solution_candidate[interaction] == 0:
                    pass
                else:
                    energy_change -= triplet.interactions[interaction]

            impact_list_values.append(abs(energy_change))
        return list(np.argsort(impact_list_values))

    @staticmethod
    def bit_flip_optimization(triplets,
                              solution_candidate,
                              impact_list,
                              reverse=True):

        """Looping over a list of triplet objects and a corresponding binary solution vector to compute if the energy
        value would improve if the binary state of a triplet (keep <-> discard) should be changed. If the energy
        decreases, the triplet state is flipped.
        :param triplets: list of triplet objects
        :param solution_candidate binary vector representing triplet states
        :param impact_list: index list of triplets
        :param reverse: True if sorting from highest impact to lowest, else from lowest to highest
        """
        if reverse:
            impact.list.reverse()
        for i, triplet in enumerate([triplets[impact_index] for impact_index in impact_list]):
            # energy change if this particular bit is flipped
            energy_change = 0

            # checking linear term
            if solution_candidate[impact_list[i]] == 0:
                energy_change += triplet.quality
            else:
                energy_change -= triplet.quality

            # Checking interactions with other triplets
            for interaction in triplet.interactions.keys():
                if solution_candidate[impact_list[i]] == 0 and solution_candidate[interaction] == 0:
                    pass
                elif solution_candidate[impact_list[i]] == 0 and solution_candidate[interaction] == 1:
                    energy_change += triplet.interactions[interaction]
                elif solution_candidate[impact_list[i]] == 1 and solution_candidate[interaction] == 0:
                    pass
                else:
                    energy_change -= triplet.interactions[interaction]

            # flip if overall energy change is negative
            if energy_change < 0:
                solution_candidate[impact_list[i]] = 1 - solution_candidate[impact_list[i]]

        return solution_candidate

    def minimal_energy_and_solution(self):
        """Calculates the minimum energy state and value.
        :return: minimum energy state, minimum energy value """
        minimum_energy_state = []
        for t in self.triplets:
            if t.is_correct_match:
                minimum_energy_state.append(1)
            else:
                minimum_energy_state.append(0)
        return minimum_energy_state, self.hamiltonian_energy(minimum_energy_state)

    def perform_tabu_search(self, tabu_search, key):
        impacts = self.make_impact_list()
        if key == "initial tabu search":
            impacts.reverse()
        # tabu_search_time_start = time.time()
        solution_candidate = copy.deepcopy(self.solution_candidate)
        solution_candidate_after_tabu, num_bit_flips = tabu_search(self.triplets, solution_candidate, impacts)
        energy_after_tabu = self.hamiltonian_energy(solution_candidate)
        if energy_after_tabu < self.energy_candidate:
            self.energy_candidate = energy_after_tabu
            self.solution_candidate = solution_candidate_after_tabu
        # tabu_search_time_end = time.time()
        # tabu_time = tabu_search_time_end - tabu_search_time_start
        # self.qubo_logging.add_tabu_log_entry("bit flips", key, num_bit_flips)
        # self.qubo_logging.add_tabu_log_entry("time tracking", key, tabu_time)
        # self.qubo_logging.add_qubo_log_entry("solution vector", key, solution_candidate_after_tabu)
        # self.qubo_logging.add_qubo_log_entry("energy", key, energy_after_tabu)
        # self.qubo_logging.add_qubo_log_entry("triplet sorting",
        #                                      key,
        #                                      compare_bits(self.triplets, self.solution_candidate))

    def hamiltonian_energy(self, binary_vector):
        """
        Calculates the energy according to a binary vector matching the triplet list
        :param binary_vector: binary solution candidate vector
        :return: energy value
        """
        hamiltonian_energy = 0
        for i, b1 in enumerate(binary_vector):
            if b1 == 1:
                hamiltonian_energy += self.triplets[i].quality
            for j in self.triplets[i].interactions.keys():
                if j < i:
                    continue
                else:
                    if binary_vector[i] == binary_vector[j] == 1:
                        hamiltonian_energy += self.triplets[i].interactions[j]
        return hamiltonian_energy

    def log_truth_energy(self):
        """Obtaining minimal energy solution and minimal energy printing information about ideal solution and energy.
        """
        min_energy_state, min_energy = self.minimal_energy_and_solution()
        self.qubo_logging.set_value("truth solution vector", min_energy_state)
        self.qubo_logging.set_value("truth minimum energy", min_energy)
        print(f"Minimal Energy of the problem: {np.around(min_energy, 3)}")

    def set_initial_vector(self):
        """Setting the initial vector according to the configuration file.
        """
        if self.config["initial binary vector"] == "zeros":
            self.solution_candidate = np.zeros(len(self.triplets))
            print("Starting with every triplet set to 0\n")
        elif self.config["initial binary vector"] == "ones":
            self.solution_candidate = np.ones(len(self.triplets))
            print("Starting with every triplet set to 1\n")
        elif self.config["initial binary vector"] is None:
            self.solution_candidate = np.random.randint(2, size=len(self.triplets))
            print("Starting with a random vector\n")

    def solve_quantum(self,
                      hamiltonian):
        """Solves the sub qubo. Returns the result in the correct order!  -> qiskit qubit ordering
        :param hamiltonian: hamiltonian in the form of a PauliSumOp list from qiskit
        :return
            result of the solving process"""
        result_quantum = self.solver.compute_minimum_eigenvalue(hamiltonian)
        if self.config["solver"]["error mitigation algorithm"] == "algebraic":
            result_quantum = self.error_mitigation.meas_filter.apply(result_quantum)

        def get_result_in_correct_order(res):
            """Returns the result with the highest probability in the correct ordering of the qubits.
            :param res: result obtained via quantum algorithm
            :return res_out: result with the highest probability in the correct ordering of the qubits"""
            highest_prob = 0
            key_best_vqe = ""
            for key, prob, in zip(res.eigenstate.keys(), res.eigenstate.values()):
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
                       index):

        result = None
        sub_qubo_time = None

        # start timer sub-qubo creation and solving
        print(f"Processing Sub-qubo {index + 1} of {num_sub_qubos}", end="\r")
        hamiltonian = Hamiltonian(triplet_slice=triplet_slice,
                                  solution_candidate=self.solution_candidate,
                                  rescaling=self.config["qubo"]["hamiltonian rescaling"])
        if self.loop_count == 0:  # Uses a lot of disk space, so just results of first iteration
            self.qubo_logging.add_qubo_log_entry("hamiltonian",
                                                 self.sub_qubo_counter,
                                                 [hamiltonian.linear_term(), hamiltonian.quadratic_term()])

        # Create hamiltonian and check if no linear and quadratic terms --> do nothing
        hamiltonian = qubo_representation_vqe()
        if hamiltonian is None:
            return

        # Overwrite ansatz for vqe if HamiltonianDriven
        if self.ansatz.layout == "HamiltonianDriven":
            triplet_slice = triplet_slice
            self.ansatz = self.ansatz.hamiltonian_driven(triplet_slice)
            self.solver.quantum_instance.ansatz = self.ansatz

        if self.config["algorithm"] == "VQE" or self.config["algorithm"] == "QAOA":

            # Timestamps and solving the SubQUBO
            start_sub_qubo = time.time()
            result = self.solve_quantum(hamiltonian=hamiltonian.qubo_representation_quantum())
            end_sub_qubo = time.time()
            sub_qubo_time = end_sub_qubo - start_sub_qubo

            # Check if comparing to Eigensolver
            if self.config["solver"]["compare to eigensolver"] and self.loop_count == 0:  # one iteration
                exact_result = self.exact_solver.solve(hamiltonian.qubo_representation()).x
                self.qubo_logging.add_qubo_log_entry("subqubo solving success",
                                                     self.sub_qubo_counter,
                                                     ([float(r) for r in result] == exact_result))

        elif self.config["algorithm"] == "Numpy Eigensolver":
            start_sub_qubo = time.time()
            result = self.exact_solver.solve(hamiltonian.qubo_representation()).x
            end_sub_qubo = time.time()
            sub_qubo_time = end_sub_qubo - start_sub_qubo

        # update solution
        for k, entry in enumerate(impact_list[sub_qubo_size * i: sub_qubo_size * (i + 1)]):
            self.solution_candidate[entry] = result[k]
        self.qubo_logging.add_qubo_log_entry("time tracking subQUBOs",
                                             self.sub_qubo_counter,
                                             sub_qubo_time)

        self.sub_qubo_counter += 1

    def configure_solver(self):
        # Set up the the chosen algorithm, VQE, QAOA or NumpyEigensolver
        # VQE --> ansatz can be chosen
        if self.config["algorithm"] == "VQE":
            self.solver.set_vqe(self.ansatz)

        # QAOA --> ansatz determined by algorithm
        elif self.config["algorithm"] == "QAOA":
            self.solver.set_qaoa()

        # Numpy Eigensolver does not need further preparation
        elif self.config["algorithm"] == "Numpy Eigensolver":
            pass
        else:
            "No valid algorithm chosen!"
