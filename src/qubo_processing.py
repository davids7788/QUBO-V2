import copy
import time
import numpy as np

from qiskit.utils import algorithm_globals
from qiskit.algorithms import NumPyMinimumEigensolver
from qiskit.optimization.algorithms import MinimumEigenOptimizer

from functions_library import TabuSearch, hms_string, compare_bits
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
        self.ansatz = ansatz
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

        # Calculate number of sub qubos:
        if len(self.triplets) % self.config["solver"]["num qubits"] <= 1:
            self.num_sub_qubos = int(len(self.triplets) / self.config["solver"]["num qubits"])
        else:
            self.num_sub_qubos = int(len(self.triplets) / sself.config["solver"]["num qubits"]) + 1
        print(f"Number of Subqubos: {num_sub_qubos}")

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

    def impact_list_solve(self):
        """Performs impact list algorithm"""

        self.qubo_process(self.make_impact_list)

    def qubo_process(self,
                     optimization_strategy):
        """Solves the QUBO with a simple impact list approach
        :param optimization_strategy: function which is used to determine triplet ordering
        """
        start_solving_process_time = time.time()

        # Here starts the sub qubo algorithm
        loop_count = 0
        pass_count = 0
        sub_qubo_counter = 0

        # global optimization iterations
        while pass_count < search_depth:
            start_loop = time.time()
            start_quantum_part = time.time()

            # # updating the impact list
            triplet_ordering = optimization_strategy

            # sub-qubos with full size
            for i in range(int(len(self.triplets) / sub_qubo_size)):
                # start timer sub-qubo creation and solving
                print(f"Processing Sub-qubo {i + 1} of {num_sub_qubos}", end="\r")
                hamiltonian = Hamiltonian(triplet_slice=[self.triplets[i] for i in
                                                         triplet_ordering[sub_qubo_size * i: sub_qubo_size * (i + 1)]],
                                          solution_candidate=self.solution_candidate,
                                          rescaling=self.config["qubo"]["hamiltonian rescaling"])
                if sub_qubo_counter < 500:   # Uses a lot of disk space, so just saving some
                    self.qubo_logging.add_qubo_log_entry("hamiltonian",
                                                         sub_qubo_counter,
                                                         [hamiltonian.linear_term(), hamiltonian.quadratic_term()])

                hamiltonian = qubo_representation_vqe()
                if hamiltonian is None:
                    continue

                if self.algorithm == "VQE":
                    if self.ansatz.layout == "HamiltonianDriven":
                        triplet_slice = [self.triplets[i] for i in
                                         impact_list[sub_qubo_size * i: sub_qubo_size * (i + 1)]]
                        self.ansatz.set_hamiltonian_driven_ansatz(triplet_slice)
                    elif self.ansatz.layout == "TwoLocal":
                        pass  # Needs to be only set once, was done in "main.py"
                    elif self.solver.ansatz.name is None:
                        self.ansatz.set_no_entanglement_circuit()
                    self.solver.set_vqe_backend

                if "QAOA" in self.mode:
                    self.solver.set_vqe_ideal_qaoa(self.ansatz)
                else:
                    self.solver.set_vqe_ideal_qasm(self.ansatz)
                    sub_qubo = SubQubo(self.solver, hamiltonian.qubo_representation_vqe(), self.ansatz)

                    start_sub_qubo = time.time()
                    result = sub_qubo.solve_vqe()
                    end_sub_qubo = time.time()

                    # compare to NumpyEigensolver
                    if compare_to_eigensolver:
                        exact_result = self.exact_solver.solve(hamiltonian.qubo_representation()).x

                        self.qubo_logging.add_qubo_log_entry("subqubo solving success",
                                                             sub_qubo_counter,
                                                             ([float(r) for r in result] == exact_result))
                        print([float(r) for r in result] == exact_result)

                    # update solution
                    for k, entry in enumerate(impact_list[sub_qubo_size * i: sub_qubo_size * (i + 1)]):
                        self.solution_candidate[entry] = result[k]
                    self.qubo_logging.add_qubo_log_entry("time tracking subqubos",
                                                         sub_qubo_counter,
                                                         (end_sub_qubo - start_sub_qubo))

                elif self.mode == "Numpy Eigensolver":
                    start_sub_qubo = time.time()
                    result = self.exact_solver.solve(hamiltonian.qubo_representation()).x
                    end_sub_qubo = time.time()

                    # update solution
                    for k, entry in enumerate(impact_list[sub_qubo_size * i: sub_qubo_size * (i + 1)]):
                        self.solution_candidate[entry] = result[k]
                    self.qubo_logging.add_qubo_log_entry("time tracking subqubos",
                                                         sub_qubo_counter,
                                                         end_sub_qubo - start_sub_qubo)
                sub_qubo_counter += 1

            # Check if more than 1 triplet is left:
            # 1) More than one triplet is left: Take the least <sub_qubo_size> triplets and perform a vqe
            # 2) One or zero are left: Just make a bit flip and compare results
            if len(self.triplets) % sub_qubo_size != 0:
                hamiltonian = Hamiltonian(triplet_slice=[self.triplets[i] for i in impact_list[- sub_qubo_size:]],
                                          solution_candidate=self.solution_candidate,
                                          normalizing=self.normalizing_hamiltonian)

                if "VQE" in self.mode:
                    hamiltonian.qubo_representation_vqe()
                    if self.ansatz.name == "hamiltonian aware":
                        triplet_slice = [self.triplets[i] for i in impact_list[- sub_qubo_size:]]
                        self.ansatz.set_hamiltonian_aware_ansatz(triplet_slice)
                    elif self.ansatz.name == "TwoLocal":
                        # already set at the main.py, because does not change during the entire program
                        pass
                    elif self.solver.ansatz.name is None:
                        self.ansatz.set_no_entanglement_circuit()

                    if "QAOA" in self.mode:
                        self.solver.set_vqe_ideal_qaoa(self.ansatz)
                    else:
                        self.solver.set_vqe_ideal_qasm(self.ansatz)
                    sub_qubo = SubQubo(self.solver, hamiltonian, self.ansatz)

                    # start_sub_qubo = time.time()
                    result = sub_qubo.solve_vqe()
                    # end_sub_qubo = time.time()

                    # compare to NumpyEigensolver
                    if compare_to_eigensolver:
                        exact_result = self.exact_solver.solve(hamiltonian.qubo_representation()).x

                        self.qubo_logging.add_qubo_log_entry("subqubo solving success",
                                                             sub_qubo_counter,
                                                             ([float(r) for r in result] == exact_result))
                        print([float(r) for r in result] == exact_result)

                    # update solution
                    for k, entry in enumerate(impact_list[- sub_qubo_size:]):
                        self.solution_candidate[entry] = result[k]
                    # self.qubo_logging.add_qubo_log_entry("time tracking subqubos",
                    #                                      sub_qubo_counter,
                    #                                      (end_sub_qubo - start_sub_qubo))

                elif self.mode == "Numpy Eigensolver":
                    # start_sub_qubo = time.time()
                    result = self.exact_solver.solve(hamiltonian.qubo_representation()).x
                    # end_sub_qubo = time.time()

                    # update solution
                    for k, entry in enumerate(impact_list[- sub_qubo_size:]):
                        self.solution_candidate[entry] = result[k]
                    # self.qubo_logging.add_qubo_log_entry("time tracking subqubos",
                    #                                      sub_qubo_counter,
                    #                                      end_sub_qubo - start_sub_qubo)

            elif len(self.triplets) % sub_qubo_size == 1:
                self.solution_candidate[impact_list[-1]] = 1 - self.solution_candidate[impact_list[-1]]
                last_bit_flip_energy = self.hamiltonian_energy(self.solution_candidate)
                if last_bit_flip_energy < self.energy_candidate:
                    self.energy_candidate = last_bit_flip_energy
                else:
                    self.solution_candidate[impact_list[-1]] = 1 - self.solution_candidate[impact_list[-1]]

            end_quantum_part = time.time()
            quantum_time = end_quantum_part - start_quantum_part
            self.qubo_logging.add_qubo_log_entry("time tracking quantum",
                                                 str(loop_count),
                                                 quantum_time)
            # update energy
            self.energy_candidate = self.hamiltonian_energy(self.solution_candidate)

            # loop wise tabu search
            if tabu_usage["loop wise"]:
                print("Performing loop-wise tabular search")
                self.perform_tabu_search(tabu_search, str(loop_count) + "_loop_wise")
            else:
                self.qubo_logging.add_qubo_log_entry("solution vector",
                                                     str(loop_count),
                                                     self.solution_candidate)
                self.qubo_logging.add_qubo_log_entry("energy",
                                                     str(loop_count),
                                                     self.energy_candidate)

            triplet_statistics = compare_bits(self.triplets, self.solution_candidate)
            self.qubo_logging.add_qubo_log_entry("triplet statistics",
                                                 str(loop_count),
                                                 triplet_statistics)

            # check if new solution candidate has a lower energy value than the one before
            if self.energy_candidate < self.qubo_logging.qubo_log["computed minimum energy"]:
                self.qubo_logging.set_qubo_log_value("computed minimum energy",
                                                     self.energy_candidate)
                self.qubo_logging.set_qubo_log_value("minimal energy solution",
                                                     self.solution_candidate)

                print(f"Minimum energy found: {np.round(float(self.energy_candidate), 3)} at pass count: {pass_count}")
                print(f"True positives: {triplet_statistics[0]}\n"
                      f"True negatives: {triplet_statistics[1]}\n"
                      f"False positives: {triplet_statistics[2]}\n"
                      f"False negatives: {triplet_statistics[3]}\n"
                      f"Accuracy: {np.around(sum(triplet_statistics[0:2]) / sum(triplet_statistics[0:4]), 4)}\n")
                pass_count = 0
            else:
                print(f"Energy after loop {pass_count}: {np.round(self.energy_candidate, 3)}              ")
                print(f"True positives: {triplet_statistics[0]}\n"
                      f"True negatives: {triplet_statistics[1]}\n"
                      f"False positives: {triplet_statistics[2]}\n"
                      f"False negatives: {triplet_statistics[3]}\n"
                      f"Accuracy: {np.around(sum(triplet_statistics[0:2]) / sum(triplet_statistics[0:4]), 4)}\n")
                pass_count += 1

            # increasing overall loop count
            loop_count += 1

            end_loop = time.time()
            loop_time = end_loop - start_loop
            self.qubo_logging.add_qubo_log_entry("time tracking qubo iteration",
                                                 str(loop_count),
                                                 loop_time)

            self.qubo_logging.save_results(self.save_folder)

        # final tabu search
        if tabu_usage["final"]:
            print("Performing final tabu search...")
            self.perform_tabu_search(tabu_search, str(loop_count) + "_final")

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
