import copy
import time
import numpy as np
import matplotlib.pyplot as plt

from qiskit.utils import algorithm_globals
from qiskit.algorithms import NumPyMinimumEigensolver
from qiskit_optimization.algorithms import MinimumEigenOptimizer

from src.solver import Solver
from src.hamiltonian import Hamiltonian
from src.ansatz import Ansatz
from src.qubo_logging import QuboLogging
from src.track import Track
from src.error_mitigation import ErrorMitigation

algorithm_globals.massive = True


class QuboProcessing:
    def __init__(self,
                 triplet_list_file: str,
                 config: dict,
                 solver: Solver,
                 ansatz: Ansatz,
                 qubo_logging: QuboLogging,
                 save_folder: str,
                 error_mitigation: ErrorMitigation):
        """Processes the Qubo and provides solving functions like a solving process via an impact list.
        :param triplet_list_file: .npy file with triplet objects
        :param config: dictionary with configuration parameters
        :param solver: solver object
        :param ansatz: ansatz circuit
        :param qubo_logging: object for handling information about the solving process of class QuboLogging
        :param save_folder: folder to which results are stored
        :param error_mitigation: class to help mitigate errors
        """
        self.triplets = np.load(triplet_list_file, allow_pickle=True)
        self.config = config
        self.solver = solver
        self.ansatz = ansatz
        if self.solver is not None:
            self.configure_solver()
        self.qubo_logging = qubo_logging
        self.error_mitigation = error_mitigation
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

        # Performing initial bit flip optimization
        if self.config["bit flip optimization"]["iterations"] > 0:
            start_bit_flip = time.time()
            for i in range(self.config["bit flip optimization"]["iterations"]):
                new_solution_candidate, energy_change = self.bit_flip_optimization(self.triplets,
                                                                                   self.solution_candidate,
                                                                                   self.make_impact_list(),
                                                                                   self.config["bit flip optimization"]
                                                                                   ["reverse"])
                self.energy_candidate += energy_change
                self.qubo_logging.add_entry("solution vector", self.iteration, self.solution_candidate)
                self.qubo_logging.add_entry("energy", self.iteration, self.energy_candidate)
                self.iteration += 1
                print(f"Energy after performing {i + 1} bit flip optimization(s): "
                      f"{np.around(self.energy_candidate, 2)}")
            if self.config["qubo"]["search depth"] == 0:
                self.create_tracks()
                self.qubo_logging.save_results(self.save_folder)
                self.write_results_to_file()
            end_bit_flip = time.time()
            print(f"Bit flip optimization needed {QuboProcessing.hms_string((end_bit_flip - start_bit_flip))}")

        # Here starts the sub qubo algorithm
        self.loop_count = 0
        self.pass_count = 0
        self.sub_qubo_counter = 0

    def impact_list_solve(self):
        """Performs impact list algorithm
        """
        self.qubo_process(optimization_strategy=self.make_impact_list)
        self.write_results_to_file()

    def qubo_process(self,
                     optimization_strategy):
        """Solves the QUBO with
        :param optimization_strategy: function which is used to determine triplet ordering
        """
        # timer
        start_solving_process_time = time.time()

        # global optimization iterations
        while self.pass_count < self.config["qubo"]["search depth"]:
            start_loop = time.time()
            start_quantum_part = time.time()

            # triplet list ordering determines also how many subqubos are built within one iteration
            # so it is possible to define a function with repeating indices resulting in a longer triplet ordering list
            # e.g. [1, 5, 3, 2, 4, 0], but also [1, 5, 1, 5, 3, 2, 4, 1, 5, ...]
            triplet_ordering = optimization_strategy()

            # Calculate number of sub qubos:
            if len(triplet_ordering) % self.config["qubo"]["num qubits"] == 0:
                num_sub_qubos = int(len(triplet_ordering) / self.config["qubo"]["num qubits"])
            else:
                num_sub_qubos = int(len(triplet_ordering) / self.config["qubo"]["num qubits"]) + 1

            # Sub-qubos with full size
            for i in range(int(len(self.triplets) / self.config["qubo"]["num qubits"])):
                print(f"Processing SubQUBO {i + 1} of {num_sub_qubos}", end="\r")
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
                print(f"Processing SubQUBO {num_sub_qubos} of {num_sub_qubos}", end="\r")
                result = self.solve_subqubos([self.triplets[i] for i in
                                              triplet_ordering[- self.config["qubo"]["num qubits"]:]])
                if result is None:
                    pass
                else:
                    for k, entry in enumerate(triplet_ordering[- self.config["qubo"]["num qubits"]:]):
                        self.solution_candidate[entry] = result[k]

            # Time tracking quantum part
            end_quantum_part = time.time()
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

            print(f"Energy after iteration {self.loop_count} at pass count {self.pass_count}: "
                  f"{np.around(self.energy_candidate, 2)}")

            # increasing overall loop count
            self.loop_count += 1

            # time tracking iteration
            end_loop = time.time()
            loop_time = end_loop - start_loop
            self.qubo_logging.add_entry("time tracking qubo iteration",
                                        str(self.loop_count),
                                        loop_time)

            self.qubo_logging.save_results(self.save_folder)

        end_solving_process_time = time.time()
        solving_process_time = end_solving_process_time - start_solving_process_time
        self.qubo_logging.add_entry("time tracking complete",
                                    "complete run",
                                    solving_process_time)
        self.create_tracks()
        self.qubo_logging.save_results(self.save_folder)
        print(f"Qubo solving process needed {QuboProcessing.hms_string(solving_process_time)}")

    def make_impact_list(self):
        """
        Creates an impact list based on how much influence on the energy a bit flip has
        :return:
            list of indices ordered from lowest to highest impact of triplets in triplet list
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
                              triplet_ordering,
                              reverse=True):

        """Looping over a list of triplet objects and a corresponding binary solution vector to compute if the energy
        value would improve if the binary state of a triplet (keep <-> discard) should be changed. If the energy
        decreases, the triplet state is flipped.
        :param triplets: list of triplet objects
        :param solution_candidate binary vector representing triplet states
        :param triplet_ordering: index list of triplets
        :param reverse: False if provided sorting order, else reversed
        """
        if reverse:
            triplet_ordering.reverse()
        energy_change_total = 0
        for i, triplet in enumerate([triplets[index] for index in triplet_ordering]):
            # energy change if this particular bit is flipped
            energy_change = 0

            # checking linear term
            if solution_candidate[triplet_ordering[i]] == 0:
                energy_change += triplet.quality
            else:
                energy_change -= triplet.quality

            # Checking interactions with other triplets
            for interaction in triplet.interactions.keys():
                if solution_candidate[triplet_ordering[i]] == 0 and solution_candidate[interaction] == 0:
                    pass
                elif solution_candidate[triplet_ordering[i]] == 0 and solution_candidate[interaction] == 1:
                    energy_change += triplet.interactions[interaction]
                elif solution_candidate[triplet_ordering[i]] == 1 and solution_candidate[interaction] == 0:
                    pass
                else:
                    energy_change -= triplet.interactions[interaction]

            # flip if overall energy change is negative
            if energy_change < 0:
                solution_candidate[triplet_ordering[i]] = 1 - solution_candidate[triplet_ordering[i]]
                energy_change_total += energy_change

        return solution_candidate, energy_change_total

    def minimal_energy_and_solution(self):
        """Calculates the minimum energy state and value.
        :return:
            minimum energy state, minimum energy value """
        minimum_energy_state = []
        for t in self.triplets:
            if t.is_correct_match:
                minimum_energy_state.append(1)
            else:
                minimum_energy_state.append(0)
        return minimum_energy_state, self.hamiltonian_energy(minimum_energy_state)

    def hamiltonian_energy(self, binary_vector):
        """
        Calculates the energy according to a binary vector matching the triplet list
        :param binary_vector: binary solution candidate vector
        :return:
            energy value
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
        result_quantum = self.solver.quantum_algorithm.compute_minimum_eigenvalue(hamiltonian)
        if self.config["qubo"]["error mitigation algorithm"] == "algebraic":
            result_quantum_vector = ErrorMitigation.dict_to_vector(dict(result_quantum.eigenstate.items()))
            result_quantum_vector_mitigated = np.dot(self.error_mitigation.meas_filter_matrix, result_quantum_vector)
            result_quantum = ErrorMitigation.vector_to_dict(result_quantum_vector_mitigated)

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
                       triplet_slice):
        """Function for solving a SubQUBO
        :param triplet_slice: slice of triplets used for the SubQUBO
        :return
            result of computation
        """
        result = None
        sub_qubo_time = None

        # start timer sub-qubo creation and solving

        hamiltonian = Hamiltonian(triplet_slice=triplet_slice,
                                  solution_candidate=self.solution_candidate,
                                  rescaling=self.config["qubo"]["hamiltonian rescaling"])
        if self.loop_count == 0:  # Uses a lot of disk space, so just results of first iteration
            self.qubo_logging.add_entry("hamiltonian",
                                        self.sub_qubo_counter,
                                        [hamiltonian.linear_term(), hamiltonian.quadratic_term()])

        # Create hamiltonian and check if no linear and quadratic terms --> do nothing
        if hamiltonian.qubo_representation_quantum() is None:
            return

        if self.config["solver"]["algorithm"] == "VQE" or self.config["solver"]["algorithm"] == "QAOA":

            # Overwrite ansatz for vqe if HamiltonianDriven
            if self.ansatz.config["ansatz"]["layout"] == "HamiltonianDriven":
                triplet_slice = triplet_slice
                self.ansatz = self.ansatz.hamiltonian_driven(triplet_slice)
                self.solver.quantum_instance.ansatz = self.ansatz

            # Timestamps and solving the SubQUBO
            start_sub_qubo = time.time()
            result = self.solve_quantum(hamiltonian=hamiltonian.qubo_representation_quantum())
            end_sub_qubo = time.time()
            sub_qubo_time = end_sub_qubo - start_sub_qubo

            # Check if comparing to Eigensolver
            if self.config["qubo"]["compare to eigensolver"] and self.loop_count == 0:  # one iteration
                exact_result = self.exact_solver.solve(hamiltonian.qubo_representation()).x
                self.qubo_logging.add_entry("compare to analytical solution",
                                            self.sub_qubo_counter,
                                            ([float(r) for r in result] == exact_result))

        elif self.config["solver"]["algorithm"] == "Numpy Eigensolver":
            start_sub_qubo = time.time()
            result = self.exact_solver.solve(hamiltonian.qubo_representation()).x
            end_sub_qubo = time.time()
            sub_qubo_time = end_sub_qubo - start_sub_qubo

        # update solution

        self.qubo_logging.add_entry("time tracking subQUBOs",
                                    self.sub_qubo_counter,
                                    sub_qubo_time)

        self.sub_qubo_counter += 1
        return result

    def configure_solver(self):
        """Configures the Solver object for VQE, QAOA or NumpyEigensolver.
        """
        # Set up the the chosen algorithm, VQE, QAOA or NumpyEigensolver
        # VQE --> ansatz can be chosen
        if self.config["solver"]["algorithm"] == "VQE":
            self.solver.set_vqe(self.ansatz)

        # QAOA --> ansatz determined by algorithm
        elif self.config["solver"]["algorithm"] == "QAOA":
            self.solver.set_qaoa()

        # Numpy Eigensolver does not need further preparation
        elif self.config["algorithm"] == "Numpy Eigensolver":
            pass
        else:
            "No valid algorithm chosen!"

    @staticmethod
    def hms_string(sec_elapsed):
        """Nicely formatted time string.
        """
        h = int(sec_elapsed / (60 * 60))
        m = int((sec_elapsed % (60 * 60)) / 60)
        s = sec_elapsed % 60
        return "{}:{:>02}:{:>05.2f}".format(h, m, s)

    def write_results_to_file(self):
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

    def create_tracks(self):
        """"Creates reco track list (for after preselection) and computed track list.
        """
        # Create all the correct tracks as reference
        truth_solution = self.qubo_logging.qubo_log["truth solution vector"]
        reco_tracks = []
        for index, result in enumerate(truth_solution):
            if result == 1:
                for key, value in zip(self.triplets[index].interactions.keys(),
                                      self.triplets[index].interactions.values()):
                    if int(key) < int(self.triplets[index].triplet_id):  # to not create tracks two times
                        if truth_solution[key] == 1:
                            new_track = Track()
                            new_track.add_triplet_to_track(self.triplets[index])
                            new_track.add_triplet_to_track(self.triplets[key])
                            reco_tracks.append(new_track)

        self.qubo_logging.set_value("max reco tracks", reco_tracks)

        solution = self.qubo_logging.qubo_log["computed solution vector"]
        found_tracks = []
        used_triplets_for_tracks = set()
        for index, result in enumerate(solution):  # looping over triplet list
            if result == 1:  # check if triplet is kept
                best_interaction_key = None  # variable for best connection if ambiguities
                best_value = 0  # variable for best connection value if ambiguities
                for key, value in zip(self.triplets[index].interactions.keys(),
                                      self.triplets[index].interactions.values()):
                    if int(key) < int(self.triplets[index].triplet_id):  # to not create tracks two times
                        if solution[key] == 1:  # only if possible partner triplet is kept, too
                            if value < best_value:
                                best_value = value
                                best_interaction_key = key

                if best_value < 0:  # if partner triplet was found
                    new_track = Track()
                    new_track.add_triplet_to_track(self.triplets[index])
                    new_track.add_triplet_to_track(self.triplets[best_interaction_key])
                    found_tracks.append(new_track)
                    used_triplets_for_tracks.add(index)
                    used_triplets_for_tracks.add(best_interaction_key)

        found_tracks_ambiguity_solved = []
        for i, track_1 in enumerate(found_tracks):
            ambiguities = False
            for j, track_2 in enumerate(found_tracks[i + 1:]):
                if track_1.is_in_conflict(track_2):
                    ambiguities = True
                    sum_interactions_track_1 = 0
                    sum_interactions_track_2 = 0
                    for triplet_1, triplet_2 in zip(track_1.triplets[:-1], track_1.triplets[1:]):
                        sum_interactions_track_1 += triplet_1.interactions[triplet_2.triplet_id]
                    for triplet_1, triplet_2 in zip(track_2.triplets[:-1], track_2.triplets[1:]):
                        sum_interactions_track_2 += triplet_1.interactions[triplet_2.triplet_id]
                    if sum_interactions_track_2 > sum_interactions_track_1:
                        found_tracks_ambiguity_solved.append(track_2)
                    else:
                        found_tracks_ambiguity_solved.append(track_1)
            if not ambiguities:
                found_tracks_ambiguity_solved.append(track_1)

            self.qubo_logging.set_value("reconstructed tracks", found_tracks_ambiguity_solved)
