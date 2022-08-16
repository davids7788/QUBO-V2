import numpy as np


class QuboLogging:
    def __init__(self):
        """Class for managing information about the solving progress of the QUBO.
        """
        self.qubo_log = {"truth solution vector": None,            # truth ground state vector
                         "truth minimum energy": None,             # truth ground state energy
                         "computed solution vector": None,         # computed vector with minimal energy
                         "computed minimum energy": None,          # computed minimal energy
                         "time tracking complete": {},             # complete time needed for program execution
                         "time tracking quantum": {},              # quantum part
                         "time tracking qubo iteration": {},       # one iteration of the global optimisation algorithm
                         "time tracking subQUBOs": {},             # time needed for solving a single subqubo
                         "time tracking bit flip search": {},      # time needed for one bit flip search
                         "compare to analytical solution": {},     # [True , False, ..., False] for each SubQUBO
                         "hamiltonian": {},                        # [linear, quadratic] for 500 SubQUBOs
                         "energy": {},                             # energy level after each iteration
                         "solution vector": {}}                    # Solution vector

    def add_entry(self,
                  sub_dict: str,
                  key,
                  value):
        """Adds an entry to the specified sub-dictionary.
        :param sub_dict : sub-dictionary of qubo log dictionary
        :param key      : key of the new entry
        :param value    : value of the new entry
        """
        self.qubo_log[sub_dict].update({key: value})

    def set_value(self,
                  key: str,
                  value: object):
        """Updates the specific value of the dictionary.
        :param key   : key of the new entry
        :param value : value of the new entry
        """
        self.qubo_log[key] = value

    def save_results(self,
                     folder: str):
        """Saves the stored information into a .npy file.
        :param folder: folder to save the results
        """
        np.save(f"{folder}/qubo_log", self.qubo_log)
