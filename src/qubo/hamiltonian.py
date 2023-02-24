import numpy as np

from qiskit_optimization import QuadraticProgram
from qiskit.opflow.primitive_ops.pauli_sum_op import PauliSumOp


class Hamiltonian:
    def __init__(self,
                 triplet_slice,
                 t_mapping,
                 solution_candidate,
                 rescaling=None,
                 only_specified_connections=None):
        """Class for handling the Hamiltonian
        :param triplet_slice: slice of triplets participating in the SubQUBO
        :param t_mapping: name of triplets consist of <hit_ID>_<hit_ID>_<hit_ID> and is mapped to its position
                          inside the triplet list
        :param solution_candidate: binary vector representing kept and discarded triplets e.g [0,1,..., 1]
        :param rescaling: "complete"    : (a_i + sum(outer_terms_bij)) / (#outer_terms_bij + 1)
                          "outer terms" :  a_i + (sum(outer_terms_bij) / #outer_terms_bij)
                          "None"        :  a_i + sum(outer_terms_bij)
        :param only_specified_connections: dictionary of triplet ids which have to be considered for accumulating
                                           relations from outside the (sub)qubo
        :param k: some mysterious value, maybe helpful, maybe not
        """
        self.only_specified_connections = only_specified_connections
        self.triplet_slice = triplet_slice
        self.t_mapping = t_mapping
        self.solution_candidate = solution_candidate
        self.triplet_ids = [triplet.triplet_id for triplet in triplet_slice]
        self.rescaling = rescaling

    def linear_term(self):
        """Pseudo-linear term, sums (and normalizes if chosen) all interaction values outside the
        current sub-QUBO and adds it to the linear term.
        :return
            list of linear terms for the Sub-QUBO
        """
        linear = np.zeros(len(self.triplet_slice))
        for i, triplet in enumerate(self.triplet_slice):
            # self-term
            linear[i] += self.triplet_slice[i].quality
            # outer-term
            lin_out_connection = []
            lin_out_conflict = []
            for interaction_key, interaction_value in triplet.interactions.items():
                if self.only_specified_connections is None:
                    if interaction_key not in self.triplet_ids:
                        if self.solution_candidate[self.t_mapping[interaction_key]] == 1:
                            if interaction_value < 0:
                                lin_out_connection.append(triplet.interactions[interaction_key])
                            else:
                                lin_out_conflict.append(triplet.interactions[interaction_key])
                else:
                    if interaction_key not in self.triplet_ids:
                        if interaction_key not in self.only_specified_connections:
                            continue
                        if self.solution_candidate[self.t_mapping[interaction_key]] == 1:
                            if interaction_value < 0:
                                lin_out_connection.append(triplet.interactions[interaction_key])
                            else:
                                lin_out_conflict.append(triplet.interactions[interaction_key])

            if self.rescaling == "complete":
                linear[i] += (sum(lin_out_connection) + sum(lin_out_conflict))
                linear[i] /= (len(lin_out_connection) + len(lin_out_conflict) + 1)
            elif self.rescaling == "outer terms":
                if lin_out_connection or lin_out_conflict:
                    linear[i] += (sum(lin_out_connection) + sum(lin_out_conflict)) / \
                                 (len(lin_out_connection) + len(lin_out_conflict) + 1)
            else:
                linear[i] += (sum(lin_out_connection) + sum(lin_out_conflict))
        return linear

    def quadratic_term(self):
        """Returns the quadratic term of the hamiltonian.
        :return
            list of quadratic terms for the Sub-QUBO
        """
        quadratic = np.zeros((len(self.triplet_slice), len(self.triplet_slice)))
        for i, t1 in enumerate(self.triplet_slice):
            for j, t2 in enumerate(self.triplet_slice[i:]):
                if t2.triplet_id in t1.interactions.keys():
                    quadratic[i, j + i] = t1.interactions[t2.triplet_id]
        return quadratic

    def qubo_representation(self):
        """Returns a qubo representation of the hamiltonian for Numpy Eigensolver.
        :return
            qubo representation for Numpy Eigensolver
        """
        qubo = QuadraticProgram()
        for i in range(len(self.triplet_slice)):
            qubo.binary_var(name='x' + str(i))
        qubo.minimize(linear=self.linear_term(), quadratic=self.quadratic_term())
        return qubo

    def qubo_representation_quantum(self):
        """Returns qubo hamiltonian representation for VQE and QAOA.
        :return
            qubo representation for VQE and QAOA
        """
        op, offset = self.qubo_representation().to_ising()
        operator_list = []
        try:
            for operator in op:
                aux_list = str(operator).split(' * ')
                t = tuple([aux_list[1], aux_list[0]])
                operator_list.append(t)
        except TypeError:
            return None
        except IndexError:
            return None
        qubo_vqe = PauliSumOp.from_list(operator_list)
        return qubo_vqe

    def error_mitigated_hamiltonian(self):
        """TODO: Implement!"""
        # returning error mitigated hamiltonian -->
        # "Measurement Error Mitigation in Quantum Computers Through Classical Bit-Flip Correction"
        pass
