import numpy as np
import matplotlib.pyplot as plt

from qiskit_optimization import QuadraticProgram
from qiskit.opflow.primitive_ops.pauli_sum_op import PauliSumOp


class Hamiltonian:
    """Class for handling the Hamiltonian
    :param triplet_slice: slice of triplets participating in the SubQUBO
    :param normalizing: "complete"    : (sum(outer_terms_bij) + a_i) / (#outer_terms_bij + 1)
                        "outer terms" : a_i + (sum(outer_terms_bij) / #outer_terms_bij)
                        "None"        : a_i + sum(outer_terms_bij)
    """
    def __init__(self,
                 triplet_slice,
                 solution_candidate,
                 rescaling=None):
        self.triplet_slice = triplet_slice
        self.solution_candidate = solution_candidate
        self.triplet_ids = [triplet.triplet_id for triplet in triplet_slice]
        self.rescaling = rescaling

    def linear_term(self):
        """Pseudo-linear term, sums (and normalizes if chosen) all interaction values outside the
        current subqubo and adds it to the linear term.
        :return linear term: List of linear terms for the SubQUBO
        """
        linear = np.zeros(len(self.triplet_slice))
        for i, triplet in enumerate(self.triplet_slice):
            # self-term
            linear[i] += self.triplet_slice[i].quality
            # outer-term
            lin_out = 0
            lin_out_counter = 0
            for interaction_key in triplet.interactions.keys():
                if interaction_key not in self.triplet_ids:
                    if self.solution_candidate[interaction_key] == 1:
                        lin_out += triplet.interactions[interaction_key]
                        lin_out_counter += 1
            if self.normalizing == "complete":
                linear[i] += lin_out
                linear[i] /= (lin_out_counter + 1)
            elif self.normalizing == "outer terms":
                if lin_out != 0:
                    linear[i] += lin_out / lin_out_counter
            else:
                linear[i] += lin_out
        return linear

    def quadratic_term(self):
        """Returns the quadratic term of the hamiltonian."""
        quadratic = np.zeros((len(self.triplet_slice), len(self.triplet_slice)))
        for i, t1 in enumerate(self.triplet_slice):
            for j, t2 in enumerate(self.triplet_slice[i:]):
                if t2.triplet_id in t1.interactions.keys():
                    quadratic[i, j + i] = t1.interactions[t2.triplet_id]
        return quadratic

    def qubo_representation(self):
        """Returns a qubo representation of the hamiltonian suited for numpy eigensolver."""
        qubo = QuadraticProgram()
        for i in range(len(self.triplet_slice)):
            qubo.binary_var(name='x' + str(i))
        qubo.minimize(linear=self.linear_term(), quadratic=self.quadratic_term())
        return qubo

    def qubo_representation_vqe(self):
        """Returns qubo hamiltonian suited for VQE"""
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
        # returning error mitigated hamiltonian --> paper from Karl's group
        pass
