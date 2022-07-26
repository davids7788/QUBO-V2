from qiskit.algorithms import VQE, QAOA
from qiskit.algorithms.optimizers import COBYLA, L_BFGS_B, SPSA, SLSQP, NFT
from qiskit.utils import QuantumInstance
from qiskit import Aer


class Solver:
    def __init__(self,
                 config):
        """Creates a solver object for solving (Sub-) QUBOs.
        :param config: dictionary with configuration parameters
        """
        self.model = None
        self.quantum_instance = None
        self.parameters = {}
        self.parameters.update(kwargs)

    def get_optimizer(self):
        """Returns the optimizer with the maxiter value created from the class attributes.
        :return
            Optimizer configured according to the set parameter in the parameter dictionary.
        """
        if self.parameters["optimizer"] == "COBYLA":
            return COBYLA(maxiter=self.parameters["maxiter"])
        if self.parameters["optimizer"] == "L_BFGS_B":
            return L_BFGS_B(maxiter=self.parameters["maxiter"])
        if self.parameters["optimizer"] == "SPSA":
            return SPSA(maxiter=self.parameters["maxiter"])
        if self.parameters["optimizer"] == "SLSQP":
            return SLSQP(maxiter=self.parameters["maxiter"])
        if self.parameters["optimizer"] == "NFT":
            return NFT(maxiter=self.parameters["maxiter"])
        if self.parameters["optimizer"] == "NELDER_MEAD":
            return NELDER_MEAD(maxiter=self.parameters["maxiter"])
        if self.parameters["optimizer"] == "ADAM":
            return ADAM(maxiter=self.parameters["maxiter"])
        if self.parameters["optimizer"] == "GSLS":
            return GSLS(maxiter=self.parameters["maxiter"])

    def set_vqe_ideal_qasm(self, ansatz):
        """
        Set vqe with a specific ansatz.
        :param ansatz: quantum circuit
        """
        self.quantum_instance = QuantumInstance(backend=Aer.get_backend('qasm_simulator'),
                                                seed_simulator=self.parameters["seed"],
                                                seed_transpiler=self.parameters["seed"],
                                                shots=self.parameters["shots"],
                                                optimization_level=self.parameters["optimization_level"])
        vqe = VQE(ansatz=ansatz.circuit,
                  optimizer=self.get_optimizer(),
                  quantum_instance=self.quantum_instance)
        self.model = vqe

    def set_vqe_ideal_qaoa(self, ansatz):
        self.quantum_instance = QuantumInstance(backend=Aer.get_backend('qasm_simulator'),
                                                seed_simulator=self.parameters["seed"],
                                                seed_transpiler=self.parameters["seed"],
                                                shots=self.parameters["shots"],
                                                optimization_level=self.parameters["optimization_level"])
        ideal_qaoa = QAOA(optimizer=self.get_optimizer(),
                          reps=ansatz.circuit_depth,
                          max_evals_grouped=1,
                          quantum_instance=self.quantum_instance)
        self.model = ideal_qaoa