from qiskit.algorithms import VQE, QAOA
from qiskit.algorithms.optimizers import COBYLA, L_BFGS_B, SPSA, SLSQP, NFT
from qiskit.utils import QuantumInstance
from qiskit import Aer
from qiskit.providers.aer import QasmSimulator
from qiskit.providers.fake_provider import FakeAthens, FakeCasablanca, FakeJakarta, FakeGuadalupe


class Solver:
    def __init__(self,
                 config):
        """Creates a solver object for solving (Sub-) QUBOs.
        :param config: dictionary with configuration parameters
        """
        self.config = config
        self.quantum_instance = None
        self.quantum_algorithm = None

        if self.config["solver"]["backend"] == "Ideal Qasm Sim":
            self.backend = Aer.get_backend('qasm_simulator')
        elif self.config["solver"]["backend"] == "FakeAthens":
            self.backend = QasmSimulator.from_backend(FakeAthens())       # 5 qubits
        elif self.config["solver"]["backend"] == "FakeCasablanca":
            self.backend = QasmSimulator.from_backend(FakeCasablanca())   # 7 qubits
        elif self.config["solver"]["backend"] == "FakeJakarta":
            self.backend = QasmSimulator.from_backend(FakeJakarta())      # 7 qubits
        elif self.config["solver"]["backend"] == "FakeGuadalupe":
            self.backend = QasmSimulator.from_backend(FakeGuadalupe())    # 16 qubits

        self.set_quantum_instance()

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

    def set_quantum_instance(self):
        """Sets the quantum instance used by the quantum algorithms with parameters of the config file.
        """
        self.quantum_instance = QuantumInstance(backend=self.backend,
                                                seed_simulator=self.config["solver"]["seed"],
                                                seed_transpiler=self.config["solver"]["seed"],
                                                shots=self.config["solver"]["shots"],
                                                optimization_level=self.config["solver"]["optimization_level"])

    def set_vqe(self, ansatz):
        """Set vqe with a specific ansatz.
        :param ansatz: quantum circuit
        """
        vqe = VQE(ansatz=ansatz.circuit,
                  optimizer=self.get_optimizer(),
                  quantum_instance=self.quantum_instance)
        self.quantum_algorithm = vqe

    def set_qaoa(self):
        """Set QAOA.
        """
        ideal_qaoa = QAOA(optimizer=self.get_optimizer(),
                          reps=ansatz.circuit_depth,
                          max_evals_grouped=1,
                          quantum_instance=self.quantum_instance)
        self.quantum_algorithm = ideal_qaoa
