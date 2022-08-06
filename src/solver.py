from qiskit.algorithms import VQE, QAOA
from qiskit.algorithms.optimizers import COBYLA, L_BFGS_B, SPSA, SLSQP, NFT, NELDER_MEAD, ADAM, GSLS
from qiskit.utils import QuantumInstance
from qiskit import Aer
from qiskit.providers.aer import QasmSimulator
from qiskit.providers.fake_provider import FakeAthens, FakeCasablanca, FakeJakarta, FakeGuadalupe


class Solver:
    def __init__(self,
                 config):
        """Class for handling a solver object for solving (Sub-) QUBOs.
        :param config: dictionary with configuration parameters
        """
        self.config = config
        self.quantum_instance = None
        self.quantum_algorithm = None
        self.set_quantum_instance()

    def get_backend(self):
        """Sets the backend set in the config file.
        """
        if self.config["solver"]["backend"] == "Qasm Ideal Sim":
            return Aer.get_backend('qasm_simulator')                # noiseless simulator
        elif self.config["solver"]["backend"] == "FakeAthens":
            return QasmSimulator.from_backend(FakeAthens())         # 5 qubits
        elif self.config["solver"]["backend"] == "FakeCasablanca":
            return QasmSimulator.from_backend(FakeCasablanca())     # 7 qubits
        elif self.config["solver"]["backend"] == "FakeJakarta":
            return QasmSimulator.from_backend(FakeJakarta())        # 7 qubits
        elif self.config["solver"]["backend"] == "FakeGuadalupe":
            return QasmSimulator.from_backend(FakeGuadalupe())      # 16 qubits

    def get_optimiser(self):
        """Returns the optimiser with the maxiter value from the config file.
        :return Optimiser: the optimiser with the maxiter value from the config file
        """
        if self.config["solver"]["optimiser"] == "COBYLA":
            return COBYLA(maxiter=self.config["solver"]["maxiter"])
        if self.config["solver"]["optimiser"] == "L_BFGS_B":
            return L_BFGS_B(maxiter=self.config["solver"]["maxiter"])
        if self.config["solver"]["optimiser"] == "SPSA":
            return SPSA(maxiter=self.config["solver"]["maxiter"])
        if self.config["solver"]["optimiser"] == "SLSQP":
            return SLSQP(maxiter=self.config["solver"]["maxiter"])
        if self.config["solver"]["optimiser"] == "NFT":
            return NFT(maxiter=self.config["solver"]["maxiter"])
        if self.config["solver"]["optimiser"] == "NELDER_MEAD":
            return NELDER_MEAD(maxiter=self.config["solver"]["maxiter"])
        if self.config["solver"]["optimiser"] == "ADAM":
            return ADAM(maxiter=self.config["solver"]["maxiter"])
        if self.config["solver"]["optimiser"] == "GSLS":
            return GSLS(maxiter=self.config["solver"]["maxiter"])

    def set_quantum_instance(self):
        """Sets the quantum instance used by the quantum algorithms with parameters of the config file.
        """
        self.quantum_instance = QuantumInstance(backend=self.get_backend(),
                                                seed_simulator=self.config["solver"]["seed"],
                                                seed_transpiler=self.config["solver"]["seed"],
                                                shots=self.config["solver"]["shots"],
                                                optimization_level=self.config["solver"]["optimisation level"])

    def set_vqe(self, ansatz):
        """Set VQE with a specific ansatz.
        :param ansatz: quantum circuit
        """
        vqe = VQE(ansatz=ansatz.circuit,
                  optimizer=self.get_optimiser(),
                  quantum_instance=self.quantum_instance)
        self.quantum_algorithm = vqe

    def set_qaoa(self):
        """Set QAOA.
        """
        ideal_qaoa = QAOA(optimizer=self.get_optimiser(),
                          reps=self.config["ansatz"]["circuit depth"],
                          quantum_instance=self.quantum_instance)
        self.quantum_algorithm = ideal_qaoa
