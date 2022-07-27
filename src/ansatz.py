from qiskit.circuit import Parameter, QuantumCircuit
from qiskit.circuit.library.n_local.two_local import TwoLocal


class Ansatz:
    def __init__(self,
                 config):
        """Class for managing the ansatz circuit.
        :param config: dictionary with configuration parameters
        """
        self.config = config
        self.circuit = None

    def set_two_local(self):
        """Configures the "TwoLocal" ansatz.
        """
        self.circuit = TwoLocal(num_qubits=self.config["qubo"]["num qubits"],
                                rotation_blocks=self.config["ansatz"]["rotation blocks"],
                                entanglement_blocks=self.config["ansatz"]["entanglement blocks"],
                                entanglement=self.config["ansatz"]["entanglement"],
                                reps=self.config["ansatz"]["circuit depth"],
                                skip_unentangled_qubits=self.config["ansatz"]["skip unentangled qubits"],
                                skip_final_rotation_layer=self.config["ansatz"]["skip final rotation layer"])

    def set_hamiltonian_driven(self,
                               triplet_list_slice):
        """Creates entanglements of the ansatz based on the hamiltonian and sets it as the circuit of the ansatz.
        Direct entanglement means if there is a direct connection, a value b_ij which would connect both.
        Currently restricted to "ry" and "cx" gates.
        :param triplet_list_slice: triplet list slice of the problem
        """
        qc = QuantumCircuit(self.config["qubo"]["num qubits"])
        # Number of repetitions
        entanglement_map = []
        for i in range(self.config["ansatz"]["circuit depth"]):
            for j in range(self.config["qubo"]["num qubits"]):
                qc.ry(Parameter(str(i) + str(j)), j)  # Parametrized ry gates on every qubit
            for k, p1 in enumerate(triplet_list_slice):
                for m, p2 in enumerate(triplet_list_slice):
                    if p2.triplet_id in p1.interactions.keys():
                        if (p1, p2) in entanglement_map or (p2, p1) in entanglement_map:
                            continue
                        else:
                            qc.cx(k, m)  # if b_ij != 0 for the qubits/triplets, set entanglement cx
                            entanglement_map.append((p1, p2))
        # final rotation layer
        if not self.config["ansatz"]["skip final rotation layer"]:
            for n in range(self.config["qubo"]["num qubits"]):
                qc.ry(Parameter(str(self.config["ansatz"]["circuit depth"]) *
                                self.config["qubo"]["num qubits"] + str(n)), n)   # Final rotation layer
        self.circuit = qc

    def set_no_entanglements(self):
        """Set quantum circuits with just rotation layers and not entanglements."""
        qc = QuantumCircuit(self.config["qubo"]["num qubits"])
        for i in range(self.config["ansatz"]["circuit depth"]):
            for j in range(self.config["qubo"]["num qubits"]):
                qc.ry(Parameter(str(i) + str(j)), j)  # Parametrized ry gates on every qubit
        self.circuit = qc
