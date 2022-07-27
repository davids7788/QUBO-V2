from qiskit.circuit import Parameter, QuantumCircuit
from qiskit.circuit.library.n_local.two_local import TwoLocal


class Ansatz:
    def __init__(self,
                 config):
        """Class for managing the ansatz circuit.
        :param config: dictionary with configuration parameters
        """
        self.layout = config["ansatz"]["layout"]
        self.num_qubits = config["qubo"]["num qubits"]
        self.circuit_depth = config["ansatz"]["circuit_depth"]
        self.skip_unentangled_qubits = config["ansatz"]["skip_unentangled_qubits"]
        self.skip_final_rotation_layer = config["ansatz"]["skip_final_rotation_layer"]
        self.circuit = None

    def set_two_local(self,
                      rotation_blocks="ry",
                      entanglement_blocks="cx",
                      entanglement="full"):
        """Configures the "TwoLocal" ansatz.
        :param rotation_blocks: Single gate or list of gates possible, e.g. "ry" or ["ry", "rz"]
        :param entanglement_blocks: Entanglement block, e.g "cx", "cxx"
        :param entanglement: "linear", "full", "circular"
        """
        self.circuit = TwoLocal(num_qubits=self.num_qubits,
                                rotation_blocks=rotation_blocks,
                                entanglement_blocks=entanglement_blocks,
                                entanglement=entanglement,
                                reps=self.circuit_depth,
                                skip_unentangled_qubits=self.skip_unentangled_qubits,
                                skip_final_rotation_layer=self.skip_final_rotation_layer)

    def set_hamiltonian_driven(self,
                               triplet_list_slice):
        """Creates entanglements of the ansatz based on the hamiltonian and sets it as the circuit of the ansatz.
        Direct entanglement means if there is a direct connection, a value b_ij which would connect both.
        Currently restricted to "ry" and "cx" gates.
        :param triplet_list_slice: triplet list slice of the problem
        """
        qc = QuantumCircuit(self.num_qubits)
        # Number of repetitions
        entanglement_map = []
        for i in range(self.circuit_depth):
            for j in range(self.num_qubits):
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
        if not self.skip_final_rotation_layer:
            for n in range(self.num_qubits):
                qc.ry(Parameter(str(self.circuit_depth) * self.num_qubits + str(n)), n)  # Final rotation layer
        self.circuit = qc

    def set_no_entanglements(self):
        """Set quantum circuits with just rotation layers and not entanglements."""
        qc = QuantumCircuit(self.num_qubits)
        for i in range(self.circuit_depth):
            for j in range(self.num_qubits):
                qc.ry(Parameter(str(i) + str(j)), j)  # Parametrized ry gates on every qubit
        # final rotation layer
        if not self.skip_final_rotation_layer:
            for n in range(self.num_qubits):
                qc.ry(Parameter(str(self.circuit_depth) * self.num_qubits + str(n)), n)

        self.circuit = qc
