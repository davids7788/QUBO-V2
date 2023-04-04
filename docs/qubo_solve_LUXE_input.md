# Configuration file for `solve_qubo.py`
The `qubo_solve.py` script takes a YAML configuration file as first input argument. 
It consists of the sections [ansatz](#ansatz), [solver](#solver), [qubo](#qubo) and 
[bit flip optimisation](#bit%20flip%20optimisation).

# ansatz
Defines the parameters of the quantum circuit used for the quantum part of the solving algorithm.
* `layout:` "TwoLocal" (see qiskit documentation), "HamiltonianDriven" (only entangling if triplets have an 
  immediate connection) or null (either when using NumpyEigensolver or when not wanting any entanglements)
* `circuit depth:` number of repetition blocks (int)
* `rotation blocks:` single string that represents a gate, or a list of strings of parametrized gates, e.g "ry" or ["ry", "rx"]
* `entanglement blocks:` multi qubit gates, such as "cx" or "cxx" 
* `entanglement:` how the entanglement structure is realised, e.g "linear", "full", "circular" 
* `skip final rotation layer:` True if rotation block should be skipped at the end of the circuit, else False
* `skip unentangled qubits:` True, if the single qubit gates are only applied to qubits that are entangled, else False


# solver
Parameters for the solver, the used algorithm and backend.
* `algorithm:` "NumpyEigensolver", "VQE" or "QAOA" or null (for using bit flip optimisation only)
* `backend:` "Qasm Ideal Sim", "FakeAthens" (5Q), "FakeCasablanca" (7Q), "FakeJakarta" (7Q) and "FakeGuadalupe" (16Q)
* `optimiser:` "ADAM", "COBYLA", "GSLS", "L_BFGS_B", NELDER_MEAD", "NFT", "SLSQP" or "SPSA" 
* `maxiter:` max number of function evaluations (int)
* `maxfev` max number of iterations of the optimising step, NFT only (int)
* `seed:` some number for setting a defined seed (int)
* `shots:` number of evaluations of the quantum circuit (int)
* `optimisation level:` 0 (no optimisation), 1 (light optimisation), 2 (heavy optimisation) or 3 
  (even heavier optimisation)  
  
# qubo
Additional parameters for solving the QUBO, such as initialisation, optimisation strategy, noise mitigation and 
comparing to analytical solution.
* `optimisation strategy:` "impact list" or "connection list" (more added in the future)
* `hamiltonian rescaling:` When using the SubQUBO approach, contributions from outside the SubQUBO are added to the 
  linear term of the participants in the SubQUBO. If "complete" then all the values, including the $a_i$ term from the SubQUBO
  are divided by the number of all the participating outer contributions + 1, if "outer term", then it is only divided by the number
  of outside contributions, else no rescaling is done
* `initial binary vector:` initial solution vector, "ones" for [1, 1, ..., 1], "zeros for" [0, 0, ..., 0], else random
* `search depth:` number of iterations the algorithm is allowed not to find a better solution in the context of a 
  minimum energy (int)
* `num qubits:` number of qubits of the SubQUBO (int)
* `error mitigation algorithm:` "algebraic" for error mitigation via matrix inverse of a calibration measurement (adding method of https://arxiv.org/pdf/2007.03663.pdf soon)
* `compare to eigensolver:` True if calculation on SubQUBO level is done with Eigensolver too

# bit flip optimisation
Parameters for the bit flip optimisation 
* `iterations:` num of iterations (int)
* `reverse:` bit flip optimisation is done in order of the impact list, if True, then of reversed impact list, 
  else False (other optimisation algorithms will be added in the future)
