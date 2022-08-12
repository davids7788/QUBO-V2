import sys
import os
import yaml
import numpy as np

from pathlib import Path

from qubo.error_mitigation import ErrorMitigation
from qubo.qubo_processing import QuboProcessing
from qubo.qubo_logging import QuboLogging
from qubo.ansatz import Ansatz
from qubo.solver import Solver


# sys argv [1]: config file
# sys argv [2]: folder containing a .npy triplet list file

with open(sys.argv[1], 'r') as f:
    config_file = yaml.safe_load(f)

# Create new folder
folder = sys.argv[2]

file_extension = ""
if config_file["solver"]["algorithm"] == "Numpy Eigensolver":
    file_extension += "_eigensolver"
elif config_file["solver"]["algorithm"] == "VQE":
    file_extension += "_vqe"
elif config_file["solver"]["algorithm"] == "QAOA":
    file_extension += "_qaoa"
else:
    file_extension += "_bit_flip_only"


new_folder = folder + "/" + str(np.random.randint(1e8, 1e9)) + file_extension
if Path(new_folder).is_dir():
    pass
else:
    print(f"Creating folder: {new_folder}")
    os.mkdir(new_folder)

# Create logger
qubo_logger = QuboLogging()

# Create ansatz object and set parameters from config file
if config_file["ansatz"]["layout"] is None:
    ansatz = None
else:
    ansatz = Ansatz(config=config_file)

# If not "hamiltonian aware" additional parameters are set
if config_file["ansatz"]["layout"] == "TwoLocal":
    ansatz.set_two_local()
elif config_file["ansatz"]["layout"] is None:
    if config_file["solver"]["algorithm"] == "QAOA":
        pass
    elif config_file["solver"]["algorithm"] is not None:
        if config_file["solver"]["algorithm"] != "Numpy Eigensolver":
            ansatz.set_no_entanglements()

# Create Solver object and set parameters from config file
if config_file["solver"]["algorithm"] != "Numpy Eigensolver" and config_file["solver"]["algorithm"] is not None:
    solver = Solver(config_file)
else:
    solver = None

if config_file["qubo"]["error mitigation algorithm"] is not None:
    error_mitigation = ErrorMitigation(solver.get_backend())
    error_mitigation.algebraic_mitigation(ansatz.circuit)
else:
    error_mitigation = None

# Create and configure solving process
qubo_processor = QuboProcessing(folder + "/triplet_list.npy",
                                config=config_file,
                                solver=solver,
                                ansatz=ansatz,
                                qubo_logging=qubo_logger,
                                save_folder=new_folder,
                                error_mitigation=error_mitigation)

# Select solving method
if config_file["qubo"]["optimisation strategy"] == "impact list":
    qubo_processor.qubo_process()
