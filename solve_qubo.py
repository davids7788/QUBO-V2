import sys
import os
import yaml
from pathlib import Path

from qubo_processing import QuboProcessing
from qubo_logging import QuboLogging
from ansatz import Ansatz
from solver import Solver


# sys argv [1]: config file
# sys argv [2]: folder with .npy file containing a triplet list

with open(sys.argv[1], 'r') as f:
    config_file = yaml.safe_load(f)

# Create new folder
folder = sys.argv[2]

file_extension = ""
if config_file["solver"]["algorithm"] == "Numpy Eigensolver":
    file_extension += "_numpy-eigensolver"
elif config_file["solver"]["algorithm"] == "VQE":
    file_extension += "_vqe"
elif config_file["solver"]["algorithm"] == "QAOA":
    file_extension += "_qaoa"


new_folder = folder + "/" + str(random.randint(1e8, 1e9)) + file_extension
if Path(new_folder).is_dir():
    pass
else:
    print(new_folder)
    os.mkdir(new_folder)

with open(new_folder + "/qubo_info.txt", "w") as f:
    f.write("Qubo solved with the following configuration: \n")
    f.write("---\n")
    for outer_key in self.configuration.keys():
        for inner_key, value in self.configuration[outer_key].items():
            f.write(f"\n{inner_key}: {value}")
        f.write("\n")
    f.write("\n\n")
    f.write("---\n")

# Create logger
qubo_logger = QuboLogging()

# Create ansatz object and set parameters from config file
if config_file["ansatz"]["name"] is None:
    ansatz = None
else:
    ansatz = Ansatz(config=config_file)

# If not "hamiltonian aware" additional parameters are set
if config_file["ansatz"]["layout"] == "TwoLocal":
    ansatz.set_two_local(config_file["ansatz"])

# Create Solver object and set parameters from config file
solver = Solver(config["solver"])


# Create and configure solving process
qubo_solver = QuboProcessing(folder + "/triplet_list.npy",
                             config=config_file,
                             solver=solver,
                             ansatz=ansatz,
                             qubo_logging=qubo_logger,
                             save_folder=new_folder + "/qubo_results")


# Select solving method
if config_file["qubo"]["optimization strategy"] == "impact list":
    qubo_solver.impact_list_solve(config)
