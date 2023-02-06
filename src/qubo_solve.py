import sys
import os
import yaml
import numpy as np

from pathlib import Path

from qubo.qubo_processing import QuboProcessing
from qubo.qubo_logging import QuboLogging
from qubo.ansatz import Ansatz
from qubo.solver import Solver
from track_reconstruction.create_reco_xplets import reco_xplets_simplified_LUXE
from track_reconstruction.track_reconstruction_efficiency import track_reconstruction_efficiency_simplified_LUXE


# sys argv [1]: config file
# sys argv [2]: folder containing a .npy triplet list file

print("\n-----------------------------------")
print("\nStarting to solve the QUBO...\n")

with open(sys.argv[1], 'r') as f:
    config_file = yaml.safe_load(f)

# Create new folder
folder = sys.argv[2]

file_extension = "_" + sys.argv[1].split("/")[-1].split(".")[0]

new_folder = folder + "/" + str(np.random.randint(1e8, 1e9)) + file_extension
if Path(new_folder).is_dir():
    pass
else:
    print(f"Creating folder: {new_folder}\n")
    os.mkdir(new_folder)

# Create logger
qubo_logger = QuboLogging()

# Create ansatz object and set parameters from config file
if config_file["ansatz"]["layout"] is None:
    ansatz = None
else:
    ansatz = Ansatz(config=config_file)

# If not "hamiltonian driven" additional parameters are set
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

# Create and configure solving process
qubo_processor = QuboProcessing(folder + "/triplet_list.npy",
                                config=config_file,
                                solver=solver,
                                ansatz=ansatz,
                                qubo_logging=qubo_logger,
                                save_folder=new_folder,
                                verbose=1)


# Select solving method
if "impact list" in config_file["qubo"]["optimisation strategy"]:
    qubo_processor.qubo_processing()
if "connection list" in config_file["qubo"]["optimisation strategy"]:
    qubo_processor.qubo_processing()
if "paired list" in config_file["qubo"]["optimisation strategy"]:
    qubo_processor.qubo_processing()
if "impact without conflicts" in config_file["qubo"]["optimisation strategy"]:
    qubo_processor.qubo_processing()

reco_xplets_simplified_LUXE(qubo_processor.get_kept_triplets(),
                            new_folder)

track_reconstruction_efficiency_simplified_LUXE(f"{new_folder}/reco_xplet_list.npy")

print("-----------------------------------\n")
print("QUBO solved successfully!\n")
