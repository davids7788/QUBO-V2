import argparse
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

parser = argparse.ArgumentParser(description='QUBO solve Muon Collider',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--config_file',
                    action='store',
                    type=str,
                    default=None,
                    help='QUBO configuration file')

parser.add_argument('--qubo_folder',
                    action='store',
                    type=str,
                    default=None,
                    help='Folder with a triplet_list.npy file')


parser_args = parser.parse_args()
config_file = parser_args.config_file
qubo_folder = parser_args.qubo_folder


print("\n-----------------------------------")
print("\nStarting to solve the QUBO...\n")

with open(config_file, 'r') as f:
    configuration = yaml.safe_load(f)

file_extension = config_file.split("/")[-1].split(".")[0]

new_folder = qubo_folder + "/" + file_extension
if Path(new_folder).is_dir():
    pass
else:
    print(f"Creating folder: {new_folder}\n")
    os.mkdir(new_folder)

# Create logger
qubo_logger = QuboLogging()

# Create ansatz object and set parameters from config file
if configuration["ansatz"]["layout"] is None:
    ansatz = None
else:
    ansatz = Ansatz(config=config_file)

# If not "hamiltonian driven" additional parameters are set
if configuration["ansatz"]["layout"] == "TwoLocal":
    ansatz.set_two_local()
elif configuration["ansatz"]["layout"] is None:
    if configuration["solver"]["algorithm"] == "QAOA":
        pass
    elif configuration["solver"]["algorithm"] is not None:
        if configuration["solver"]["algorithm"] != "Numpy Eigensolver":
            ansatz.set_no_entanglements()

# Create Solver object and set parameters from config file
if configuration["solver"]["algorithm"] != "Numpy Eigensolver" and configuration["solver"]["algorithm"] is not None:
    solver = Solver(configuration)
else:
    solver = None

# Create and configure solving process
qubo_processor = QuboProcessing(qubo_folder + "/triplet_list.npy",
                                config=configuration,
                                solver=solver,
                                ansatz=ansatz,
                                qubo_logging=qubo_logger,
                                save_folder=new_folder,
                                verbose=1)
qubo_processor.qubo_processing()


# reco_xplets_simplified_LUXE(qubo_processor.get_kept_triplets(),
#                             new_folder)

# track_reconstruction_efficiency_simplified_LUXE(f"{new_folder}/reco_xplet_list.npy")

print("-----------------------------------\n")
print("QUBO solved successfully!\n")
