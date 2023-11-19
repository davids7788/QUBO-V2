import yaml
import argparse
import os

from pathlib import Path

from pattern_building.pattern_builder import PatternBuilder
from pattern_building.segment_manager import LUXESegmentManager
from pattern_building.qubo_coefficients import QuboCoefficients

from qubo.qubo_processing import QuboProcessing
from qubo.qubo_logging import QuboLogging
from qubo.ansatz import Ansatz
from qubo.solver import Solver

from track_building.gen_multiplets import GenMultiplet
from track_building.reco_multiplets import make_reco_multiplets
from track_building.efficiency import get_efficiency

parser = argparse.ArgumentParser(description='QUBO pattern_building Simplified LUXE',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--steering_file',
                    action='store',
                    type=str,
                    default=None,
                    help='Preselection configuration file')

parser.add_argument('--tracking_data',
                    action='store',
                    type=str,
                    default=None,
                    help='Tracking data csv file')

# get input
parser_args = parser.parse_args()
steering_file = parser_args.steering_file
tracking_data = parser_args.tracking_data

# loading steering_file, creating folder
with open(steering_file, 'r') as f:
    configuration = yaml.safe_load(f)

# information about the detector setting, the given data format of the tracking file and the sample composition
geometry_file = configuration['tracking data']['detector geometry']
tracking_data_format = configuration['tracking data']['tracking data format']
sample_composition = configuration['tracking data']['sample composition']

# make a new folder which is named <tracking_data>_qubo in the directory where the tracking file is located
qubo_preparation_folder = f'{tracking_data}_qubo'
if Path(qubo_preparation_folder).is_dir():
    pass
else:
    os.mkdir(qubo_preparation_folder)

# start of program
print('\n-----------------------------------')
print('Starting pattern recognition...\n')

# setting up segmentation manager which is used to faster map hits from consecutive layers
s_manager = LUXESegmentManager(configuration, geometry_file)
s_manager.create_LUXE_segments()
s_manager.segment_mapping_LUXE()

# setting up pattern builder object and loading data from specified file
pattern_builder = PatternBuilder(configuration)
pattern_builder.load_tracking_data(tracking_data, s_manager, tracking_data_format)
pattern_builder.create_multiplets(s_manager)

if sample_composition == 'blinded':
    print('Truth information about particle tracks is not available in blinded sample!\n')
    print('-----------------------------------\n')
    print('Building multiplets is not possible from blinded sample!')
else:
    pattern_builder.information_about_particle_tracks(z_position_layers=s_manager.z_position_to_layer,
                                                      setup=s_manager.setup,
                                                      sample_composition=sample_composition,
                                                      min_track_length=configuration['track']['minimum track length'])
    print('Building multiplets on generator level with truth information...\n')
    gen_multiplets = GenMultiplet(tracking_data_file=tracking_data,
                                  min_track_length=configuration['track']['minimum track length'])

    gen_multiplets.make_gen_multiplets(tracking_data_format)
    gen_multiplets.save_multiplets(save_to_folder=qubo_preparation_folder)


# set qubo parameters
print('Calculate triplet coefficients a_i and b_ij...')
qubo_coefficients = QuboCoefficients(configuration, qubo_preparation_folder)
qubo_coefficients.set_triplet_coefficients(s_manager)

# rescale parameters
qubo_coefficients.coefficient_rescaling()

print('-----------------------------------\n')
print('QUBO preparation finished successfully!\n')


print('-----------------------------------\n')

file_extension = f'{configuration["qubo"]["optimisation strategy"].replace(" ", "_")}_' \
                 f'{configuration["solver"]["algorithm"].replace(" ","_")}_' \
                 f'{configuration["qubo"]["num qubits"]}'
qubo_folder = qubo_preparation_folder + "/" + file_extension

if Path(qubo_folder).is_dir():
    pass
else:
    os.mkdir(qubo_folder)
                 
# Create logger
qubo_logger = QuboLogging()

# Create ansatz object and set parameters from config file
if configuration['ansatz']['layout'] is None:
    ansatz = None
else:
    ansatz = Ansatz(config=steering_file)

# If not "hamiltonian driven" additional parameters are set
if configuration['ansatz']['layout'] == 'TwoLocal':
    ansatz.set_two_local()
elif configuration['ansatz']['layout'] is None:
    if configuration['solver']['algorithm'] == 'QAOA':
        pass
    elif configuration['solver']['algorithm'] is not None:
        if configuration['solver']['algorithm'] != 'Numpy Eigensolver':
            ansatz.set_no_entanglements()

# Create Solver object and set parameters from config file
if configuration['solver']['algorithm'] != 'Numpy Eigensolver' and configuration['solver']['algorithm'] is not None:
    solver = Solver(configuration)
else:
    solver = None

# Create and configure solving process
qubo_processor = QuboProcessing(qubo_preparation_folder + '/triplet_list.npy',
                                config=configuration,
                                solver=solver,
                                ansatz=ansatz,
                                qubo_logging=qubo_logger,
                                save_folder=qubo_folder,
                                verbose=1)
qubo_processor.qubo_processing()        
make_reco_multiplets(qubo_processor.get_kept_triplets(),
                     qubo_folder,
                     tracking_data)


gen_prefix = f'{qubo_folder}/reco_xplet_list.npy'.split('/')[-3].split('-')[0]
gen_xplet = f'{qubo_preparation_folder}/{tracking_data.split("/")[-1].split(".")[0]}_gen_xplet_list'

try:
    get_efficiency(f'{qubo_folder}/reco_xplet_list_ambiguity_solved.npy',
                   f'{gen_xplet}.npy')
except FileNotFoundError:
    print('\nNo track reconstruction efficiency available in blinded samples!')

print('-----------------------------------\n')
print('Pattern recognition finished successfully!\n')
