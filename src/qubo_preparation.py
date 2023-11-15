import yaml
import argparse
import os

from pathlib import Path
from track_reconstruction.gen_multiplets import GenMultiplet

from pattern_building.pattern_builder import PatternBuilder
from pattern_building.segment_manager import SegmentManager
from pattern_building.qubo_coefficients import QuboCoefficients

parser = argparse.ArgumentParser(description='QUBO pattern_building Simplified LUXE',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--config_file',
                    action='store',
                    type=str,
                    default=None,
                    help='Preselection configuration file')

parser.add_argument('--tracking_data',
                    action='store',
                    type=str,
                    default=None,
                    help='Tracking data csv file')

# set variables from inputs
parser_args = parser.parse_args()
config_file = parser_args.config_file
tracking_data = parser_args.tracking_data

# loading arguments, creating folder
with open(config_file, 'r') as f:
    configuration = yaml.safe_load(f)

geometry_file = configuration['detector geometry']
tracking_data_format = configuration['tracking data format']
sample_composition = configuration['sample composition']


# creating a new folder which is named as the config_file, folder is stored inside the specified target_folder
new_folder = f'{tracking_data}_qubo'
print(new_folder)

if Path(new_folder).is_dir():
    pass
else:
    os.mkdir(new_folder)

# start of program
print("\n-----------------------------------")
print("\nStarting QUBO creation...\n")

# setting up segmentation manager and algorithm to reduce computation time of pattern building
s_manager = SegmentManager(configuration,
                           geometry_file)
s_manager.create_LUXE_segments()
s_manager.segment_mapping_LUXE()

# setting up pattern builder object and loading data from specified file
pattern_builder = PatternBuilder(configuration)
pattern_builder.load_tracking_data(tracking_data, 
                                   s_manager, 
                                   tracking_data_format)

if sample_composition == 'blinded':
    print("Truth information about particle tracks is not available in blinded sample!\n")
else:
    pattern_builder.information_about_particle_tracks(z_position_layers=s_manager.z_position_to_layer,
                                                      setup=s_manager.setup,
                                                      sample_composition=sample_composition,
                                                      min_track_length=configuration['track']['minimum track length'])
pattern_builder.create_multiplets(s_manager)

# pattern creation
print("-----------------------------------\n")
if sample_composition == 'blinded':
    print("Creating truth multiplets is not possible from blinded sample!\n")
else:
    print("Create multiplets on generator level with truth information...\n")
    gen_multiplets = GenMultiplet(tracking_data_file=tracking_data,
                                  min_track_length=configuration['track']['minimum track length'])

    gen_multiplets.make_gen_multiplets(tracking_data_format)
    gen_multiplets.information_about_tracking_data()
    gen_multiplets.save_multiplets(save_to_folder=new_folder)


# set qubo parameters
print("Calculate triplet coefficients a_i and b_ij...")
qubo_coefficients = QuboCoefficients(configuration, new_folder)
qubo_coefficients.set_triplet_coefficients(s_manager)

# rescale parameters
qubo_coefficients.coefficient_rescaling()

print("-----------------------------------\n")
print("QUBO preparation finished successfully!\n")
