import yaml
import argparse

from pathlib import Path
from track_reconstruction.create_gen_xplets import *

from pattern_building.pattern_builder import PatternBuilder
from pattern_building.segment_manager import SegmentManager
from pattern_building.qubo_coefficients import QuboCoefficients
from pattern_building.plot_statistics import plot_coefficients_statistics

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

parser.add_argument('--geometry_file',
                    action='store',
                    type=str,
                    default=None,
                    help='LUXE geometry csv file')

parser.add_argument('--target_folder',
                    action='store',
                    type=str,
                    default=None,
                    help='Folder to store results')

parser_args = parser.parse_args()
config_file = parser_args.config_file
tracking_data = parser_args.tracking_data
geometry_file = parser_args.geometry_file
target_folder = parser_args.target_folder


# loading arguments, creating folder
with open(config_file, 'r') as f:
    configuration = yaml.safe_load(f)


# Creating a new folder which is named as the config file and stored inside the folder according to a set folder
new_folder = target_folder + "/" + tracking_data.split("/")[-1].split(".csv")[0] + "-" + \
             ".".join(config_file.split("/")[-1].split(".")[0:-1])

if Path(new_folder).is_dir():
    pass
else:
    os.mkdir(new_folder)


# Program
print("\n-----------------------------------")
print("\nStarting QUBO creation...\n")

# Segmentation algorithm --> reduce combinatorial tasks
s_manager = SegmentManager(configuration,
                           geometry_file)
s_manager.create_LUXE_segments()
s_manager.segment_mapping_LUXE()

# Triplet creation
pattern_builder = PatternBuilder(configuration, new_folder)
pattern_builder.load_tracking_data(tracking_data, s_manager)
pattern_builder.create_x_plets_LUXE(s_manager)
pattern_builder.write_info_file()

# Create truth Xplets
gen_xplets_simplified_LUXE(tracking_data, "/".join(new_folder.split("/")[0:-1]))


# Set and rescale parameters plot statistics
qubo_coefficients = QuboCoefficients(configuration, new_folder)
qubo_coefficients.set_triplet_coefficients(s_manager)
qubo_coefficients.coefficient_rescaling()

# Just for visualising coefficient distribution
plot_coefficients_statistics(pattern_builder.num_particles,
                             qubo_coefficients)

print("-----------------------------------\n")
print("QUBO preparation finished successfully!\n")
