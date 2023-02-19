import sys
import yaml

from pathlib import Path
from track_reconstruction.create_gen_xplets import *

from preselection.triplet_creator_LUXE import TripletCreatorLUXE
from preselection.segment_manager import SegmentManager
from preselection.qubo_coefficients import QuboCoefficients

# sys argv [1]: config file
# sys argv [2]: tracking data, csv file
# sys argv [3]: geometry file
# sys argv [4]: folder to store result

# loading arguments, creating folder
with open(sys.argv[1], 'r') as f:
    config_file = yaml.safe_load(f)
tracking_data = sys.argv[2]

geometry_file = sys.argv[3]

# Creating a new folder which is named as the config file and stored inside the folder according to a set folder
new_folder = sys.argv[4] + "/" + tracking_data.split("/")[-1].split(".csv")[0] + "-" + \
             ".".join(sys.argv[1].split("/")[-1].split(".")[0:-1])

if Path(new_folder).is_dir():
    pass
else:
    os.mkdir(new_folder)


# Program

print("\n-----------------------------------")
print("\nStarting QUBO creation...\n")

# Segmentation algorithm --> reduce combinatorial tasks
s_manager = SegmentManager([config_file['binning']['num bins x'],
                            config_file['binning']['num bins y']],
                           config_file['doublet'],
                           geometry_file)
s_manager.create_LUXE_segments()
s_manager.segment_mapping_LUXE()

# Triplet creation
triplet_creator = TripletCreatorLUXE(config_file, new_folder)
triplet_creator.load_tracking_data(tracking_data, s_manager)
triplet_creator.create_x_plets_simplified_LUXE(s_manager)
triplet_creator.write_info_file()


# Create truth Xplets
gen_xplets_simplified_LUXE(tracking_data, "/".join(new_folder.split("/")[0:-1]))
exit()

# Set and rescale parameters plot statistics
qubo_coefficients = QuboCoefficients(config_file, new_folder)
qubo_coefficients.set_triplet_coefficients(s_manager)
qubo_coefficients.filling_lists_for_statistics()
qubo_coefficients.parameter_rescaling()

# Just for visualising coefficients
qubo_coefficients.plot_and_save_statistics(triplet_creator.num_complete_tracks,
                                           triplet_creator.preselection_statistic_dx_x0,
                                           triplet_creator.preselection_statistic_scattering)

print("-----------------------------------\n")
print("QUBO preparation finished successfully!\n")
