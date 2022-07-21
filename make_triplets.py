import sys
import yaml
import os

from pathlib import Path

from x_plets import XpletCreatorLUXE
from segment_manager import SegmentManager


# sys argv [1]: config file
# sys argv [2]: tracking data, csv file
# sys argv [3]: geometry file
# sys argv [4]: folder to store result

# loading arguments, creating folder
with open(sys.argv[1], 'r') as f:
    config_file = yaml.safe_load(f)
tracking_data = sys.argv[2]

# Creating a new folder which is named as the config file and stored inside the folder according to a set folder
new_folder = sys.argv[4] + "/" + tracking_data.split("/")[-1].split(".csv")[0] + "-" + \
             ".".join(sys.argv[1].split("/")[-1].split(".")[0:-1])

if Path(new_folder).is_dir():
    pass
else:
    os.mkdir(new_folder)

# Program
s_manager = SegmentManager(config_file, sys.argv[3])

x_plet_creator = XpletCreatorLUXE(config_file, new_folder)
x_plet_creator.load_tracking_data(tracking_data, s_manager)
x_plet_creator.make_x_plet_list_simplified_setup(s_manager)
x_plet_creator.write_info_file()

qubo_coefficients = QuboCoefficients(config_file, new_folder)
qubo_coefficients.set_triplet_coefficients(s_manager)
qubo_coefficients.filling_lists_for_statistics()
qubo_coefficients.parameter_rescaling()
qubo_coefficient.plot_and_save_statistics()
