import yaml
import argparse
import time

from simplified_simulation.ptarmigan import PtargmiganSimData
from simplified_simulation.experimental_results_MC_toy import ExperimentalResults
from simplified_simulation.toy_experiment import MCToyExperiment
from simplified_simulation.detector_plane import DetectorPlane
from simplified_simulation.dipole_magnet import DipoleMagnet
from simplified_simulation.convert_to_csv import *
from utility.time_tracking import hms_string

parser = argparse.ArgumentParser(description='Simplified Simulation',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--config_file',
                    action='store',
                    type=str,
                    default=None,
                    help='Simplified simulation configuration file')

parser.add_argument('--ptarmigan_file',
                    action='store',
                    type=str,
                    default=None,
                    help='Ptarmigan event file')

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
ptarmigan = parser_args.ptarmigan_file
geometry = parser_args.geometry_file
target_folder = parser_args.target_folder

outfile_name = ".".join(ptarmigan.split("/")[-1].split(".")[0:-1])
print("\nStarting simulation...")
print(f"Loading file {ptarmigan}")

with open(config_file) as f:
    input_file = yaml.safe_load(f)

layer_inf = []
with open(geometry, 'r') as file:
    csv_reader_geometry = csv.reader(file)
    csv_header_geometry = next(csv_reader_geometry)
    x_start = csv_header_geometry.index("x_start")
    x_end = csv_header_geometry.index("x_end")
    y_start = csv_header_geometry.index("y_start")
    y_end = csv_header_geometry.index("y_end")
    z = csv_header_geometry.index("z")
    for row in csv_reader_geometry:
        layer_inf.append((float(row[x_start]),
                          float(row[x_end]),
                          float(row[y_start]),
                          float(row[y_end]),
                          float(row[z])))

# In the simplified case, there is one big chip instead of 18 per layer. The pixel scale accounts for that.
if "sl" in geometry:
    pixel_scale = 18
    outfile_appendix = "_sl"
elif "fl" in geometry:
    pixel_scale = 1
    outfile_appendix = "_fl"
else:
    print("No valid geometry chosen! Exiting...")
    exit()

detector_plane_list = []
for i in range(len(layer_inf)):
    detector_plane_list.append(DetectorPlane((layer_inf[i][0], layer_inf[i][1]),
                                             (layer_inf[i][2], layer_inf[i][3]),
                                             layer_inf[i][4],
                                             (input_file["detector"]["num pixel x"] * pixel_scale,
                                              input_file["detector"]["num pixel y"]),
                                             input_file["detector"]["thickness"],
                                             input_file["detector"]["radiation length"],
                                             input_file["detector"]["rho"] *
                                             input_file["settings"]["scattering"]))

dipole_magnet = DipoleMagnet(input_file["dipole magnet"]["dipole start"],
                             input_file["dipole magnet"]["dipole end"],
                             input_file["dipole magnet"]["dipole field"])


source = PtargmiganSimData(ptarmigan)
result = ExperimentalResults()


Experiment_1 = MCToyExperiment(source,
                               'positron',
                               detector_plane_list,
                               dipole_magnet,
                               result,
                               scattering=input_file["settings"]["scattering"])

time_simulation_start = time.time()
Experiment_1.start_experiment()
time_simulation_end = time.time()
print(f"Time to run simulation: {hms_string(time_simulation_end - time_simulation_start)}")

result.save_results(target_folder + "/" + outfile_name + outfile_appendix)
print("Simulation finished!\n")

time_converting_start = time.time()
print("Converting to .csv files:")
print("Processing true hits...")
convert_to_csv_true(target_folder + "/" + outfile_name + outfile_appendix)

print("Smearing hits...")
convert_to_csv_smeared(target_folder + "/" + outfile_name + outfile_appendix)

# print("Calculating mid of fired pixels...\n")
# convert_to_csv_pixel(target_folder + "/" + outfile_name + outfile_appendix)

print("Converting finished successfully!")
time_converting_end = time.time()
print(f"Time to convert results to .csv files: {hms_string(time_converting_end - time_converting_start)}")
