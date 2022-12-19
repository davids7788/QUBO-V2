import csv
import sys
import yaml

from simplified_simulation.ptarmigan import PtargmiganSimData
from simplified_simulation.experimental_results_MC_toy import ExperimentalResults
from simplified_simulation.toy_experiment import MCToyExperiment
from simplified_simulation.detector_plane import DetectorPlane
from simplified_simulation.dipole_magnet import DipoleMagnet
from simplified_simulation.convert_to_csv import *

config_file = sys.argv[1]
ptarmigan = sys.argv[2]
geometry = sys.argv[3]
save_location = sys.argv[4]

outfile_name = ".".join(ptarmigan.split("/")[-1].split(".")[0:-1])

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

Experiment_1.start_experiment()

result.save_results(save_location + "/" + outfile_name + outfile_appendix)
convert_to_csv_true(save_location + "/" + outfile_name + outfile_appendix)
convert_to_csv_smeared(save_location + "/" + outfile_name + outfile_appendix)
# convert_to_csv_pixel(save_location + "/" + outfile_name + outfile_appendix)
