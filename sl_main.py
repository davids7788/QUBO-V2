from src.ptarmigan import PtargmiganSimData
from src.experimental_results_MC_toy import ExperimentalResults
from src.toy_experiment import MCToyExperiment
from src.detector_plane import DetectorPlane
from src.dipole_magnet import DipoleMagnet

import sys
import yaml

config_file = sys.argv[1]
ptarmigan = sys.argv[2]
save_location = sys.argv[3]

outfile_name = ".".join(ptarmigan.split("/")[-1].split(".")[0:-1])


with open(config_file) as f:
    input_file = yaml.safe_load(f)

detector_plane_list = []
for i in range(4):
    detector_plane_list.append(DetectorPlane((input_file["detector"]["layer positions x start"],
                                              input_file["detector"]["layer positions x end"]),
                                             (input_file["detector"]["layer positions y start"],
                                              input_file["detector"]["layer positions y end"]),
                                             input_file["detector"]["layer positions z"][i],
                                             (input_file["detector"]["num pixel x"] * 18,          # simple setup
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


Experiment_1 = MCToyExperiment('SL',
                               source,
                               'positron',
                               detector_plane_list,
                               dipole_magnet,
                               result,
                               scattering=input_file["settings"]["scattering"])

Experiment_1.start_experiment()

result.save_results(save_location + "/" + outfile_name + '_sl')
