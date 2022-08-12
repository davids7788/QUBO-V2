from simplified_simulation.ptarmigan import PtargmiganSimData
from simplified_simulation.experimental_results_MC_toy import ExperimentalResults
from simplified_simulation.toy_experiment import MCToyExperiment
from simplified_simulation.detector_plane import DetectorPlane
from simplified_simulation.dipole_magnet import DipoleMagnet

import sys
import yaml

config_file = sys.argv[1]
ptarmigan = sys.argv[2]
save_location = sys.argv[3]

outfile_name = ".".join(ptarmigan.split("/")[-1].split(".")[0:-1])

with open(config_file) as f:
    input_file = yaml.safe_load(f)

detector_plane_list = []
for i, z12 in enumerate(zip(input_file["detector"]["layers right z positions"],
                            input_file["detector"]["layers left z positions"])):
    for j, chip1 in enumerate(input_file["detector"]["layers right x start"]):
        detector_plane_list.append(DetectorPlane((input_file["detector"]["layers right x start"][j],
                                                  input_file["detector"]["layers right x end"][j]),
                                                 (input_file["detector"]["layers positions y start"],
                                                  input_file["detector"]["layers positions y end"]),
                                                 z12[0],
                                                 (input_file["detector"]["num pixel x"],
                                                  input_file["detector"]["num pixel y"]),
                                                 input_file["detector"]["thickness"],
                                                 input_file["detector"]["radiation length"],
                                                 input_file["detector"]["rho"] *
                                                 input_file["settings"]["scattering"]))
    for k, chip2 in enumerate(input_file["detector"]["layers left x start"]):
        detector_plane_list.append(DetectorPlane((input_file["detector"]["layers left x start"][k],
                                                  input_file["detector"]["layers left x end"][k]),
                                                 (input_file["detector"]["layers positions y start"],
                                                  input_file["detector"]["layers positions y end"]),
                                                 z12[1],
                                                 (input_file["detector"]["num pixel x"],
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

Experiment_1 = MCToyExperiment('FL',
                               source,
                               'positron',
                               detector_plane_list,
                               dipole_magnet,
                               result,
                               scattering=input_file["settings"]["scattering"])

Experiment_1.start_experiment()

result.save_results(save_location + "/" + outfile_name + '_fl')
