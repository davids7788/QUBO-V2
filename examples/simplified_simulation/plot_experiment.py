import numpy as np
import os
import sys
sys.path.append('../../src'),

from simplified_simulation.detector_plane import DetectorPlane
from simplified_simulation.visualisation import Visualisation
from simplified_simulation.ptarmigan import PtargmiganSimData

# .h5 ptarmigan source file, make sure it exists
source = PtargmiganSimData("../../ptarmigan_files/xi_4/e0gpc_4.0_0000_particles.h5")
# folder with Simplified Simulation result,

result = np.load('../../simplified_simulation_files/xi_4/e0gpc_4.0_0000_particles_sl.npy', allow_pickle=True)
os.chdir('../../simplified_simulation_files/xi_4')

visualisation = Visualisation(source=source,
                              result=result)

visualisation.beam_source_visualisation(),
visualisation.plot_experiment_with_results([-0.015, 0.015],       # frame limits in 3D,
                                           [0.025, 0.575],
                                           [3.8, 4.3],
                                           number_of_displayed_tracks=50)
