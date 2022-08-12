import numpy as np


class ExperimentalResults:
    def __init__(self):
        self.true_detector_hits_dictionary = {}
        self.particle_momentum_history = {}
        self.particle_position_history = {}
        self.list_of_planes = {}

    def save_results(self, file_name):
        """Saves the result as a numpy file."""
        np.save(file_name, np.array({'True detector hits': self.true_detector_hits_dictionary,
                                     'Particle momentum history log': self.particle_momentum_history,
                                     'Particle position history log': self.particle_position_history,
                                     'Plane list': self.list_of_planes}))