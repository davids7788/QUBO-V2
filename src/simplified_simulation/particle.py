import numpy as np


class Particle:
    """Class for creating Particles used for a Toy MC experiment
    :param:
        position    : position of the particle (4-vector)
        momentum    : momentum of the particle (4-vector) in
        particle_ID : ID of a particle
    """
    def __init__(self,
                 position,
                 momentum,
                 particle_id):

        self.particle_ID = particle_id
        self.position = position
        self.momentum = momentum
        self.time = 0

    def move_along_beam(self, distance):
        """Propagates x, y and z value for moving along the beam.
        :param distance: distance in z direction
        """
        [px, py, pz] = self.momentum
        self.position = [self.position[0] + distance / pz * px,
                         self.position[1] + distance / pz * py,
                         self.position[2] + distance]

    def update_time(self, distance):
        c = 299792458   # m/s
        e_total = np.sqrt(momentum[0]**2 + momentum[1]**2 + momentum[2]**2)   # MeV
        e_kin = e_total - 0.511   #
        velocity = c * np.sqrt(1 - (1 / (1 + e_kin/0.511))**2)
        self.time += distance / velocity
