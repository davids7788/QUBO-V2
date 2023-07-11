import numpy as np


class Particle:
    """Class for creating Particles used for a Toy MC experiment
    :param position    : position of the particle (4-vector)
    :param momentum    : momentum of the particle (4-vector) in
    :param particle_ID : ID of a particle
    """
    def __init__(self,
                 position,
                 momentum,
                 particle_id):

        self.particle_ID = particle_id
        self.position = position
        self.momentum = momentum
        self.time = 0

    def move_along_beam(self, dz):
        """Propagates x, y and z value for moving along the beam.
        :param dz: distance in z direction
        """
        [px, py, pz] = self.momentum
        self.position = [self.position[0] + dz / pz * px,
                         self.position[1] + dz / pz * py,
                         self.position[2] + dz]

    def update_time(self, distance):
        """Updates the time attribute for a particle which traveled the given distance
        :param distance: distance [m]
        """
        c = 299792458   # m/s
        e_total = np.sqrt(self.momentum[0]**2 + self.momentum[1]**2 + self.momentum[2]**2)   # MeV
        e_kin = e_total - 0.511   #
        velocity = c * np.sqrt(1 - (1 / (1 + e_kin/0.511))**2)
        self.time += distance / velocity
