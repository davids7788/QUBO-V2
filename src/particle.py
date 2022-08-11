class Particle:
    """Class for creating Particles used for a Toy MC experiment
    :param:
        position    : position of the particle (4-vector)
        momentum    : momentum of the particle (4-vector)
        particle_ID : ID of a particle
    """
    def __init__(self,
                 position,
                 momentum,
                 particle_id):

        self.particle_ID = particle_id
        self.position = position
        self.momentum = momentum

    def move_along_beam(self, distance):
        """Propagates x, y and z value for moving along the beam.
        :param
            distance: distance in z direction
        """
        [px, py, pz] = self.momentum
        self.position = [self.position[0] + distance / pz * px,
                         self.position[1] + distance / pz * py,
                         self.position[2] + distance]
