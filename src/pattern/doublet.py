import numpy as np


class Doublet:
    def __init__(self,
                 hit_1_particle_key: int,
                 hit_2_particle_key: int,
                 hit_1_position: list[float],
                 hit_2_position: list[float],
                 hit_1_id: int,
                 hit_2_id: int,
                 energy_1: float,
                 energy_2: float,
                 time_1: float = 0,
                 time_2: float = 0):
        """Class for doublet objects, consisting of two hits on different detector layers.
        :param
            hit_1_particle_key: particle true number, from simulation
            hit_2_particle_key: particle true number, from simulation
            hit_1_position: particle position [x, y, z] on detector layer of hit 1
            hit_2_position: particle position [x, y, z] or detector layer of hit 2
            hit_1_id: unique integer hit id on detector of hit 1
            hit_2_id: unique integer hit id on detector of hit 2
            energy_1: energy value of particle from hit 1
            energy_2: energy value of particle from hit 2
            time_1: absolute time of the hit in s
            time_2: absolute time of the hit in s
        """
        self.hit_1_position = hit_1_position
        self.hit_1_particle_key = hit_1_particle_key
        self.hit_2_position = hit_2_position
        self.hit_2_particle_key = hit_2_particle_key
        self.hit_1_id = hit_1_id
        self.hit_2_id = hit_2_id
        self.energy_1 = energy_1
        self.energy_2 = energy_2
        self.time_1 = time_1
        self.time_2 = time_2

    def is_correct_match(self):
        """Checks if doublet hits stem from the same particle.
        :return
            True if created from same particle, else False.
        """
        if self.hit_1_particle_key == self.hit_2_particle_key:
            return True
        return False

    def xz_angle(self):
        """Returns the angle in the xz-plane with respect to beam axis in z direction.
        :return
            angle in the xz plane.
        """
        return np.arctan2((self.hit_2_position[0] - self.hit_1_position[0]),
                          (self.hit_2_position[2] - self.hit_1_position[2]))

    def yz_angle(self):
        """Returns the angle in the xz-plane with respect to beam axis in z direction.
        :return
            angle in the yz plane.
        """
        return np.arctan2((self.hit_2_position[1] - self.hit_1_position[1]),
                          (self.hit_2_position[2] - self.hit_1_position[2]))
