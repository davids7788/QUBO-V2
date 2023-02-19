import numpy as np


class Doublet:
    def __init__(self,
                 hit_1_particle_key: str = '',
                 hit_2_particle_key: str = '',
                 hit_1_position: tuple[float] = (0.0, 0.0, 0.0),
                 hit_2_position: tuple[float] = (0.0, 0.0, 0.0),
                 hit_1_id: str = '',
                 hit_2_id: str = '',
                 energy_1: float = -1.0,
                 energy_2: float = -1.0,
                 time_1: float = -1e10,
                 time_2: float = -1e10):
        """Class for doublet objects, consisting of two hits on consecutive detector layers.
        :param hit_1_particle_key: particle true number (simulation only), default = ''
        :param hit_2_particle_key: particle true number (simulation only), default = ''
        :param hit_1_position: particle position (x, y, z) of hit 1, default: (0.0, 0.0, 0.0)
        :param hit_2_position: particle position (x, y, z) of hit 2, default: (0.0, 0.0, 0.0)
        :param hit_1_id: unique integer hit id on detector of hit 1, default: ''
        :param hit_2_id: unique integer hit id on detector of hit 2, default: ''
        :param energy_1: energy value of particle from hit 1, default: -1.0
        :param energy_2: energy value of particle from hit 2, default: -1.0
        :param time_1: absolute time of the hit 1 [s], default: -1e10
        :param time_2: absolute time of the hit 2 [s], default: -1e10
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

    def is_correct_match(self) -> bool:
        """Checks if doublet hits stem from the same particle.
        :return
            True if created from same particle, else False. If both particles have
        """
        if self.hit_1_particle_key == self.hit_2_particle_key:
            return True
        return False

    def xz_angle(self) -> tuple[float]:
        """Calculating the angle in the xz-plane with respect to beam axis in z direction.
        :return
            angle in the xz plane.
        """
        return np.arctan2((self.hit_2_position[0] - self.hit_1_position[0]),
                          (self.hit_2_position[2] - self.hit_1_position[2]))

    def yz_angle(self) -> tuple[float]:
        """Calculating the angle in the xz-plane with respect to beam axis in z direction.
        :return
            angle in the yz plane.
        """
        return np.arctan2((self.hit_2_position[1] - self.hit_1_position[1]),
                          (self.hit_2_position[2] - self.hit_1_position[2]))
