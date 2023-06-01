from math_functions.geometry import xyz_angle
from pattern.detector_hit import DetectorHit


class Triplet:
    def __init__(self,
                 hit_1: DetectorHit,
                 hit_2: DetectorHit,
                 hit_3: DetectorHit):
        """Class for triplet objects consisting of three hits on consecutive layers.
        :param hit_1: detector hit 1
        :param hit_2: detector hit 2
        :param hit_3: detector hit 3
        """
        self.hit_1 = hit_1
        self.hit_2 = hit_2
        self.hit_3 = hit_3

        # hit IDs are used to generate the triplet id -> unique identifier
        self.triplet_id = "_".join([self.hit_1.hit_id,
                                    self.hit_2.hit_id,
                                    self.hit_3.hit_id])

        self.interactions = {}   # Interactions with other triplets, {<other triplet id> : value}
        self.quality = 0.0   # Describing how well it fits the expected particle trajectory

    def angles_between_doublets(self) -> tuple[float, float]:
        """Calculating the angles between the doublets in xz and yz
        :return
             (angle xz of doublets, angle yz of doublets)
        """
        xz_12 = xyz_angle(self.hit_1.x,
                          self.hit_2.x,
                          self.hit_1.z,
                          self.hit_2.z)

        xz_23 = xyz_angle(self.hit_2.x,
                          self.hit_3.x,
                          self.hit_2.z,
                          self.hit_3.z)

        yz_12 = xyz_angle(self.hit_1.y,
                          self.hit_2.y,
                          self.hit_1.z,
                          self.hit_2.z)

        yz_23 = xyz_angle(self.hit_2.y,
                          self.hit_3.y,
                          self.hit_2.z,
                          self.hit_3.z)
        return xz_23 - xz_12, yz_23 - yz_12

    def is_correct_match(self) -> bool:
        """Checks if all hits of the triplet stem from the same particle.
        :return
            True if triplet originates from one single particle, else False.
        """
        if self.hit_1.particle_id == self.hit_2.particle_id == self.hit_3.particle_id:
            return True
        return False
