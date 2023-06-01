from math_functions.geometry import xyz_angle
from pattern.detector_hit import DetectorHit


class Doublet:
    def __init__(self,
                 hit_1: DetectorHit,
                 hit_2: DetectorHit):
        """Class for doublet objects, consisting of two hits on consecutive detector layers.
        :param hit_1: detector hit 1
        :param hit_2: detector hit 2

        """
        self.hit_1 = hit_1
        self.hit_2 = hit_2

    def is_correct_match(self) -> bool:
        """Checks if doublet hits stem from the same particle.
        :return
            True if created from same particle, else False.
        """
        if self.hit_1.particle_id == self.hit_2.particle_id:
            return True
        return False

    def xz_angle(self) -> float:
        """Calculating the angle in the xz-plane with respect to beam axis in z direction.
        :return
            angle in the xz plane.
        """
        return xyz_angle(self.hit_1.x,
                         self.hit_2.x,
                         self.hit_1.z,
                         self.hit_2.z)

    def yz_angle(self) -> float:
        """Calculating the angle in the xz-plane with respect to beam axis in z direction.
        :return
            angle in the yz plane.
        """
        return xyz_angle(self.hit_1.y,
                         self.hit_2.y,
                         self.hit_1.z,
                         self.hit_2.z)
