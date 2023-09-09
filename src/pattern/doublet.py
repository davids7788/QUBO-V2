from math_functions.geometry import xyz_angle
from pattern.detector_hit import DetectorHit


class Doublet:
    """Class for combining DetectorHits to a doublet structure. It provides methods to analyse the doublet.
    """
    def __init__(self,
                 hit_1: DetectorHit,
                 hit_2: DetectorHit):
        """Set fields according to the provided DetectorHit objects.
        :param hit_1: detector hit 1
        :param hit_2: detector hit 2
        """
        self.hit_1 = hit_1
        self.hit_2 = hit_2

    def from_same_particle(self) -> bool:
        """Checks if doublet hits stem from the same particle.
        :return
            True if created from same signal particle, else False.
        """
        # check if particles stem from signal
        if self.hit_1.particle_info != 'signal':
            return False
        if self.hit_2.particle_info != 'signal':
            return False

        # check if same particle
        if self.hit_1.particle_id == self.hit_2.particle_id:
            return True
        return False

    def xz_angle(self) -> float:
        """Returns angle in the xz-plane with respect to the beam axis in z direction.
        :return
            angle in the xz-plane.
        """
        return xyz_angle(self.hit_1.x,
                         self.hit_2.x,
                         self.hit_1.z,
                         self.hit_2.z)

    def yz_angle(self) -> float:
        """Return the angle in the xz-plane with respect to the beam axis in z direction.
        :return
            angle in the yz-plane.
        """
        return xyz_angle(self.hit_1.y,
                         self.hit_2.y,
                         self.hit_1.z,
                         self.hit_2.z)
