from math_functions.geometry import xyz_angle
from pattern.detector_hit import DetectorHit


class Triplet:
    """Class for combining DetectorHits to a triplet structure. It provides methods to analyse the triplet.
    """
    def __init__(self,
                 hit_1: DetectorHit,
                 hit_2: DetectorHit,
                 hit_3: DetectorHit):
        """Set fields according to the provided DetectorHit objects.

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

        # interactions with other triplets, {<other triplet id> : value}
        self.interactions = {}

        # quality of the triplet, e.g. how well it fits the expected particle trajectory
        self.quality = 0.0

    def from_same_particle(self) -> bool:
        """Check if all hits of the triplet stem from the same particle.

        :return
            True if triplet originates from one single particle, else False.
        """
        # check if particles stem from signal
        if not self.hit_1.is_signal:
            return False
        if not self.hit_2.is_signal:
            return False
        if not self.hit_3.is_signal:
            return False

        # check if same particle
        if self.hit_1.particle_id == self.hit_2.particle_id == self.hit_3.particle_id:
            return True
        return False
