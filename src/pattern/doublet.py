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
        if not self.hit_1.is_signal:
            return False
        if not self.hit_2.is_signal:
            return False

        # check if same particle
        if self.hit_1.particle_id == self.hit_2.particle_id:
            return True
        return False
