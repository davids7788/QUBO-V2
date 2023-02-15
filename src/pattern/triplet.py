from pattern.doublet import Doublet


class Triplet:
    def __init__(self,
                 doublet_1: Doublet,
                 doublet_2: Doublet):
        """Class for triplet objects consisting of two connected Doublet objects. The two Doublet objects need
        to have a common hit. This has to be the second hit of one Doublet and the first hit of the other Doublet.
        :param doublet_1: doublet object 1
        :param doublet_2: doublet object 2
        """
        self.doublet_1 = doublet_1
        self.doublet_2 = doublet_2
        if self.doublet_1.hit_2_id != self.doublet_2.hit_1_id:
            print("Doublets are not forming a triplet. Please check the triplet creation procedure!")

        # hit IDs are used to generate the triplet id -> unique identifier
        self.triplet_id = "_".join([self.doublet_1.hit_1_id,
                                    self.doublet_1.hit_2_id,
                                    self.doublet_2.hit_2_id])
        self.interactions = {}   # Interactions with other triplets, {<other triplet id> : value}
        self.quality = 0.0   # Describing how well it fits the expected particle trajectory

    def angles_between_doublets(self) -> tuple[float, float]:
        """Calculating the angles between the doublets in xz and yz
        :return
             (angle xz of doublets, angle yz of doublets)
        """
        return self.doublet_2.xz_angle() - self.doublet_1.xz_angle(), \
            self.doublet_2.yz_angle() - self.doublet_1.yz_angle()

    def is_correct_match(self) -> bool:
        """Checks if all hits of the triplets stem from the same particle.
        :return
            True if triplet originates from one single particle, else False.
        """
        if self.doublet_1.hit_1_particle_key == self.doublet_1.hit_2_particle_key == \
                self.doublet_2.hit_1_particle_key == self.doublet_2.hit_2_particle_key:
            return True
        return False
