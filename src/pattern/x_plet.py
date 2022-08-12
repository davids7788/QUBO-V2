from pattern.triplet import Triplet


class Xplet:
    def __init__(self):
        self.hit_ids = {}
        self.particle_ids = {}
        self.coordinates = {}
        self.triplet_ids = []
        self.is_empty = True
    """This class takes a series of hits in the form of connectable triplets and combines them into a x-plet pattern.
    """
    def add_triplet(self,
                    triplet: Triplet):
        self.triplet_ids.append(triplet.triplet_id)
        if self.is_empty:
            self.hit_ids.update({0: triplet.doublet_1.hit_1_id})
            self.hit_ids.update({1: triplet.doublet_1.hit_2_id})
            self.hit_ids.update({2: triplet.doublet_2.hit_2_id})

            self.particle_ids.update({0: triplet.doublet_1.hit_1_particle_key})
            self.particle_ids.update({1: triplet.doublet_1.hit_2_particle_key})
            self.particle_ids.update({2: triplet.doublet_2.hit_2_particle_key})

            self.coordinates.update({0: triplet.doublet_1.hit_1_position})
            self.coordinates.update({1: triplet.doublet_1.hit_2_position})
            self.coordinates.update({2: triplet.doublet_2.hit_2_position})
        else:
            insert_position = len(self.triplet_ids) + 2
            self.hit_ids.update({insert_position: triplet.doublet_2.hit_2_id})
            self.particle_ids.update({insert_position: triplet.doublet_2.hit_2_particle_key})
            self.coordinates.update({insert_position: triplet.doublet_2.hit_2_position})

        self.is_empty = False
