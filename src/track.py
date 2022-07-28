from src.doublet import Doublet
from src.triplet import Triplet


class Track:
    def __init__(self):
        """Class for handling track objects created from triplets.
        """
        self.triplets = []
        self.is_matched_track = True

    def check_match(self):
        """Checks if track has only entries stemming from a single particle and sets the corresponding attribute
        accordingly.
        """
        for t in self.triplets:
            if not t.is_correct_match:
                self.is_matched_track = False

    def add_triplet_to_track(self,
                             triplet: Triplet):
        """Adds a hit to the track and checks if it's still a matched track --> stemming from one particle.
        :param triplet: another Triplet"""
        self.triplets.append(triplet)
        self.check_match()

    def is_in_conflict(self,
                       other_track):
        """Checks if there's some ambiguity and returns the result.
        :param other_track: other track object
        :return
            True if conflict, else False"""
        if not set([t.triplet_id for t in self.triplets]).isdisjoint([t.triplet_id for t in other_track.triplets]):
            return True
        return False
