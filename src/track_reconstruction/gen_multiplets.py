import csv
import numpy as np

from pattern.detector_hit import DetectorHit
from pattern.multiplet import Multiplet
from utility.data_format_handler import load_data


class GenMultiplet:
    """Class for multiplet creation on generator level information
    """
    def __init__(self,
                 tracking_data_file: str,
                 min_track_length: int) -> None:
        """Set fields
        :param tracking_data_file: file with tracking data information
        :param min_track_length: minimum required length of multiplet to be considered a track candidate
        """
        self.multiplets_dict = {}
        self.gen_multiplets = []
        self.tracking_data_file = tracking_data_file
        self.output_name = ".".join(tracking_data_file.split("/")[-1].split(".")[0:-1])

        self.num_signal_hits = 0
        self.num_background_hits = 0
        self.min_track_length = min_track_length

    def save_multiplets(self,
                        save_to_folder: str):
        """Saves the multiplets in the specified folder.
        :param save_to_folder: folder to store the multiplets
        """
        np.save(f"{save_to_folder}/{self.output_name}_gen_xplet_list", self.gen_multiplets)
        print(f"Number of generated multiplets: {len(self.gen_multiplets)}")
        print(f"Number of tracks with at least {self.min_track_length} hits: "
              f"{len([g for g in self.gen_multiplets if len(g.hit_id) >= self.min_track_length])}\n")
        print("-----------------------------------\n")

    def create_multiplets(self):
        """Creates the multiplets from the sorted particle dictionary.
        """
        for particle_entry in self.multiplets_dict.values():
            particle_entry.sort(key=lambda entry: entry.z)

            truth_multiplet = Multiplet()
            for value in particle_entry:
                truth_multiplet.add_hit(value)

            self.gen_multiplets.append(truth_multiplet)

    def add_hit_to_multiplet_dict(self,
                                  hit: DetectorHit) -> None:
        """Adds a hit to the multiplet dictionary.
        :param hit: DetectorHit object
        """
        # background hit
        if not hit.is_signal:
            self.num_background_hits += 1

        # signal hit
        else:
            self.num_signal_hits += 1
            if hit.particle_id in self.multiplets_dict.keys():
                self.multiplets_dict[hit.particle_id].append(hit)
            else:
                self.multiplets_dict.update({hit.particle_id: [hit]})

    def make_gen_multiplets(self,
                            tracking_data_format) -> None:
        """Creates all signal multiplets from the specified simplified simulation tracking data file.
        :param tracking_data_format: format of the tracking data
                                     --> 'key4hep slcio', 'key4hep csv', 'simplified simulation csv'
        """
        hits = load_data(self.tracking_data_file, tracking_data_format)
        for hit in hits:
            # add to multiplet dictionary
            self.add_hit_to_multiplet_dict(hit)

        # create multiplets from sorted particle dictionary
        self.create_multiplets()

    def information_about_tracking_data(self):
        """Prints information about signal and background hits.
        """
        print(f'Number of signal hits: {self.num_signal_hits}')
        print(f'Number of background hits: {self.num_background_hits}\n')
