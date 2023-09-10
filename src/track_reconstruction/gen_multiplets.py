import csv
import numpy as np
import os

from pattern.detector_hit import DetectorHit
from pattern.multiplet import Multiplet


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

    def gen_multiplets_key4hep_csv(self) -> None:
        """Creates all signal multiplets from the specified key4hep simulation tracking data file
        """
        fieldnames_key4hep_csv = ['particle_id',
                                  'geometry_id',
                                  'system_id',
                                  'layer_id',
                                  'side_id',
                                  'module_id',
                                  'sensor_id',
                                  'tx',
                                  'ty',
                                  'tz',
                                  'tt',
                                  'tpx',
                                  'tpy',
                                  'tpz',
                                  'te',
                                  'deltapx',
                                  'deltapy',
                                  'deltapz',
                                  'deltae',
                                  'index']

        def get_cell_id(row_for_cell_id: list[str]) -> str:
            """Returns the cell id, used for ACTS Kalman Filter for LUXE inside key4hep environment for.
            :param row_for_cell_id: list of str with information about particle hit in detector

            :return
                concatenated string of module and stave to identify region where hit stems from"""
            return f'{bin(int(row_for_cell_id[fieldnames_key4hep_csv.index("module_id")])).zfill(3)[2:]}' \
                   f'{bin(int(row_for_cell_id[fieldnames_key4hep_csv.index("layer_id")])).zfill(3)[2:]}10'

        with open(self.tracking_data_file, 'r') as file:
            csv_reader = csv.reader(file)
            _ = next(csv_reader)  # access header, csv files should consist of one line of header

            for row in csv_reader:
                # information needs to be trimmed to adjust for the DetectorHit class model
                row_trimmed = [row[fieldnames_key4hep_csv.index('index')],                           # hit_id
                               np.round(1e-3 * float(row[fieldnames_key4hep_csv.index('tx')]), 3),   # x-coordinate
                               np.round(1e-3 * float(row[fieldnames_key4hep_csv.index('ty')]), 3),   # y-coordinate
                               np.round(1e-3 * float(row[fieldnames_key4hep_csv.index('tz')]), 3),   # z-coordinate
                               row[fieldnames_key4hep_csv.index('layer_id')],        # layer_id
                               row[fieldnames_key4hep_csv.index('module_id')],       # module_id
                               int(get_cell_id(row), base=2),                        # cell_id for ACTS tracking
                               row[fieldnames_key4hep_csv.index('particle_id')],     # particle_id
                               True,                                                 # not provided in key4hep csv
                               row[fieldnames_key4hep_csv.index('te')],              # energy of particle
                               row[fieldnames_key4hep_csv.index('tt')]]              # time of particle
                hit = DetectorHit(row_trimmed)

                # add to multiplet dictionary
                self.add_hit_to_multiplet_dict(hit)

        # create multiplets from sorted particle dictionary
        self.create_multiplets()

    def gen_multiplets_simplified_LUXE(self) -> None:
        """Creates all signal multiplets from the specified simplified simulation tracking data file.
        Assumes, that the structure of the csv file from the simplified simulation matches the DetectorHit class model.
        """

        with open(self.tracking_data_file, 'r') as file:
            csv_reader_tracking_file = csv.reader(file)
            _ = next(csv_reader_tracking_file)

            for row in csv_reader_tracking_file:
                # Create new DetectorHit object
                hit = DetectorHit(row)

                # add to multiplet dictionary
                self.add_hit_to_multiplet_dict(hit)

        # create multiplets from sorted particle dictionary
        self.create_multiplets()

    def information_about_tracking_data(self):
        """Prints information about signal and background hits.
        """
        print(f'Number of signal hits: {self.num_signal_hits}')
        print(f'Number of background hits: {self.num_background_hits}\n')
