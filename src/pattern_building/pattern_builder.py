import time

import numpy as np

from math_functions.checks import is_valid_doublet, is_valid_triplet, x0_at_z_ref

from utility.time_tracking import hms_string
from utility.data_format_handler import load_data
from pattern.detector_hit import DetectorHit
from pattern.doublet import Doublet
from pattern.triplet import Triplet
from pattern_building.segment_manager import LUXESegmentManager


class PatternBuilder:
    """Class for pattern building."""
    def __init__(self,
                 configuration: dict):
        """Set fields.
        :param configuration: dictionary with pattern building configuration
        """
        self.configuration = configuration
        self.doublet_creation_time = None
        self.triplet_creation_time = None

        # signal particles
        self.particle_dict_signal = {}

        # signal hits
        self.signal_hits = 0

        # background particle
        self.particle_dict_background = {}

        # background hits
        self.background_hits = 0

        # undecided (blind sample)
        self.particle_dict_blinded = {}

        # undecided particles (blind sample)
        self.blinded_hits = 0

        # some values to check if computation successful
        self.all_truth_doublets = set()
        self.all_truth_triplets = set()
        self.found_correct_doublets = 0
        self.found_correct_triplets = 0
        self.found_doublets = 0
        self.found_triplets = 0

        # track level
        self.num_signal_tracks = 0

    def load_tracking_data(self,
                           tracking_data_file: str,
                           segment_manager: LUXESegmentManager,
                           tracking_data_format: str) -> None:
        """Loads tracking data from a file and stores and places the data in the corresponding segments.
        :param tracking_data_file: LUXE tracking data file
        :param segment_manager: SegmentManager object
        :param tracking_data_format: format of tracking data file
        """
        print('\n-----------------------------------\n')
        print(f"Using tracking data file {tracking_data_file.split('/')[-1]}\n"
              f"Placing data in segments ...\n")

        list_of_hits = load_data(tracking_data_file, tracking_data_format)

        for hit in list_of_hits:
            if hit.is_signal:
                self.fill_signal_particle_dictionary(hit)
                self.signal_hits += 1
            elif hit.is_signal is None:
                self.fill_blinded_particle_dictionary(hit)
                self.blinded_hits += 1
            else:
                self.fill_background_particle_dictionary(hit)
                self.background_hits += 1

            # adding particle to a segment, used to reduce combinatorial candidates
            segment = segment_manager.get_segment_at_known_xyz_value(hit)
            if segment:
                segment.data.append(hit)
            else:
                print('Segment not found. Please check location of hits or detector geometry file!')

        print(f'Number of signal hits found: {self.signal_hits}')
        print(f'Number of background hits found: {self.background_hits}')
        print(f'Number of blinded hits found: {self.blinded_hits}')

    def fill_signal_particle_dictionary(self,
                                        hit: DetectorHit) -> None:
        """Fill signal particle dict.
        :param hit: detector hit object
        """
        if hit.particle_id in self.particle_dict_signal.keys():
            self.particle_dict_signal[hit.particle_id].append(hit)
        else:
            self.particle_dict_signal.update({hit.particle_id: [hit]})

    def fill_blinded_particle_dictionary(self,
                                         hit: DetectorHit) -> None:
        """Fill unknown particle dict.
        :param hit: detector hit object
        """
        if hit.particle_id in self.particle_dict_blinded.keys():
            self.particle_dict_blinded[hit.particle_id].append(hit)
        else:
            self.particle_dict_blinded.update({hit.particle_id: [hit]})

    def fill_background_particle_dictionary(self,
                                            hit: DetectorHit):
        """Fill unknown particle dict.
        :param hit: detector hit object
        """
        if hit.particle_id in self.particle_dict_background.keys():
            self.particle_dict_background[hit.particle_id].append(hit)
        else:
            self.particle_dict_background.update({hit.particle_id: [hit]})

    def information_about_particle_tracks(self,
                                          z_position_layers: list[float],
                                          setup: str,
                                          sample_composition: str,
                                          min_track_length) -> None:
        """Prints information about how many particles interact at least once with a detector chip
        and how many complete tracks can be reconstructed.
        :param z_position_layers: z position of detector layers:
        :param setup 'full' or 'simplified'
        :param sample_composition: information about the tracking sample, signal, signal+background or blinded
        :param min_track_length: minimum number of hits required from a signal positron to be counted as track
        """

        max_layer_dist = None

        for key, values in self.particle_dict_signal.items():
            if len(values) >= min_track_length:
                self.num_signal_tracks += 1
            if setup == 'simplified':
                max_layer_dist = 1
            if setup == 'full':
                max_layer_dist = 2

            doublet_list = []
            for i, hit_1 in enumerate(values[0:-1]):
                for hit_2 in values[i + 1:]:
                    if 1 <= z_position_layers.index(hit_2.z) - z_position_layers.index(hit_1.z) <= max_layer_dist:
                        doublet_list.append(Doublet(hit_1, hit_2))
            for d in doublet_list:
                self.all_truth_doublets.add(d)

            for m, d1 in enumerate(doublet_list):
                for d2 in doublet_list[m + 1:]:
                    if d1.hit_2.hit_id == d2.hit_1.hit_id:
                        self.all_truth_triplets.add(Triplet(d1.hit_1, d2.hit_1, d2.hit_2))

        if sample_composition == 'blinded':
            print('Truth information about particles cannot be accessed in blinded example!')
        else:
            print(f"Number of signal particles with at least one hit: {len(self.particle_dict_signal.keys())}")
            print(f"Number of complete signal tracks: {self.num_signal_tracks}\n")

    def create_doublet_list(self,
                            segment_manager: LUXESegmentManager) -> None:
        """Creates the doublets. For LUXE detector model only.
        :param segment_manager: SegmentManager object with already set segments and mapping
        """
        num_layers = len(segment_manager.segment_storage.keys())
        for layer in range(num_layers - 1):
            for segment in segment_manager.segment_storage[layer]:
                if segment.name not in segment_manager.segment_mapping.keys():
                    continue
                next_segments = segment_manager.target_segments(segment.name)  # target segments

                for first_hit in segment.data:
                    for target_segment in next_segments:
                        for second_hit in target_segment.data:
                            if not is_valid_doublet(first_hit,
                                                    second_hit,
                                                    segment_manager.z_position_to_layer[0],
                                                    self.configuration):
                                continue
                            doublet = Doublet(first_hit, second_hit)
                            segment.doublet_data.append(doublet)
                            self.found_doublets += 1

                            if doublet.from_same_particle():
                                self.found_correct_doublets += 1

    def create_triplet_list(self,
                            segment_manager: LUXESegmentManager) -> None:
        """Creates the doublets. For LUXE detector model only.
        :param segment_manager: SegmentManager object with already set segments and mapping
        """
        num_layers = len(segment_manager.segment_storage.keys())
        for layer in range(num_layers - 2):
            for segment in segment_manager.segment_storage[layer]:
                if segment.name not in segment_manager.segment_mapping.keys():
                    continue
                next_segments = segment_manager.target_segments(segment.name)  # target segment list
                for target_segment in next_segments:
                    for first_doublet in segment.doublet_data:
                        for second_doublet in target_segment.doublet_data:
                            if first_doublet.hit_2.hit_id != second_doublet.hit_1.hit_id:  # check if match
                                continue
                            if not is_valid_triplet(first_doublet.hit_1,
                                                    first_doublet.hit_2,
                                                    second_doublet.hit_2,
                                                    self.configuration["triplet"]["max scattering"]):
                                continue
                            triplet = Triplet(first_doublet.hit_1,
                                              first_doublet.hit_2,
                                              second_doublet.hit_2)
                            self.found_triplets += 1
                            segment.triplet_data.append(triplet)

                            # filling lists for statistical purposes
                            if triplet.from_same_particle():
                                self.found_correct_triplets += 1

    def create_multiplets(self,
                          segment_manager: LUXESegmentManager) -> None:
        """Creates multiplets. For LUXE detector model only.
        :param segment_manager: SegmentManager object with already set segments and mapping
        """
        print("-----------------------------------\n")
        print("Forming doublets ...\n")

        doublet_list_start = time.process_time()
        self.create_doublet_list(segment_manager)
        doublet_list_end = time.process_time()

        self.doublet_creation_time = hms_string(doublet_list_end - doublet_list_start)
        print(f"Time elapsed for forming doublets: "
              f"{self.doublet_creation_time}")
        print(f"Number of doublets found: {self.found_doublets}")
        if len(self.all_truth_doublets) > 0:
            print(f"Doublet selection efficiency: "
                  f"{np.around(100 * self.found_correct_doublets / len(self.all_truth_doublets), 2)} %\n")
        else:
            print(f"Doublet selection efficiency cannot be calculated from blinded sample\n")

        print("-----------------------------------\n")
        print("Forming triplets ...\n")
        list_triplet_start = time.process_time()
        self.create_triplet_list(segment_manager)
        list_triplet_end = time.process_time()
        self.triplet_creation_time = hms_string(list_triplet_end - list_triplet_start)

        print(f"Time elapsed for forming triplets: "
              f"{self.triplet_creation_time}")
        print(f"Number of triplets found: {self.found_triplets}")
        if len(self.all_truth_doublets) > 0:
            print(f"Triplets selection efficiency: "
                  f"{np.around(100 * self.found_correct_triplets / len(self.all_truth_triplets), 2)} %\n")
        else:
            print(f"Triplets selection efficiency cannot be calculated from blinded sample\n")
