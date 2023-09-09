import csv
import time

import numpy as np

from math_functions.checks import jit_is_valid_triplet, dxy_x0_check
from math_functions.geometry import x0_at_z_ref, xyz_angle
from utility.time_tracking import hms_string
from pattern.detector_hit import DetectorHit
from pattern.doublet import Doublet
from pattern.triplet import Triplet
from pattern_building.segment_manager import SegmentManager


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

        # only signal particles
        self.particle_dict_signal = {}

        # background particle
        self.particle_dict_background = {}

        # undecided (blind example)
        self.particle_dict_blinded = {}

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
                           segment_manager: SegmentManager,
                           simulation_tool: str) -> None:
        """Loads tracking data from a file and stores and places the data in the corresponding segments.
        :param tracking_data_file: LUXE tracking data file
        :param segment_manager: SegmentManager object
        :param simulation_tool: simplified_simulation_csv or key4hep_csv
        """
        print('\n-----------------------------------\n')
        print(f"Using tracking data file {tracking_data_file.split('/')[-1]}\n"
              f"Placing data in segments ...\n")

        list_of_hits = None   # Initialising variable for list of detector hits
        if simulation_tool == 'simplified_simulation':
            list_of_hits = PatternBuilder.load_tracking_data_from_simplified_simulation_csv(tracking_data_file)
        elif simulation_tool == 'key4hep':
            if '.csv' in tracking_data_file:
                list_of_hits = PatternBuilder.load_tracking_data_from_key4hep_csv(tracking_data_file)
            else:
                pass   # slcio implementation
        else:
            print('Loading data from given file format not implemented yet!')
            print('Exiting...')
            exit()

        signal_hits = 0
        background_hits = 0
        blinded_hits = 0

        for hit in list_of_hits:
            if hit.is_signal:
                self.fill_signal_particle_dictionary(hit)
                signal_hits += 1
            elif hit.is_signal is None:
                self.fill_blinded_particle_dictionary(hit)
                blinded_hits += 1
            else:
                self.fill_background_particle_dictionary(hit)
                background_hits += 1

            # adding particle to a segment, used to reduce combinatorial candidates
            segment = segment_manager.get_segment_at_known_xyz_value(hit.x,
                                                                     hit.y,
                                                                     hit.z)
            if segment:
                segment.data.append(hit)

        print(f'Number of signal hits found: {signal_hits}')
        print(f'Number of background hits found: {background_hits}')
        print(f'Number of blinded hits found: {blinded_hits}')

    @staticmethod
    def load_tracking_data_from_simplified_simulation_csv(tracking_data_file: str) -> list[DetectorHit]:
        """Loads tracking data from .csv file and returns a list of DetectorHit objects created from the file entries.
        :param tracking_data_file: .csv tracking data file

        :return:
            list of DetectorHit objects
        """
        detector_hits = []
        with open(tracking_data_file, 'r') as file:
            csv_reader = csv.reader(file)
            _ = next(csv_reader)  # access header, csv files should consist of one line of header

            for row in csv_reader:
                detector_hits.append(DetectorHit(row))

        return detector_hits

    @staticmethod
    def load_tracking_data_from_key4hep_csv(tracking_data_file: str) -> list[DetectorHit]:
        """Loads tracking data from .csv file and returns a list of DetectorHit objects created from the file entries.
        :param tracking_data_file: .csv tracking data file

        :return:
            list of DetectorHit objects
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

        detector_hits = []
        with open(tracking_data_file, 'r') as file:
            csv_reader = csv.reader(file)
            _ = next(csv_reader)  # access header, csv files should consist of one line of header

            for row in csv_reader:
                row_trimmed = [row[fieldnames_key4hep_csv.index('index')],              # hit_id
                               np.round(1e-3 * float(row[fieldnames_key4hep_csv.index('tx')]), 3),   # x-coordinate
                               np.round(1e-3 * float(row[fieldnames_key4hep_csv.index('ty')]), 3),   # y-coordinate
                               np.round(1e-3 * float(row[fieldnames_key4hep_csv.index('tz')]), 3),   # z-coordinate
                               row[fieldnames_key4hep_csv.index('layer_id')],           # layer_id
                               row[fieldnames_key4hep_csv.index('module_id')],          # module_id
                               int(get_cell_id(row), base=2),                           # cell_id for ACTS tracking
                               row[fieldnames_key4hep_csv.index('particle_id')],        # particle_id
                               True,                                                    # not provided in key4hep csv
                               row[fieldnames_key4hep_csv.index('te')],                 # energy of particle
                               row[fieldnames_key4hep_csv.index('tt')]]                 # time of particle
                detector_hits.append(DetectorHit(row_trimmed))

        return detector_hits

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
            print('Blinded example, no truth information about signal or background are provided!')
        else:
            print(f"Number of signal particles with at least one hit: {len(self.particle_dict_signal.keys())}")
            print(f"Number of complete signal tracks: {self.num_signal_tracks}\n")

    def is_valid_doublet_candidate(self,
                                   first_hit: DetectorHit,
                                   second_hit: DetectorHit,
                                   z_ref: float) -> bool:
        """Checks if two hits are actually a doublet candidate.
        :param first_hit: hit with lower z-value
        :param second_hit: hit with higher z-value
        :param z_ref: reference layer z-value

        :return
            True if valid doublet candidate, else False
        """

        # calculate x0
        x0 = x0_at_z_ref(second_hit.x,
                         first_hit.x,
                         second_hit.z,
                         first_hit.z,
                         z_ref)

        # check dx / x0 criteria
        if not dxy_x0_check(first_hit.x,
                            second_hit.x,
                            first_hit.z,
                            second_hit.z,
                            x0,
                            criteria_mean=self.configuration["doublet"]["dx/x0"],
                            criteria_eps=self.configuration["doublet"]["dx/x0 eps"]):
            return False

        # check dy / x0 criteria
        if not dxy_x0_check(first_hit.y,
                            second_hit.y,
                            first_hit.z,
                            second_hit.z,
                            x0,
                            criteria_mean=self.configuration["doublet"]["dy/x0"],
                            criteria_eps=self.configuration["doublet"]["dy/x0 eps"]):
            return False

        return True

    def create_doublet_list(self,
                            segment_manager: SegmentManager) -> None:
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
                            is_valid_doublet = self.is_valid_doublet_candidate(first_hit,
                                                                               second_hit,
                                                                               segment_manager.z_position_to_layer[0])
                            if not is_valid_doublet:
                                continue
                            doublet = Doublet(first_hit, second_hit)
                            segment.doublet_data.append(doublet)
                            self.found_doublets += 1

                            if doublet.from_same_particle():
                                self.found_correct_doublets += 1

    def create_triplet_list(self,
                            segment_manager: SegmentManager) -> None:
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
                            if self.is_valid_triplet(first_doublet.hit_1, first_doublet.hit_2, second_doublet.hit_2):
                                triplet = Triplet(first_doublet.hit_1, first_doublet.hit_2, second_doublet.hit_2)
                                self.found_triplets += 1
                                segment.triplet_data.append(triplet)

                                # filling lists for statistical purposes
                                if triplet.from_same_particle():
                                    self.found_correct_triplets += 1

    def create_x_plets_LUXE(self,
                            segment_manager: SegmentManager) -> None:
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
        print(f"Doublet selection efficiency: "
              f"{np.around(100 * self.found_correct_doublets / len(self.all_truth_doublets), 2)} %\n")

        print("-----------------------------------\n")
        print("Forming triplets ...\n")
        list_triplet_start = time.process_time()
        self.create_triplet_list(segment_manager)
        list_triplet_end = time.process_time()
        self.triplet_creation_time = hms_string(list_triplet_end - list_triplet_start)

        print(f"Time elapsed for forming triplets: "
              f"{self.triplet_creation_time}")
        print(f"Number of triplets found: {self.found_triplets}")
        print(f"Triplet selection efficiency: "
              f"{np.around(100 * self.found_correct_triplets / len(self.all_truth_triplets), 2)} %\n")
        print("-----------------------------------\n")

    def is_valid_triplet(self,
                         hit_1: DetectorHit,
                         hit_2: DetectorHit,
                         hit_3: DetectorHit) -> bool:
        """Checks if doublets may be combined to a triplet, depending on the doublet angles -> scattering
        :param hit_1: first hit of a triplets
        :param hit_2: second hit of a triplets
        :param hit_3: third hit of a triplets

        :return:
            True if criteria applies, else False
        """
        xz_12 = xyz_angle(hit_1.x, hit_2.x, hit_1.z, hit_2.z)

        xz_23 = xyz_angle(hit_2.x, hit_3.x, hit_2.z, hit_3.z)

        yz_12 = xyz_angle(hit_1.y, hit_2.y, hit_1.z, hit_2.z)

        yz_23 = xyz_angle(hit_2.y, hit_3.y, hit_2.z, hit_3.z)

        return jit_is_valid_triplet(xz_12,
                                    xz_23,
                                    yz_12,
                                    yz_23,
                                    self.configuration["triplet"]["max scattering"])

    def write_info_file(self,
                        target_folder) -> None:
        """Writes information about the Preselection parameters and some statistics into
        'preselection_info.txt' which is stored inside the output folder.
        """
        with open(target_folder + "/preselection_info.txt", "w") as f:
            f.write("Preselection performed with the following configuration: \n")
            for outer_key in self.configuration.keys():
                for inner_key, value in self.configuration[outer_key].items():
                    f.write(f"\n{inner_key}: {value}")
                f.write("\n")
            f.write("\n\n")
            f.write("---\n")
            f.write(f"Number of particles with at least one hit on the detector: "
                    f"{len(self.particle_dict_signal.keys())}\n")
            f.write(f"Number of generated tracks: {self.num_signal_tracks}\n")
            f.write("\n\n")

            f.write(f"Time elapsed for forming doublets: "
                    f"{self.doublet_creation_time}\n")
            f.write(f"Number of doublets found: {self.found_doublets}\n")
            f.write(f"Doublet selection efficiency: "
                    f"{np.around(100 * self.found_correct_doublets / len(self.all_truth_doublets), 3)} %\n")
            f.write("\n\n")

            f.write(f"Time elapsed for creating triplets: "
                    f"{self.triplet_creation_time}\n")
            f.write(f"Number of triplets found: {self.found_triplets}\n")
            f.write(f"Triplet selection efficiency: "
                    f"{np.around(100 * self.found_correct_triplets / len(self.all_truth_triplets), 3)} %\n")
