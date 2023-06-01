import csv
import time

import numpy as np

from math_functions.checks import w_is_valid_triplet, dxy_x0_check
from math_functions.geometry import x0_at_z_ref, xyz_angle
from utility.time_tracking import hms_string
from pattern.detector_hit import DetectorHit
from pattern.doublet import Doublet
from pattern.triplet import Triplet
from pattern_building.segment_manager import SegmentManager


class PatternBuilder:
    def __init__(self,
                 configuration: dict,
                 save_to_folder: str):
        """Class for creating x_plets from LUXE detector hits
        :param configuration: information needed for the creating segmentation list,
              provided as a nested dictionary
        :param save_to_folder : folder in which results are stored
        """
        self.configuration = configuration
        self.save_to_folder = save_to_folder
        self.doublet_creation_time = None
        self.triplet_creation_time = None
        self.num_particles = 0
        self.particle_dict = {}

        # some values to check if computation successful
        self.num_complete_tracks = 0
        self.num_all_doublets = 0
        self.num_all_triplets = 0
        self.found_correct_doublets = 0
        self.found_correct_triplets = 0
        self.found_doublets = 0
        self.found_triplets = 0

    def load_tracking_data(self,
                           tracking_data_file: str,
                           segment_manager: SegmentManager) -> None:
        """Loads data from a .csv file and stores it into a 2-dim array. Also converts values which are used for
        calculations to float from string. The file structure is displayed in the class description.
        :param tracking_data_file: Luxe tracking data file
        :param segment_manager: SegmentManager object with already set segments and mapping
        """
        print(f"Using tracking data file {tracking_data_file.split('/')[-1]}\n"
              f"Distributing data into segments ...\n")
        with open(tracking_data_file, 'r') as file:
            csv_reader = csv.reader(file)
            _ = next(csv_reader)  # access header, csv files should consist of one line of header

            for row in csv_reader:
                # Create new DetectorHit object
                hit = DetectorHit(row)

                # check if particle ID was already seen
                if hit.particle_id in self.particle_dict.keys():
                    self.particle_dict[hit.particle_id].append(hit)
                else:
                    self.particle_dict.update({hit.particle_id: [hit]})

                # adding particle to a segment, used to reduce combinatorial candidates
                segment = segment_manager.get_segment_at_known_xyz_value(hit.x,
                                                                         hit.y,
                                                                         hit.z)
                if segment:
                    segment.data.append(hit)

        if segment_manager.setup == 'full':
            self.information_about_particle_tracks(first_layer=[segment_manager.z_position_to_layer[0],
                                                                segment_manager.z_position_to_layer[1]],
                                                   last_layer=[segment_manager.z_position_to_layer[-2],
                                                               segment_manager.z_position_to_layer[-1]])
        elif segment_manager.setup == 'simplified':
            self.information_about_particle_tracks(first_layer=[segment_manager.z_position_to_layer[0]],
                                                   last_layer=[segment_manager.z_position_to_layer[-1]])

    def information_about_particle_tracks(self,
                                          first_layer: list[float],
                                          last_layer: list[float]) -> None:
        """Prints information about how many particles interact at least once with a detector chip
        and how many complete tracks can be reconstructed.
        :param first_layer: depending on the setup the first layer consists of two overlapping chips  or one
        :param last_layer: depending on the setup the last layer consists of two overlapping chips or one
        """
        for key, value in self.particle_dict.items():
            last = False
            first = False
            for entry in value:
                if entry.z in last_layer:
                    last = True
                if entry.z in first_layer:
                    first = True
            self.num_particles += 1
            if last and first:
                self.num_complete_tracks += 1
            self.num_all_doublets += len(value) - 1
            self.num_all_triplets += len(value) - 2

        print(f"Number of particles with at least one hit: {self.num_particles}")
        print(f"Number of complete tracks: {self.num_complete_tracks}\n")

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
        x0 = x0_at_z_ref(first_hit.x,
                         second_hit.x,
                         first_hit.z,
                         second_hit.z,
                         z_ref)

        # check dx / x0 criteria
        if not dxy_x0_check(first_hit.x,
                            second_hit.x,
                            x0,
                            criteria_mean=self.configuration["doublet"]["dx/x0"],
                            criteria_eps=self.configuration["doublet"]["dx/x0 eps"]):
            return False

        # check dy / x0 criteria
        if not dxy_x0_check(first_hit.y,
                            second_hit.y,
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

                            if doublet.is_correct_match():
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
                                if triplet.is_correct_match():
                                    self.found_correct_triplets += 1
                segment.doublet_data.clear()   # --> lower memory usage, num doublets are not needed anymore

    def create_x_plets_LUXE(self,
                            segment_manager: SegmentManager) -> None:
        """Creates xplets. For LUXE detector model only.
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
              f"{np.around(100 * self.found_correct_doublets / self.num_all_doublets, 2)} %\n")

        print("-----------------------------------\n")
        print("Forming triplets ...\n")
        list_triplet_start = time.process_time()
        self.create_triplet_list(segment_manager)
        list_triplet_end = time.process_time()
        self.triplet_creation_time = hms_string(list_triplet_end - list_triplet_start)

        print(f"Time elapsed for  forming triplets: "
              f"{self.triplet_creation_time}")
        print(f"Number of triplets found: {self.found_triplets}")
        print(f"Triplet selection efficiency: "
              f"{np.around(100 * self.found_correct_triplets / self.num_all_triplets, 2)} %\n")
        print("-----------------------------------\n")

    def is_valid_triplet(self,
                         hit_1: DetectorHit,
                         hit_2: DetectorHit,
                         hit_3: DetectorHit):
        """Checks if doublets may be combined to a triplet, depending on the doublet angles -> scattering
        :param hit_1: first hit of a triplets
        :param hit_2: second hit of a triplets
        :param hit_3: third hit of a triplets
            True if criteria applies, else False
        """
        xz_12 = xyz_angle(hit_1.x, hit_2.x, hit_1.z, hit_2.z)

        xz_23 = xyz_angle(hit_2.x, hit_3.x, hit_2.z, hit_3.z)

        yz_12 = xyz_angle(hit_1.y, hit_2.y, hit_1.z, hit_2.z)

        yz_23 = xyz_angle(hit_2.y, hit_3.y, hit_2.z, hit_3.z)

        return w_is_valid_triplet(xz_12,
                                  xz_23,
                                  yz_12,
                                  yz_23,
                                  self.configuration["triplet"]["max scattering"])

    def write_info_file(self) -> None:
        """Writes information about the Preselection parameters and some statistics into
        'preselection_info.txt' which is stored inside the output folder.
        """
        with open(self.save_to_folder + "/preselection_info.txt", "w") as f:
            f.write("Preselection performed with the following configuration: \n")
            for outer_key in self.configuration.keys():
                for inner_key, value in self.configuration[outer_key].items():
                    f.write(f"\n{inner_key}: {value}")
                f.write("\n")
            f.write("\n\n")
            f.write("---\n")
            f.write(f"Number of particles with at least one hit on the detector:: {self.num_particles}\n")
            f.write(f"Number of generated tracks: {self.num_complete_tracks}\n")
            f.write("\n\n")

            f.write(f"Time elapsed for forming doublets: "
                    f"{self.doublet_creation_time}\n")
            f.write(f"Number of doublets found: {self.found_doublets}\n")
            f.write(f"Doublet selection efficiency: "
                    f"{np.around(100 * self.found_correct_doublets / self.num_all_doublets,  3)} %\n")
            f.write("\n\n")

            f.write(f"Time elapsed for creating triplets: "
                    f"{self.triplet_creation_time}\n")
            f.write(f"Number of triplets found: {self.found_triplets}\n")
            f.write(f"Triplet selection efficiency: "
                    f"{np.around(100 * self.found_correct_triplets / self.num_all_triplets, 3)} %\n")
