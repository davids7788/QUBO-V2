import numpy as np
import time

from pattern.triplet import Triplet
from pattern.doublet import Doublet
from utility.time_tracking import hms_string
from math_functions.geometry import angle_based_measure
from pattern_building.segment_manager import SegmentManager


class QuboCoefficients:

    def __init__(self,
                 configuration: dict,
                 save_to_folder: str):
        """Class for handling and setting QUBO coefficients
        :param configuration: information needed for the segments, delivered by loading yaml file as a nested
        :param save_to_folder: folder to store results
        """
        self.configuration = configuration
        self.save_to_folder = save_to_folder

        # storing triplets
        self.triplet_list = set()

        # Dictionary of quality and conflict functions
        self.match_mode = None
        self.conflict_mode = None
        self.quality_mode = None

        self.match = None
        self.conflict = None
        self.quality = None
        self.parameter_setting()

    def parameter_setting(self) -> None:
        """Reads in the parameter setting and sets the fields accordingly.
        """
        b_ij_conflict = self.configuration["qubo parameters"]["b_ij conflict"]
        b_ij_match = self.configuration["qubo parameters"]["b_ij match"]
        a_i_quality = self.configuration["qubo parameters"]["a_i"]

        if b_ij_match == 'default':
            self.match = angle_based_measure
            self.match_mode = 'default'
        else:
            self.match = b_ij_match
            self.match_mode = 'constant'

        if b_ij_conflict == "default":
            self.conflict = angle_based_measure
            self.conflict_mode = 'default'
        else:
            self.conflict = b_ij_conflict
            self.conflict_mode = 'constant'

        if a_i_quality == 'default':
            self.quality = angle_based_measure
            self.quality_mode = 'default'
        else:
            self.quality = a_i_quality
            self.quality_mode = 'constant'

    def set_triplet_coefficients(self,
                                 segment_manager: SegmentManager) -> None:
        """Sets the triplet coefficients according to the configuration files. If a (re-)normalization was
        set it is also applied. If the process is successful a message and the target folder location containing the
        triplet list is displayed.
        :param segment_manager: segment manager object
        """
        print("Calculate triplet coefficients a_i and b_ij...")
        set_triplet_coefficients_start = time.process_time()

        num_layers = len(segment_manager.segment_storage.keys())

        for layer in range(num_layers - 2):
            for segment in segment_manager.segment_storage[layer]:
                if segment.name not in segment_manager.segment_mapping.keys():
                    continue
                next_segments = segment_manager.target_segments(segment.name)  # target segments
                for t1 in segment.triplet_data:
                    if self.quality_mode == 'constant':
                        t1.quality = self.quality
                    elif self.quality_mode == "default":
                        t1.quality = angle_based_measure([[t1.hit_1.x, t1.hit_2.x, t1.hit_3.x],
                                                          [t1.hit_1.y, t1.hit_2.y, t1.hit_3.y],
                                                          [t1.hit_1.z, t1.hit_2.z, t1.hit_3.z]])

                    for target_segment in next_segments + [segment]:   # checking all combinations with other triplets
                        for t2 in target_segment.triplet_data:
                            interaction_value = self.triplet_interaction(t1, t2)
                            # Only interactions != 0 are treated
                            if interaction_value == 0:
                                continue
                            t1.interactions.update({t2.triplet_id: interaction_value})
                            t2.interactions.update({t1.triplet_id: interaction_value})

        # filling list structure
        for layer in range(num_layers - 2):
            for segment in segment_manager.segment_storage[layer]:
                for triplet in segment.triplet_data:
                    self.triplet_list.add(triplet)

        set_triplet_coefficients_end = time.process_time()
        apply_coefficients_time = hms_string(set_triplet_coefficients_end - set_triplet_coefficients_start)
        self.triplet_list = list(self.triplet_list)
        print(f"Time elapsed for setting triplet coefficients: "
              f"{apply_coefficients_time}\n")

    def coefficient_rescaling(self) -> None:
        """Rescaling parameters according to the config file.
        """
        print("Starting rescaling of a_i and b_ij parameters...")
        # additional processing of qubo parameter

        rescale_triplet_coefficients_start = time.process_time()

        if self.configuration["scale range parameters"]["z_scores"]:
            mu = np.mean([t.quality for t in self.triplet_list])
            sigma = np.std([t.quality for t in self.triplet_list])
            for triplet in self.triplet_list:
                triplet.quality = (triplet.quality - mu) / sigma

        # Scaling in the following way to [a, b] : X' = a + (X - X_min) (b - a) / (X_max - X_min)
        if self.configuration["scale range parameters"]["quality"] is not None:

            a = self.configuration["scale range parameters"]["quality"][0]
            b = self.configuration["scale range parameters"]["quality"][1]
            min_quality = min([t.quality for t in self.triplet_list])
            max_quality = max([t.quality for t in self.triplet_list])
            for triplet in self.triplet_list:
                triplet.quality = a + (triplet.quality - min_quality) * (b - a) / (max_quality - min_quality)

        # scaling connectivity
        if self.configuration["scale range parameters"]["interaction"] is not None:
            conflict_term = float(self.configuration["qubo parameters"]["b_ij conflict"])
            # excluding conflict terms
            connectivity_values = [interaction
                                   for t in self.triplet_list
                                   for interaction in t.interactions.values() if interaction <= 0]

            min_connectivity = min(connectivity_values)
            max_connectivity = max(connectivity_values)
            range_connectivity = max_connectivity - min_connectivity

            a = self.configuration["scale range parameters"]["interaction"][0]
            b = self.configuration["scale range parameters"]["interaction"][1]

            for triplet in self.triplet_list:
                for key in triplet.interactions.keys():
                    if triplet.interactions[key] == conflict_term:
                        continue
                    else:
                        triplet.interactions[key] = a + (triplet.interactions[key] - min_connectivity) * (b - a) / \
                                                    range_connectivity

        rescale_triplet_coefficients_end = time.process_time()
        rescale_coefficients_time = hms_string(rescale_triplet_coefficients_end - rescale_triplet_coefficients_start)

        print(f"Time elapsed for rescaling triplet coefficients: "
              f"{rescale_coefficients_time}\n")

        print("Finished setting and rescaling of parameters.\n")
        print(f"Save triplet list...")
        np.save(f"{self.save_to_folder}/triplet_list", self.triplet_list)

    def triplet_interaction(self,
                            t1: Triplet,
                            t2: Triplet) -> float:
        """Compares two triplets and  how they match.
        :param t1: first triplet
        :param t2: second triplet
        :return
            value based on connectivity/conflict and chosen set of parameters
        """
        # checking number of shared hits

        t1_set = {t1.hit_1.hit_id,
                  t1.hit_2.hit_id,
                  t1.hit_3.hit_id}
        t2_set = {t2.hit_1.hit_id,
                  t2.hit_2.hit_id,
                  t2.hit_3.hit_id}
        intersection = len(t1_set.intersection(t2_set))

        # same and not interacting triplets get a zero as a coefficient
        if intersection in [0, 3]:
            return 0

        if intersection == 1:
            if t1.hit_3 == t2.hit_1:
                return - 1 + self.match([[t1.hit_1.x, t1.hit_2.x, t1.hit_3.x, t2.hit_2.x, t2.hit_3.x],
                                         [t1.hit_1.y, t1.hit_2.y, t1.hit_3.y, t2.hit_2.y, t2.hit_3.x],
                                         [t1.hit_1.z, t1.hit_2.z, t1.hit_3.z, t2.hit_2.z, t2.hit_3.x]])
            if self.conflict_mode == 'constant':
                return self.conflict
            else:
                return 1 + self.conflict([[t1.hit_1.x, t1.hit_2.x, t1.hit_3.x, t2.hit_1.x, t2.hit_2.x, t2.hit_3.x],
                                          [t1.hit_1.y, t1.hit_2.y, t1.hit_3.y, t2.hit_1.y, t2.hit_2.y, t2.hit_3.y],
                                          [t1.hit_1.z, t1.hit_2.z, t1.hit_3.z, t2.hit_1.z, t2.hit_2.z, t2.hit_3.z]])

        if intersection == 2:
            # triplets on same layer always have a conflict
            if t1.hit_2 == t2.hit_1 and t1.hit_3 == t2.hit_2:
                if self.match_mode == 'constant':
                    return self.match
                else:
                    return - 1 + self.match([[t1.hit_1.x, t1.hit_2.x, t1.hit_3.x, t2.hit_3.x],
                                             [t1.hit_1.y, t1.hit_2.y, t1.hit_3.y, t2.hit_3.y],
                                             [t1.hit_1.z, t1.hit_2.z, t1.hit_3.z, t2.hit_3.z]])

            else:
                if self.conflict_mode == 'constant':
                    return self.conflict
                else:
                    return 1 + self.conflict([[t1.hit_1.x, t1.hit_2.x, t1.hit_3.x, t2.hit_3.x],
                                              [t1.hit_1.y, t1.hit_2.y, t1.hit_3.y, t2.hit_3.y],
                                              [t1.hit_1.z, t1.hit_2.z, t1.hit_3.z, t2.hit_3.z]])
