import numpy as np
import matplotlib.pyplot as plt

from preselection.connection_metrics import *
from preselection.segment_manager import SegmentManager


class QuboCoefficients:

    def __init__(self,
                 configuration: dict,
                 save_to_folder: str):
        """Class for handling and setting QUBO coefficients
        :param configuration: information needed for the segments, delivered by loading yaml file as a nested
               python dictionary:
                {
                doublet: {dx/x0: float,
                          eps:   float>,
                          dy/y0: float},
                triplet: {max scattering: float}
                binning: {num bins x: int,
                          num bins y: int}
                qubo parameters: {b_ij conflict: float
                                  b_ij match: <name of an implemented function> or float,
                                  a_i: <name of an implemented function> or float}
                scale range parameters: {z_scores: True or False,
                                         quality: null (= None) or [float, float],
                                         interaction: null (= None) or [float, float]}
                }
        :param save_to_folder: folder to store results
        """
        self.configuration = configuration
        self.mode = self.recognize_parameter_setup()
        self.save_to_folder = save_to_folder

        # storing triplets
        self.triplet_list = set()

        # Keeping track of coefficients and saving them for statistic purposes
        self.quality_wrong_match_list = []
        self.connectivity_wrong_match_list = []
        self.quality_correct_match_list = []
        self.connectivity_correct_match_list = []

        # Dictionary of quality and conflict functions
        self.quality_functions = {}
        self.conflict_functions = {}

        # Adding a very simple set of quality and conflict functions
        self.add_quality_function("two norm angle standard deviation", two_norm_std_angle)
        self.add_conflict_function("two norm angle standard deviation", two_norm_std_angle)

    def set_triplet_coefficients(self,
                                 segment_manager: SegmentManager) -> None:
        """Sets the triplet coefficients according to the configuration files. If a (re-)normalization was
        set it is also applied. If the process is successful a message and the target folder location containing the
        triplet list is displayed.
        :param segment_manager: segment manager object
        """
        # self checking configuration
        print("Calculate triplet coefficients a_i and b_ij...")
        quality = None
        try:
            quality = float(self.configuration["qubo parameters"]["a_i"])
            quality_mode = "constant"
        except ValueError:
            quality_mode = "angle function"

        num_layers = len(segment_manager.segment_storage.keys())
        for layer in range(num_layers - 2):
            for segment in segment_manager.segment_storage[layer]:
                if segment.name not in segment_manager.segment_mapping.keys():
                    continue
                next_segments = segment_manager.target_segments(segment.name)  # target segments
                for t1 in segment.triplet_data:
                    if quality_mode == "constant":
                        t1.quality = quality
                    elif quality_mode == "angle function":
                        triplet_angles_xz, triplet_angles_yz = t1.angles_between_doublets()
                        t1.quality = np.sqrt(triplet_angles_xz ** 2 + triplet_angles_yz ** 2)

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
            segment.triplet_data.clear()

    def collecting_qubo_parameters(self) -> None:
        """Function for collecting and storing information about quality and interaction values,
        as well as truth information about the triplets.
        """
        self.triplet_list = list(self.triplet_list)
        t_mapping = {key: value for key, value in zip([t.triplet_id for t in self.triplet_list],
                                                      np.arange(len(self.triplet_list)))}

        print("Collect statistics about a_i and b_ij...")
        for i, t1 in enumerate(self.triplet_list):
            if t1.is_correct_match():
                self.quality_correct_match_list.append(t1.quality)
            else:
                self.quality_wrong_match_list.append(t1.quality)

            for i_key, i_value in t1.interactions.items():
                if t_mapping[i_key] > i:
                    continue
                t2 = self.triplet_list[t_mapping[i_key]]
                if t1.is_correct_match() and t2.is_correct_match():   # need to have 2 overlap hits, so this is enough
                    self.connectivity_correct_match_list.append(i_value)
                else:
                    self.connectivity_wrong_match_list.append(i_value)

    def coefficient_rescaling(self) -> None:
        """Rescaling parameters according to the config file.
        """
        print("Starting rescaling of a_i and b_ij parameters...")
        # additional processing of qubo parameter

        if self.configuration["scale range parameters"]["z_scores"]:
            mu = np.mean(self.quality_correct_match_list + self.quality_wrong_match_list)
            sigma = np.std(self.quality_correct_match_list + self.quality_wrong_match_list)
            for triplet in self.triplet_list:
                triplet.quality = (triplet.quality - mu) / sigma
            self.quality_correct_match_list = [(quality - mu) / sigma for quality in self.quality_correct_match_list]
            self.quality_wrong_match_list = [(quality - mu) / sigma for quality in self.quality_wrong_match_list]

        # Scaling in the following way to [a, b] : X' = a + (X - X_min) (b - a) / (X_max - X_min)
        if self.configuration["scale range parameters"]["quality"] is not None:

            a = self.configuration["scale range parameters"]["quality"][0]
            b = self.configuration["scale range parameters"]["quality"][1]
            min_quality = min(min(self.quality_correct_match_list), min(self.quality_wrong_match_list))
            max_quality = max(min(self.quality_correct_match_list), max(self.quality_wrong_match_list))
            for triplet in self.triplet_list:
                triplet.quality = a + (triplet.quality - min_quality) * (b - a) / (max_quality - min_quality)

            # rewriting a_i lists
            self.quality_correct_match_list = [a + (quality - min_quality) * (b - a) / (max_quality - min_quality)
                                               for quality in self.quality_correct_match_list]
            self.quality_wrong_match_list = [a + (quality - min_quality) * (b - a) / (max_quality - min_quality)
                                             for quality in self.quality_wrong_match_list]

        # scaling connectivity
        if self.configuration["scale range parameters"]["interaction"] is not None:
            conflict_term = float(self.configuration["qubo parameters"]["b_ij conflict"])
            # excluding conflict terms
            connectivity_values = [con for con in self.connectivity_correct_match_list if con != conflict_term] + \
                                  [con for con in self.connectivity_wrong_match_list if con != conflict_term]
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

            # rewriting connectivity lists
            self.connectivity_wrong_match_list = [a + (v - min_connectivity) * (b - a) / range_connectivity
                                                  for v in self.connectivity_wrong_match_list if v != conflict_term]
            self.connectivity_correct_match_list = [a + (v - min_connectivity) * (b - a) / range_connectivity
                                                    for v in self.connectivity_correct_match_list if v != conflict_term]

        print("\nFinished setting and rescaling of parameters.")
        print(f"\nSaving triplet list...")
        np.save(f"{self.save_to_folder}/triplet_list", self.triplet_list)

    def triplet_interaction(self,
                            triplet,
                            other_triplet):
        """Compares two triplets and  how they match.
        :param triplet: first triplet
        :param other_triplet: triplet to compare with the first one
        :return
            value based on connectivity/conflict and chosen set of parameters
        """
        # checking number of shared hits
        z_position_set = {triplet.doublet_1.hit_1_position[2],
                          triplet.doublet_1.hit_2_position[2],
                          triplet.doublet_2.hit_2_position[2],
                          other_triplet.doublet_1.hit_1_position[2],
                          other_triplet.doublet_1.hit_2_position[2],
                          other_triplet.doublet_2.hit_2_position[2]}
        if len(z_position_set) == 3:
            same_layer = True
        else:
            same_layer = False

        hit_list = {triplet.doublet_1.hit_1_position,
                    triplet.doublet_1.hit_2_position,
                    triplet.doublet_2.hit_2_position,
                    other_triplet.doublet_1.hit_1_position,
                    other_triplet.doublet_1.hit_2_position,
                    other_triplet.doublet_2.hit_2_position}
        intersection = 6 - len(hit_list)
        # same and not interacting triplets get a zero as a coefficient
        if intersection in [0, 3]:
            return 0

        # different mode types
        if self.mode[0] == "constant" and self.mode[1] == "constant":
            if intersection == 1:
                return float(self.configuration["qubo parameters"]["b_ij conflict"])
            elif intersection == 2 and not same_layer:
                return float(self.configuration["qubo parameters"]["b_ij match"])
            else:
                return float(self.configuration["qubo parameters"]["b_ij conflict"])

        # constant for conflict, quality for function
        elif self.mode[0] == "constant" and self.mode[1] != "constant":
            if intersection == 1:
                return float(self.configuration["qubo parameters"]["b_ij conflict"])
            elif intersection == 2 and not same_layer:
                return - 1 + self.quality_functions[self.configuration["qubo parameters"]["b_ij match"]](
                    triplet.doublet_1, triplet.doublet_2, other_triplet.doublet_2)
            else:
                return float(self.configuration["qubo parameters"]["b_ij conflict"])

        # function for conflict, constant for quality
        elif self.mode[0] != "constant" and self.mode[1] == "constant":
            if intersection == 1:
                return 1 + self.conflict_functions[self.configuration["qubo parameters"]["b_ij conflict"]](
                    triplet.doublet_1, triplet.doublet_2, other_triplet.doublet_2)
            elif intersection == 2 and not same_layer:
                return float(self.configuration["qubo parameters"]["b_ij match"])
            else:
                return 1 + self.conflict_functions[self.configuration["qubo parameters"]["b_ij conflict"]](
                    triplet.doublet_1, triplet.doublet_2, other_triplet.doublet_2)

        # only functions for conflicts / quality
        elif self.mode[0] != "constant" and self.mode[1] != "constant":
            if intersection == 1:
                return 1 + self.conflict_functions[self.configuration["qubo parameters"]["b_ij conflict"]](
                    triplet.doublet_1, triplet.doublet_2, other_triplet.doublet_2)
            elif intersection == 2 and not same_layer:
                return - 1 + self.quality_functions[self.configuration["qubo parameters"]["b_ij match"]](
                    triplet.doublet_1, triplet.doublet_2, other_triplet.doublet_2)
            else:
                return 1 + self.conflict_functions[self.configuration["qubo parameters"]["b_ij conflict"]](
                    triplet.doublet_1, triplet.doublet_2, other_triplet.doublet_2)

    def recognize_parameter_setup(self):
        """Recognizes the chosen parameter setup for the QUBO
        :return:
            tuple of strings with setup information
        """
        try:
            float(self.configuration["qubo parameters"]["b_ij conflict"])
            conflict = "constant"
        except ValueError:
            conflict = self.configuration["qubo parameters"]["b_ij conflict"]
        try:
            float(self.configuration["qubo parameters"]["b_ij match"])
            match = "constant"
        except ValueError:
            match = self.configuration["qubo parameters"]["b_ij match"]

        return conflict, match

    def add_quality_function(self,
                             quality_function_name,
                             quality_function_object):
        """Adds a function to the quality function dictionary
        :param quality_function_name: name of the function
        :param quality_function_object: a function object
        """
        self.quality_functions.update({quality_function_name: quality_function_object})

    def add_conflict_function(self,
                              conflict_function_name,
                              conflict_function_object):
        """Adds a function to the conflict function dictionary
        :param conflict_function_name: name of the function
        :param conflict_function_object: a function object
        """
        self.conflict_functions.update({conflict_function_name: conflict_function_object})

