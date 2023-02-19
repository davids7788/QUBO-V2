import numpy as np
import matplotlib.pyplot as plt

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
        self.add_quality_function("two norm angle standard deviation", QuboCoefficients.two_norm_std_angle)
        self.add_conflict_function("two norm angle standard deviation", QuboCoefficients.two_norm_std_angle)

    def set_triplet_coefficients(self,
                                 segment_manager: SegmentManager) -> None:
        """Sets the triplet coefficients according to the configuration files. If a (re-)normalization was
        set it is also applied. If the process is successful a message and the target folder location containing the
        triplet list is displayed.
        :param segment_manager: segment manager object
        """
        # self checking configuration
        print("Setting triplet coefficients a_i and b_ij...")
        quality = None
        try:
            quality = float(self.configuration["qubo parameters"]["a_i"])
            quality_mode = "constant"
        except ValueError:
            quality_mode = "angle function"

        for segment in segment_manager.segment_list:
            if segment.layer > len(segment_manager.detector_layers) - 3:
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
        for segment in segment_manager.segment_list:
            for triplet in segment.triplet_data:
                self.triplet_list.add(triplet)
            segment.triplet_data.clear()

    def filling_lists_for_statistics(self):
        """Function for collecting and storing information about quality and interaction values,
        as well as truth information about the triplets.
        """
        self.triplet_list = list(self.triplet_list)
        t_mapping = {key: value for key, value in zip([t.triplet_id for t in self.triplet_list],
                                                      np.arange(len(self.triplet_list)))}

        print("Collecting statistics about a_i and b_ij...")
        for i, t1 in enumerate(self.triplet_list):
            if t1.is_correct_match():
                self.quality_correct_match_list.append(t1.quality)
            else:
                self.quality_wrong_match_list.append(t1.quality)

            for i_key, i_value in t1.interactions.items():
                if t_mapping[i_key] < i:
                    t2 = self.triplet_list[t_mapping[i_key]]
                    if t1.doublet_1.hit_1_particle_key == t1.doublet_1.hit_2_particle_key == \
                       t1.doublet_2.hit_1_particle_key == t1.doublet_2.hit_2_particle_key == \
                       t2.doublet_1.hit_1_particle_key == t2.doublet_1.hit_2_particle_key == \
                       t2.doublet_2.hit_1_particle_key == t2.doublet_2.hit_2_particle_key:
                        self.connectivity_correct_match_list.append(i_value)
                    else:
                        self.connectivity_wrong_match_list.append(i_value)

    def parameter_rescaling(self):
        """Rescaling parameters according to the config file.
        """
        print("Starting rescaling of a_i and b_ij parameters...")
        # additional processing of qubo parameters
        if self.configuration["scale range parameters"]["z_scores"]:
            quality_values = self.quality_correct_match_list + self.quality_wrong_match_list
            mu = np.mean(quality_values)
            sigma = np.std(quality_values)
            for triplet in self.triplet_list:
                triplet.quality = (triplet.quality - mu) / sigma
            for i, quality in enumerate(self.quality_correct_match_list):
                self.quality_correct_match_list[i] = (quality - mu) / sigma
            for i, quality in enumerate(self.quality_wrong_match_list):
                self.quality_wrong_match_list[i] = (quality - mu) / sigma

                # Scaling in the following way to [a, b] : X' = a + (X - X_min) (b - a) / (X_max - X_min)
        if self.configuration["scale range parameters"]["quality"] is not None:
            a = self.configuration["scale range parameters"]["quality"][0]
            b = self.configuration["scale range parameters"]["quality"][1]
            quality_values = self.quality_correct_match_list + self.quality_wrong_match_list
            min_quality = min(quality_values)
            max_quality = max(quality_values)
            for triplet in self.triplet_list:
                triplet.quality = a + (triplet.quality - min_quality) * (b - a) / (max_quality - min_quality)

            # rewriting a_i lists
            for i, quality in enumerate(self.quality_correct_match_list):
                self.quality_correct_match_list[i] = a + (quality - min_quality) * (b - a) / (max_quality - min_quality)
            for i, quality in enumerate(self.quality_wrong_match_list):
                self.quality_wrong_match_list[i] = a + (quality - min_quality) * (b - a) / (max_quality - min_quality)

        # scaling connectivity
        if self.configuration["scale range parameters"]["interaction"] is not None:
            # excluding conflict terms
            connectivity_values = [value for value in self.connectivity_correct_match_list if value !=
                                   float(self.configuration["qubo parameters"]["b_ij conflict"])] + \
                                  [value for value in self.connectivity_wrong_match_list if value !=
                                   float(self.configuration["qubo parameters"]["b_ij conflict"])]
            min_connectivity = min(connectivity_values)
            max_connectivity = max(connectivity_values)

            a = self.configuration["scale range parameters"]["interaction"][0]
            b = self.configuration["scale range parameters"]["interaction"][1]
            for triplet in self.triplet_list:
                for key in triplet.interactions.keys():
                    if triplet.interactions[key] == self.configuration["qubo parameters"]["b_ij conflict"]:
                        continue
                    else:
                        triplet.interactions[key] = a + (triplet.interactions[key] - min_connectivity) * (b - a) / \
                                                    (max_connectivity - min_connectivity)

            # rewriting connectivity lists
            for i in range(len(self.connectivity_wrong_match_list)):
                if self.connectivity_wrong_match_list[i] == float(self.configuration["qubo parameters"]
                                                                  ["b_ij conflict"]):
                    continue
                else:
                    self.connectivity_wrong_match_list[i] = a + (self.connectivity_wrong_match_list[i] -
                                                                 min_connectivity) * (b - a) / \
                                                                (max_connectivity - min_connectivity)
            for i in range(len(self.connectivity_correct_match_list)):
                if self.connectivity_correct_match_list[i] == float(self.configuration["qubo parameters"]
                                                                    ["b_ij conflict"]):
                    continue
                else:
                    self.connectivity_correct_match_list[i] = a + (self.connectivity_correct_match_list[i] -
                                                                   min_connectivity) * (b - a) / \
                                                                  (max_connectivity - min_connectivity)

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
        intersection = 0

        t1 = [triplet.doublet_1.hit_1_position,
              triplet.doublet_1.hit_2_position,
              triplet.doublet_2.hit_2_position]
        t2 = [other_triplet.doublet_1.hit_1_position,
              other_triplet.doublet_1.hit_2_position,
              other_triplet.doublet_2.hit_2_position]
        for value_t1 in t1:
            for value_t2 in t2:
                if value_t1 == value_t2:
                    intersection += 1

        # check if triplet origin from same layer
        if t1[0][2] == t2[0][2] or t1[2][2] == t2[0][2]:
            same_layer = True
        else:
            same_layer = False

        # same and not interacting triplets get a zero as a coefficient
        if intersection == 0 or intersection == 3:
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

