import csv
import time
import numpy as np
import sys

from doublet import Doublet
from triplet import Triplet
from segment import Segment
from segment_manager import SegmentManager


class XpletCreatorLUXE:
    def __init__(self,
                 luxe_tracking_data_file_name,
                 detector_geometry,
                 segment_manager,
                 configuration,
                 save_to_folder):
        """Class for creating x_plets from detector hits
        :param luxe_tracking_data_file_name: .csv file name with tracking data from the LUXE positron detector:
            hit_ID          : unique hit ID
            x               : x value [m]
            y               : y value [m]
            z               : z value [m]
            layer_ID        : layer ID, starting from 0
            particle_ID     : number to identify particles
            particle_energy : particle energy [MeV]
        :param detector_geometry: simplified LUXE (sl) or full LUXE (fl)
        :param segment_manager : managing segments with hits
        :param configuration  : dictionary, configuration for detector setup and xplet selection
            {
            doublet: {dx/x0: <value>,
                       eps:   <value>},
            triplet: {angle diff x: <value>,
                      angle diff y: <value>},
            binning: {num bins x: <value>},
            qubo parameters: {b_ij conflict: <value>,
                              b_ij match: value,
                              a_i: <value>}
            scale range parameters: {z_scores: <value>,
                                     quality: <value>,
                                     interaction: <value>}
            }
        :param save_to_folder : folder in which results are stored
        """
        self.luxe_tracking_data_file_name = luxe_tracking_data_file_name
        self.detector_geometry = detector_geometry
        self.configuration = configuration
        self.save_to_folder = save_to_folder
        self.segment_manager = segment_manager
        self.mode = self.recognize_parameter_setup()

        # Simplified model only has 4 detector planes
        if self.detector_geometry == "sl":
            self.num_layers = 4
        elif self.detector_geometry == "fl":
            self.num_layers = 8
        else:
            print("No valid detector geometry chosen.")
            print("Exit program...")
            exit()

        # some values to check if computation successful
        self.num_particles = 0
        self.found_correct_doublets = 0
        self.found_correct_triplets = 0
        self.found_doublets = 0
        self.found_triplets = 0

        # indices from .csv file, code gets not messed up if some information might be added in the future to them
        self.x = None
        self.y = None
        self.z = None
        self.hit_id = None
        self.particle_id = None
        self.layer_id = None
        self.particle_energy = None

        # loading data and setting indices for .csv file
        self.num_segments = self.configuration["binning"]["num bins"]
        if self.detector_geometry == "sl":
            self.create_segments_simplified_model()
            self.segment_manager.map_segments_simple_model()
        elif self.detector_geometry == "fl":
            pass  # to be implemented --> FullLUXE

        self.load_tracking_data()
        print(f"Number of particles found: {self.num_particles}")

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
        self.add_quality_function("two norm angle standard deviation", XpletCreatorLUXE.two_norm_std_angle)
        self.add_conflict_function("two norm angle standard deviation", XpletCreatorLUXE.two_norm_std_angle)

    def recognize_parameter_setup(self):
        """Recognizes the chosen parameter setup for the QUBO
        :return: tuple of strings with setup information"""
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

    def load_tracking_data(self):
        """Loads data from a .csv file and stores it into a 2-dim array. Also converts values which are used for
        calculations to float from string. The file structure is displayed in the class description.
        """
        particle_numbers = set()  # counting particle ID
        with open(self.luxe_tracking_data_file_name, 'r') as file:
            csv_reader = csv.reader(file)
            csv_header = next(csv_reader)  # access header, csv files should consist of one line of header
            self.x = csv_header.index("x")
            self.y = csv_header.index("y")
            self.z = csv_header.index("z")
            self.hit_id = csv_header.index("hit_ID")
            self.particle_id = csv_header.index("particle_ID")
            self.layer_id = csv_header.index("layer_ID")
            self.particle_energy = csv_header.index("particle_energy")
            for row in csv_reader:
                row_converted = []
                for i in range(len(row)):
                    if i in [self.x, self.y, self.z]:
                        row_converted.append(float(row[i]))  # convert coordinates from string to float
                    else:
                        row_converted.append(row[i])
                    self.segment_manager.segment_list[
                        self.segment_manager.get_segment_at_known_xz_value(row[self.x],
                                                                         row[self.z])].data.append(row_converted)
                particle_numbers.add(row[self.particle_id])
            self.num_particles = len(particle_numbers)

    def make_x_plet_list_simplified_setup(self):
        """Creates doublets and triplets. For the simplified model only."""

        print("\nCreating doublet lists...\n")
        doublet_list_start = time.time()  # doublet list timer
        for segment in self.segment_manager.segment_list:
            if segment.layer > 2:  # no doublets start from last layer
                continue
            next_segments = self.segment_manager.target_segments(segment.name)   # target segments
            for first_hit in segment.data:
                for target_segment in next_segments:
                    for second_hit in target_segment.data:
                        if self.doublet_criteria_check(first_hit[self.x],
                                                       second_hit[self.x],
                                                       first_hit[self.z],
                                                       second_hit[self.z]):
                            doublet = Doublet(first_hit[self.particle_id],
                                              second_hit[self.particle_id],
                                              (first_hit[self.x],
                                               first_hit[self.y],
                                               first_hit[self.z]),
                                              (second_hit[self.x],
                                               second_hit[self.y],
                                               second_hit[self.z]),
                                              first_hit[self.hit_id],
                                              second_hit[self.hit_id])
                            if doublet.is_correct_match:
                                self.found_correct_doublets += 1
                            self.found_doublets += 1
                            segment.doublet_data.append(doublet)

        doublet_list_end = time.time()  # doublet list timer
        print(f"Time elapsed for creating list of doublets: "
              f"{XpletCreatorLUXE.hms_string(doublet_list_end - doublet_list_start)}")
        print(f"Number of tracks approximately possible to reconstruct at doublet level: "
              f"{int(self.found_correct_doublets / 3)}")
        print(f"Number of doublets found: {self.found_doublets}\n")

        list_triplet_start = time.time()
        print("\nCreating triplet lists...\n")
        for segment in self.segment_manager.segment_list:
            if segment.layer > 1:   # triplets start only at first two layers
                continue
            next_segments = self.segment_manager.target_segments(segment.name)  # target segments
            for target_segment in next_segments:
                for first_doublet in segment.doublet_data:
                    for second_doublet in target_segment.doublet_data:
                        if first_doublet.hit_2_position != second_doublet.hit_1_position:  # check if match
                            continue
                        if self.triplet_criteria_check(first_doublet, second_doublet):
                            triplet = Triplet(first_doublet, second_doublet, self.found_triplets)
                            self.found_triplets += 1
                            segment.triplet_data.append(triplet)
                            if triplet.is_correct_match:
                                self.found_correct_triplets += 1
            segment.doublet_data.clear()   # --> lower memory usage, num doublets are >> num triplets

        list_triplet_end = time.time()

        print(f"Time elapsed for creating list of triplets: "
              f"{XpletCreatorLUXE.hms_string(list_triplet_end - list_triplet_start)}")
        print(f"Number of tracks approximately possible to reconstruct at triplet level: "
              f"{int(self.found_correct_triplets / 2)}")
        print(f"Number of triplets found: {self.found_triplets}")

    def triplet_criteria_check(self, doublet1, doublet2):
        """
        Checks if doublets may be combined to a triplet, depending on the doublet angles -> scattering
        :param doublet1: first doublet, nearer to IP
        :param doublet2: second doublet, further from IP
        :return: True if criteria applies, else False
        """
        if abs(doublet2.xz_angle() - doublet1.xz_angle()) < self.configuration["triplet"]["angle diff x"]:
            if abs(doublet2.yz_angle() - doublet1.yz_angle()) < self.configuration["triplet"]["angle diff y"]:
                return True
        return False

    def doublet_criteria_check(self, x1, x2, z1, z2):
        """
        Checks if hits may combined to doublets, applying dx/x0 criterion
        :param x1: x value first hit
        :param x2: x value second hit
        :param z1: z value first hit
        :param z2: z value second hit
        :return: True if criteria applies, else False
        """
        if abs(((x2 - x1) / self.z_at_x0(x1, x2, z1, z2) - self.configuration["doublet"]["dx/x0"])) > \
                self.configuration["doublet"]["eps"]:
            return False
        return True

    def z_at_x0(self, x_end, x_start, z_end, z_start):
        """
        Help function for calculation x position of doublet at a z-reference value, usually the first detector layer
        counted from the IP
        :param x_end: x-position of target segment
        :param x_start: x-position of segment
        :param z_end: z-position of target segment
        :param z_start: z-position of segment
        :return: x_position at the reference layer
        """
        dx = x_end - x_start
        dz = z_end - z_start
        return x_end - dx * abs(z_end - self.segment_manager.reference_layer_z) / dz

    def add_quality_function(self,
                             quality_function_name,
                             quality_function_object):
        """Adds a function to the quality function dictionary
        :param quality_function_name: name of the function
        :param quality_function_object:  function object
        """
        self.quality_functions.update({quality_function_name: quality_function_object})

    def add_conflict_function(self,
                              conflict_function_name,
                              conflict_function_object):
        """Adds a function to the conflict function dictionary
        :param conflict_function_name: name of the function
        :param conflict_function_object:  function object
        """
        self.conflict_functions.update({conflict_function_name: conflict_function_object})

    def set_triplet_coefficients(self):
        """Sets the triplet coefficients according to the configuration files. If a (re-)normalization was
        set it is also applied. If successful a message and the target folder location containing the triplet
        list is displayed.
        """
        quality = None
        try:
            quality = float(self.configuration["qubo parameters"]["a_i"])
            quality_mode = "constant"
        except ValueError:
            quality_mode = "angle function"

        for segment in self.segment_manager.segment_list:
            if segment.layer > len(self.segment_manager.detector_layers) - 3:
                continue
            next_segments = self.segment_manager.target_segments(segment.name)  # target segments
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
        for segment in self.segment_manager.segment_list:
            for triplet in segment.triplet_data:
                self.triplet_list.add(triplet)
            segment.triplet_data.clear()
        self.triplet_list = list(self.triplet_list)
        self.triplet_list.sort(key=lambda t: t.triplet_id)   # triplet id = index in triplet list

        for t1 in self.triplet_list:
            # keeping track of statistics
            if t1.is_correct_match:
                self.quality_correct_match_list.append(t1.quality)
            else:
                self.quality_wrong_match_list.append(t1.quality)

            for i_key, i_value in t1.interactions.items():
                if i_key < t1.triplet_id:
                    t2 = self.triplet_list[i_key]
                    if t1.doublet_1.hit_1_particle_key == t1.doublet_1.hit_2_particle_key == \
                       t1.doublet_2.hit_1_particle_key == t1.doublet_2.hit_2_particle_key == \
                       t2.doublet_1.hit_1_particle_key == t2.doublet_1.hit_2_particle_key == \
                       t2.doublet_2.hit_1_particle_key == t2.doublet_2.hit_2_particle_key:
                        self.connectivity_correct_match_list.append(i_value)
                    else:
                        self.connectivity_wrong_match_list.append(i_value)

        # additional processing of qubo parameters
        if self.configuration["scale range parameters"]["z_scores"] is not None:
            quality_values = self.quality_correct_match_list + self.quality_wrong_match_list
            mu = np.mean(quality_values)
            sigma = np.std(quality_values)
            max_z_value = max([abs((value - mu) / sigma) for value in quality_values])
            for triplet in self.triplet_list:
                triplet.quality = (triplet.quality - mu) / sigma / max_z_value

        # Scaling in the following way to [a, b] : X' = a + (X - X_min) (b - a) / (X_max - X_min)
        if self.configuration["scale range parameters"]["quality"] is not None:
            a = self.configuration["settings"]["quality scale range"][0]
            b = self.configuration["settings"]["quality scale range"][1]
            quality_values = self.quality_correct_match_list + self.quality_wrong_match_list
            min_quality = min(quality_values)
            max_quality = max(quality_values)
            for triplet in self.triplet_list:
                triplet.quality = a + (triplet.quality - min_quality) * (b - a) / (max_quality - min_quality)

            # rewriting a_i lists
            for i in range(len(self.quality_correct_match_list)):
                self.quality_correct_match_list[i] = \
                    a + (self.quality_correct_match_list[i] - min_quality) * (b - a) / (max_quality - min_quality)
            for i in range(len(self.quality_wrong_match_list)):
                self.quality_wrong_match_list[i] = \
                    a + (self.quality_wrong_match_list[i] - min_quality) * (b - a) / (max_quality - min_quality)

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

        print("\nCoefficients set successfully")
        print(f"\nSaving triplet list to folder: {self.save_to_folder}")
        np.save(f"{self.save_to_folder}/triplet_list", self.triplet_list)

    def triplet_interaction(self,
                            triplet,
                            other_triplet):
        """Compares two triplets and  how they match.
        :param triplet: first triplet
        :param other_triplet: triplet to compare with the first one
        :return value based on connectivity/conflict and chosen set of parameters
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
        if t1[0][2] == t2[0][2]:
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
        elif self.configuration["settings"]["qubo mode"] == "f_con f_ct":
            if intersection == 1:
                return 1 + self.conflict_functions[self.configuration["qubo parameters"]["b_ij conflict"]](
                    triplet.doublet_1, triplet.doublet_2, other_triplet.doublet_2)
            elif intersection == 2 and not same_layer:
                return - 1 + self.quality_functions[self.configuration["qubo parameters"]["b_ij match"]](
                    triplet.doublet_1, triplet.doublet_2, other_triplet.doublet_2)
            else:
                return 1 + self.conflict_functions[self.configuration["qubo parameters"]["b_ij conflict"]](
                    triplet.doublet_1, triplet.doublet_2, other_triplet.doublet_2)

    @staticmethod
    def hms_string(sec_elapsed):
        """Nicely formatted time string."""
        h = int(sec_elapsed / (60 * 60))
        m = int((sec_elapsed % (60 * 60)) / 60)
        s = sec_elapsed % 60
        return "{}:{:>02}:{:>05.2f}".format(h, m, s)

    @staticmethod
    def two_norm_std_angle(doublet_1, doublet_2, doublet_3):
        """Returns 2-norm of angle difference in xz and yz.
        :param
            doublet_1 : doublet from hit 1 + 2
            doublet_2 : doublet from hit 2 + 3
            doublet_3 : doublet from hit 3 + 4
        :return
            2-norm of angle difference in xz and yz."""
        angle_xz_doublet_1 = doublet_1.xz_angle()
        angle_yz_doublet_1 = doublet_1.yz_angle()

        angle_xz_doublet_2 = doublet_2.xz_angle()
        angle_yz_doublet_2 = doublet_2.yz_angle()

        angle_xz_doublet_3 = doublet_3.xz_angle()
        angle_yz_doublet_3 = doublet_3.yz_angle()

        return np.sqrt(np.std([angle_xz_doublet_1, angle_xz_doublet_2, angle_xz_doublet_3]) ** 2 +
                       np.std([angle_yz_doublet_1, angle_yz_doublet_2, angle_yz_doublet_3]) ** 2)

    def write_info_file(self):
        """        Writes information about the Preselection parameters and some statistics into
        'Preselection.txt' which is stored inside the output folder.
        """
        with open(self.save_to_folder + "/preselection_info.txt", "w") as f:
            f.write("Preselection performed with the following set of parameters: \n")
            f.write("---\n")
            for outer_key in self.configuration.keys():
                for inner_key, value in self.configuration[outer_key]:
                    f.write(f"\n{innter_key}: {value}")
                f.write("\n")

            f.write("\n\n---\n")
            f.write("Statistics:\n\n")
            f.write(f"Number of particles hitting at leas one detector layer: {self.num_particles}\n\n")
            f.write(f"Number of doublets found: {self.found_doublets}\n")
            f.write(f"Number of tracks approximately possible to reconstruct at doublet level: "
                    f"{int(self.found_correct_doublets / 3)}\n\n")
            f.write(f"Number of triplets found: {self.found_triplets}\n")
            f.write(f"Number of tracks approximately possible to reconstruct at triplet level: "
                    f"{int(self.found_correct_triplets / 2)}\n")

    def plot_and_save_statistics(self):
        """This functions plots and saves various statistics to the same folder where the triplet list is saved to.
        The results are saved as 'triplet_coefficients_statistics.pdf' and 'triplet_interactions.pdf' into the
        target folder.
        """
        # Number of interactions with other triplets
        interactions_list = []
        for t in self.triplet_list:
            interactions_list.append(len(list(t.interactions.keys())))
        n, bins, _ = plt.hist(interactions_list, bins=max(interactions_list))
        plt.figure(figsize=(12, 9))
        plt.hist(np.array(interactions_list),
                 weights=1 / sum(n) * np.ones(len(interactions_list)),
                 bins=max(interactions_list),
                 edgecolor="firebrick",
                 linewidth=3,
                 histtype='step',
                 label=f"Number of particles: {self.num_particles}\n"
                       f"Number of triplets: {len(self.triplet_list)}")
        plt.yscale("log")
        plt.legend(loc="best", fontsize=20)

        plt.xlabel("Number of interactions with other triplets", fontsize=20)
        plt.ylabel("Fraction of counts", fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.savefig(f"{self.save_to_folder}/triplet_interactions.pdf")

        # distribution of coefficients
        fig, (ax1, ax2) = plt.subplots(2, figsize=(12, 12))

        ax1.hist([self.quality_correct_match_list, self.quality_wrong_match_list],
                 bins=50,
                 label=[r"quality correct match", r"quality wrong match"],
                 edgecolor='k',
                 color=["goldenrod", "royalblue"],
                 align="left")
        ax1.set_yscale('log')
        ax1.legend(loc='best', fontsize=20)
        ax1.tick_params(axis='both', labelsize=20)
        ax1.set_ylabel("counts", fontsize=20)
        ax1.set_xlabel("[a.u]", fontsize=20)

        n2, bins2, patches2 = ax2.hist([self.connectivity_correct_match_list, self.connectivity_wrong_match_list],
                                       bins=50,
                                       label=[r"connectivity correct match",
                                              r"connectivity wrong match"],
                                       edgecolor='k',
                                       color=["goldenrod", "royalblue"],
                                       align="left",
                                       rwidth=1)

        width = patches2[1][-1].get_width()
        height = patches2[1][-1].get_height()

        patches2[1][-1].remove()
        ax2.bar(1, align="center", height=height, width=width, color='red', edgecolor="k", label="conflicts")

        ax2.set_yscale('log')
        ax2.legend(loc='best', fontsize=20)
        ax2.set_ylabel("counts", fontsize=20)
        ax2.tick_params(axis='both', labelsize=20)
        ax2.set_xlabel("[a.u]", fontsize=20)
        plt.savefig(f"{self.save_to_folder}/triplet_coefficients_statistics.pdf")