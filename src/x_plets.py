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
        self.doublet_creation_time = None
        self.triplet_creation_time = None

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
        print(f"Time elapsed for creating doublets: "
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

        print(f"Time elapsed for creating triplets: "
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
            f.write(f"Number of particles hitting at least one detector layer: {self.num_particles}\n\n")
            f.write(f"Number of doublets found: {self.found_doublets}\n")
            f.write(f"Time needed to create doublet candidates: "
                    f"{XpletCreatorLUXE.hms_string(self.doublet_creation_time)}\n")
            f.write(f"Number of tracks approximately possible to reconstruct at doublet level: "
                    f"{int(self.found_correct_doublets / 3)}\n\n")

            f.write(f"Number of triplets found: {self.found_triplets}\n")
            f.write(f"Time needed to create triplet candidates: "
                    f"{XpletCreatorLUXE.hms_string(self.triplet_creation_time)}\n")
            f.write(f"Number of tracks approximately possible to reconstruct at triplet level: "
                    f"{int(self.found_correct_triplets / 2)}\n")
