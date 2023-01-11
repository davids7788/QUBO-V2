import csv
import time
import numpy as np

from pattern.doublet import Doublet
from pattern.triplet import Triplet
from preselection.segment_manager import SegmentManager


class TripletCreatorLUXE:
    def __init__(self,
                 configuration,
                 save_to_folder):
        """Class for creating x_plets from detector hits
        :param configuration  : dictionary, configuration for detector setup and xplet selection
            {
            doublet: {dx/x0: <value>,
                      eps:   <value>,
                      dy: <value>},
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
        self.configuration = configuration
        self.save_to_folder = save_to_folder
        self.doublet_creation_time = None
        self.triplet_creation_time = None
        self.num_particles = 0

        self.preselection_statistic_dx_x0 = []
        self.preselection_statistic_scattering = []

        # some values to check if computation successful
        self.num_complete_tracks = 0
        self.found_correct_doublets = 0
        self.found_correct_triplets = 0
        self.found_doublets = 0
        self.found_triplets = 0
        self.xplet_numbers = set()

        # indices from .csv file, code gets not messed up if some information might be added in the future to them
        self.x_index = None
        self.y_index = None
        self.z_index = None
        self.hit_id_index = None
        self.particle_id_index = None
        self.layer_id_index = None
        self.particle_energy_index = None

    def load_tracking_data(self,
                           tracking_data_file: str,
                           segment_manager: SegmentManager,
                           geometry_file: str):
        """Loads data from a .csv file and stores it into a 2-dim array. Also converts values which are used for
        calculations to float from string. The file structure is displayed in the class description.
        :param tracking_data_file: Luxe tracking data file
        :param segment_manager: SegmentManager object with already set segments and mapping
        :param geometry_file: geometry .csv file
        """
        with open(geometry_file, 'r') as file:
            csv_reader_geometry_file = csv.reader(file)
            csv_header_geometry_file = next(csv_reader_geometry_file)
            z = csv_header_geometry_file.index("z")
            z_values = []
            for row in csv_reader_geometry_file:
                z_values.append(float(row[z]))
            z_values_unique = list(set(z_values))
            z_values_unique.sort()
            if len(z_values_unique) == 4:
                # Simplified LUXE setup consists of 4 layers, assuming consisting one big chip each
                last_layer = [z_values_unique[-1], z_values_unique[-1]]
                first_layer = [z_values_unique[0], z_values_unique[0]]
            else:
                # Full LUXE setup consists of 4 layers, each layer has 2 rows of chips, the rows overlap in the middle
                last_layer = [z_values_unique[-2], z_values_unique[-1]]
                first_layer = [z_values_unique[0], z_values_unique[1]]

        particle_numbers = {}  # counting particle ID

        print(f"Loading tracking data: {tracking_data_file}\n"
              f"Distributing data into segments ...\n")

        with open(tracking_data_file, 'r') as file:
            csv_reader = csv.reader(file)
            csv_header = next(csv_reader)  # access header, csv files should consist of one line of header
            self.x_index = csv_header.index("x")
            self.y_index = csv_header.index("y")
            self.z_index = csv_header.index("z")
            self.hit_id_index = csv_header.index("hit_ID")
            self.particle_id_index = csv_header.index("particle_ID")
            self.layer_id_index = csv_header.index("layer_ID")
            self.particle_energy_index = csv_header.index("particle_energy")
            for row in csv_reader:
                row_converted = []
                for i in range(len(row)):
                    if i in [self.x_index, self.y_index, self.z_index, self.particle_energy_index]:
                        row_converted.append(float(row[i]))  # convert coordinates from string to float
                    else:
                        row_converted.append(row[i])
                # check if particle ID was already seen
                if row[self.particle_id_index] in particle_numbers:
                    particle_numbers[row[self.particle_id_index]].update(
                                        {z_values_unique.index(row_converted[self.z_index]):
                                         row_converted[self.z_index]})
                else:
                    particle_numbers.update({row[self.particle_id_index]:
                                            {z_values_unique.index(row_converted[self.z_index]):
                                             row_converted[self.z_index]}})

                # adding particle to a segment, used to reduce combinatorial candidates
                segment_index_for_entry = segment_manager.get_segment_at_known_xyz_value(row_converted[self.x_index],
                                                                                         row_converted[self.y_index],
                                                                                         row_converted[self.z_index])
                segment_manager.segment_list[segment_index_for_entry].data.append(row_converted)

            for key, outer_value in particle_numbers.items():
                last = False
                first = False
                for inner_value in outer_value.values():
                    if inner_value in last_layer:
                        last = True
                    if inner_value in first_layer:
                        first = True
                self.num_particles += 1
                if last and first:
                    self.num_complete_tracks += 1
                    self.xplet_numbers.add(key)

            print(f"Number of particles with at least one hit: {self.num_particles}")
            print(f"Number of complete tracks: {self.num_complete_tracks}\n")

    def create_x_plets_simplified_LUXE(self,
                                       segment_manager: SegmentManager):
        """Creates doublets and triplets. For the simplified model only.
        :param segment_manager: SegmentManager object with already set segments and mapping
        """
        print("-----------------------------------")
        print("Creating doublet lists...\n")
        doublet_list_start = time.time()  # doublet list timer
        for segment in segment_manager.segment_list:
            if segment.layer > len(segment_manager.detector_layers) - 2:  # no doublets start from last layer
                continue
            next_segments = segment_manager.target_segments(segment.name)   # target segments

            for first_hit in segment.data:
                for target_segment in next_segments:
                    for second_hit in target_segment.data:
                        x0 = self.x0_at_z_ref(first_hit[self.x_index],
                                              second_hit[self.x_index],
                                              first_hit[self.z_index],
                                              second_hit[self.z_index],
                                              segment_manager.get_z_reference_layer_LUXE())
                        if abs(second_hit[self.y_index] - first_hit[self.y_index]) / x0 > \
                                self.configuration["doublet"]["dy/x0"]:
                            continue
                        if self.doublet_criteria_check(first_hit[self.x_index],
                                                       second_hit[self.x_index],
                                                       first_hit[self.z_index],
                                                       second_hit[self.z_index],
                                                       segment_manager.get_z_reference_layer_LUXE()):
                            doublet = Doublet(first_hit[self.particle_id_index],
                                              second_hit[self.particle_id_index],
                                              (first_hit[self.x_index],
                                               first_hit[self.y_index],
                                               first_hit[self.z_index]),
                                              (second_hit[self.x_index],
                                               second_hit[self.y_index],
                                               second_hit[self.z_index]),
                                              first_hit[self.hit_id_index],
                                              second_hit[self.hit_id_index],
                                              first_hit[self.particle_energy_index],
                                              second_hit[self.particle_energy_index])
                            if doublet.is_correct_match() and doublet.hit_1_particle_key in self.xplet_numbers:
                                self.found_correct_doublets += 1
                                self.preselection_statistic_dx_x0.append((doublet.hit_2_position[0] -
                                                                          doublet.hit_1_position[0]) /
                                                                         self.x0_at_z_ref(doublet.hit_2_position[0],
                                                                                          doublet.hit_1_position[0],
                                                                                          doublet.hit_2_position[2],
                                                                                          doublet.hit_1_position[2],
                                                                                          segment_manager.
                                                                                          get_z_reference_layer_LUXE()))
                            self.found_doublets += 1
                            segment.doublet_data.append(doublet)
        doublet_list_end = time.time()  # doublet list timer
        self.doublet_creation_time = TripletCreatorLUXE.hms_string(doublet_list_end - doublet_list_start)
        print(f"Time elapsed for creating doublets: "
              f"{self.doublet_creation_time}")
        print(f"Number of doublets found: {self.found_doublets}\n")
        print(f"Number of tracks approximately possible to reconstruct with set doublet preparation parameters: "
              f"{int(self.found_correct_doublets / 3)}\n")
        print(f"Doublet selection efficiency d_matched / d_max: "
              f"{np.around(100 * self.found_correct_doublets / (3 * self.num_complete_tracks), 3)} %\n"
              f"d_max: maximum possible correctly matched doublets\n"
              f"d_matched: number of doublets stemming from only one particle")

        list_triplet_start = time.time()
        print("\nCreating triplet lists...\n")
        for segment in segment_manager.segment_list:
            if segment.layer > len(segment_manager.detector_layers) - 3:   # triplets start only at first two layers
                continue
            next_segments = segment_manager.target_segments(segment.name)  # target segments
            for target_segment in next_segments:
                for first_doublet in segment.doublet_data:
                    for second_doublet in target_segment.doublet_data:
                        if first_doublet.hit_2_position != second_doublet.hit_1_position:  # check if match
                            continue
                        if self.triplet_criteria_check(first_doublet, second_doublet):
                            triplet = Triplet(first_doublet, second_doublet)
                            self.found_triplets += 1
                            segment.triplet_data.append(triplet)

                            # filling lists for statistical purposes
                            if triplet.is_correct_match() and triplet.doublet_1.hit_1_particle_key \
                                    in self.xplet_numbers:
                                self.preselection_statistic_scattering.append(
                                    np.sqrt(triplet.angles_between_doublets()[0]**2 +
                                            triplet.angles_between_doublets()[1]**2))
                                self.found_correct_triplets += 1
            segment.doublet_data.clear()   # --> lower memory usage, num doublets are >> num triplets

        list_triplet_end = time.time()
        self.triplet_creation_time = TripletCreatorLUXE.hms_string(list_triplet_end - list_triplet_start)

        print(f"Time elapsed for creating triplets: "
              f"{self.triplet_creation_time}")
        print(f"Number of triplets found: {self.found_triplets}")
        print(f"Number of tracks approximately possible to reconstruct with set triplet preparation parameters: "
              f"{int(self.found_correct_triplets / 2)}\n")
        print(f"Triplet selection efficiency t_matched / t_max: "
              f"{np.around(100 * self.found_correct_triplets / (2 * self.num_complete_tracks), 3)} %\n"
              f"d_max: maximum possible correctly matched triplets\n"
              f"d_matched: number of triplets stemming from only one particle")

    def triplet_criteria_check(self,
                               doublet1: Doublet,
                               doublet2: Doublet):
        """Checks if doublets may be combined to a triplet, depending on the doublet angles -> scattering
        :param doublet1: first doublet, nearer to IP
        :param doublet2: second doublet, further from IP
        :return:
            True if criteria applies, else False
        """
        if np.sqrt((doublet2.xz_angle() - doublet1.xz_angle())**2 +
                   (doublet2.yz_angle() - doublet1.yz_angle())**2) < self.configuration["triplet"]["max scattering"]:
            return True
        return False

    def doublet_criteria_check(self,
                               x1: float,
                               x2: float,
                               z1: float,
                               z2: float,
                               z_ref: float):
        """Checks if hits may be combined to doublets, applying dx/x0 criterion
        :param x1: x value first hit
        :param x2: x value second hit
        :param z1: z value first hit
        :param z2: z value second hit
        :param z_ref: z-value reference
        :return:
            True if criteria applies, else False
        """
        if abs(((x2 - x1) / TripletCreatorLUXE.x0_at_z_ref(x1, x2, z1, z2, z_ref) -
                self.configuration["doublet"]["dx/x0"])) > \
                self.configuration["doublet"]["eps"]:
            return False
        return True

    @staticmethod
    def x0_at_z_ref(x_end: float,
                    x_start: float,
                    z_end: float,
                    z_start: float,
                    z_ref: float):
        """Function for calculation x position of a doublet at a z-reference value.
        :param x_end: x-position of second hit
        :param x_start: x-position of first hit
        :param z_end: z-position of second hit
        :param z_start: z-position of first hit
        :param z_ref: z-value reference layer
        :return:
            x_position at the reference layer
        """
        dx = x_end - x_start
        dz = z_end - z_start
        return x_end - dx * abs(z_end - z_ref) / dz

    @staticmethod
    def hms_string(sec_elapsed):
        """Nicely formatted time string.
        :param sec_elapsed time in ms
        :return
            hh:mm:ss.msms
        """
        h = int(sec_elapsed / (60 * 60))
        m = int((sec_elapsed % (60 * 60)) / 60)
        s = sec_elapsed % 60
        return "{}:{:>02}:{:>05.2f}".format(h, m, s)

    def write_info_file(self):
        """Writes information about the Preselection parameters and some statistics into
        'preselection_info.txt' which is stored inside the output folder.
        """
        with open(self.save_to_folder + "/preselection_info.txt", "w") as f:
            f.write("Preselection performed with the following set of parameters: \n")
            f.write("---\n")
            for outer_key in self.configuration.keys():
                for inner_key, value in self.configuration[outer_key].items():
                    f.write(f"\n{inner_key}: {value}")
                f.write("\n")
            f.write("\n\n")
            f.write("---\n")
            f.write(f"Number of particles with at least one hit on the detector:: {self.num_particles}\n")
            f.write(f"Number of generated tracks: {self.num_complete_tracks}\n")

            f.write(f"Time elapsed for creating doublets: "
                    f"{self.doublet_creation_time}\n")
            f.write(f"Number of doublets found: {self.found_doublets}\n")
            f.write(f"Number of tracks approximately possible to reconstruct with set doublet preparation parameters: "
                    f"{int(self.found_correct_doublets / 3)}\n\n")
            f.write(f"Doublet selection efficiency d_matched / d_max: "
                    f"{np.around(100 * self.found_correct_doublets / (3 * self.num_complete_tracks), 3)} %\n"
                    f"d_max: maximum possible correctly matched doublets, known from truth number of particles\n"
                    f"d_matched: number of doublets stemming from only one particle\n\n")

            f.write(f"Time elapsed for creating triplets: "
                    f"{self.triplet_creation_time}\n")
            f.write(f"Number of triplets found: {self.found_triplets}\n")
            f.write(f"Number of tracks approximately possible to reconstruct with set triplet preparation parameters: "
                    f"{int(self.found_correct_triplets / 2)}\n\n")
            f.write(f"Triplet selection efficiency t_matched / t_max: "
                    f"{np.around(100 * self.found_correct_triplets / (2 * self.num_complete_tracks), 3)} %\n"
                    f"t_max: maximum possible correctly matched triplets, known from truth number of particles\n"
                    f"t_matched: number of triplets stemming from only one particle")
