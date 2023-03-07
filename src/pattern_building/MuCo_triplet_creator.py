import csv
import time

import numpy as np

from utility.time_tracking import hms_string
from pattern.doublet import Doublet
from pattern.triplet import Triplet


class MuCoTripletCreator:
    def __init__(self):
        """Class for creating triplets from Muon Collider detector hits.
        """
        self.doublet_creation_time = None
        self.triplet_creation_time = None
        self.num_signal_events = 0
        self.num_all_doublets = 0
        self.num_all_triplets = 0

        self.hit_id_index = None
        self.x_index = None
        self.y_index = None
        self.z_index = None
        self.mc_particle_id_index = None
        self.time_index = None
        self.pdg_index = None

        self.data = []
        self.doublets = []
        self.triplets = []

    def load_tracking_data(self,
                           tracking_data_file: str) -> None:
        """Loads data from a .csv file and stores it into a 2-dim array. Also converts values which are used for
        calculations to float from string. The file structure is displayed in the class description.
        :param tracking_data_file: Luxe tracking data file
        """
        print(f"Using tracking data file {tracking_data_file.split('/')[-1]}\n")

        with open(tracking_data_file, 'r') as file:
            csv_reader = csv.reader(file)
            csv_header = next(csv_reader)  # access header, csv files should consist of one line of header
            self.x_index = csv_header.index("x")
            self.y_index = csv_header.index("y")
            self.z_index = csv_header.index("z")
            self.hit_id_index = csv_header.index("hit_ID")
            self.mc_particle_id_index = csv_header.index("MC_particle_ID")
            self.pdg_index = csv_header.index("PDG")
            self.time_index = csv_header.index("time")

            for row in csv_reader:
                self.data.append([float(r) if row.index(r) not in [self.time_index, self.hit_id_index]
                                  else r for r in row])

    def create_doublet(self,
                       first_hit: list[float],
                       second_hit: list[float]) -> type(Doublet):
        """Creates a doublet from two hits.
        :param first_hit: first hit of doublet
        :param second_hit: second hit of doublet
        :return
            Doublet object
        """
        doublet = Doublet(hit_1_particle_key=first_hit[self.mc_particle_id_index],
                          hit_2_particle_key=second_hit[self.mc_particle_id_index],
                          hit_1_position=(first_hit[self.x_index],
                                          first_hit[self.y_index],
                                          first_hit[self.z_index]),
                          hit_2_position=(second_hit[self.x_index],
                                          second_hit[self.y_index],
                                          second_hit[self.z_index]),
                          hit_1_id=first_hit[self.hit_id_index],
                          hit_2_id=second_hit[self.hit_id_index],
                          time_1=first_hit[self.time_index],
                          time_2=second_hit[self.time_index])

        return doublet

    def distance_to_IP(self,
                       hit):
        """"""
        return np.sqrt(hit[self.x_index]**2 + hit[self.y_index]**2 + hit[self.z_index]**2)

    def create_xplet_list(self) -> None:
        """Creates the doublets. Only distance to IP check, pure combinatorial.
        """
        print("-----------------------------------\n")
        print("Forming doublets ...\n")
        doublet_list_start = time.process_time()

        for i, hit_1 in enumerate(self.data):
            for j, hit_2 in enumerate(self.data[i+1:]):
                if self.distance_to_IP(hit_1) > self.distance_to_IP(hit_2):
                    continue
                if abs(hit_1[self.x_index] - hit_2[self.x_index]) > 15:
                    continue
                if abs(hit_1[self.y_index] - hit_2[self.y_index]) > 15:
                    continue
                self.doublets.append(self.create_doublet(hit_1,
                                                         hit_2))

        doublet_list_end = time.process_time()
        doublet_creation_time = hms_string(doublet_list_end - doublet_list_start)

        print(f"Time elapsed for forming doublets: "
              f"{doublet_creation_time}")
        print(f"Number of doublets found: {len(self.doublets)}")

        print("-----------------------------------\n")
        print("Forming triplets ...\n")
        list_triplet_start = time.process_time()

        for i, d1 in enumerate(self.doublets):
            for j, d2 in enumerate(self.doublets):
                if d1.hit_2_id != d2.hit_1_id:
                    continue
                self.triplets.append(Triplet(d1, d2))

        list_triplet_end = time.process_time()
        triplet_creation_time = hms_string(list_triplet_end - list_triplet_start)
        print(f"Time elapsed for  forming triplets: "
              f"{triplet_creation_time}")
        print(f"Number of triplets found: {len(self.triplets)}")

    def set_qubo_coefficients(self,
                              target_folder: str) -> None:

        connection = -1
        conflict = 1

        for i, triplet in enumerate(self.triplets):
            print(i, end="\r")
            for other_triplet in self.triplets[i + 1:]:

                t1_set = {triplet.doublet_1.hit_1_id,
                          triplet.doublet_1.hit_2_id,
                          triplet.doublet_2.hit_2_id}
                t2_set = {other_triplet.doublet_1.hit_1_id,
                          other_triplet.doublet_1.hit_2_id,
                          other_triplet.doublet_2.hit_2_id}
                intersection = len(t1_set.intersection(t2_set))

                if intersection in [0, 3]:
                    continue

                elif intersection == 1:
                    triplet.interactions.update({other_triplet.triplet_id: conflict})
                    other_triplet.interactions.update({triplet.triplet_id: conflict})

                else:
                    # triplets on same layer always have a conflict
                    if triplet.doublet_1.hit_1_position[2] == other_triplet.doublet_1.hit_1_position[2]:
                        triplet.interactions.update({other_triplet.triplet_id: conflict})
                        other_triplet.interactions.update({triplet.triplet_id: conflict})
                    else:
                        triplet.interactions.update({other_triplet.triplet_id: connection})
                        other_triplet.interactions.update({triplet.triplet_id: connection})

        print(f"Save triplet list...")
        np.save(f"{target_folder}/triplet_list", self.triplets)





