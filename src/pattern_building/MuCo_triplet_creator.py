import csv
import time

import numpy as np

from utility.time_tracking import hms_string
from pattern.doublet import Doublet
from pattern.triplet import Triplet
import matplotlib.pyplot as plt
from random import randint


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

        self.doublets_dict = {"L-01": [],
                              "L-12": [],
                              "L-23": [],
                              "L-34": [],
                              "L-45": [],
                              "L-56": [],
                              "L-67": []}
        self.triplets_dict = {"L-012": [],
                              "L-123": [],
                              "L-234": [],
                              "L-345": [],
                              "L-456": [],
                              "L-567": []}
        self.muon_hits = 0
        
        # Dictionary for storying hits with respect to their radial distance to the IP
        # L-0: 30.0 < d < 31.4
        # L-1: 32.1 < d < 33.4
        # L-2: 51.0 < d < 52.9
        # L-3: 53.1 < d < 54.9
        # L-4: 74.0 < d < 75.5
        # L-5: 76.1 < d < 77.5
        # L-6: 102.0 < d < 103.2
        # L-7: 104.1 < d < 105.2

        self.vxd_hits_dict = {"L-0": [],
                              "L-1": [],
                              "L-2": [],
                              "L-3": [],
                              "L-4": [],
                              "L-5": [],
                              "L-6": [],
                              "L-7": []}

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
                hit_converted = ([float(r) if row.index(r) not in [self.time_index, self.hit_id_index]
                                  else r for r in row])
                if hit_converted[self.pdg_index] == 13:
                    self.muon_hits += 1
                d = self.radial_distance_to_IP(hit_converted)
                if 30.0 < d < 31.4:
                    self.vxd_hits_dict['L-0'].append(hit_converted)
                elif 32.1 < d < 33.4:
                    self.vxd_hits_dict['L-1'].append(hit_converted)
                elif 51.0 < d < 53.0:
                    self.vxd_hits_dict['L-2'].append(hit_converted)
                elif 53.1 < d < 54.9:
                    self.vxd_hits_dict['L-3'].append(hit_converted)
                elif 74.0 < d < 75.5:
                    self.vxd_hits_dict['L-4'].append(hit_converted)
                elif 76.1 < d < 77.5:
                    self.vxd_hits_dict['L-5'].append(hit_converted)
                elif 102.0 < d < 103.2:
                    self.vxd_hits_dict['L-6'].append(hit_converted)
                elif 104.1 < d < 105.2:
                    self.vxd_hits_dict['L-7'].append(hit_converted)
                else:
                    print(f"Hit not in VXD Tracker! d = {d}")
                
            print(f"{self.muon_hits} muon hits found")
            print(f"Number of expected signal doublets: {self.muon_hits / 8 * 7}")
            print(f"Number of expected signal triplets: {self.muon_hits / 8 * 6}")

    def create_doublet(self,
                       first_hit: list[float],
                       second_hit: list[float]) -> type(Doublet):
        """Creates a doublet from two hits.
        :param first_hit: first hit of doublet
        :param second_hit: second hit of doublet
        :return
            Doublet object
        """
        first_hit_particle_key = first_hit[self.pdg_index]
        if first_hit_particle_key != 13:
            first_hit_particle_key = randint(14, 1e16)
        second_hit_particle_key = second_hit[self.pdg_index]
        if second_hit_particle_key != 13:
            second_hit_particle_key  = randint(14, 1e16)
            
        doublet = Doublet(hit_1_particle_key=first_hit_particle_key,
                          hit_2_particle_key=second_hit_particle_key,
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

    def radial_distance_to_IP(self,
                              hit):
        """"""
        return np.sqrt(hit[self.x_index]**2 + hit[self.y_index]**2)

    def create_xplet_list(self) -> None:
        """Creates the doublets. Only distance to IP check, pure combinatorial.
        """
        print("-----------------------------------\n")
        print("Forming doublets ...\n")
        doublet_list_start = time.process_time()
        
        rejected_doublet_because_of_time = 0
        found_correct_doublets = 0
        
        doublets_plotting = []

        for i in range(len(self.vxd_hits_dict.keys()) - 1):
            for hit_1 in self.vxd_hits_dict[f'L-{i}']:
                for hit_2 in self.vxd_hits_dict[f'L-{i + 1}']:
                    muon_match=False
                    if abs(np.arctan2(hit_1[self.y_index], hit_1[self.x_index]) - \
                           np.arctan2(hit_2[self.y_index], hit_2[self.x_index])) * 180 / np.pi > 5:
                        continue
                    if abs(np.arctan2(hit_1[self.y_index], hit_1[self.z_index]) - \
                           np.arctan2(hit_2[self.y_index], hit_2[self.z_index])) * 180 / np.pi > 5:
                        continue
                    if hit_1[self.pdg_index] == hit_2[self.pdg_index] == 13: 
                        muon_match = True
                    if abs(float(hit_2[self.time_index]) - float(hit_1[self.time_index])) > 0.09:
                        rejected_doublet_because_of_time += 1
                        continue
                    if muon_match:
                        found_correct_doublets += 1

                    self.doublets_dict[f'L-{i}{i + 1}'].append(self.create_doublet(hit_1,
                                                                                   hit_2))
                    doublets_plotting.append(([hit_1[self.x_index], hit_2[self.x_index]], 
                                              [hit_1[self.y_index], hit_2[self.y_index]])) 
        plt.figure(figsize=(12,12))
        for d in doublets_plotting:
            plt.plot(d[0], d[1], marker="o", c="blue", mfc='red', markersize=12)
                    
        if self.muon_hits > 0:
            print(f"Found correct doublets: {found_correct_doublets} --> "
                  f"{100 * np.around(found_correct_doublets / (self.muon_hits / 8 * 7), 2)} %")
        else:
            print(f"Found correct doublets: {found_correct_doublets} --> {0.0} %")
        print(f"Rejected doublets because of time: {rejected_doublet_because_of_time}")
                        
        doublet_list_end = time.process_time()
        doublet_creation_time = hms_string(doublet_list_end - doublet_list_start)

        print(f"Time elapsed for forming doublets: "
              f"{doublet_creation_time}")
        print(f"Number of doublets found: {len([d for layer in self.doublets_dict.values() for d in layer])}")

        print("-----------------------------------\n")
        print("Forming triplets ...\n")
        list_triplet_start = time.process_time()
        
        found_correct_triplets = 0
        
        for i in range(len(self.vxd_hits_dict.keys()) - 2):
            for d1 in (self.doublets_dict[f'L-{i}{i + 1}']):
                for d2 in (self.doublets_dict[f'L-{i + 1}{i + 2}']):
                    if d1.hit_2_id != d2.hit_1_id:
                        continue
                    if d1.hit_1_particle_key == d1.hit_2_particle_key == d2.hit_1_particle_key == d2.hit_2_particle_key == 13:
                        found_correct_triplets += 1
                    new_triplet = Triplet(d1, d2)
                    new_triplet.quality = 0.1
                    self.triplets_dict[f'L-{i}{i + 1}{i + 2}'].append(new_triplet)
        
        list_triplet_end = time.process_time()
        triplet_creation_time = hms_string(list_triplet_end - list_triplet_start)
        if self.muon_hits > 0:          
            print(f"Found correct triplets: {found_correct_triplets} --> "
                  f"{100 * np.around(found_correct_triplets / (self.muon_hits / 8 * 6), 2)} %")
        else:
             print(f"Found correct triplets: {found_correct_triplets} --> {0.0} %")     
        print(f"Time elapsed for  forming triplets: "
              f"{triplet_creation_time}")
        print(f"Number of triplets found: {len([t for layer in self.triplets_dict.values() for t in layer])}\n")
        

    def set_qubo_coefficients(self,
                              target_folder: str) -> None:

        connection = -1
        conflict = 2
        for i in range(len(self.vxd_hits_dict.keys()) - 3):
            print(f'Setting coeffciecints for Layer L-{i}{i + 1}{i + 2} --> L-{i + 1}{i + 2}{i + 3}')
            for t1 in self.triplets_dict[f'L-{i}{i + 1}{i + 2}']:
                for t2 in self.triplets_dict[f'L-{i}{i + 1}{i + 2}'] + self.triplets_dict[f'L-{i + 1}{i + 2}{i + 3}']:

                    t1_set = {t1.doublet_1.hit_1_id,
                              t1.doublet_1.hit_2_id,
                              t1.doublet_2.hit_2_id}
                    t2_set = {t2.doublet_1.hit_1_id,
                              t2.doublet_1.hit_2_id,
                              t2.doublet_2.hit_2_id}
                    intersection = len(t1_set.intersection(t2_set))

                    if intersection in [0, 3]:
                        continue

                    elif intersection == 1:
                        t1.interactions.update({t2.triplet_id: conflict})
                        t2.interactions.update({t1.triplet_id: conflict})

                    elif intersection == 2:
                        t1.interactions.update({t2.triplet_id: connection})
                        t2.interactions.update({t1.triplet_id: connection})

        print(f"Save triplet list...\n")
        np.save(f"{target_folder}/triplet_list", [t for layer in self.triplets_dict.values() for t in layer])


