import csv
import time
import os

import numpy as np
import matplotlib.pyplot as plt

from math_functions.geometry import w_angle_diff
from utility.time_tracking import hms_string
from pattern.doublet import Doublet
from pattern.triplet import Triplet

from random import randint


class MuCoTripletCreator:
    def __init__(self):
        """Class for creating triplets from Muon Collider detector hits.
        """
        self.doublet_creation_time = 0
        self.triplet_creation_time = 0
        self.num_signal_events = 0
        self.num_all_doublets = 0
        self.num_all_triplets = 0
        self.muon_tracks = {}
        self.muon_hits = 0

        # fieldnames according to .csv naming
        self.fieldnames = ['hit_ID', 'x', 'y', 'z', 'layer', 'MC_particle_ID', 'px', 'py', 'pz', 'time', 'PDG', 'event']

    def load_tracking_data(self,
                           event_folder: str,
                           dlfilter: bool,
                           s_manager,
                           no_endcaps=True) -> None:
        """Loads data from the event stored in various .csv files inside the specified folder.
        Also converts values which are used for calculations to float from string.
        :param event_folder: folder containing a muon collider event
        :param dlfilter: if True, then the double layer filtered hits are used instead of the whole VXDTracker sample
        :param s_manager: segment manager
        :param no_endcaps: if True no data from endcaps are used
        """
        print(f'Processing folder: {event_folder}')
        event_files = [f for f in os.listdir(event_folder) if '.csv' in f]

        if not dlfilter:
            for f in event_files:
                if "_DLFiltered_" in f:
                    event_files.remove(f)
        else:
            for f in event_files:
                if "_VXD_" in f:
                    event_files.remove(f)

        if no_endcaps:
            event_files = [f for f in event_files if 'Endcap' not in f]

        arguments_to_convert = [self.fieldnames.index('x'),
                                self.fieldnames.index('y'),
                                self.fieldnames.index('z'),
                                self.fieldnames.index('px'),
                                self.fieldnames.index('py'),
                                self.fieldnames.index('pz'),
                                self.fieldnames.index('time')]

        for e_file in event_files:
            if '.csv' not in e_file:
                continue
            print(f'\nLoading data from file: {e_file}')
            with open(f'{event_folder}/{e_file}', 'r') as file:
                csv_reader = csv.reader(file)
                next(csv_reader)  # skip header
                for i, row in enumerate(csv_reader):
                    print(f'Processing detector hit : {i}', end="\r")
                    hit_converted = ([float(r) if row.index(r) in arguments_to_convert else r for r in row])
                    target_segment = s_manager.get_segment_at_known_xyz_value(hit_converted[self.fieldnames.index('x')],
                                                                              hit_converted[self.fieldnames.index('y')],
                                                                              hit_converted[self.fieldnames.index('z')],
                                                                              e_file,
                                                                              str(hit_converted[
                                                                                  self.fieldnames.index('layer')]))
                    target_segment.data.append(hit_converted)

                    if int(hit_converted[self.fieldnames.index('PDG')]) == 13:
                        self.muon_hits += 1
                        mc_particle_id = hit_converted[self.fieldnames.index('MC_particle_ID')]
                        if mc_particle_id in self.muon_tracks.keys():
                            self.muon_tracks[mc_particle_id].append(hit_converted)
                        else:
                            self.muon_tracks.update({mc_particle_id: [hit_converted]})
        print(f'\n{len(list(self.muon_tracks.keys()))} muon(s) caused {self.muon_hits} hits in the detector together')

    def create_doublet(self,
                       first_hit: list[float],
                       second_hit: list[float]) -> type(Doublet):
        """Creates a doublet from two hits.
        :param first_hit: first hit of doublet
        :param second_hit: second hit of doublet
        :return
            Doublet object
        """
        first_hit_particle_key = first_hit[self.fieldnames.index('PDG')]
        if int(first_hit_particle_key) != 13:
            first_hit_particle_key = randint(14, 1e16)
        second_hit_particle_key = second_hit[self.fieldnames.index('PDG')]
        if int(second_hit_particle_key) != 13:
            second_hit_particle_key = randint(14, 1e16)
            
        doublet = Doublet(hit_1_particle_key=first_hit_particle_key,
                          hit_2_particle_key=second_hit_particle_key,
                          hit_1_position=(first_hit[self.fieldnames.index('x')],
                                          first_hit[self.fieldnames.index('y')],
                                          first_hit[self.fieldnames.index('z')]),
                          hit_2_position=(second_hit[self.fieldnames.index('x')],
                                          second_hit[self.fieldnames.index('y')],
                                          second_hit[self.fieldnames.index('z')]),
                          hit_1_id=first_hit[self.fieldnames.index('MC_particle_ID')],
                          hit_2_id=second_hit[self.fieldnames.index('MC_particle_ID')],
                          time_1=first_hit[self.fieldnames.index('time')],
                          time_2=second_hit[self.fieldnames.index('time')])

        return doublet

    def create_xplet_list(self, s_manager) -> None:
        """Creates the doublets. Only distance to IP check, pure combinatorial.
        :param s_manager: segment manager object
        """
        print("-----------------------------------\n")
        print("Forming doublets ...\n")
        doublet_list_start = time.process_time()

        doublet_hits = []
        
        rejected_doublet_because_of_time = 0
        found_correct_doublets = 0

        num_segments = len(list(s_manager.segment_mapping_key.keys()))
        segment_process_counter = 0

        for name, segment in s_manager.segment_mapping_key.items():
            print(f'Processing segment {segment_process_counter} of {num_segments}', end='\r')
            for target_segment in s_manager.segment_mapping[name]:

                for hit_1 in segment.data:
                    phi_1 = np.arctan2(hit_1[self.fieldnames.index('y')],
                                       hit_1[self.fieldnames.index('x')])
                    for hit_2 in target_segment.data:
                        phi_2 = np.arctan2(hit_2[self.fieldnames.index('y')],
                                           hit_2[self.fieldnames.index('x')])
                        if abs(phi_2 - phi_1) > 0.1:
                            continue
                        if w_angle_diff(hit_1[self.fieldnames.index('x')],
                                        hit_2[self.fieldnames.index('x')],
                                        hit_1[self.fieldnames.index('y')],
                                        hit_2[self.fieldnames.index('y')],
                                        hit_1[self.fieldnames.index('z')],
                                        hit_2[self.fieldnames.index('z')]) > 0.01:
                            continue
                        self.create_doublet(hit_1, hit_2)
                        if "VXDTracker_7" in segment.name and "ITracker_0" in target_segment.name and \
                                hit_1[self.fieldnames.index('PDG')] == hit_2[self.fieldnames.index('PDG')] == 13:
                            print('found')
                        if hit_1[self.fieldnames.index('MC_particle_ID')] == \
                           hit_2[self.fieldnames.index('MC_particle_ID')]:
                            PDG_1 = hit_1[self.fieldnames.index('PDG')]
                            PDG_2 = hit_2[self.fieldnames.index('PDG')]
                            if int(PDG_1) == int(PDG_2) == 13:
                                doublet_hits.append([hit_1, hit_2])
            segment_process_counter += 1
        plt.figure(dpi=500)
        for item in doublet_hits:
            r_1 = np.sqrt(item[0][self.fieldnames.index('x')]**2 + item[0][self.fieldnames.index('y')]**2)
            r_2 = np.sqrt(item[1][self.fieldnames.index('x')]**2 + item[1][self.fieldnames.index('y')]**2)
            plt.plot([item[0][self.fieldnames.index('z')], item[1][self.fieldnames.index('z')]],
                     [r_1, r_2],
                     marker='o')
        # plt.xscale('log')
        # plt.yscale('log')
        plt.show()
        plt.savefig("tracking.pdf")

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
                    found_correct_doublets += 1
                    if muon_match:
                        pass

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


