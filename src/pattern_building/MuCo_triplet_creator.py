import csv
import time
import os

import numpy as np
import matplotlib.pyplot as plt

from math_functions.geometry import w_rz_angle_diff, w_phi_angle_diff
from utility.time_tracking import hms_string
from pattern.doublet import Doublet
from pattern.triplet import Triplet

from random import randint


class MuCoTripletCreator:
    def __init__(self, event_folder):
        """Class for creating triplets from Muon Collider detector hits.
        :param event_folder folder containing .csv files from a sinlge muon collider event
                """
        self.event_folder = event_folder

        self.doublet_creation_time = 0
        self.triplet_creation_time = 0
        self.num_signal_events = 0
        self.num_all_doublets = 0
        self.correct_doublets_tracker = set()
        self.correct_triplets_tracker = set()
        self.num_all_triplets = 0

        self.num_hits_vxd_barrel = 0
        self.num_hits_inner_tracker_barrel = 0
        self.num_hits_outer_tracker_barrel = 0

        self.num_hits_vxd_endcap = 0
        self.num_hits_inner_tracker_endcap = 0
        self.num_hits_outer_tracker_endcap = 0

        self.muon_hits = 0

        self.triplet_list = []

        # fieldnames according to .csv naming
        self.fieldnames = ['hit_ID', 'x', 'y', 'z', 'layer', 'MC_particle_ID', 'px', 'py', 'pz', 'time', 'PDG', 'event']

    def load_tracking_data(self,
                           dlfilter: bool,
                           s_manager,
                           no_endcaps=True) -> None:
        """Loads data from the event stored in various .csv files inside the specified folder.
        Also converts values which are used for calculations to float from string.

        :param dlfilter: if True, then the double layer filtered hits are used instead of the whole VXDTracker sample
        :param s_manager: segment manager
        :param no_endcaps: if True no data from endcaps are used
        """
        print(f'Processing folder: {self.event_folder}')
        event_files = [f for f in os.listdir(self.event_folder) if '.csv' in f]

        # Check if using Double Layer Filtered hits or all hits from the Vertex detector
        if not dlfilter:
            event_files = [f for f in event_files if "_DLFiltered_" not in f]
        else:
            event_files = [f for f in event_files if "_VXD_" not in f]

        # Check if endcap region is included
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
            print(f'\nLoading data from file: {e_file}')
            num_detector_hits = 0
            with open(f'{self.event_folder}/{e_file}', 'r') as file:
                csv_reader = csv.reader(file)
                next(csv_reader)  # skip header
                for i, row in enumerate(csv_reader):
                    num_detector_hits += 1
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
            if 'DLFiltered' or '_VXDTracker_' in e_file:
                self.num_hits_vxd_barrel += num_detector_hits
            if '_VXDTrackerEndcap_' in e_file:
                self.num_hits_vxd_endcap += num_detector_hits
            if '_ITracker_' in e_file:
                self.num_hits_inner_tracker_barrel += num_detector_hits
            if '_ITrackerEndcap_' in e_file:
                self.num_hits_inner_tracker_endcap += num_detector_hits
            if '_OTracker_' in e_file:
                self.num_hits_outer_tracker_barrel += num_detector_hits
            if '_OTrackerEndcap_' in e_file:
                self.num_hits_outer_tracker_endcap += num_detector_hits

        print(f'\n{self.muon_hits} muon hit(s) found in the detector!')

    def create_doublet(self,
                       first_hit: list[float],
                       second_hit: list[float]) -> type(Doublet):
        """Creates a doublet from two hits.
        :param first_hit: first hit of doublet
        :param second_hit: second hit of doublet
        :return
            Doublet object
        """
        hit_1_id = '_'.join([first_hit[self.fieldnames.index('layer')],
                             first_hit[self.fieldnames.index('hit_ID')]])
        hit_2_id = '_'.join([second_hit[self.fieldnames.index('layer')],
                             second_hit[self.fieldnames.index('hit_ID')]])

        doublet = Doublet(hit_1_particle_key=first_hit[self.fieldnames.index('MC_particle_ID')],
                          hit_2_particle_key=second_hit[self.fieldnames.index('MC_particle_ID')],
                          hit_1_position=(first_hit[self.fieldnames.index('x')],
                                          first_hit[self.fieldnames.index('y')],
                                          first_hit[self.fieldnames.index('z')]),
                          hit_2_position=(second_hit[self.fieldnames.index('x')],
                                          second_hit[self.fieldnames.index('y')],
                                          second_hit[self.fieldnames.index('z')]),
                          hit_1_id=hit_1_id,
                          hit_2_id=hit_2_id,
                          time_1=first_hit[self.fieldnames.index('time')],
                          time_2=second_hit[self.fieldnames.index('time')])

        return doublet

    def create_xplet_list(self,
                          s_manager,
                          configuration) -> None:
        """Creates the doublets. Only distance to IP check, pure combinatorial.
        :param s_manager: segment manager object
        :param configuration: Muon Collider qubo preparation configuration file
        """
        sigma_time_vxd = configuration['doublet']['sigma time VXD']
        sigma_time_inner_tracker = configuration['doublet']['sigma time inner tracker']
        sigma_time_outer_tracker = configuration['doublet']['sigma time outer tracker']
        rz_doublet = configuration['doublet']['rz']
        phi_doublet = configuration['doublet']['phi']

        print("-----------------------------------\n")
        print("Forming doublets ...\n")
        doublet_list_start = time.process_time()

        rejected_doublet_because_of_time = 0

        num_segments = len(list(s_manager.segment_mapping_key.keys()))
        segment_process_counter = 0

        hit_1_time_check = None
        hit_2_time_check = None

        for name, segment in s_manager.segment_mapping_key.items():
            print(f'Processing segment {segment_process_counter} of {num_segments}', end='\r')

            # Set time check information for first hit according to the layer
            if "VXD" in segment.name:
                if sigma_time_vxd:
                    hit_1_time_check = sigma_time_sigma_time_vxd
            if "ITracker" in segment.name:
                if sigma_time_inner_tracker:
                    hit_1_time_check = sigma_time_inner_tracker
            if "OTracker" in segment.name:
                if sigma_time_outer_tracker:
                    hit_1_time_check = sigma_time_outer_tracker

            for target_segment in s_manager.segment_mapping[name]:

                # Set time check information for first hit according to the layer
                if "VXD" in target_segment.name:
                    if sigma_time_vxd:
                        hit_2_time_check = sigma_time_vxd
                if "ITracker" in target_segment.name:
                    if sigma_time_inner_tracker:
                        hit_2_time_check = sigma_time_inner_tracker
                if "OTracker" in target_segment.name:
                    if sigma_time_outer_tracker:
                        hit_2_time_check = sigma_time_outer_tracker

                for hit_1 in segment.data:
                    if hit_1_time_check:
                        if abs(hit_1[self.fieldnames['time']]) > hit_1_time_check:
                            continue
                    for hit_2 in target_segment.data:
                        if hit_2_time_check:
                            if abs(hit_2[self.fieldnames['time']]) > hit_2_time_check:
                                continue

                        # check difference in phi
                        if w_phi_angle_diff(hit_1[self.fieldnames.index('x')],
                                            hit_2[self.fieldnames.index('x')],
                                            hit_1[self.fieldnames.index('y')],
                                            hit_2[self.fieldnames.index('y')]) > phi_doublet:
                            continue

                        # check difference in rz
                        if w_rz_angle_diff(hit_1[self.fieldnames.index('x')],
                                           hit_2[self.fieldnames.index('x')],
                                           hit_1[self.fieldnames.index('y')],
                                           hit_2[self.fieldnames.index('y')],
                                           hit_1[self.fieldnames.index('z')],
                                           hit_2[self.fieldnames.index('z')]) > rz_doublet:
                            continue
                        # create doublet object
                        doublet = self.create_doublet(hit_1, hit_2)
                        segment.doublet_data.append(doublet)
                        self.num_all_doublets += 1
                        pdg_1 = hit_1[self.fieldnames.index('PDG')]
                        pdg_2 = hit_2[self.fieldnames.index('PDG')]

                        # correct_doublet = doublet stemming from two muon hits
                        if int(pdg_1) == int(pdg_2) == 13:
                            self.correct_doublets_tracker.add('_'.join(segment.name.split('_')[0:2]))
            segment_process_counter += 1

        # Give estimate if doublet building procedure was effective
        print(f"Found correct doublets: {len(self.correct_doublets_tracker)} --> "
              f"{100 * np.around(len(self.correct_doublets_tracker) / 13, 2)} %")
        print(f"Rejected doublets because of time: {rejected_doublet_because_of_time}")
                        
        doublet_list_end = time.process_time()
        self.doublet_creation_time = hms_string(doublet_list_end - doublet_list_start)

        print(f"Time elapsed for forming doublets: "
              f"{self.doublet_creation_time}")
        print(f"Number of doublets found: {self.num_all_doublets}")

        print("-----------------------------------\n")
        print("Forming triplets ...\n")
        list_triplet_start = time.process_time()

        for name, segment in s_manager.segment_mapping_key.items():
            print(f'Processing segment {segment_process_counter} of {num_segments}', end='\r')
            for target_segment in s_manager.segment_mapping[name]:
                if not target_segment.doublet_data:
                    continue
                for d1 in segment.doublet_data:
                    for d2 in target_segment.doublet_data:
                        if d1.hit_2_id == d2.hit_1_id:
                            triplet = Triplet(d1, d2)
                            segment.triplet_data.append(triplet)
                            self.num_all_triplets += 1
                            if d1.hit_1_particle_key != '-1000' and \
                                    d1.hit_2_particle_key != '-1000' and\
                                    d2.hit_2_particle_key != '-1000':
                                self.correct_triplets_tracker.add('_'.join(segment.name.split('_')[0:2]))

        list_triplet_end = time.process_time()
        self.triplet_creation_time = hms_string(list_triplet_end - list_triplet_start)
        print(f"Found correct triplets: {len(self.correct_triplets_tracker)} --> "
              f"{100 * np.around(len(self.correct_triplets_tracker) / 12, 2)} %")
        print(f"Time elapsed for  forming triplets: "
              f"{self.triplet_creation_time}")
        print(f"Number of triplets found: {self.num_all_triplets}\n")

    def set_qubo_coefficients(self,
                              s_manager,
                              configuration: dict,
                              save_to_folder: str) -> None:
        """Sets the QUBO coefficients.
        :param s_manager: segment manager
        :param configuration: configuration of pattern building
        :param save_to_folder: folder to save triplet list
        """
        print("-----------------------------------\n")
        print("Setting QUBO coefficients...\n")
        connection = configuration['qubo parameters']['b_ij match']
        conflict = configuration['qubo parameters']['b_ij conflict']

        num_segments = len(list(s_manager.segment_mapping_key.keys()))
        segment_process_counter = 0

        for name, segment in s_manager.segment_mapping_key.items():
            print(f'Processing segment {segment_process_counter} of {num_segments}', end='\r')
            for target_segment in s_manager.segment_mapping[name] + [segment]:
                if not target_segment.triplet_data:
                    continue
                for t1 in segment.triplet_data:
                    for t2 in target_segment.triplet_data:
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
                            if segment.name == target_segment.name:
                                t1.interactions.update({t2.triplet_id: conflict})
                                t2.interactions.update({t1.triplet_id: conflict})
                            else:
                                t1.interactions.update({t2.triplet_id: connection})
                                t2.interactions.update({t1.triplet_id: connection})

        for segment in s_manager.segment_mapping_key.values():
            if not segment.triplet_data:
                continue
            for triplet in segment.triplet_data:
                self.triplet_list.append(triplet)

        print("-----------------------------------\n")
        print(f"Save triplet list...\n")
        np.save(f"{save_to_folder}/triplet_list", self.triplet_list)

    def write_info_file(self,
                        save_to_folder: str,
                        configuration: dict) -> None:
        """Writes information about the Preselection parameters and some statistics into
        'preselection_info.txt' which is stored inside the output folder.
        :param save_to_folder: folder to write the information into
        :param configuration: configuration of the pattern building procedure
        """
        with open(save_to_folder + "/preselection_info.txt", "w") as f:
            f.write("Preselection performed with the following configuration: \n")
            for outer_key in configuration.keys():
                for inner_key, value in configuration[outer_key].items():
                    f.write(f"\n{inner_key}: {value}")
                f.write("\n")
            f.write("\n\n")
            f.write("---\n")
            f.write(f"Number of hits in VXDTrackerBarrel: {self.num_hits_vxd_barrel}\n")
            f.write(f"Number of hits in VXDEndcap: {self.num_hits_vxd_endcap}\n")
            f.write(f"Number of hits InnerTrackerBarrel: {self.num_hits_inner_tracker_barrel}\n")
            f.write(f"Number of hits InnerTrackerEndcap: {self.num_hits_inner_tracker_endcap}\n")
            f.write(f"Number of hits OuterTrackerBarrel: {self.num_hits_outer_tracker_barrel}\n")
            f.write(f"Number of hits OuterTrackerEndcap: {self.num_hits_outer_tracker_endcap}\n")

            f.write("\n\n")

            f.write(f"Time elapsed for forming doublets: "
                    f"{self.doublet_creation_time}\n")
            f.write(f"Number of doublets found: {self.num_all_doublets}\n")
            f.write(f"Doublet selection efficiency: "
                    f"{np.around(100 * len(self.correct_doublets_tracker) / 13,  3)} %\n")
            f.write("\n\n")

            f.write(f"Time elapsed for creating triplets: "
                    f"{self.triplet_creation_time}\n")
            f.write(f"Number of triplets found: {self.num_all_triplets}\n")
            f.write(f"Triplet selection efficiency: "
                    f"{np.around(100 * len(self.correct_triplets_tracker) / 12, 3)} %\n")
