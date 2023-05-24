import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import circle_fit as cf
from math import ceil

directory = '../../muon_collider_files'
preselection = 'MuCo_example.yaml'

folders = os.listdir(directory)

reconstructed_muon = []
pt = []
reconstructed_correct_doublets = []
reconstructed_correct_triplets = []


layer_info = {0: [30.1, 31.5],
              1: [32.1, 33.5],
              2: [51.1, 53.05],
              3: [53.1, 55.0],
              4: [74.1, 75.8],
              5: [76.1, 77.8],
              6: [102.1, 103.2],
              7: [104.1, 105.2],
              8: [126.0, 131.5],
              9: [338.0, 343.5],
              10: [553.5, 557.0],
              11: [817.0, 824.0],
              12: [1150.0, 1158.0],
              13: [1486.0, 1490]}


def radius(t_object):
    return np.sqrt(t_object.doublet_1.hit_1_position[0]**2 + t_object.doublet_1.hit_1_position[1]**2)


def pT(t_obj):
    pt_values = []
    if t_obj[0].doublet_1.hit_1_pdg == '13':
        pt_values.append(np.sqrt(t_obj[0].doublet_1.hit_1_momentum[0]**2 + t_obj[0].doublet_1.hit_1_momentum[1]**2))
    if t_obj[0].doublet_1.hit_2_pdg == '13':
        pt_values.append(np.sqrt(t_obj[0].doublet_1.hit_2_momentum[0]**2 + t_obj[0].doublet_1.hit_2_momentum[1]**2))
    if t_obj[0].doublet_2.hit_2_pdg == '13':
        pt_values.append(np.sqrt(t_obj[0].doublet_2.hit_2_momentum[0]**2 + t_obj[0].doublet_2.hit_2_momentum[1]**2))
    for t in t_obj[1:]:
        if t.doublet_2.hit_2_pdg == '13':
            pt_values.append(np.sqrt(t.doublet_2.hit_2_momentum[0] ** 2 + t.doublet_2.hit_2_momentum[1] ** 2))

    return sum(pt_values) / len(pt_values)


def build_tracks(t_sorted):
    tracks = []
    for start_triplet in t_sorted[0]:

        track_candidates = [[start_triplet]]
        for j in range(1, len(list(layer_info))):
            for t_next in t_sorted[j]:
                for part_track in track_candidates:
                    if part_track[-1].doublet_2 == t_next.doublet_1:
                        new_part_track = [tx for tx in part_track]
                        new_part_track.append(t_next)
                        track_candidates.append(new_part_track)
        for t_candidate in track_candidates:
            tracks.append(t_candidate)

    max_length_tracks = 0
    for tr in tracks:
        if len(tr) > max_length_tracks:
            max_length_tracks = len(tr)

    for tr in tracks:
        if len(tr) == max_length_tracks:
            return tr    # rework if ambiguity solving is enabled!!!!!!!


for folder in folders:
    if os.path.isdir(f'{directory}/{folder}'):
        triplets = np.load(f'{directory}/{folder}/{preselection}/triplet_list.npy', allow_pickle=True)
        correct = []
        for t in triplets:
            if t.is_correct_match():
                correct.append((t.doublet_1.hit_1_position[0], t.doublet_1.hit_1_position[1]))

        xc, yc, r, _ = cf.least_squares_circle(correct)
        print(r * )
        exit()
        result_folder = os.listdir(f'{directory}/{folder}/{preselection}')[0]
        qubo_log = np.load(f'{directory}/{folder}/{preselection}/{result_folder}/qubo_log.npy', allow_pickle=True)[()]

        kept_triplets = {0: [],
                         1: [],
                         2: [],
                         3: [],
                         4: [],
                         5: [],
                         6: [],
                         7: [],
                         8: [],
                         9: [],
                         10: [],
                         11: [],
                         12: [],
                         13: []}

        for t, result in zip(triplets, qubo_log['computed solution vector']):
            if result == 1:
                r = radius(t)
                for key, value in layer_info.items():
                    if value[0] <= r <= value[1]:
                        kept_triplets[key].append(t)

        tr = build_tracks(kept_triplets)
        if tr is not None:
            if len(tr) > ceil(len(list(layer_info.keys())) / 2):
                num_correct = 0

                if tr[0].doublet_1.hit_1_pdg == '13':
                    num_correct += 1
                if tr[0].doublet_1.hit_2_pdg == '13':
                    num_correct += 1
                if tr[0].doublet_2.hit_2_pdg == '13':
                    num_correct += 1
                for element in tr[1:]:
                    if element.doublet_2.hit_2_pdg == '13':
                        num_correct += 1
                if num_correct >= ceil(len(tr) / 2):
                    p_t = pT(tr)
                    if p_t:
                        print(f'Reconstructed track with {np.around(p_t, 2)} GeV')
                else:
                    print(f'Track candidate mostly/purely combinatorial: length of track = {len(tr)} hits, '
                          f'num_correct_hits = {num_correct}')
            else:
                print('No track candidate with at least 7 hits found')


