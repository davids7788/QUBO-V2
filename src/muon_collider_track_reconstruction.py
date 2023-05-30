import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import circle_fit as cf
from math import ceil


energy = 50
track_length_requirement = 8

qubo = 'eigensolver_10q_impact_list_reverse'

directory = f'/nfs/dust/luxe/user/spatarod/Muon-Collider/MuonColliderCompleteStudy/signal_plus_full_background_test_set/{energy}GeV'
preselection = 'MuCo_example.yaml'

folders = os.listdir(directory)

pt_truth = {}
pt_reconstructed = {}
signal_or_background = {}

reconstruction_status = []

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
    if not t_obj:
        return -10
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
        if len(tr) + 2 >= max_length_tracks:
            max_length_tracks = len(tr) + 2

    best_choice = [None, 0]
    for tr in tracks:
        if len(tr) + 2 == max_length_tracks:
            connection_sum = 0
            for t1, t2 in zip(tr[0:-1], tr[1:]):
                connection_sum += t2.interactions[t1.triplet_id]
            if connection_sum < best_choice[1]:
                best_choice = [tr, connection_sum]
    return best_choice[0]


for folder in folders:
    if os.path.isdir(f'{directory}/{folder}'):
        print(f'\nProcessing event: {folder.split("-")[-1]}')
        try:
            triplets = np.load(f'{directory}/{folder}/{preselection}/triplet_list.npy', allow_pickle=True)
        except FileNotFoundError:
            continue
        pT_from_truth_triplets = []
        for t in triplets:
            if t.is_correct_match():
                pT_from_truth_triplets.append(t)

        pT_true = pT(pT_from_truth_triplets)
        print(f'Signal track truth p_T = {np.around(pT_true, 2)} GeV')
        pt_truth.update({folder: pT_true})

        try:
            qubo_log = np.load(f'{directory}/{folder}/{preselection}/{qubo}/qubo_log.npy', allow_pickle=True)[()]
        except FileNotFoundError:
            continue
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
            if len(tr) >= track_length_requirement:
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
                # check if purely combinatorial or stemming from signal
                if num_correct >= ceil(len(tr) / 2):
                    signal_or_background.update({folder: True})
                else:
                    signal_or_background.update({folder: False})
                pt_reconstruction = []

                # for t in tr:
                #     pt_reconstruction.append((t.doublet_1.hit_1_position[0], t.doublet_1.hit_1_position[1]))
                # xc, yc, r, _ = cf.least_squares_circle(pt_reconstruction)

                # p_t = (1 / 1000 * 0.3 * 3.57 * r)
                # pt_reconstructed.update({folder: p_t})

                # print(f'Reconstructed track with p_T = {np.around(p_t, 2)} GeV')
                print(f'Length of track = {len(tr) + 2} hits, num_correct_hits = {num_correct}')
                reconstruction_status.append((pT_true, 1))
            else:
                reconstruction_status.append((pT_true, 0))
                print(f'No reconstructed track candidate with at least {track_length_requirement} hits found')



pt_failed = []
pt_succeeded = []

for result in reconstruction_status:
    if result[0] != -10:
        if result[1] == 1:
            pt_succeeded.append(result[0])
        else:
            pt_failed.append(result[0])
                
plt.figure(figsize=(9,6), dpi=400)
plt.hist(pt_succeeded, label='succeeded', linewidth=2.5, bins=20, histtype='step')
plt.hist(pt_failed, label='failed', linewidth=2.5, bins=20, histtype='step')
plt.title(f'Reconstruction Efficiency = {np.around(100 * len(pt_succeeded) / (len(pt_succeeded) + len(pt_failed)), 1)} %', fontsize=18)
plt.xlabel('Energy [GeV]', fontsize=18)
plt.ylabel('# Reconstrcuted tracks', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(loc='best', fontsize=18)
plt.savefig(f'Reconstrcution_effciciency_{preselection}_{energy}.pdf')











