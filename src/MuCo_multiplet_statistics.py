import numpy as np
import matplotlib.pyplot as plt
import os
import sys


track_length_requirement = 8

directory = '../../muon_collider_files'
preselection = 'MuCo_example.yaml'

folders = os.listdir(directory)

triplet_time_diff_signal = []
triplet_time_diff_background = []

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

doublets_used = []

for folder in folders:
    if os.path.isdir(f'{directory}/{folder}'):
        print(f'\nProcessing event: {folder.split("-")[-1]}')
        triplets = np.load(f'{directory}/{folder}/{preselection}/triplet_list.npy', allow_pickle=True)

        for t in triplets:
            t_0, t_1, t_2 = t.doublet_1.hit_1_time, t.doublet_1.hit_2_time, t.doublet_2.hit_2_time
            max_t = max([t_0, t_1, t_2])
            min_t = min([t_0, t_1, t_2])
            if t.is_correct_match():
                triplet_time_diff_signal.append(1000 * (max_t - min_t))
            else:
                triplet_time_diff_background.append(1000 * (max_t - min_t))

plt.hist(triplet_time_diff_signal, label='S', bins=25, histtype='step', linewidth=2.5)
plt.hist(triplet_time_diff_background, label='BG', bins=25, histtype='step', linewidth=2.5)
plt.xlabel('diff time [ps]')
plt.ylabel('#triplets')
plt.yscale('log')
plt.legend(loc='best')
plt.show()

