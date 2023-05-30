import numpy as np
import matplotlib.pyplot as plt
import os
import sys


max_energy = 1000

directory = f'/nfs/dust/luxe/user/spatarod/Muon-Collider/MuonColliderCompleteStudy/signal_plus_full_background/{max_energy}GeV'
preselection = 'MuCo_example.yaml'

folders = os.listdir(directory)

triplet_time_diff_signal = set()
triplet_time_diff_background = set()

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


pt = []
num_triplets = []

num_samples = 0
for folder in folders:
    if os.path.isdir(f'{directory}/{folder}'):
        print(f'\nProcessing event: {folder.split("-")[-1]}')
        
        try:
            triplets = np.load(f'{directory}/{folder}/{preselection}/triplet_list.npy', allow_pickle=True)
            pt_temp = 0
            num_triplets.append(len(triplets))
            num_samples += 1

            for t in triplets:
                t_0, t_1, t_2 = t.doublet_1.hit_1_time, t.doublet_1.hit_2_time, t.doublet_2.hit_2_time
                max_t = max([t_0, t_1, t_2])
                min_t = min([t_0, t_1, t_2])
                if t.is_correct_match():
                    triplet_time_diff_signal.add(1000 * (max_t - min_t))
                    if pt_temp == 0:
                        p_t = np.sqrt(t.doublet_1.hit_1_momentum[0]**2 + t.doublet_1.hit_1_momentum[1]**2)
                        pt_temp = p_t
                else:
                    triplet_time_diff_background.add(1000 * (max_t - min_t))
            pt.append(pt_temp)
        except:
            pass

plt.figure(figsize=(9,6), dpi=400)
plt.title(f'Muon energy uniformly distributed from 0 - {max_energy} GeV', fontsize=18)
plt.hist(list(triplet_time_diff_signal), label='S', bins=25, histtype='step', linewidth=2.5)
plt.hist(list(triplet_time_diff_background), label='BG', bins=25, histtype='step', linewidth=2.5)
plt.xlabel('triplet time spread [ps]', fontsize=16, loc='right')
plt.ylabel('#triplets', fontsize=16, loc='top')
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.yscale('log')
plt.legend(loc='best', fontsize=16)
plt.savefig(f'MuCo_triplet_time_{max_energy}GeV.pdf')
plt.close()

plt.figure(figsize=(9,6), dpi=400)
plt.title(f'Muon energy uniformly distributed from 0 - {max_energy} GeV', fontsize=18)
plt.plot(pt, num_triplets, label=f'N={num_samples} samples',linestyle='', marker='o')
plt.xlabel('truth $p_T$ [GeV]', fontsize=16, loc='right')
plt.ylabel('#triplets', fontsize=16, loc='top')
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.yscale('log')
plt.legend(loc='best', fontsize=16)
plt.savefig(f'MuCo_triplet_numbers_{max_energy}GeV.pdf')
plt.close()