import os
import sys
import numpy as np
sys.path.append("../src")
from utility.time_tracking import hms_string
from pattern.doublet import Doublet
from pattern.triplet import Triplet
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

energy = '1000GeV'

tracks = []
combinatorial_triplets = set()


t_complete = []

def radial_distance(triplet):
    return np.sqrt(triplet.doublet_1.hit_1_position[0]**2 + triplet.doublet_1.hit_1_position[1]**2)


for i in range(250):
    prefix = f'/nfs/dust/luxe/user/spatarod/Muon-Collider/DPG2023-update/{energy}/{i}'
    triplet_list = np.load(f'{prefix}/triplet_list.npy', allow_pickle=True)
    for t_i in triplet_list:
        t_complete.append(t_i)
    
    folder_entry = os.listdir(prefix)
    folder = None
    for file_or_folder in folder_entry:
        if 'eigensolver_7q_connection_list' in file_or_folder:
            folder = file_or_folder     

    qubo_result = np.load(f"{prefix}/{folder}/qubo_log.npy", allow_pickle=True)[()]
    print(f'Load {prefix}/triplet_list.npy')
    print(f'Processing {prefix}/{folder}/qubo_log.npy \n')
    
    triplets_dict = {"L-012": [],
                     "L-123": [],
                     "L-234": [],
                     "L-345": [],
                     "L-456": [],
                     "L-567": []}
   
    for triplet, result in zip(triplet_list, qubo_result["computed solution vector"]):
        if result == 1:
            d = radial_distance(triplet)
            if 30.0 < d < 31.4:
                triplets_dict['L-012'].append(triplet)
            elif 32.1 < d < 33.4:
                triplets_dict['L-123'].append(triplet)
            elif 51.0 < d < 53.0:
                triplets_dict['L-234'].append(triplet)
            elif 53.1 < d < 54.9:
                triplets_dict['L-345'].append(triplet)
            elif 74.0 < d < 75.5:
                triplets_dict['L-456'].append(triplet)
            elif 76.1 < d < 77.5:
                triplets_dict['L-567'].append(triplet)

    full_tracks = []
    for t_start in triplets_dict['L-012']:
        f_track = [t_start]
        for j in range(5):
            for t_next in triplets_dict[f'L-{j + 1}{j + 2}{j + 3}']:
                if t_next.doublet_1 == f_track[-1].doublet_2:
                    f_track.append(t_next)
        if len(f_track) == 6:
            full_tracks.append(f_track)
            
    for full_track in full_tracks:
        count = 0
        for triplet in full_track:
            if triplet.doublet_1.hit_1_particle_key != 13:
                count += 1
        if full_track[-1].doublet_2.hit_1_particle_key != 13:
            count += 1
        if full_track[-1].doublet_2.hit_2_particle_key != 13:
            count += 1
        if count > 3:
            for triplet in full_track:
                combinatorial_triplets.add(triplet)
                   
        else:
            tracks.append(full_track)    
        
    for triplet, result in zip(triplet_list, qubo_result["computed solution vector"]):
        if result == 1:
            for t in full_tracks:
                if triplet not in t:
                    combinatorial_triplets.add(triplet)

fig = plt.figure(figsize=(12,12), dpi=500)
ax = fig.add_subplot()
 
for triplet in combinatorial_triplets:
    ax.plot([triplet.doublet_1.hit_1_position[0],
             triplet.doublet_1.hit_2_position[0],
             triplet.doublet_2.hit_2_position[0]], 
            [triplet.doublet_1.hit_1_position[1],
             triplet.doublet_1.hit_2_position[1],
             triplet.doublet_2.hit_2_position[1]], marker="o", c="k", mfc='k', markersize=12, alpha=0.7)
    
for i, track in enumerate(tracks):
    if i > 20:
        break
    for triplet in track:
        ax.plot([triplet.doublet_1.hit_1_position[0],
                 triplet.doublet_1.hit_2_position[0],
                 triplet.doublet_2.hit_2_position[0]], 
                [triplet.doublet_1.hit_1_position[1],
                 triplet.doublet_1.hit_2_position[1],
                 triplet.doublet_2.hit_2_position[1]], marker="o", c="blue", mfc='red', markersize=12)

# for triplet in t_complete:
#     ax.plot([triplet.doublet_1.hit_1_position[0],
#              triplet.doublet_1.hit_2_position[0],
#              triplet.doublet_2.hit_2_position[0]], 
#             [triplet.doublet_1.hit_1_position[1],
#              triplet.doublet_1.hit_2_position[1],
#              triplet.doublet_2.hit_2_position[1]], marker="o", c="blue", mfc='red', markersize=12)


legend_elements = [Line2D([0], [0], marker='o', mfc='r', color='w', label='tracks', markersize=12),
                   Line2D([0], [0], marker='o', mfc='k', color='w', label='background', markersize=12)]    
                           
ax.tick_params(axis="x", labelsize=20) 
ax.tick_params(axis="y", labelsize=20) 
ax.legend(handles=legend_elements, loc='best', fontsize=20)                          
 
plt.xlabel("x [mm]", fontsize=20)
plt.ylabel("y [mm]", fontsize=20)
plt.xlim(-135, 135)
plt.ylim(-135, 135)
plt.title(f'$E_{{muon}}$ = {energy}, track reconstruction efficiency for 250 events: {100 * np.around(len(tracks) / 250, 3)} %', fontsize=20, loc='left') 

plt.savefig(f"reconstructed_tracks_{energy}.pdf")
plt.savefig(f"reconstructed_tracks_{energy}.jpg")