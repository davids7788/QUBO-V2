import numpy as np
import csv
import os
import sys
from numpy.linalg import norm

folder = sys.argv[1]
sys.path.append("..")

os.chdir(folder)
if not os.path.isdir("true"):
    os.mkdir("true")
if not os.path.isdir("smeared"):
    os.mkdir("smeared")


for file in os.listdir():
    if ".npy" not in file:
        continue
    data = np.load(file, allow_pickle=True)[()]

    hit_id = 0
    with open('true/true_hits_' + '.'.join(file.split('.')[0:-1]) + '.csv', 'w', newline='') as csvfile:
        fieldnames = ['hit_ID', 'x', 'y', 'z', 'layer_ID', 'particle_ID', 'particle_energy']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for key, plane in zip(data['Plane list'].keys(), data['Plane list'].values()):
            for particle, hit in zip(plane.true_hits_dictionary.keys(), plane.true_hits_dictionary.values()):
                writer.writerow({'hit_ID': hit_id,
                                 'x': hit[0],
                                 'y': hit[1],
                                 'z': hit[2],
                                 'layer_ID': key,
                                 'particle_ID': particle,
                                 'particle_energy': norm(data["Particle momentum history log"][particle][key])})
                hit_id += 1
        
    hit_id = 0
    with open('smeared/smeared_hits_' + '.'.join(file.split('.')[0:-1]) + '.csv', 'w', newline='') as csvfile:
        fieldnames = ['hit_ID', 'x', 'y', 'z', 'layer_ID', 'particle_ID', 'particle_energy']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for key, plane in zip(data['Plane list'].keys(), data['Plane list'].values()):
            sigma_x_resolution = (plane.limits_x[1] - plane.limits_x[0]) / (plane.num_bins[0]) / np.sqrt(12)
            sigma_y_resolution = (plane.limits_y[1] - plane.limits_y[0]) / (plane.num_bins[1]) / np.sqrt(12)
            for particle, hit in zip(plane.true_hits_dictionary.keys(), plane.true_hits_dictionary.values()):
                smeared_x = np.random.normal(hit[0], sigma_x_resolution)
                smeared_y = np.random.normal(hit[0], sigma_y_resolution)
                writer.writerow({'hit_ID': hit_id,
                                 'x': smeared_x,
                                 'y': smeared_y,
                                 'z': hit[2],
                                 'layer_ID': key,
                                 'particle_ID': particle,
                                 'particle_energy': norm(data["Particle momentum history log"][particle][key])})
                hit_id += 1
