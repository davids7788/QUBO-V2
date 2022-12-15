import os
import sys
from src.detector_plane import DetectorPlane
import matplotlib.pyplot as plt
import numpy as np


def acceptance_plots_simplified_LUXE(file):
    """Takes a file and makes acceptance plots for that file.
    :param file: .npy file with simplified tracking data.
    """
    data = np.load(file, allow_pickle=True)[()]
    num_particles = len(list(data["Particle momentum history log"].keys()))
    num_of_hits = {}
    full_acceptance_energy = []
    partial_acceptance_energy = []
    no_acceptance_energy = []
    for particle in list(data["Particle momentum history log"].keys()):
        counter = 0
        for key in data["Particle momentum history log"][particle].keys():
            if "Plane" in key:
                counter += 1
        num_of_hits.update({particle: counter})
        if counter == 0:
            no_acceptance_energy.append(np.linalg.norm(data["Particle momentum history log"][particle]['IP']))
        elif 0 < counter < 4:
            partial_acceptance_energy.append(np.linalg.norm(data["Particle momentum history log"][particle]['IP']))
        else:
            full_acceptance_energy.append(np.linalg.norm(data["Particle momentum history log"][particle]['IP']))

    plt.figure(figsize=(8, 6), dpi=250)
    plt.hist(num_of_hits.values(),
             weights=np.ones(num_particles) / num_particles,
             bins=8,
             range=(0, 8),
             align='left',
             histtype='step',
             linewidth=2.5,
             label=f"N = {num_particles}\n" +
             f"a(>0) = {np.round(sum([1 for value in num_of_hits.values() if value > 0]) / num_particles, 3)}\n" 
             f"a(>3) = {np.round(sum([1 for value in num_of_hits.values() if value > 3]) / num_particles, 3)}")
    plt.yscale('log')
    plt.xticks(fontsize=14)
    plt.xlabel("Number of detector interactions", fontsize=14)
    plt.ylabel("Fraction of particles", fontsize=14)
    plt.yticks(fontsize=14)
    plt.title("Geometrical detector acceptance", fontsize=16)
    plt.legend(loc='upper left', fontsize=14)
    plt.grid(True, which="both", ls="-")
    plt.savefig('acceptance/' + '.'.join(file.split('.')[0:-1]) + '_acceptance.pdf')
    plt.close()

    plt.figure(figsize=(8, 6), dpi=250)
    plt.hist(np.array(full_acceptance_energy) / 1000,
             weights=np.ones(len(full_acceptance_energy)) / num_particles,
             bins=70,
             range=(0, 14),
             align='left',
             histtype='step',
             linewidth=2.5,
             label="> 3 interactions")
    plt.hist(np.array(partial_acceptance_energy) / 1000,
             weights=np.ones(len(partial_acceptance_energy)) / num_particles,
             bins=70,
             range=(0, 14),
             align='left',
             histtype='step',
             linewidth=2.5,
             label="0 < i < 4 interactions")
    plt.hist(np.array(no_acceptance_energy) / 1000,
             weights=np.ones(len(no_acceptance_energy)) / num_particles,
             bins=70,
             range=(0, 14),
             align='left',
             histtype='step',
             linewidth=2.5,
             label="0 interactions")
    plt.hist([],
             label=f"N = {num_particles}",
             histtype="step",
             linewidth=0)

    plt.yscale('log')
    plt.xticks(fontsize=14)
    plt.xlabel("Energy [GeV]", fontsize=14)
    plt.ylabel("Fraction of particles / 0.5 GeV", fontsize=14)
    plt.yticks(fontsize=14)
    plt.title("Geometrical acceptance vs. energy", fontsize=16)
    plt.legend(loc='upper right', fontsize=14)
    plt.grid(True, which="both", ls="-")
    plt.savefig('acceptance/' + '.'.join(file.split('.')[0:-1]) + '_acceptance_vs_energy.pdf')
    plt.close()
