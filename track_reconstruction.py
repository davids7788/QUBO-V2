import sys
import csv
from typing import List, Any, Union

import matplotlib.pyplot as plt
import numpy as np

from src.triplet import Triplet
from src.doublet import Doublet
from src.track import Track

folder = sys.argv[1]
xi = folder.split("_")[2]

data = np.load(folder + "/qubo_log.npy", allow_pickle=True)

generated_tracks = np.load("/".join(folder.split("/")[0:-1]) + "/track_list.npy", allow_pickle=True)
reconstructed_tracks = data[()]["reconstructed tracks"]

matched_tracks = [track for track in reconstructed_tracks if track.is_matched_track]
fake_tracks = [track for track in reconstructed_tracks if not track.is_matched_track]


def get_track_energy(track_for_energy):
    """Sums over triplets in a track and calculates the mean energy.
     :return
        energy of track
    """
    return track_for_energy.triplets[0].doublet_1.energy_1


def get_track_angle(track_for_angle):
    """Sums over triplets in a track and calculates the mean energy.
     :return
        angle of track
    """
    return track_for_angle.triplets[0].doublet_1.xz_angle()


def get_standard_deviation(k, n):
    if k == 0:
        return 0
    if k > n:
        k = n
    return np.sqrt(((k + 1) * (k + 2)) / ((n + 2) * (n + 3)) - ((k + 1)**2 / (n + 2)**2))


angle_distribution_reconstructed = [get_track_angle(track) for track in reconstructed_tracks]
angle_distribution_generated = [get_track_angle(track) for track in generated_tracks]
angle_distribution_matched = [get_track_angle(track) for track in matched_tracks]
angle_distribution_fake = [get_track_angle(track) for track in fake_tracks]

energy_distribution_reconstructed = [get_track_energy(track) for track in reconstructed_tracks]
energy_distribution_generated = [get_track_energy(track) for track in generated_tracks]
energy_distribution_matched = [get_track_energy(track) for track in matched_tracks]
energy_distribution_fake = [get_track_energy(track) for track in fake_tracks]


# energy dependent part
n_matched_energy, bins_matched_energy, patches_matched_energy = plt.hist(energy_distribution_matched,
                                                                         label="matched",
                                                                         alpha=0.4,
                                                                         bins=25,
                                                                         range=(1000, 12000))
n_generated_energy, bins_generated_energy, patches_generated_energy = plt.hist(energy_distribution_generated,
                                                                               label="gen",
                                                                               alpha=0.4,
                                                                               bins=25,
                                                                               range=(1000, 12000))
n_fake_energy, bins_fake_energy, patches_fake_energy = plt.hist(energy_distribution_fake,
                                                                label="fake",
                                                                alpha=0.4,
                                                                bins=25,
                                                                range=(1000, 12000))

n_reco_energy, bins_reco_energy, patches_reco_energy = plt.hist(energy_distribution_reconstructed,
                                                                label="reco",
                                                                alpha=0.4,
                                                                bins=25,
                                                                range=(1000, 12000))
plt.close()

# efficiency
efficiency_e = [min(k / n, 1) for k, n in zip(n_matched_energy, n_generated_energy)]
eff_standard_deviation_e_lower = []
for k, n in zip(n_matched_energy, n_generated_energy):
    eff_standard_deviation_e_lower.append(
        get_standard_deviation(k, n) if k / n - get_standard_deviation(k, n) >= 0 else k / n)
eff_standard_deviation_e_upper = []
for k, n in zip(n_matched_energy, n_generated_energy):
    eff_standard_deviation_e_upper.append(
        get_standard_deviation(k, n) if min(k / n, 1) + get_standard_deviation(k, n) <= 1 else 1 - min(k / n, 1))

# fake rate
fake_rate_e = [k / n for k, n in zip(n_fake_energy, n_reco_energy)]
f_r_standard_deviation_lower = []
for k, n in zip(n_fake_energy, n_reco_energy):
    f_r_standard_deviation_lower.append(
        get_standard_deviation(k, n) if k / n - get_standard_deviation(k, n) >= 0 else k / n)

f_r_standard_deviation_upper = []
for k, n in zip(n_fake_energy, n_reco_energy):
    f_r_standard_deviation_upper.append(
        get_standard_deviation(k, n) if min(k / n, 1) + get_standard_deviation(k, n) <= 1 else 1 - min(k / n, 1))


fig, ax = plt.subplots(dpi=300, figsize=(8, 6))
ax2 = ax.twinx()
ax.set_title(rf"$\xi$={xi}, {len(generated_tracks)} generated tracks",
             loc="left",
             fontsize=14)
ax.errorbar(bins_matched_energy[:-1],
            np.array(efficiency_e) * 100,
            yerr=(np.array(eff_standard_deviation_e_lower) * 100, np.array(eff_standard_deviation_e_upper) * 100),
            linestyle="",
            marker="o",
            color="blue",
            markersize=6,
            ecolor="k",
            elinewidth=1.5,
            label="efficiency")
ax.errorbar(bins_matched_energy[:-1],
            np.array(fake_rate_e) * 100,
            yerr=(np.array(f_r_standard_deviation_lower) * 100, np.array(f_r_standard_deviation_upper) * 100),
            linestyle="",
            marker="o",
            color="red",
            markersize=6,
            ecolor="k",
            elinewidth=1.5,
            label="fake rate")
ax2.hist(energy_distribution_generated,
         color="grey",
         label="gen tracks",
         alpha=0.8,
         bins=25,
         histtype="step",
         linewidth=1.5,
         align="left",
         range=(1000, 11000))
ax.tick_params(labelsize=14)
ax.set_xlabel("True particle energy [MeV]", fontsize=14)
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_ylabel("[%]", fontsize=14)
ax2.tick_params(labelsize=14)
ax2.legend(loc=[0.75, 0.6], fontsize=12)
ax.legend(loc=[0.75, 0.4], fontsize=12)
ax2.set_ylabel("counts / 400 MeV", fontsize=14)
fig.savefig(f"{folder}/efficiency_{xi}.pdf")
ax.grid()
