import sys
import csv
from typing import List, Any, Union

import matplotlib.pyplot as plt
import numpy as np

from src.triplet import Triplet
from src.doublet import Doublet
from src.track import Track

folder = sys.argv[1]
print(folder.split("_"))
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
print(n_fake_energy)
print(n_reco_energy)
print(fake_rate_e)
f_r_standard_deviation_lower = []
for k, n in zip(n_fake_energy, n_reco_energy):
    f_r_standard_deviation_lower.append(
        get_standard_deviation(k, n) if k / n - get_standard_deviation(k, n) >= 0 else k / n)

f_r_standard_deviation_upper = []
for k, n in zip(n_fake_energy, n_reco_energy):
    f_r_standard_deviation_upper.append(
        get_standard_deviation(k, n) if min(k / n, 1) + get_standard_deviation(k, n) <= 1 else 1 - min(k / n, 1))


plt.figure(dpi=200, figsize=(8, 6))
plt.title(f"Efficiency and fake rate for xi = {xi}, {len(generated_tracks)} generated tracks", loc="left")
plt.errorbar(bins_matched_energy[:-1],
             efficiency_e,
             yerr=(eff_standard_deviation_e_lower, eff_standard_deviation_e_upper),
             linestyle="",
             marker="^",
             color="blue",
             markersize=8,
             ecolor="k",
             elinewidth=2.5,
             label="efficiency")
plt.errorbar(bins_matched_energy[:-1],
             fake_rate_e,
             yerr=(f_r_standard_deviation_lower, f_r_standard_deviation_upper),
             linestyle="",
             marker="^",
             color="red",
             markersize=8,
             ecolor="k",
             elinewidth=2.5,
             label="efficiency")
plt.xticks(fontsize=16)
plt.show()

