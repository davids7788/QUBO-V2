import sys

sys.path.append("../src")

import matplotlib.pyplot as plt
import numpy as np

folder = sys.argv[1]
gen = "/".join(folder.split("/")[0:-1]).split("-")[0]
xi = folder.split("e0gpc_")[1].split("_")[0]

generated_xplets = np.load(gen + "_gen_xplet_list.npy", allow_pickle=True)
reconstructed_xplets = np.load(f"{folder}/reco_xplet_list.npy", allow_pickle=True)

def get_efficiency(reco, gen, num_hits_from_same_particle):
    remove_shared = set()
    shared_hits = []
    reco = np.load(reco, allow_pickle=True)
    gen = np.load(gen, allow_pickle=True)
    fake = 0
    matched = set()
    for xplet in reco:
        matched_xplet = False
        intersection = set(xplet.hit_ids.values()).intersection(remove_shared)
        if remove_shared == set():
            pass
        else:
            if intersection != set():
                shared_hits.append(intersection)
                continue
        for hit in xplet.hit_ids.values():
            remove_shared.add(hit)
        for particle in set(xplet.particle_ids.values()):
            if list(xplet.particle_ids.values()).count(particle) >= num_hits_from_same_particle:
                matched_xplet = True
        if matched_xplet:
            matched.add(particle)
        else:
            fake += 1
    # print(shared_hits)
    return len(matched) / len(gen), fake / len(reco)

matched_xplets, fake_xplets = get_efficiency(f"{folder}/reco_xplet_list.npy", gen + "_gen_xplet_list.npy", 3)
print(matched_xplets, fake_xplets)


def get_xplet_x(xplet_for_x):
    """Returns the x-value on the first layer of the xplet.
     :return
        x of xplet
    """
    return xplet_for_x.coordinates[0][0]


def get_standard_deviation(k, n):
    if k == 0:
        return 0
    if k > n:
        k = n
    return np.sqrt(((k + 1) * (k + 2)) / ((n + 2) * (n + 3)) - ((k + 1)**2 / (n + 2)**2))


x_distribution_reconstructed = [get_xplet_x(xplet) for xplet in reconstructed_xplets]
x_distribution_generated = [get_xplet_x(xplet) for xplet in generated_xplets]
x_distribution_matched = [get_xplet_x(xplet) for xplet in matched_xplets]
x_distribution_fake = [get_xplet_x(xplet) for xplet in fake_xplets]

print(np.around(len(matched_xplets) / len(generated_xplets), 3))

# x dependent part
n_matched_x, bins_matched_x, patches_matched_x = plt.hist(x_distribution_matched,
                                                          label="matched",
                                                          alpha=0.4,
                                                          bins=25,
                                                          range=(0.02, 0.58))
n_generated_x, bins_generated_x, patches_generated_x = plt.hist(x_distribution_generated,
                                                                label="gen",
                                                                alpha=0.4,
                                                                bins=25,
                                                                range=(0.02, 0.58))
n_fake_x, bins_fake_x, patches_fake_x = plt.hist(x_distribution_fake,
                                                 label="fake",
                                                 alpha=0.4,
                                                 bins=25,
                                                 range=(0.02, 0.58))
n_reco_x, bins_reco_x, patches_reco_x = plt.hist(x_distribution_reconstructed,
                                                 label="reco",
                                                 alpha=0.4,
                                                 bins=25,
                                                 range=(0.02, 0.58))
plt.close()

# efficiency
efficiency_e = [min(k / n, 1) for k, n in zip(n_matched_x, n_generated_x)]
eff_standard_deviation_e_lower = []
for k, n in zip(n_matched_x, n_generated_x):
    eff_standard_deviation_e_lower.append(
        get_standard_deviation(k, n) if k / n - get_standard_deviation(k, n) >= 0 else k / n)
eff_standard_deviation_e_upper = []
for k, n in zip(n_matched_x, n_generated_x):
    eff_standard_deviation_e_upper.append(
        get_standard_deviation(k, n) if min(k / n, 1) + get_standard_deviation(k, n) <= 1 else 1 - min(k / n, 1))

# fake rate
fake_rate_e = [k / n for k, n in zip(n_fake_x, n_reco_x)]
f_r_standard_deviation_lower = []
for k, n in zip(n_fake_x, n_reco_x):
    f_r_standard_deviation_lower.append(
        get_standard_deviation(k, n) if k / n - get_standard_deviation(k, n) >= 0 else k / n)

f_r_standard_deviation_upper = []
for k, n in zip(n_fake_x, n_reco_x):
    f_r_standard_deviation_upper.append(
        get_standard_deviation(k, n) if min(k / n, 1) + get_standard_deviation(k, n) <= 1 else 1 - min(k / n, 1))


fig, ax = plt.subplots(dpi=300, figsize=(8, 6))
ax2 = ax.twinx()
ax.set_title(rf"$\xi$={xi}, {len(generated_xplets)} generated xplets",
             loc="left",
             fontsize=14)
ax.errorbar(bins_matched_x[:-1],
            np.array(efficiency_e) * 100,
            yerr=(np.array(eff_standard_deviation_e_lower) * 100, np.array(eff_standard_deviation_e_upper) * 100),
            linestyle="",
            marker="o",
            color="blue",
            markersize=6,
            ecolor="k",
            elinewidth=1.5,
            label="efficiency")
ax.errorbar(bins_matched_x[:-1],
            np.array(fake_rate_e) * 100,
            yerr=(np.array(f_r_standard_deviation_lower) * 100, np.array(f_r_standard_deviation_upper) * 100),
            linestyle="",
            marker="o",
            color="red",
            markersize=6,
            ecolor="k",
            elinewidth=1.5,
            label="fake rate")
ax2.hist(x_distribution_generated,
         color="grey",
         label="gen xplets",
         alpha=0.8,
         bins=25,
         histtype="step",
         linewidth=1.5,
         align="left",
         range=(0.02, 0.58))
ax.tick_params(labelsize=14)
ax.set_xlabel("xplet position on first detector layer [m]", fontsize=14)
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_ylabel("[%]", fontsize=14)
ax2.tick_params(labelsize=14)
ax2.legend(loc=[0.75, 0.6], fontsize=12)
ax.legend(loc="best", fontsize=12)
ax2.set_ylabel("counts", fontsize=14)
fig.savefig(f"{folder}/efficiency_{xi}.pdf")
ax.grid()
