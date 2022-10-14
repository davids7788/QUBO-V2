import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../src")


data = np.load("../qubo_files/e0gpc_7.0_0000_sl-c_4/triplet_list.npy", allow_pickle=True)

connection_strength_matched = []
connection_strength_fake = []

quality_matched = []
quality_fake = []

for triplet in data:
    if triplet.is_correct_match:
        quality_matched.append(triplet.quality)
        for value in triplet.interactions.values():
            if value < 0:
                connection_strength_matched.append(value)
    else:
        quality_fake.append(triplet.quality)
        for value in triplet.interactions.values():
            if value < 0:
                connection_strength_fake.append(value)

fig, ax = plt.subplots()
ax.hist([connection_strength_matched, connection_strength_fake],
        range=(-1, -0.9),
        bins=20,
        label=["matched", "fake"],
        color=["teal", "firebrick"],
        align="left")
ax.text(-1.005, 34e4, r"phase-0, $\xi=7.0$", fontsize=16)
ax.set_xlabel("Connection strength", fontsize=16)
ax.set_ylabel("Number of triplets connections", fontsize=16)
ax.tick_params(axis='both', labelsize=16)
plt.legend(loc="best", fontsize=16)
plt.savefig("triplet_connections.pdf", bbox_inches='tight')

fig2, ax = plt.subplots()
ax.hist([quality_matched, quality_fake],
        range=(-1, 1),
        bins=20,
        label=["matched", "fake"],
        color=["teal", "firebrick"],
        align="left")
ax.text(-1.0, 6e4, r"phase-0, $\xi=7.0$", fontsize=16)
ax.set_xlabel("Triplet quality", fontsize=16)
ax.set_ylabel("Number of triplets", fontsize=16)
ax.tick_params(axis='both', labelsize=16)
plt.legend(loc="best", fontsize=16)
plt.savefig("triplet_qualities.pdf", bbox_inches='tight')