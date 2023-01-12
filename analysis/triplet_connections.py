import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../src")


data = np.load("../../QUBO/e0gpc_7.0_0000_sl-c_4/triplet_list.npy", allow_pickle=True)

connection_strength_matched = []
connection_strength_fake = []

quality_matched = []
quality_fake = []
quality_matched_more_0_8 = 0
quality_fake_more_0_8 = 0

matched_more_than_0_92 = 0
fake_more_than_0_92 = 0


for triplet in data:
    if triplet.is_correct_match:
        quality_matched.append(triplet.quality)
        if quality_matched[-1] > 0.8:
            quality_matched_more_0_8 +=1
        for value in triplet.interactions.values():
            if value < 0:
                connection_strength_matched.append(value)
                if value > -0.92:
                    matched_more_than_0_92 += 1
    else:
        quality_fake.append(triplet.quality)
        if quality_fake[-1] > 0.8:
            quality_fake_more_0_8 +=1
        for value in triplet.interactions.values():
            if value < 0:
                connection_strength_fake.append(value)
                if value > -0.92:
                    fake_more_than_0_92 += 1

print(quality_matched_more_0_8/ (len(quality_matched)))
print(quality_fake_more_0_8/ (len(quality_fake)))


fig, ax = plt.subplots(figsize=(8,6), dpi=200)
ax.hist([connection_strength_matched, connection_strength_fake],
        range=(-1, -0.9),
        bins=50,
        label=["matched", "fake"],
        color=["teal", "firebrick"],
        linewidth=2.5,
        histtype="step",)
ax.text(-1.0005, 1e5, r"phase-0, $\xi=e0gpc_7.0$", fontsize=16)
ax.set_xlabel("Connection strength", fontsize=16)
ax.set_ylabel("Number of connections", fontsize=16)
ax.tick_params(axis='both', labelsize=16)
plt.grid()
plt.legend(loc="best", fontsize=16)
plt.savefig("triplet_connections.pdf", bbox_inches='tight')
plt.savefig("triplet_connections.png", bbox_inches='tight')

fig2, ax = plt.subplots(figsize=(8, 6), dpi=200)
ax.hist([quality_matched, quality_fake],
        range=(-1, 1),
        bins=50,
        label=["matched", "fake"],
        color=["teal", "firebrick"],
        linewidth=2.5,
        histtype="step",)
ax.text(-1.0, 25e3, r"phase-0, $\xi=e0gpc_7.0$", fontsize=16)
ax.set_xlabel("Triplet quality", fontsize=16)
ax.set_ylabel("Number of triplets", fontsize=16)
ax.tick_params(axis='both', labelsize=16)
plt.grid()
plt.legend(loc="upper left", fontsize=16)
plt.savefig("triplet_qualities.pdf", bbox_inches='tight')
plt.savefig("triplet_qualities.png", bbox_inches='tight')