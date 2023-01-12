import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../src")


data = np.load("../../QUBO/e0gpc_7.0_0000_sl-c_4/triplet_list.npy", allow_pickle=True)

number_of_conflicts_matched = []
number_of_conflicts_fake = []
quality_conflicts_more_than_150 = 0
fake_conflicts_more_than_150 = 0

for triplet in data:
    if triplet.is_correct_match:
        number_of_conflicts_matched.append(list(triplet.interactions.values()).count(1))
        if number_of_conflicts_matched[-1] > 135:
            quality_conflicts_more_than_150 += 1
    else:
        number_of_conflicts_fake.append(list(triplet.interactions.values()).count(1))
        if number_of_conflicts_fake[-1] > 135:
            fake_conflicts_more_than_150 += 1

# print(quality_conflicts_more_than_150 / (quality_conflicts_more_than_150 + len(number_of_conflicts_matched)))
# print(fake_conflicts_more_than_150 / (fake_conflicts_more_than_150 + len(number_of_conflicts_fake)))

fig2, ax = plt.subplots(figsize=(8,6), dpi=200)
ax.hist([number_of_conflicts_matched, number_of_conflicts_fake],
        range=(0, 250),
        bins=250,
        histtype="step",
        linewidth=2.5,
        # cumulative=True,
        # density=True,
        label=["matched triplets", "fake triplets"],
        color=["teal", "firebrick"])
ax.text(150, 8e3, r"phase-0, $\xi=e0gpc_7.0$", fontsize=16)
ax.set_xlabel("Number of conflicts", fontsize=16)
ax.set_ylabel("Number of triplets", fontsize=16)
ax.tick_params(axis='both', labelsize=16)
plt.legend(loc="best", fontsize=16)
plt.grid()
plt.savefig("number_conflicts.pdf", bbox_inches='tight')
plt.savefig("number_conflicts.png", bbox_inches='tight')
