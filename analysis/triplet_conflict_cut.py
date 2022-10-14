import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../src")


data = np.load("../qubo_files/e0gpc_7.0_0000_sl-c_4/triplet_list.npy", allow_pickle=True)

number_of_conflicts_matched = []
number_of_conflicts_fake = []

for triplet in data:
    if triplet.is_correct_match:
        number_of_conflicts_matched.append(list(triplet.interactions.values()).count(1))
    else:
        number_of_conflicts_fake.append(list(triplet.interactions.values()).count(1))

fig2, ax = plt.subplots()
ax.hist([number_of_conflicts_matched, number_of_conflicts_fake],
        range=(0, 250),
        bins=25,
        label=["matched triplets", "fake triplets"],
        color=["teal", "firebrick"],
        align="left")
ax.text(150, 8e4, r"phase-0, $\xi=7.0$", fontsize=16)
ax.set_xlabel("Number of conflicts", fontsize=16)
ax.set_ylabel("Number of triplets", fontsize=16)
ax.tick_params(axis='both', labelsize=16)
plt.legend(loc="best", fontsize=16)
plt.savefig("number_conflicts.pdf", bbox_inches='tight')
