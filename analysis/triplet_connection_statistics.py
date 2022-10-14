import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../src")


data = np.load("../qubo_files/e0gpc_7.0_0000_sl-c_4/triplet_list.npy", allow_pickle=True)

all_connections = []
correct_connections = []

for triplet in data:
    connections = []
    correct = None
    for partner, connection in zip(triplet.interactions.keys(), triplet.interactions.values()):
        if connection < 0:
            if triplet.is_correct_match and data[partner].is_correct_match:
                correct = connection
            connections.append(connection)
    connections.sort()
    if correct is not None:
        all_connections.append(len(connections))
        correct_connections.append(1 + connections.index(correct))

fig, ax = plt.subplots()
ax.hist(correct_connections,
        bins=4,
        range=(1, 5),
        label=r"phase-0, $\xi=7.0$",
        align="left",
        color="teal",
        rwidth=0.5)
ax.set_xlabel("Position of correct match in interaction list", fontsize=16)
ax.set_ylabel("Number of triplets", fontsize=16)
ax.tick_params(axis='both', labelsize=16)

ax.legend(loc="best", fontsize=16)
plt.savefig("position_of_correct_match.pdf", bbox_inches='tight')

fig2, ax = plt.subplots()
ax.hist(all_connections,
        bins=8,
        range=(1, 9),
        label=r"phase-0, $\xi=7.0$",
        align="left",
        color="darkkhaki",
        rwidth=0.5)
ax.set_xlabel("Number of partner triplets to form a quadruplet", fontsize=16)
ax.set_ylabel("Number of triplets", fontsize=16)
ax.set_xticks([1, 2, 3, 4, 5, 6, 7, 8])
ax.tick_params(axis='both', labelsize=16)

ax.legend(loc="best", fontsize=16)
plt.savefig("number_connections.pdf", bbox_inches='tight')
