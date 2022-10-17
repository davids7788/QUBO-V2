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
        correct_connections.append(1 + connections.index(correct))
    all_connections.append(len(connections))

print(correct_connections.count(4) / len(correct_connections))
print((all_connections.count(4) +
       all_connections.count(5) +
       all_connections.count(6) +
       all_connections.count(7) +
       all_connections.count(8)) / len(all_connections))
fig, ax = plt.subplots(figsize=(8,6), dpi=200)
ax.hist(correct_connections,
        bins=4,
        range=(1, 5),
        label=r"phase-0, $\xi=7.0$",
        histtype="bar",
        color="teal",
        align="left",
        rwidth=0.5)
ax.set_xticks([1, 2, 3, 4])
ax.set_xlabel("Position of correct match in interaction list", fontsize=16)
ax.set_ylabel("Number of triplets", fontsize=16)
ax.tick_params(axis='both', labelsize=16)

ax.legend(loc="best", fontsize=16)
plt.savefig("position_of_correct_match.pdf", bbox_inches='tight')
plt.savefig("position_of_correct_match.png", bbox_inches='tight')

fig2, ax = plt.subplots(figsize=(8,6), dpi=200)
ax.hist(all_connections,
        bins=8,
        range=(1, 9),
        label=r"phase-0, $\xi=7.0$",
        histtype="bar",
        color="firebrick",
        align="left",
        rwidth=0.5)
ax.set_xlabel("Number of partner triplets to form a quadruplet", fontsize=16)
ax.set_ylabel("Number of triplets", fontsize=16)
ax.set_xticks([1, 2, 3, 4, 5, 6, 7, 8])
ax.tick_params(axis='both', labelsize=16)

ax.legend(loc="best", fontsize=16)
plt.savefig("number_connections.pdf", bbox_inches='tight')
plt.savefig("number_connections.png", bbox_inches='tight')
