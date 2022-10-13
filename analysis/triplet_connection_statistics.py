import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../src")


data = np.load("../test_files_qubo/e0gpc_7.0_0000_sl-c_4/triplet_list.npy", allow_pickle=True)

all_connections = []
correct_connections = []

for triplet in data:
    if int(triplet.triplet_id) % 1000:
        print(int(triplet.triplet_id))
    connections = []
    correct = None
    for partner, connection in zip(triplet.interactions.keys(), triplet.interactions.values()):
        if triplet.is_correct_match and data[partner].is_correct_match:
            correct = connection
            connections.append(connection)
        if connection < 0:
            connections.append(connection)
    connections.sort()
    if correct is not None:
        all_connections.append(len(connections))
        correct_connections.append(1 + connections.index(correct))


plt.figure()
plt.hist(all_connections, bins=30, label="number of connections")
plt.hist(correct_connections, bins=30, label="correct connection position")
plt.legend(loc="best")
plt.show()
