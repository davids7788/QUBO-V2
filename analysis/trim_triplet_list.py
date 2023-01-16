import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../src")

data = np.load("../../QUBO/e0gpc_7.0_0000_sl-c_6/253254035_bit_flip/qubo_log.npy", allow_pickle=True)[()]
data2 = np.load("../../QUBO/e0gpc_7.0_0000_sl-c_4/968446171_bit_flip/qubo_log.npy", allow_pickle=True)[()]
print(data["truth minimum energy"], data["computed minimum energy"])
print(data2["truth minimum energy"], data2["computed minimum energy"])


correct = 0
false = 0
for computed, truth in zip(data["truth solution vector"], data["computed solution vector"]):
    if truth == computed:
        correct += 1
    else:
        false += 1
print(correct, false)
exit()


data = np.load("../../QUBO/e0gpc_7.0_0000_sl-c_4/triplet_list.npy", allow_pickle=True)

print(len(data))

# remove triplet if > 135 conflicts
# remove triplet if quality > 0.8

# remove connection if > -0.92
# remove position >3 in interaction list


t_removed_list = set()

for i, triplet in enumerate(data):
    if i % 10000 == 0:
        print(f"{i} of {len(data)}")
    if triplet.quality > 0.8:
        triplet.quality = 1
        for key, value in zip(triplet.interactions.keys(), triplet.interactions.values()):
            triplet.interactions.update({key: 1})
            data[key].interactions.update({triplet.triplet_id: 1})

        continue
    if list(triplet.interactions.values()).count(1) > 135:
        for key, value in zip(triplet.interactions.keys(), triplet.interactions.values()):
            triplet.interactions.update({key: 1})
            data[key].interactions.update({triplet.triplet_id: 1})
        triplet.quality = 1
        continue
    sorted_interactions = list(triplet.interactions.values())
    for key, value in zip(triplet.interactions.keys(), triplet.interactions.values()):
        if value > -0.92:
            triplet.interactions.update({key: 1})
            data[key].interactions.update({triplet.triplet_id: 1})
        elif len(sorted_interactions) > 2:
            if sorted_interactions[2] < 0:
                if value > sorted_interactions[2]:
                    triplet.interactions.update({key: 1})
                    data[key].interactions.update({triplet.triplet_id: 1})

print(len(list(t_removed_list)))

np.save("../../QUBO/e0gpc_7.0_0000_sl-c_6/triplet_list.npy", data)