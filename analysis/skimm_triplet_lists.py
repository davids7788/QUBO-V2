import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../src")
for c in range(5, 6):
    for j in range(10):
        data = np.load(f"/nfs/dust/luxe/user/spatarod/towards_paper/e-laser/phase-0/gpc/7.0/smeared/e0gpc_7.0_000{j}_sl-c_{c}/triplet_list.npy", allow_pickle=True)
        print(f"file e0gpc_7.0_000{j}_sl-c_{c}")

        # remove triplet if > 135 conflicts
        # remove triplet if quality > 0.8

        # remove connection if > -0.92
        # remove position >3 in interaction list

        for i, triplet in enumerate(data):
            # if i % 100000 == 0:
            #     print(f"{i} of {len(data)}")
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

        np.save(f"/nfs/dust/luxe/user/spatarod/towards_paper/e-laser/phase-0/gpc/7.0/smeared/e0gpc_7.0_000{j}_sl-c{c}_v2/triplet_list.npy", data)