import numpy as np
import os
import sys
import matplotlib.pyplot as plt
sys.path.append("../src")

triplets = np.load("C:/Users/david/PycharmProjects/QUBO_advanced/qubo_files/xi_5/"
                   "e0gpc_5.0_0000_sl-example/triplet_list.npy", allow_pickle=True)

t_mapping = {key: value for key, value in zip([t.triplet_id for t in triplets], np.arange(len(triplets)))}


def hamiltonian_energy(binary_vector, triplet_subset=None):
    """Calculates the energy according to a binary vector matching the triplet list
    :param binary_vector: binary solution candidate vector
    :param triplet_subset: subset of triplet, represented by indices
    :return:
        energy value
    """
    hamiltonian_energy = 0
    if triplet_subset is None:
        for i, b1 in enumerate(binary_vector):
            if b1 == 1:
                hamiltonian_energy += triplets[i].quality
            for j in triplets[i].interactions.keys():
                if t_mapping[j] < i:
                    continue
                else:
                    pass
                    if binary_vector[i] == binary_vector[t_mapping[j]] == 1:
                        hamiltonian_energy += triplets[i].interactions[j]
        return hamiltonian_energy

    else:
        for i, b1 in enumerate(binary_vector):
            if b1 == 1:
                hamiltonian_energy += triplets[triplet_subset[i]].quality

            for j in triplets[triplet_subset[i]].interactions.keys():
                if j < triplet_subset[i]:
                    continue
                if j not in triplet_subset:
                    continue
                position = triplet_subset.index(j)
                if binary_vector[i] == binary_vector[position] == 1:
                    hamiltonian_energy += triplets[triplet_subset[i]].interactions[j]
        return hamiltonian_energy


def make_impact_list(triplet_list,
                     t_mapping,
                     solution_candidate,
                     local=False):
    """Creates an impact list based on how much influence on the energy a bit flip has
    :param triplet_list: list of triplet objects
    :param t_mapping: name of triplets consist of <hit_ID>_<hit_ID>_<hit_ID> and is mapped to its position
                      inside the triplet list
    :param solution_candidate: binary vector representing kept and discarded triplets e.g [0,1,..., 1]
    :param local: if True, then triplets outside the given triplet subset are not considered
    :return:
        list of indices ordered from lowest to highest impact of triplets in triplet list
    """
    impact_list_values = []
    t_indices = []
    for triplet in triplet_list:
        t_indices.append(triplet.triplet_id)

    for triplet, t_i in zip(triplet_list, solution_candidate):
        energy_change = 0
        if t_i == 0:
            energy_change += triplet.quality
        else:
            energy_change -= triplet.quality

        for interaction in triplet.interactions.keys():

            if local:
                if interaction not in t_indices:
                    continue
            if t_i == 0 and solution_candidate[t_mapping[interaction]] == 0:
                pass
            elif t_i == 0 and solution_candidate[t_mapping[interaction]] == 1:
                energy_change += triplet.interactions[interaction]
            elif t_i == 1 and solution_candidate[t_mapping[interaction]] == 0:
                pass
            else:
                energy_change -= triplet.interactions[interaction]

        impact_list_values.append(abs(energy_change))
    return list(np.argsort(impact_list_values))
def make_connection_list(triplet_list,
                         t_mapping):
    """Creates a preferred connections list and returns a list of indices for the triplets
    :param triplet_list: list of triplet objects
    :param t_mapping: triplet mapping
    :return:
        list of indices ordered from  highest connection values to lowest connection values
    """
    def sort_by_connection(connections_map_entry):
        return connections_map_entry[0]

    connections_map = []   # connection_value, triplet_1, triplet_2

    for triplet in triplet_list:
        for connection, value in zip(triplet.interactions.keys(), triplet.interactions.values()):
            if value < 0:
                connections_map.append([value, triplet.triplet_id, connection])
        if len(list(triplet.interactions.values())) == 0:
            connections_map.append([0, triplet.triplet_id, None])
        if len(list(triplet.interactions.values())) > 0:
            if min(list(triplet.interactions.values())) > 0:
                connections_map.append([1, triplet.triplet_id, None])

    connections_map.sort(key=sort_by_connection)

    triplet_used = set()
    triplet_ordering = []

    for entry in connections_map:
        for i in [1, 2]:
            if entry[i] is None:
                continue
            if entry[i] not in triplet_used:
                triplet_ordering.append(t_mapping[entry[i]])
                triplet_used.add(entry[i])
                for connection, value in zip(triplet_list[t_mapping[entry[i]]].interactions.keys(),
                                             triplet_list[t_mapping[entry[i]]].interactions.values()):
                    if value < 0:
                        if connection not in triplet_used:
                            triplet_used.add(connection)
                            triplet_ordering.append(t_mapping[connection])

    return triplet_ordering


impact_list = make_impact_list(triplets, t_mapping, np.ones(len(triplets)))
connection_list = make_connection_list(triplets, t_mapping)
reverse_connection_map = {value: key for key, value in zip([i for i in range(len(triplets))], connection_list)}
reverse_impact_map = {value: key for key, value in zip([i for i in range(len(triplets))], impact_list)}


connection_measure_impact = []
connection_measure_connection = []


for position, index in enumerate(impact_list):
    for connection, value in triplets[index].interactions.items():
        if value < 0:
            connection_measure_impact.append(abs(position - reverse_impact_map[t_mapping[connection]]))

for position, index in enumerate(connection_list):
    for connection, value in triplets[index].interactions.items():
        if value < 0:
            connection_measure_connection.append(abs(position - reverse_connection_map[t_mapping[connection]]))


print(sum(i <= 3 for i in connection_measure_connection) / len(connection_measure_connection))
print(sum(i <= 3 for i in connection_measure_impact) / len(connection_measure_impact))
print("-----------")
print(sum(i <= 7 for i in connection_measure_connection) / len(connection_measure_connection))
print(sum(i <= 7 for i in connection_measure_impact) / len(connection_measure_impact))
print("-----------")
print(sum(i <= 12 for i in connection_measure_connection) / len(connection_measure_connection))
print(sum(i <= 12 for i in connection_measure_impact) / len(connection_measure_impact))
print("-----------")
print(sum(i <= 16 for i in connection_measure_connection) / len(connection_measure_connection))
print(sum(i <= 16 for i in connection_measure_impact) / len(connection_measure_impact))
print("-----------")
print(sum(i <= 100 for i in connection_measure_connection) / len(connection_measure_connection))
print(sum(i <= 100 for i in connection_measure_impact) / len(connection_measure_impact))

plt.figure(figsize=(8, 6), dpi=400)
plt.hist(connection_measure_impact,
         bins=100,
         histtype="step",
         label="impact list",
         color="royalblue",
         linewidth=2)
plt.hist(connection_measure_connection,
         bins=100,
         histtype="step",
         label="connection list",
         color="black",
         linewidth=2)
plt.legend(loc="best", fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.ylabel("Number of connections", loc="top", fontsize=16)
plt.xlabel("distance of entries in sorted list        ", loc="right", fontsize=16)
plt.title("xi=7, 1 representative BX", loc="left", fontsize=18)
plt.yscale("log")
plt.show()

