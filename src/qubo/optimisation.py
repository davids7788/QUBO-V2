import numpy as np


def bit_flip_optimisation(triplets,
                          solution_candidate,
                          triplet_ordering,
                          reverse=True):

    """Looping over a list of triplet objects and a corresponding binary solution vector to compute if the energy
    value would improve if the binary state of a triplet (keep <-> discard) should be changed. If the energy
    decreases, the triplet state is flipped.
    :param
        triplets: list of triplet objects
        t_mapping: mapping of names to positions
        solution_candidate binary vector representing triplet states
        reverse: False if provided sorting order, else reversed
    """
    if reverse:
        triplet_ordering.reverse()
    energy_change_total = 0
    for i, triplet in enumerate([triplets[index] for index in triplet_ordering]):
        # energy change if this particular bit is flipped
        energy_change = 0

        # checking linear term
        if solution_candidate[triplet_ordering[i]] == 0:
            energy_change += triplet.quality
        else:
            energy_change -= triplet.quality

        # Checking interactions with other triplets
        for interaction in triplet.interactions.keys():
            if solution_candidate[triplet_ordering[i]] == 0 and solution_candidate[interaction] == 0:
                pass
            elif solution_candidate[triplet_ordering[i]] == 0 and solution_candidate[interaction] == 1:
                energy_change += triplet.interactions[interaction]
            elif solution_candidate[triplet_ordering[i]] == 1 and solution_candidate[interaction] == 0:
                pass
            else:
                energy_change -= triplet.interactions[interaction]

        # flip if overall energy change is negative
        if energy_change < 0:
            solution_candidate[triplet_ordering[i]] = 1 - solution_candidate[triplet_ordering[i]]
            energy_change_total += energy_change

    return solution_candidate, energy_change_total


def make_impact_list(triplet_list,
                     t_mapping,
                     solution_candidate,
                     local=False):
    """Creates an impact list based on how much influence on the energy a bit flip has
    :param triplet_list: list of triplet objects
    :param t_mapping: triplet mapping
    :param solution_candidate: binary solution vector, representing triplet state
    :param local: if True, then triplets outside of the given triplet subset are not considered
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
    """Creates an preferred connections list and returns a list of indices for the triplets
    :param triplet_list: list of triplet objects
    :param t_mapping: triplet mapping
    :return:
        list of indices ordered from lowest to highest impact of triplets in triplet list
    """
    def sort_by_connection(connections_map_entry):
        return connections_map_entry[0]

    connections_map = []   # connection_value, triplet_1, triplet_2
    for triplet in triplet_list:
        for connection, value in zip(triplet.interactions.keys(), triplet.interactions.values()):
            if value < 0:
                connections_map.append([value, triplet.triplet_id, connection])

    connections_map.sort(key=sort_by_connection)

    triplet_used = set()
    triplet_ordering = []

    for entry in connections_map:
        for i in [1, 2]:
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
