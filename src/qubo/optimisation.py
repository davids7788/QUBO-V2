import numpy as np


def bit_flip_optimisation(triplets,
                          solution_candidate,
                          t_mapping,
                          triplet_ordering,
                          reverse=True):

    """Looping over a list of triplet objects and a corresponding binary solution vector to compute if the energy
    value would improve if the binary state of a triplet (keep <-> discard) should be changed. If the energy
    decreases, the triplet state is flipped.
    :param triplets: list of triplet objects
    :param solution_candidate: binary vector representing kept and discarded triplets e.g [0,1,..., 1]
    :param t_mapping: name of triplets consist of <hit_ID>_<hit_ID>_<hit_ID> and is mapped to its position
                      inside the triplet list
    :param triplet_ordering: order in which the triplets are sorted, e.g by impact, connectivity,...
    :param reverse: False if provided sorting order, else reversed
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
            if solution_candidate[triplet_ordering[i]] == 0 and solution_candidate[t_mapping[interaction]] == 0:
                pass
            elif solution_candidate[triplet_ordering[i]] == 0 and solution_candidate[t_mapping[interaction]] == 1:
                energy_change += triplet.interactions[interaction]
            elif solution_candidate[triplet_ordering[i]] == 1 and solution_candidate[t_mapping[interaction]] == 0:
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
        list of indices ordered from the highest connection value to the lowest connection value
    """
    # function to pass to .sorted()
    def sort_by_connection(connections_map_entry):
        return connections_map_entry[0]

    connections_map = []   # connection_value, triplet_1, triplet_2
                           # -0.998          , t1       , t2
                           # -0.993          , t2       , t234
                           # ...
    # To not run into value errors and making sure the whole list is filled, one has to take care of certain cases
    for triplet in triplet_list:
        # there are connections -> easy
        for connection, value in zip(triplet.interactions.keys(), triplet.interactions.values()):
            if value < 0:
                connections_map.append([value, triplet.triplet_id, connection])
        # no connection, no conflict, "single triplet"
        if len(list(triplet.interactions.values())) == 0:
            connections_map.append([0, triplet.triplet_id, None])
        # no connection, but conflicts (maybe remove them before!?)
        if len(list(triplet.interactions.values())) > 0:
            if min(list(triplet.interactions.values())) > 0:
                connections_map.append([1, triplet.triplet_id, None])

    connections_map.sort(key=sort_by_connection)

    # check if something is in a set is faster, also jsut using a string instead of an object
    triplet_used = set()
    triplet_ordering = []

    for entry in connections_map:
        # check if there are connections, otherwise skip
        for i in [1, 2]:
            if entry[i] is None:
                continue
            # don"t put the same connection two times in the list
            if entry[i] not in triplet_used:
                triplet_ordering.append(t_mapping[entry[i]])
                triplet_used.add(entry[i])
                # add connection, triplet 1 and triplet 2
                for connection, value in zip(triplet_list[t_mapping[entry[i]]].interactions.keys(),
                                             triplet_list[t_mapping[entry[i]]].interactions.values()):
                    if value < 0:
                        if connection not in triplet_used:
                            triplet_used.add(connection)
                            triplet_ordering.append(t_mapping[connection])

    return triplet_ordering


def make_paired_list(triplet_list,
                     t_mapping):
    """Creates a preferred connections list and returns a list of indices for the triplets
    :param triplet_list: list of triplet objects
    :param t_mapping: triplet mapping
    :return:
        list of connected pairs
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

    triplet_ordering = []

    for entry in connections_map:
        triplet_ordering.append(t_mapping[entry[1]])
        if entry[2]:
            triplet_ordering.append(t_mapping[entry[2]])

    return triplet_ordering


def impact_without_conflicts(triplet_list,
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
        list of indices ordered from lowest to highest impact of triplets in triplet list, but no
        conflicts are taken into account!
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

        for interaction, value in triplet.interactions.items():
            if value > 0:
                continue
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
