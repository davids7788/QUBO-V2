def bit_flip_optimisation(triplets,
                          solution_candidate,
                          triplet_ordering,
                          reverse=True):

    """Looping over a list of triplet objects and a corresponding binary solution vector to compute if the energy
    value would improve if the binary state of a triplet (keep <-> discard) should be changed. If the energy
    decreases, the triplet state is flipped.
    :param triplets: list of triplet objects
    :param solution_candidate binary vector representing triplet states
    :param triplet_ordering: index list of triplets
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
