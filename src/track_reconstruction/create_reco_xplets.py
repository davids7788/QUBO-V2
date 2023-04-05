import numpy as np

from pattern.x_plet import Xplet
from pattern.triplet import Triplet


def reco_xplets_simplified_LUXE(triplets,
                                save_folder,
                                fit="chi squared lin track"):
    """Creates xplets from the kept triplets.
    :param fit: 
            "chi squared lin track: assuming a linear track and fitting  a chi squared assuming f(z) = x_i, f(z) = y_j
    :param triplets: list of kept triplets
    :param save_folder: folder to store reco Xplet list 
    """
    print("\n-----------------------------------\n")
    print("Creating reco Xplets and fitting resulting tracks...\n")

    def triplet_start(t_obj: Triplet):
        """Returns z value of first triplet hit
        """
        return t_obj.doublet_1.hit_1_position[2]

    triplet_start_value = set()
    for triplet in triplets:
        triplet_start_value.add(triplet.doublet_1.hit_1_position[2])
    triplets.sort(key=triplet_start)

    # check for
    sorted_z = sorted(triplet_start_value)
    if len(triplet_start_value) == 2:
        t_possible_start = [0]
    else:
        t_possible_start = [0, 1]

    triplets_ordered_by_start_value = {}
    for v in sorted_z:
        triplets_ordered_by_start_value.update({sorted_z.index(v): []})
    for t in triplets:
        triplets_ordered_by_start_value[sorted_z.index(triplet_start(t))].append(t)

    reco_x_plets_candidates = []

    for i in range(len(sorted_z)):
        # start of xplets
        if i in t_possible_start:
            reco_x_plets_candidates_update = []
            for t in triplets_ordered_by_start_value[i]:
                reco_x_plets_candidates_update.append([t])
            reco_x_plets_candidates += reco_x_plets_candidates_update

        else:
            reco_x_plets_candidates_update = []
            for t_root in reco_x_plets_candidates:
                found_next = False
                for t_next in triplets_ordered_by_start_value[i]:
                    if t_next.doublet_1 == t_root[-1].doublet_2:
                        growing_xplet = [t_r for t_r in t_root]
                        growing_xplet.append(t_next)
                        reco_x_plets_candidates_update.append(growing_xplet)
                        found_next = True
                if not found_next:
                    reco_x_plets_candidates_update.append(t_root)
            reco_x_plets_candidates = reco_x_plets_candidates_update

    reco_xplets_subsets_removed = []
    for xplet_candidate in reco_x_plets_candidates:
        is_subset = False
        for other_xplet_candidate in reco_x_plets_candidates:
            if xplet_candidate != other_xplet_candidate:
                if set(xplet_candidate).issubset(set(other_xplet_candidate)):
                    is_subset = True
        if not is_subset:
            reco_xplets_subsets_removed.append(xplet_candidate)

    reco_x_plets_candidates = []
    for xplet_candidate in reco_xplets_subsets_removed:
        reco_pattern = Xplet()
        for triplet in xplet_candidate:
            reco_pattern.add_triplet(triplet)
        if fit == "chi squared lin track":
            reco_pattern.fit_lin_track()
        reco_x_plets_candidates.append(reco_pattern)

    print(f"Number of combinatorial reco Xplets: {len(reco_x_plets_candidates)}\n"
          f"Saving combinatorial Xplets to file {save_folder}/reco_xplet_list")
    np.save(f"{save_folder}/reco_xplet_list", reco_x_plets_candidates)

    hit_to_xplet_map = {}
    for xplet in reco_x_plets_candidates:
        for hit_id in xplet.hit_ids.values():
            if hit_id not in hit_to_xplet_map.keys():
                hit_to_xplet_map.update({hit_id: 1})
            else:
                ambiguity_number = hit_to_xplet_map[hit_id]
                hit_to_xplet_map.update({hit_id: ambiguity_number + 1})

    list_shared = list(hit_to_xplet_map.values())
    num_shared_hits = 0
    for value in list_shared:
        if value > 1:
            num_shared_hits += 1

    print(f"{num_shared_hits} shared hits found... \n")

    print("Solving ambiguities...\n")

    def select_triplets_with_max_overlap(reco_track, max_overlaps):
        """Returns the number of overlaps of a track with other tracks.
        :param reco_track: Xplet object
        :param max_overlaps: max allowed overlaps
        :return:
            True if more than max overlaps, else False"""
        overlaps = 0
        for h_id in reco_track.hit_ids.values():
            if hit_to_xplet_map[h_id] > 1:
                overlaps += 1
        if overlaps > max_overlaps:
            return True
        return False

    reco_ambiguity_solved = None
    for i in range(3):
        selected_tracks = []
        for track in reco_x_plets_candidates:
            if not select_triplets_with_max_overlap(track, 4 - i):
                selected_tracks.append(track)
        reco_ambiguity_solved = selected_tracks

    print(f"Number of reco Xplets after ambiguity solving: {len(reco_ambiguity_solved)}\n"
          f"Saving ambiguity solved Xplets to file {save_folder}/reco_xplet_list_ambiguity_solved.npy")

    np.save(f"{save_folder}/reco_xplet_list_ambiguity_solved", reco_ambiguity_solved)
