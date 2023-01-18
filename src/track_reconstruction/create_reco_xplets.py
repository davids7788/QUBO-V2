import numpy as np

from pattern.x_plet import Xplet


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
    print("Creating reco Xplets and fitting resulting track...\n")

    def triplet_start(t):
        """Returns z value of first triplet hit
        """
        return t.doublet_1.hit_1_position[2]

    triplet_start_value = set()
    for triplet in triplets:
        triplet_start_value.add(triplet.doublet_1.hit_1_position[2])
    triplets.sort(key=triplet_start)
    xplet_start = min(list(triplet_start_value))

    reco_x_plets = []
    triplets_used = set()
    for i, t1 in enumerate(triplets):
        if t1.doublet_1.hit_1_position[2] == xplet_start:
            t_list = [t1]
            for t2 in triplets[i + 1:]:
                if t_list[-1].doublet_2 == t2.doublet_1:
                    t_list.append(t2)

            reco_pattern = Xplet()
            for triplet in t_list:
                triplets_used.add(triplet.triplet_id)
                reco_pattern.add_triplet(triplet)
            if fit == "chi squared lin track":
                reco_pattern.fit_lin_track()
            reco_x_plets.append(reco_pattern)

    # tracks starting from the second layer
    for t3 in triplets:
        if t3.triplet_id not in triplets_used:
            reco_pattern = Xplet()
            reco_pattern.add_triplet(t3)
            if fit == "chi squared lin track":
                reco_pattern.fit_lin_track()
            reco_x_plets.append(reco_pattern)

    hit_to_xplet_map = {}
    for i, xplet in enumerate(reco_x_plets):
        for hit_id in xplet.hit_ids.values():
            if hit_id not in hit_to_xplet_map.values():
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
        for track in reco_x_plets:
            if not select_triplets_with_max_overlap(track, 4 - i):
                selected_tracks.append(track)
        reco_ambiguity_solved = selected_tracks

    print(f"Number of combinatorial reco Xplets: {len(reco_x_plets)}\n")
    print(f"Number of reco Xplets after ambiguity solving: {len(reco_ambiguity_solved)}\n")

    np.save(f"{save_folder}/reco_xplet_list", reco_ambiguity_solved)
