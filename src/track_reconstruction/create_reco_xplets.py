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
    for i, t1 in enumerate(triplets):
        if t1.doublet_1.hit_1_position[2] == xplet_start:
            t_list = [t1]
            for t2 in triplets[i + 1:]:
                if t_list[-1].doublet_2 == t2.doublet_1:
                    t_list.append(t2)

            reco_pattern = Xplet()
            for triplet in t_list:
                reco_pattern.add_triplet(triplet)
            if fit == "chi squared lin track":
                reco_pattern.fit_lin_track()
            reco_x_plets.append(reco_pattern)

    hit_to_xplet_map = {}
    ambiguity_candidates = []
    for i, xplet in enumerate(reco_x_plets):
        for hit_id in xplet.hit_ids.values():
            if hit_id not in hit_to_xplet_map.values():
                hit_to_xplet_map.update({hit_id: i})
            else:
                ambiguity_candidates.append(i)
                ambiguity_candidates.append(hit_to_xplet_map[hit_id])

    print("Solving ambiguities...\n")
    remove_triplets_indices = []
    for i, xplet_one in enumerate(ambiguity_candidates):
        for j, xplet_two in enumerate(ambiguity_candidates[i + 1:]):
            for hit_id_one in reco_x_plets[xplet_one].hit_ids:
                for hit_id_two in reco_x_plets[xplet_two].hit_ids:
                    if hit_id_one == hit_id_two:
                        if reco_x_plets[xplet_one].p_value > reco_x_plets[xplet_two].p_value:
                            remove_triplets_indices.append(i + j)
                        else:
                            remove_triplets_indices.append(j)

    if len(ambiguity_candidates) > 0:
        reco_ambiguity_solved = []
        for i, xplet in reco_x_plets:
            if i not in remove_triplets_indices:
                reco_ambiguity_solved.append(xplet)
    else:
        reco_ambiguity_solved = reco_x_plets

    print(f"Number of combinatorial reco Xplets: {len(reco_ambiguity_solved)}\n")
    print(f"Number of reco Xplets after ambiguity solving: {len(reco_ambiguity_solved)}\n")

    np.save(f"{save_folder}/reco_xplet_list", reco_ambiguity_solved)
