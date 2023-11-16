import numpy as np

from pattern.multiplet import Multiplet
from pattern.triplet import Triplet

from utility.data_format_handler import write_tracks_to_slcio_file


def make_reco_multiplets(triplets,
                         save_folder,
                         tracking_data,
                         fit="chi squared lin track"):
    """Creates xplets from the kept triplets.
    :param fit: 
            "chi squared lin track: assuming a linear track and fitting  a chi squared assuming f(z) = x_i, f(z) = y_j
    :param triplets: list of kept triplets
    :param save_folder: folder to store reco xplet list
    """
    print("\n-----------------------------------\n")
    print("Creating reco Xplets and fitting resulting tracks...\n")

    def triplet_start(t_obj: Triplet):
        """Returns z value of first triplet hit
        """
        return t_obj.hit_1.z

    # set for storing possible start values in z of triplets --> start of tracks in z
    triplet_start_value = set()
    for triplet in triplets:
        triplet_start_value.add(triplet.hit_1.z)
    triplets.sort(key=triplet_start)

    # check for simple or full LUXE setup, len==2 --> simple,  else full setup
    sorted_z = sorted(triplet_start_value)
    if len(triplet_start_value) == 2:
        t_possible_start = [sorted_z[0]]
    else:
        t_possible_start = list(sorted_z)[0:2]

    # order triplets by their lowest z coordinate
    triplets_ordered_by_start_value = {}
    # create empty lists per start z value
    for z_start in sorted_z:
        triplets_ordered_by_start_value.update({z_start: []})
    # fill lists with triplets
    for t in triplets:
        triplets_ordered_by_start_value[t.hit_1.z].append(t)

    reco_track_candidates = []

    for start_z in t_possible_start:
        # start of xplets
        current_z = start_z
        for t_start in triplets_ordered_by_start_value[start_z]:
            reco_multiplets = [[t_start]]
            for next_z in triplets_ordered_by_start_value.keys():
                if next_z > current_z:
                    for t_next in triplets_ordered_by_start_value[next_z]:
                        for multiplets in reco_multiplets:
                            if multiplets[-1].hit_2 == t_next.hit_1 and multiplets[-1].hit_3 == t_next.hit_2:
                                new_multiplet = [t for t in multiplets] + [t_next]
                                reco_multiplets.append(new_multiplet)
            skimmed_reco_multiplets = []
            for i, r_m in enumerate(reco_multiplets):
                is_part_of = False
                for r_m_amb in reco_multiplets[i + 1:]:
                    if len(r_m_amb) >= len(r_m):
                        if any([t1 == t2 for t1, t2 in zip(r_m, r_m_amb[:len(r_m)])]):
                            is_part_of = True
                if not is_part_of:
                    skimmed_reco_multiplets.append(r_m)
            for s_k in skimmed_reco_multiplets:
                new_multiplet = Multiplet()
                for t_s_k in s_k:
                    new_multiplet.add_hit(t_s_k.hit_1)
                new_multiplet.add_hit(s_k[-1].hit_2)
                new_multiplet.add_hit(s_k[-1].hit_3)
                reco_track_candidates.append(new_multiplet)

    print(f"Number of combinatorial reco Xplets: {len(reco_track_candidates)}\n"
          f"Saving combinatorial Xplets to file {save_folder}/reco_xplet_list")

    if '.slcio' in tracking_data:
        write_tracks_to_slcio_file(tracking_data,
                                   reco_track_candidates,
                                   save_folder)

    np.save(f"{save_folder}/reco_xplet_list", reco_track_candidates)

    hit_to_xplet_map = {}
    for multiplet in reco_track_candidates:
        for hit_id in multiplet.hit_id:
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
        for h_id in reco_track.hit_id:
            if hit_to_xplet_map[h_id] > 1:
                overlaps += 1
        if overlaps > max_overlaps:
            return True
        return False

    reco_ambiguity_solved = None
    for i in range(3):
        selected_tracks = []
        for track in reco_track_candidates:
            if not select_triplets_with_max_overlap(track, 4 - i):
                selected_tracks.append(track)
        reco_ambiguity_solved = [t for t in selected_tracks if len(t.hit_id) > 3]

    print(f"Number of reco multiplets after ambiguity solving: {len(reco_ambiguity_solved)}\n"
          f"Saving multiplets")

    np.save(f"{save_folder}/reco_xplet_list_ambiguity_solved", reco_ambiguity_solved)
