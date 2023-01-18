import os

import numpy as np

from pattern.x_plet import Xplet


def track_reconstruction_efficiency_simplified_LUXE(reco_xplet_file: str):
    """Takes a reco Xplet file and looks for the corresponding gen_Xplet file and computes the
    track reconstruction efficiency. The metrics are:
            efficiency: #matched tracks / #generated tracks
            fake rate: # fake tracks/ #reconstructed tracks
    with a matched track being a track in which the majority of hits belong to a single particle. A fake track is
    defined as not a matched track.
    :param reco_xplet_file: string
    """
    reco_xplets = np.load(reco_xplet_file, allow_pickle=True)
    gen_prefix = reco_xplet_file.split("/")[-3].split("-")[0]
    gen_xplet = np.load("/".join(reco_xplet_file.split("/")[0:-3]) + "/" + gen_prefix + "_gen_xplet_list.npy",
                        allow_pickle=True)

    matched_tracks = 0
    fake_tracks = 0

    matched_tracks_bookkeeping = set()

    for track in reco_xplets:
        matched = False
        p_id = None
        ids = set(track.particle_ids.values())
        for test_id in ids:
            count = 0
            for particle_id in track.particle_ids.values():
                if test_id == particle_id:
                    count += 1
            if count >= 3:
                p_id = test_id
                matched = True
        if matched:
            if p_id not in matched_tracks_bookkeeping:
                matched_tracks += 1
                matched_tracks_bookkeeping.add(p_id)
            else:
                pass
        else:
            fake_tracks += 1

    print("\n-----------------------------------\n")
    print(f"Track reconstruction statistics:\n"
          f"Efficiency: {100 * np.round(matched_tracks / len(gen_xplet), 3)} %\n"
          f"Fake Rate: {100 * np.round(fake_tracks / len(reco_xplets), 3)} %")
