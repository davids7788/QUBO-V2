import numpy as np

from math import floor


def track_reconstruction_efficiency_simplified_LUXE(reco_xplet_file: str,
                                                    gen_xplet_file: str):
    """Takes a reco Xplet file and looks for the corresponding gen_Xplet file and computes the
    track reconstruction efficiency. The metrics are:
            efficiency: #matched tracks / #generated tracks
            fake rate: #fake tracks/ #reconstructed tracks
    with a matched track being a track in which the majority of hits belong to a single particle. A fake track is
    defined as not a matched track.
    :param reco_xplet_file: string
    :param gen_xplet_file: string
    """
    reco_xplets = np.load(reco_xplet_file, allow_pickle=True)
    gen_xplet = [g for g in np.load(gen_xplet_file, allow_pickle=True) if len(g.hit_id) >= 4]

    matched_tracks = 0
    fake_tracks = 0

    matched_tracks_bookkeeping = set()

    for track in reco_xplets:
        matched = False
        p_id = None
        # ids  corresponds to the particle id's of each hit
        for test_id in track.particle_id:
            # more than half of the hits from same particle
            if track.particle_id.count(test_id) > len(track.particle_id) / 2:
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
          f"Efficiency: {np.round(100 * matched_tracks / len(gen_xplet), 2)} %\n"
          f"Fake Rate: {np.round(100 * fake_tracks / len(reco_xplets), 2)} %")
