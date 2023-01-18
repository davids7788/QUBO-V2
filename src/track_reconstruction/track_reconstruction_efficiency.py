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
    print(gen_prefix)
    print("/".join(reco_xplet_file.split("/")[0:-3]))
    gen_xplet = np.load("/".join(reco_xplet_file.split("/")[0:-3]) + "/" + gen_prefix + "_gen_xplet_list.npy",
                        allow_pickle=True)

    matched_tracks = 0
    fake_tracks = 0

    for track in reco_xplets:
        if len(set(track.particle_ids.values())) < 0.5 * len(list(track.particle_ids.values())):
            matched_tracks += 1
        else:
            fake_tracks += 1

    print("\n-----------------------------------\n")
    print(f"Track reconstruction statistics:"
          f"Efficiency: {100 * np.around(matched_tracks / len(gen_xplet), 3)} %\n"
          f"Fake Rate: {100 * np.around(fake_tracks / len(reco_xplets), 3)} %")
