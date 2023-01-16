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
    print("Creating reco Xplets and fitting resulting track...")

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

    print(f"Number of generated Xplets: {len(reco_x_plets)}\n")

    np.save(f"{save_folder}/reco_xplet_list", reco_x_plets)
