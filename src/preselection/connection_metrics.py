import numpy as np


def two_norm_std_angle(doublet_1,
                       doublet_2,
                       doublet_3):
    """Returns 2-norm of angle difference in xz and yz.
    :param doublet_1 : doublet from hit 1 + 2
    :param doublet_2 : doublet from hit 2 + 3
    :param doublet_3 : doublet from hit 3 + 4
    :return
        2-norm of angle difference in xz and yz.
    """
    angle_xz_doublet_1 = doublet_1.xz_angle()
    angle_yz_doublet_1 = doublet_1.yz_angle()

    angle_xz_doublet_2 = doublet_2.xz_angle()
    angle_yz_doublet_2 = doublet_2.yz_angle()

    angle_xz_doublet_3 = doublet_3.xz_angle()
    angle_yz_doublet_3 = doublet_3.yz_angle()

    return np.sqrt(np.std([angle_xz_doublet_1, angle_xz_doublet_2, angle_xz_doublet_3]) ** 2 +
                   np.std([angle_yz_doublet_1, angle_yz_doublet_2, angle_yz_doublet_3]) ** 2)
