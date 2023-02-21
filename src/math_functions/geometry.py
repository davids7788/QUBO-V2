import numpy as np

from numba import jit
from math import atan2, sqrt


@jit(nopython=True)
def x0_at_z_ref(x_end: float,
                x_start: float,
                z_end: float,
                z_start: float,
                z_ref: float) -> float:
    """Function for calculation x position of a doublet at a z-reference value.
    :param x_end: x-position of second hit
    :param x_start: x-position of first hit
    :param z_end: z-position of second hit
    :param z_start: z-position of first hit
    :param z_ref: z-value reference layer
    :return:
        x_position at the reference layer
    """
    dx = x_end - x_start
    dz = z_end - z_start
    return x_end - dx * abs(z_end - z_ref) / dz


@jit(nopython=True)
def xyz_angle(xy_1: float,
              xy_2: float,
              z_2: float,
              z_1: float) -> float:
    """Returns the angle in xz or yz with respect to the given coordinates.
    :param xy_1: x or y coordinate of the second hit
    :param xy_2: x or y coordinate of the first hit
    :param z_1: z coordinate of the second hit
    :param z_2: z coordinate of the first hit

    :return:
        angle in rad
    """
    return atan2((xy_2 - xy_1), (z_2 - z_1))


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


def two_norm_std_angle(angle_xz_1,
                       angle_xz_2,
                       angle_xz_3,
                       angle_yz_1,
                       angle_yz_2,
                       angle_yz_3):
    """Returns 2-norm of angle difference in xz and yz.
    :param angle_xz_1 : xz angle of first doublet
    :param angle_xz_2 : xz angle of second doublet
    :param angle_xz_3 : xz angle of third doublet
    :param angle_yz_1 : yz angle of first doublet
    :param angle_yz_2 : yz angle of second doublet
    :param angle_yz_3 : yz angle of third doublet
    :return
        2-norm of angle difference in xz and yz.
    """
    return sqrt(np.std(np.array([angle_xz_1, angle_xz_2, angle_xz_3]))**2
                + np.std(np.array([angle_yz_1, angle_yz_2, angle_yz_3])) ** 2)


