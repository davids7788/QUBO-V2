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
              z_1: float,
              z_2: float) -> float:
    """Returns the angle in xz or yz with respect to the given coordinates.
    :param xy_1: x or y coordinate of the second hit
    :param xy_2: x or y coordinate of the first hit
    :param z_1: z coordinate of the first hit
    :param z_2: z coordinate of the second hit

    :return:
        angle in rad
    """
    return atan2((xy_2 - xy_1), (z_2 - z_1))


@jit(nopython=True)
def _default_angle_based_interaction(angle_xz_1,
                                     angle_xz_2,
                                     angle_xz_3,
                                     angle_yz_1,
                                     angle_yz_2,
                                     angle_yz_3) -> float:
    """Returns angle based in xz and yz.
    :param angle_xz_1 : xz angle of first doublet
    :param angle_xz_2 : xz angle of second doublet
    :param angle_xz_3 : xz angle of third doublet
    :param angle_yz_1 : yz angle of first doublet
    :param angle_yz_2 : yz angle of second doublet
    :param angle_yz_3 : yz angle of third doublet
    :return
        2-norm of the standard deviation of the angle difference in xz and yz.
    """
    return sqrt(np.std(np.array([angle_xz_1, angle_xz_2, angle_xz_3]))**2
                + np.std(np.array([angle_yz_1, angle_yz_2, angle_yz_3]))**2)


@jit(nopython=True)
def _default_angle_based_quality(angle_between_doublets_xz: float,
                                 angle_between_doublets_yz: float) -> float:
    """Returns angle based in xz and yz.
    :param angle_between_doublets_xz : xz angle between doublets
    :param angle_between_doublets_yz : xz angle between doublets

    :return
        two norm of angles in xy and xz between doublets
    """
    return sqrt(angle_between_doublets_xz**2 + angle_between_doublets_yz**2)
