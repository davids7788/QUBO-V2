import numpy as np

from numba import jit
from math import atan2


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


def angle_based_measure(detector_hits: list[list[float, float, float]]) -> float:
    """Returns a measure of how well the triplets match using angle information, assuming a linear track
    :return
        max_diff_xz angle * max_diff_yz_angle
    """
    x = detector_hits[0]
    y = detector_hits[1]
    z = detector_hits[2]

    xz = [xyz_angle(x[i], x[i + 1], z[i], z[i + 1]) for i in range(len(detector_hits[0]) - 1)]
    yz = [xyz_angle(y[i], y[i + 1], z[i], z[i + 1]) for i in range(len(detector_hits[1]) - 1)]

    return np.sqrt(np.std(xz)**2 + np.std(yz)**2)
