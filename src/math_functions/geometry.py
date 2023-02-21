import numpy as np
from numba import jit
from math import sqrt, atan2


@jit(nopython=True)
def triplet_criteria_check(xz_angle_1: float,
                           xz_angle_2: float,
                           yz_angle_1: float,
                           yz_angle_2: float,
                           max_angle) -> bool:
    """Checks if doublets may be combined to a triplet, depending on the doublet angles -> scattering
    :param xz_angle_1: angle xz doublet 1
    :param xz_angle_2: angle xz doublet 2
    :param yz_angle_1: angle yz doublet 1
    :param yz_angle_2: angle yz doublet 2
    :param max_angle: max angle criteria [rad]
    :return:
        True if criteria applies, else False
    """
    if sqrt((xz_angle_2 - xz_angle_1) ** 2 + (yz_angle_2 - yz_angle_1) ** 2) < max_angle:
        return True
    return False


@jit(nopython=True)
def dy_x0_check(y1: float,
                y2: float,
                x0: float,
                criteria: float) -> bool:
    """Auxiliary function  for checking dy_x0 criteria.
    :param y1: y value of first hit
    :param y2: y value of second hit
    :param x0: extrapolated x value on reference layer
    :param criteria to decide if doublet is kept or discarded
    :return
        True if condition is satisfied, else False
    """
    if abs(y2 - y1) / x0 > criteria:
        return False
    return True


@jit(nopython=True)
def dx_x0_check(x1: float,
                x2: float,
                x0: float,
                criteria_mean,
                criteria_eps) -> bool:
    """Checks if hits may be combined to doublets, applying dx/x0 criterion
    :param x1: x value first hit
    :param x2: x value second hit
    :param x0: extrapolated x value on reference layer
    :param criteria_mean: mean value of expected criteria
    :param criteria_eps: epsilon range of expected criteria
    :return:
        True if criteria applies, else False
    """
    if abs((x2 - x1) / x0 - criteria_mean) > criteria_eps:
        return False
    return True


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
