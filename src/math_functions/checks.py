from numba import jit
from math import sqrt


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


@jit(npython=True)
def _is_in_segment(x: float,
                   y: float,
                   z: float,
                   x_start: float,
                   x_end: float,
                   y_start: float,
                   y_end: float,
                   z_position: float) -> float:
    """Checks if position x is inside the segment
    :param x: position in x [m]
    :param y: position in y [m]
    :param z: position in z [m]
    :param x_start: start of segment in x [m]
    :param x_end: end of segment in x [m]
    :param y_start: start of segment in   [m]
    :param y_end: end of segment in y [m]
    :param z_position position of layer in z [m]
    :return:
        True if coordinate in segment, else False
    """
    if z != z_position or x < x_start or x > x_end or y < y_start or y > y_end:
        return False
    return True


@jit(nopython=True)
def _is_valid_triplet(xz_angle_1: float,
                      xz_angle_2: float,
                      yz_angle_1: float,
                      yz_angle_2: float,
                      max_angle: float) -> bool:
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
