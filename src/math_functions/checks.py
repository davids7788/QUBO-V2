from numba import jit
from math import sqrt

from pattern.detector_hit import DetectorHit
from math_functions.geometry import xyz_angle, x0_at_z_ref


@jit(nopython=True)
def dxy_x0_check(xy1: float,
                 xy2: float,
                 z1: float,
                 z2: float,
                 x0: float,
                 criteria_mean: float = 0,
                 criteria_eps: float = 0) -> bool:
    """Checks doublet dx/x0 and dy/x0 criteria.

    :param xy1: x or y value og the first hit
    :param xy2: x or y value of the second hit
    :param z1: z value of first hit
    :param z2: z value of second hit
    :param x0: extrapolated x value on reference layer
    :param criteria_mean: mean value of dx/x0 criteria
    :param criteria_eps: allowed epsilon range of dx/x0 criteria

    :return
        True if criteria is fulfilled, else False
    """
    if abs((xy2 - xy1) / x0 / abs(z2 - z1) - criteria_mean) > criteria_eps:
        return False
    return True


@jit(nopython=True)
def jit_is_valid_triplet(xz_angle_1: float,
                         xz_angle_2: float,
                         yz_angle_1: float,
                         yz_angle_2: float,
                         max_angle: float) -> bool:
    """jit-enhanced function for "is_valid_triplet" function.

    :param xz_angle_1: angle xz doublet 1
    :param xz_angle_2: angle xz doublet 2
    :param yz_angle_1: angle yz doublet 1
    :param yz_angle_2: angle yz doublet 2
    :param max_angle: max angle criteria

    :return:
        True if criteria applies, else False
    """
    if sqrt((xz_angle_2 - xz_angle_1) ** 2 + (yz_angle_2 - yz_angle_1) ** 2) < max_angle:
        return True
    return False


def is_valid_triplet(hit_1: DetectorHit,
                     hit_2: DetectorHit,
                     hit_3: DetectorHit,
                     max_scattering: float) -> bool:
    """Checks if doublets may be combined to a triplet, depending on the doublet angles -> scattering
    :param hit_1: first hit of a triplets
    :param hit_2: second hit of a triplets
    :param hit_3: third hit of a triplets
    :param max_scattering: maximum allowed scattering

    :return:
        True if criteria applies, else False
    """
    xz_12 = xyz_angle(hit_1.x, hit_2.x, hit_1.z, hit_2.z)
    xz_23 = xyz_angle(hit_2.x, hit_3.x, hit_2.z, hit_3.z)

    yz_12 = xyz_angle(hit_1.y, hit_2.y, hit_1.z, hit_2.z)
    yz_23 = xyz_angle(hit_2.y, hit_3.y, hit_2.z, hit_3.z)

    return jit_is_valid_triplet(xz_12,
                                xz_23,
                                yz_12,
                                yz_23,
                                max_scattering)


def is_valid_doublet(first_hit: DetectorHit,
                     second_hit: DetectorHit,
                     z_ref: float,
                     configuration) -> bool:
    """Checks if two hits are actually a doublet candidate.
    :param first_hit: hit with lower z-value
    :param second_hit: hit with higher z-value
    :param z_ref: reference layer z-value
    :param configuration: configuration file

    :return
        True if valid doublet candidate, else False
    """
    # calculate x0
    x0 = x0_at_z_ref(second_hit.x,
                     first_hit.x,
                     second_hit.z,
                     first_hit.z,
                     z_ref)

    # check dx / x0 criteria
    if not dxy_x0_check(first_hit.x,
                        second_hit.x,
                        first_hit.z,
                        second_hit.z,
                        x0,
                        criteria_mean=configuration["doublet"]["dx/x0"],
                        criteria_eps=configuration["doublet"]["dx/x0 eps"]):
        return False

    # check dy / x0 criteria
    if not dxy_x0_check(first_hit.y,
                        second_hit.y,
                        first_hit.z,
                        second_hit.z,
                        x0,
                        criteria_mean=configuration["doublet"]["dy/x0"],
                        criteria_eps=configuration["doublet"]["dy/x0 eps"]):
        return False
    return True
