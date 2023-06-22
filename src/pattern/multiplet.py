from pattern.triplet import DetectorHit
from scipy.optimize import curve_fit
from scipy.stats import chisquare, chi2
from numba import jit


class Multiplet:
    def __init__(self):
        self.hit_ids = []
        self.particle_ids = []
        self.x = []
        self.y = []
        self.z = []
        self.particle_energy = []
        self.time = []

        self.chi_squared = None
        self.p_value = None
        """Class for multiplet objects. Multiplets store information about hits from consecutive detector layers.
        The multiplet can be fitted to a track with a linear function. 
        """

    def add_hit(self,
                hit: DetectorHit) -> None:
        """Adds a triplet to the X-plet structure
        :param hit: detector hit
        """
        self.hit_ids.append(hit.hit_id)
        self.particle_ids.append(hit.particle_id)
        self.x.append(hit.x)
        self.y.append(hit.y)
        self.z.append(hit.z)
        self.particle_energy.append(hit.particle_energy)
        if hit.time:
            self.time.append(hit.time)

    @staticmethod
    @jit(nopython=True)
    def lin_func(x: float,
                 a: float,
                 b: float) -> float:
        """Linear function with slope a and bias b
        :param x: data points
        :param a: slope
        :param b: bias
        """
        return a * x + b

    def fit_lin_track(self,
                      detector_resolution: float = 5e-6) -> None:
        """Linear fit of the track in xz and yz direction. Chi squared values are averaged. Chi squared and p-value
        attributes are set for the x-plet.
        :param detector_resolution: resolution of the detector to calculate the reduced chi-squared [m]
        """

        # fit x(z) and y(z)
        popt_xz, _ = curve_fit(Multiplet.lin_func, self.z, self.x)
        popt_yz, _ = curve_fit(Multiplet.lin_func, self.z, self.y)

        f_exp = [popt_xz[0] * z_i + popt_xz[1] for z_i in self.z]
        chi_xz = sum([((x_i - e) / detector_resolution) ** 2 for x_i, e in zip(self.x, f_exp)])

        f_exp = [popt_yz[0] * z_i + popt_yz[1] for z_i in self.z]
        chi_yz = sum([((y_i - e) / detector_resolution) ** 2 for y_i, e in zip(self.y, f_exp)])

        self.chi_squared = (chi_xz + chi_yz) / (len(self.x) + len(self.y) - 4)  # is never zero, len(x) + len(y) >= 6
        self.p_value = chi2.sf(0.5 * (chi_xz + chi_yz), df=len(self.x) - 1 - (len(self.x) - 2))
