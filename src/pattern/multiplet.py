from pattern.triplet import DetectorHit
from scipy.optimize import curve_fit
from scipy.stats import chi2
from numba import jit


class Multiplet:
    """Class for storing information from DetectorHit objects in a multiplet object.
    """
    def __init__(self):
        """Preparing fields.
        """
        self.hit_id = []
        self.x = []
        self.y = []
        self.z = []
        self.layer_id = []
        self.module_id = []
        self.cell_id = []
        self.particle_id = []
        self.is_signal = []
        self.particle_energy = []

        self.chi_squared = None
        self.p_value = None

    def add_hit(self,
                hit: DetectorHit) -> None:
        """Adds a hit to the multiplet structure.
        :param hit: detector hit
        """
        self.hit_id.append(hit.hit_id)
        self.x.append(hit.x)
        self.y.append(hit.y)
        self.z.append(hit.z)
        self.layer_id.append(hit.layer_id)
        self.cell_id.append(hit.cell_id)
        self.module_id.append(hit.module_id)
        self.particle_id.append(hit.particle_id)
        self.is_signal.append(hit.is_signal)
        self.particle_energy.append(hit.particle_energy)

    @staticmethod
    @jit(nopython=True)
    def lin_func(x: float,
                 a: float,
                 b: float) -> float:
        """Linear function with slope a and bias b
        :param x: data points
        :param a: slope
        :param b: bias

        :return
            f(x) = a * x + b
        """
        return a * x + b

    def fit_lin_track(self,
                      detector_resolution: float = 5e-6) -> None:
        """Sets a chi squared and p-value for the track assuming a linear fit in xz and yz direction.
        :param detector_resolution: resolution of the detector to calculate the reduced chi-squared [m]
        """

        # fit x(z) and y(z)
        popt_xz, _ = curve_fit(Multiplet.lin_func, self.z, self.x)
        popt_yz, _ = curve_fit(Multiplet.lin_func, self.z, self.y)

        # chi squared x(z)
        f_exp_xz = [popt_xz[0] * z_i + popt_xz[1] for z_i in self.z]
        chi_xz = sum([((x_i - e) / detector_resolution) ** 2 for x_i, e in zip(self.x, f_exp_xz)])

        # chi squared y(z)
        f_exp_yz = [popt_yz[0] * z_i + popt_yz[1] for z_i in self.z]
        chi_yz = sum([((y_i - e) / detector_resolution) ** 2 for y_i, e in zip(self.y, f_exp_yz)])

        # set chi squared and p-value
        self.chi_squared = (chi_xz + chi_yz) / (len(self.x) + len(self.y) - 4)  # is never zero, len(x) + len(y) >= 6
        self.p_value = chi2.sf(0.5 * (chi_xz + chi_yz), df=len(self.x) - 1 - (len(self.x) - 2))
