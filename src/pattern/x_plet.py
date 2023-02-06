from pattern.triplet import Triplet
from scipy.optimize import curve_fit
from scipy.stats import chisquare, chi2


class Xplet:
    def __init__(self):
        self.hit_ids = {}
        self.particle_ids = {}
        self.coordinates = {}
        self.energy = {}
        self.triplet_ids = []
        self.time = {}
        self.is_empty = True
        self.chi_squared = None
        self.p_value = None
        """Class for combining connectable triplet structures into a x-plet pattern. Inside the Xplet pattern structure
           information about the single hits are stored. The chi_squared and p-value can be set via a linear fit. 
           Additional options will be added in the future. 
        """

    def add_triplet(self,
                    triplet: Triplet):
        """Adds a triplet to the Xplet structure
        :param triplet: Triplet object
        """
        if self.is_empty:
            self.hit_ids.update({0: triplet.doublet_1.hit_1_id})
            self.hit_ids.update({1: triplet.doublet_1.hit_2_id})
            self.hit_ids.update({2: triplet.doublet_2.hit_2_id})

            self.particle_ids.update({0: triplet.doublet_1.hit_1_particle_key})
            self.particle_ids.update({1: triplet.doublet_1.hit_2_particle_key})
            self.particle_ids.update({2: triplet.doublet_2.hit_2_particle_key})

            self.coordinates.update({0: triplet.doublet_1.hit_1_position})
            self.coordinates.update({1: triplet.doublet_1.hit_2_position})
            self.coordinates.update({2: triplet.doublet_2.hit_2_position})

            self.energy.update({0: triplet.doublet_1.energy_1})
            self.energy.update({1: triplet.doublet_1.energy_2})
            self.energy.update({2: triplet.doublet_2.energy_2})

            self.time.update({0: triplet.doublet_1.time_1})
            self.time.update({1: triplet.doublet_1.time_2})
            self.time.update({2: triplet.doublet_2.time_2})
            self.is_empty = False

        else:
            insert_position = len(self.triplet_ids) + 2
            self.hit_ids.update({insert_position: triplet.doublet_2.hit_2_id})
            self.particle_ids.update({insert_position: triplet.doublet_2.hit_2_particle_key})
            self.coordinates.update({insert_position: triplet.doublet_2.hit_2_position})
            self.energy.update({insert_position: triplet.doublet_2.energy_2})
            self.time.update({insert_position: triplet.doublet_2.time_2})

        self.triplet_ids.append(triplet.triplet_id)

    @staticmethod
    def lin_func(x, a, b):
        """Linear function with slope a and bias b
        :param x: data points
        :param a: slope
        :param b: bias
        """
        return a * x + b

    def get_coordinates(self):
        """Auxiliary function to extract x,y and z coordinates.
        :return
            (x, y, z)
        """
        return ([value[0] for value in list(self.coordinates.values())],
                [value[1] for value in list(self.coordinates.values())],
                [value[2] for value in list(self.coordinates.values())])

    def fit_lin_track(self):
        """Linear fit of the track in xz and yz direction. Chi squared values are averaged. Chi squared and p-value
        attributes are set for the xplet.
        """
        x = [value[0] for value in self.coordinates.values()]
        y = [value[1] for value in self.coordinates.values()]
        z = [value[2] for value in self.coordinates.values()]

        popt_xz, _ = curve_fit(Xplet.lin_func, z, x)
        popt_yz, _ = curve_fit(Xplet.lin_func, z, y)

        f_exp = [popt_xz[0] * z_i + popt_xz[1] for z_i in z]
        chi_xz = sum([((x_i - e) / 5e-6) ** 2 for x_i, e in zip(x, f_exp)]) / (len(x) - 2)

        f_exp = [popt_yz[0] * z_i + popt_yz[1] for z_i in z]
        chi_yz = sum([((y_i - e) / 5e-6) ** 2 for y_i, e in zip(y, f_exp)]) / (len(y) - 2)

        self.chi_squared = 0.5 * (chi_xz + chi_yz)
        self.p_value = chi2.sf(0.5 * (chi_xz + chi_yz), df=len(x) - 1 - (len(x) - 2))
