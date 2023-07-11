import numpy as np


class DetectorPlane:
    """Class for simulating a detector plane for a Toy MC experiment.
    :param limits_x                   : (x_start, x_end) [m]
    :param limits_y                   : (y_start, y_end) [m]
    :param z_position:                : z_position [m]
    :param num_bins                   : (bins_x, bins_y) -> 2D-histogram
    :param effective_thickness        : the detector consists of more than one material, to simplify it is assumed a
                                        combined length (=interaction length), where particles are passing through [m]
    :param effective_radiation_length : length on which particles are interacting with the detector in [kg / m^2]
    :param detector_rho               : effective density of the detector material [kg / m³]
    """
    def __init__(self,
                 limits_x,
                 limits_y,
                 z_position,
                 num_bins,
                 effective_thickness,
                 effective_radiation_length,
                 detector_rho):
        self.limits_x = limits_x
        self.limits_y = limits_y
        self.num_bins = num_bins
        self.z_position = z_position
        self.effective_thickness = effective_thickness
        self.effective_radiation_length = effective_radiation_length
        self.detector_rho = detector_rho

        self.true_hits_dictionary = {}

    def bin_width_in_x(self):
        """Returns bin width in x calculated from number of bins and length of detector
        :return bin width in x [m]
        """
        return (self.limits_x[1] - self.limits_x[0]) / self.num_bins[0]
    
    def bin_width_in_y(self):
        """Returns bin width in y calculated from number of bins and length of detector
        :return bin width in y [m]
        """
        return (self.limits_y[1] - self.limits_y[0]) / self.num_bins[1]
    
    def bin_edges_for_hit(self, hit):
        """Returns bin edges in x and y for a particle position (hit) on the detector
        :return (x_low, x_up), (y_low, y_up) [m]
        """
        x_steps_upper_limit = int((self.limits_x[1] - hit[0]) / self.bin_width_in_x())
        y_steps_upper_limit = int((self.limits_y[1] - hit[1]) / self.bin_width_in_y())
        return (self.limits_x[1] - (x_steps_upper_limit + 1) * self.bin_width_in_x(),
                self.limits_x[1] - x_steps_upper_limit * self.bin_width_in_x()), \
               (self.limits_y[1] - (y_steps_upper_limit + 1) * self.bin_width_in_y(),
                self.limits_y[1] - y_steps_upper_limit * self.bin_width_in_y())

    def gaussian_scattering(self, particle_momentum):
        """Returns a gaussian distributed value from a sample with mu and sigma as parameters
        :param: particle_momentum : vector of momenta
        """
        # angle in rad
        angle_xz = np.arctan2(particle_momentum[0], particle_momentum[2])
        angle_yz = np.arctan2(particle_momentum[1], particle_momentum[2])

        # energy
        p_abs = np.sqrt(particle_momentum[0]**2 + particle_momentum[1]**2 + particle_momentum[2]**2)

        # using rossi formula
        def calculate_rossi_theta(d_rho, thickness, x_0, momentum):
            """Auxiliary function for the calculation of the scattering angle
            :param d_rho           : rho of detector in [kg / m³]
            :param thickness       : thickness of the detector layer
            :param x_0             : effective radiation length
            :param momentum        : momentum of the particle
            :return rossi theta
            """
            beta = 1  # very fast particles
            c = 1  # set c = 1
            return 14 / (beta * c * momentum) * np.sqrt(d_rho * thickness / x_0)

        # 'scattering' existing angles with gaussian distribution around existing angle with sigma from rossi formula
        angle_xz_after_scattering = np.random.normal(angle_xz,
                                                     calculate_rossi_theta(self.detector_rho,
                                                                           self.effective_thickness /
                                                                           np.cos(angle_xz),
                                                                           self.effective_radiation_length,
                                                                           p_abs),
                                                     size=1)[0]
        angle_yz_after_scattering = np.random.normal(angle_yz,
                                                     calculate_rossi_theta(self.detector_rho,
                                                                           self.effective_thickness /
                                                                           np.cos(angle_yz),
                                                                           self.effective_radiation_length,
                                                                           p_abs),
                                                     size=1)[0]

        # new momenta after scattering
        new_momentum_x = np.tan(angle_xz_after_scattering) * particle_momentum[2]
        new_momentum_y = np.tan(angle_yz_after_scattering) * particle_momentum[2]
        new_momentum_z = np.sqrt(p_abs**2 - new_momentum_x**2 - new_momentum_y**2)

        # energy loss per layer
        new_p_abs = p_abs * np.exp(- (np.sqrt((self.effective_thickness / np.cos(angle_xz)) ** 2 +
                                              (self.effective_thickness / np.cos(angle_yz)) ** 2)) /
                                   self.effective_radiation_length)

        # assuming to be the same correction factor for px, py, pz
        p_eff_correction_factor = abs(new_p_abs / p_abs)

        return [new_momentum_x * p_eff_correction_factor,
                new_momentum_y * p_eff_correction_factor,
                new_momentum_z * p_eff_correction_factor]
