import numpy as np


class DipoleMagnet:
    """Specifies a dipole magnet and provides function for calculating particle trajectories through it.
    :param:
        dipole start    : start of dipole location [m]
        dipole end      : end of dipole location [m]
        dipole_strength : B-Field [T]
    """

    def __init__(self,
                 dipole_start,
                 dipole_end,
                 dipole_strength):
        self.dipole_start = dipole_start
        self.dipole_end = dipole_end
        self.dipole_length = self.dipole_end - self.dipole_start
        self.dipole_strength = dipole_strength

    @staticmethod
    def xz_rotation_matrix(angle):
        """Returns a rotational matrix for the xz_plane
        Parameters:
            angle   : angle in rad to rotate
        """
        return np.array([[np.cos(angle), 0, - np.sin(angle)],
                         [0, 1, 0],
                         [np.sin(angle), 0, np.cos(angle)]])
    
    @staticmethod
    def compute_bending_radius(particle_momentum,
                               dipole_field_strength):
        """Returns the bending radius for a particle with given momentum.
        :param
            particle_momentum: momentum of a particle"""
        return 10 / 2.99792458 * 1e-3 * particle_momentum / dipole_field_strength

    def traversing_through_dipole_magnet(self,
                                         particle):
        """Computes the particles momenta and positions after traversing through the magnetic dipole.
           Therefore a transfer matrix is used. Since the magnetic field conserves the energy,
           pz (beam is in z-direction) is adjusted to account for that.
        :param
            particle_momentum: [px, py, pz]
            particle_position: [x, y, z]
        """
        theta_i = np.arctan2(particle.momentum[0], particle.momentum[2])
        for dipole_strength_value in self.dipole_strength:
            # radius
            r = abs(DipoleMagnet.compute_bending_radius(np.sqrt(particle.momentum[0]**2 + particle.momentum[2]**2),
                                                        dipole_strength_value))
            # print(R, "radius")
            
            # theta_i
            dz = self.dipole_length / len(self.dipole_strength)
            
            # theta 
            theta = np.arcsin(np.sin(theta_i) + dz / r) - theta_i

            dx = dz * np.tan(theta_i + 0.5 * theta)
            theta_i += theta
            
            dy = r * theta / np.sqrt(particle.momentum[2]**2 + particle.momentum[0]**2) * particle.momentum[1]
            
            particle.position += np.array([dx, dy, dz])
            particle.momentum = DipoleMagnet.xz_rotation_matrix(- theta).dot(particle.momentum) 
