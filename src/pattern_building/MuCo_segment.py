class MuCoDetectorSegment:
    def __init__(self,
                 name: str,
                 phi_start: float,
                 phi_end: float,
                 theta_start: float,
                 theta_end: float,
                 dist_start: float,
                 dist_end: float):
        """Manages a segment of the Muon Collider detector layers.
        :param name: unique identifier for the segment
        :param phi_start: start of segment in phi [rad]
        :param phi_end: end of segment in phi [rad]
        :param theta_start: start of segment in theta[rad]
        :param theta_end: end of segment in theta [rad]
        :param dist_start: start of segment in distance from IP [m]
        :param dist_end: end of segment in distance from IP [m]
        """
        self.name = name
        self.phi_start = phi_start
        self.phi_end = phi_end
        self.theta_start = theta_start
        self.theta_end = theta_end
        self.dist_start = dist_start
        self.dist_end = dist_end
        self.data = []   # stored info from .csv file
        self.doublet_data = []   # doublet list, first hit of the doublet is considered to be inside the segment
        self.triplet_data = []   # triplet list, first hit of the triplet is considered to be inside the segment

