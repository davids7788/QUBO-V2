class MuCoDetectorSegment:
    def __init__(self,
                 name: str,
                 phi_start: float,
                 phi_end: float,
                 z_start: float,
                 z_end: float,
                 r_start: float,
                 r_end: float):
        """Manages a segment of the Muon Collider detector layers.
        :param name: unique identifier for the segment
        :param phi_start: start of segment in phi [rad]
        :param phi_end: end of segment in phi [rad]
        :param z_start: start of segment in z [mm]
        :param z_end: end of segment in z [mm]
        :param r_start: start of segment in distance from IP [mm]
        :param r_end: end of segment in distance from IP [mm]
        """
        self.name = name
        self.phi_start = phi_start
        self.phi_end = phi_end
        self.z_start = z_start
        self.z_end = z_end
        self.r_start = r_start
        self.r_end = r_end

        self.data = []   # stored info from .csv file
        self.doublet_data = []   # doublet list, first hit of the doublet is considered to be inside the segment
        self.triplet_data = []   # triplet list, first hit of the triplet is considered to be inside the segment
