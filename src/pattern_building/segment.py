class DetectorSegment:
    """Class for handling parts (segments) of the detector.
    """
    def __init__(self,
                 name: str,
                 layer: int,
                 x_start: float,
                 x_end: float,
                 y_start: float,
                 y_end: float,
                 z_start: float,
                 z_end: float):
        """Setting field according to the given geometry values
        :param name: unique identifier for the segment
        :param layer: using layer numbering 0-3 for simplified LUXE and 0-7 for full setup
        :param x_start: start of segment in x [m]
        :param x_end: end of segment in x [m]
        :param y_start: start of segment in y [m]
        :param y_end: end of segment in y [m]
        :param z_start: start of segment in z [m]
        :param z_end: end of segment in z [m]
        """
        self.name = name
        self.layer = layer
        self.x_start = x_start
        self.x_end = x_end
        self.y_start = y_start
        self.y_end = y_end
        self.z_start = z_start
        self.z_end = z_end

        # stored info from .csv file
        self.data = []

        # doublet list, first hit of the doublet is considered to be inside the segment
        self.doublet_data = []

        # triplet list, first hit of the triplet is considered to be inside the segment
        self.triplet_data = []
