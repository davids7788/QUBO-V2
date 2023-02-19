from numba import jit


class LUXEDetectorSegment:
    def __init__(self,
                 name: str,
                 layer: int,
                 x_start: float,
                 x_end: float,
                 y_start: float,
                 y_end: float,
                 z_position: float):
        """Manages segments of the detector layers
        :param name: unique identifier for the segment
        :param layer: using layer numbering 0-3 for simplified LUXE and 0-7 for full setup
        :param x_start: start of segment in x [m]
        :param x_end: end of segment in x [m]
        :param y_start: start of segment in y [m]
        :param y_end: end of segment in y [m]
        :param z_position: position in z [m]
        """
        self.name = name
        self.layer = layer
        self.x_start = x_start
        self.x_end = x_end
        self.y_start = y_start
        self.y_end = y_end
        self.z_position = z_position
        self.data = []   # stored info from .csv file
        self.doublet_data = []   # doublet list, first hit of the doublet is considered to be inside the segment
        self.triplet_data = []   # triplet list, first hit of the triplet is considered to be inside the segment

    @jit(nopython=True)
    def is_in_segment(self,
                      x: float,
                      y: float,
                      z: float) -> bool:
        """Checks if position x is inside the segment
        :param x: position in x [m]
        :param y: position in y [m]
        :param z: position in z [m]
        :return:
            True if coordinate in segment, else False
        """
        if self.z_position == z:
            if self.x_start <= x <= self.x_end:
                if self.y_start <= y <= self.y_end:
                    return True
        return False
