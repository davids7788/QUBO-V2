class Segment:
    def __init__(self,
                 name: str,
                 layer: int,
                 x_start: float,
                 x_end: float,
                 z_position: float):
        """Manages segments of the detector layers
        :param name: string created from integer
        :param layer: 0-3 for simplified model, 0-7 for full setup
        :param x_start: start value in x of the segment [m]
        :param x_end: end value of in of the segment in x [m]
        :param z_position: position in z [m]
        """
        self.name = name
        self.layer = layer
        self.x_start = x_start
        self.x_end = x_end
        self.z_position = z_position
        self.data = []   # stored info from .csv file
        self.doublet_data = []   # doublet list, first hit of the doublet is considered to be inside the segment
        self.triplet_data = []   # triplet list, first hit of the triplet is considered to be inside the segment

    def is_in_segment(self,
                      x: str):
        """
        Checks if position x is inside the segment
        :param x: position in x [m]
        :return: True if x in segment, else False
        """
        if self.x_start <= x <= self.x_end:
            return True
        return False
