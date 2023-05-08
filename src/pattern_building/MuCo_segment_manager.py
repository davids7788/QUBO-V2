from pattern_building.MuCo_segment import MuCoDetectorSegment

class MuCoSegmentManager:
    def __init__(self,
                 configuration):
        """Class for handling segements for doublet and triplet creation for Muon Collider data.
        :param configuration: pattern building configuration file
        """
        self.mapping_criteria = configuration['doublet']
        self.binning = [configuration['binning']['phi'], configuration['binning']['theta']]
        self.segment_storage = {}  # organising segment objects, subdicts ordered by construction in the following way:
        # <layer> (key):
        # [<segment_phi_0_theta0>, <segment_phi_0_theta_1>, ..., <segment_phi0_theta_n>
        #  <segment_phi_1_theta0>, <segment_phi_1_theta_1>, ...,
        #  ..., ...                                            , <segment_phi_m_theta_n>]
        self.segment_mapping = {}  # <segment> (key): [<target_segment_0>, <target_segment_1>, ...] (value)

    def create_MuCo_segments(self) -> None:
        """Segments are created according to their phi, theta and distance from IP coordinates.
        The name of the segments gives information about their position and layer.
        """
        for angle_phi in range(self.binning[0]):
            for angle_theta in range(self.binning[1]):
                new_segment = MuCoDetectorSegment(f'Sphi{angle_phi}_Stheta{angle_theta}',
                                                  0 + angle_phi * 360 / self.binning[0],
                                                  0 + (angle_phi + 1) * 360 / self.binning[0],
                                                  0 + angle_theta * 360 / self.binning[1],
                                                  0 + (angle_theta + 1) * 360 / self.binning[1],
                                                  cd



