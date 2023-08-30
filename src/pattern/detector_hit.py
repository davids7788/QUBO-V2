class DetectorHit:
    """Class for storing information about particle-detector interactions.
    """
    def __init__(self,
                 detector_hit: list[str]):
        """Set fieldnames according to the information provided by the list of values
        :param detector_hit: list of values matching the definitions of the fieldnames provided in the
        fieldnames_matching method.
        """
        self.hit_id = DetectorHit.fieldnames_matching(detector_hit, 'hit_ID')
        self.x = DetectorHit.fieldnames_matching(detector_hit, 'x')
        self.y = DetectorHit.fieldnames_matching(detector_hit, 'y')
        self.z = DetectorHit.fieldnames_matching(detector_hit, 'z')
        self.layer_id = DetectorHit.fieldnames_matching(detector_hit, 'layer_ID')
        self.cell_id = DetectorHit.fieldnames_matching(detector_hit, 'cell_ID')
        self.is_signal_hit = DetectorHit.fieldnames_matching(detector_hit, 'is_signal')
        self.particle_id = DetectorHit.fieldnames_matching(detector_hit, 'particle_ID')
        self.particle_energy = DetectorHit.fieldnames_matching(detector_hit, 'particle_energy')
        self.time = DetectorHit.fieldnames_matching(detector_hit, 'time')

    @staticmethod
    def fieldnames_matching(detector_hit, value):
        """Returns the variable in the correct format
        """
        fieldnames = ['hit_ID',
                      'x',
                      'y',
                      'z',
                      'layer_ID',
                      'cell_id',
                      'particle_ID',
                      'is_signal',
                      'particle_energy',
                      'time']

        if value not in fieldnames:
            return None
        if value in ['hit_ID', 'layer_ID', 'is_signal']:
            return detector_hit[fieldnames.index(value)]
        elif value in ['particle_ID', 'cell_ID']:
            return int(detector_hit[fieldnames.index(value)])
        else:
            return float(detector_hit[fieldnames.index(value)])
