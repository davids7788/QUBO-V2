class DetectorHit:
    def __init__(self,
                 detector_hit: list[str]):
        """Stores information about a detector hit
        :param detector_hit: information from a <tracking>.csv file
        """
        self.hit_id = DetectorHit.fieldnames_matching(detector_hit, 'hit_ID')
        self.x = DetectorHit.fieldnames_matching(detector_hit, 'x')
        self.y = DetectorHit.fieldnames_matching(detector_hit, 'y')
        self.z = DetectorHit.fieldnames_matching(detector_hit, 'z')
        self.layer_id = DetectorHit.fieldnames_matching(detector_hit, 'layer_ID')
        self.particle_id = DetectorHit.fieldnames_matching(detector_hit, 'particle_ID')
        self.particle_energy = DetectorHit.fieldnames_matching(detector_hit, 'particle_energy')
        try:
            self.time = DetectorHit.fieldnames_matching(detector_hit, 'time')
        except IndexError:
            pass

    @staticmethod
    def fieldnames_matching(detector_hit, value):
        """Returns the variable in the correct format
        """
        fieldnames = ['hit_ID', 'x', 'y', 'z', 'layer_ID', 'particle_ID', 'particle_energy', 'time']
        if value in ['hit_ID', 'layer_ID']:
            return detector_hit[fieldnames.index(value)]
        elif value == 'particle_ID':
            return int(detector_hit[fieldnames.index(value)])
        else:
            return float(detector_hit[fieldnames.index(value)])
