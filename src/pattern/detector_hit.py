class DetectorHit:
    """Class for storing information about particle-detector interactions.
    """
    fieldnames = ['hit_ID',
                  'x',
                  'y',
                  'z',
                  'layer_ID',
                  'module_id',
                  'cell_ID',
                  'particle_ID',
                  'particle_info',
                  'particle_energy',
                  'time']

    def __init__(self,
                 detector_hit: list[str]):
        """Set fields according to the information provided by the list of values
        :param detector_hit: list of values matching the definitions of the fieldnames provided in the
        fieldnames_matching method.
        """
        self.hit_id = DetectorHit.fieldnames_matching(detector_hit, 'hit_ID')
        self.x = DetectorHit.fieldnames_matching(detector_hit, 'x')
        self.y = DetectorHit.fieldnames_matching(detector_hit, 'y')
        self.z = DetectorHit.fieldnames_matching(detector_hit, 'z')
        self.layer_id = DetectorHit.fieldnames_matching(detector_hit, 'layer_ID')
        self.module_id = DetectorHit.fieldnames_matching(detector_hit, 'module_ID')
        self.cell_id = DetectorHit.fieldnames_matching(detector_hit, 'cell_ID')
        self.particle_id = DetectorHit.fieldnames_matching(detector_hit, 'particle_ID')
        self.particle_info = DetectorHit.fieldnames_matching(detector_hit, 'particle_info')
        self.particle_energy = DetectorHit.fieldnames_matching(detector_hit, 'particle_energy')
        try:
            self.time = DetectorHit.fieldnames_matching(detector_hit, 'time')
        except IndexError:
            self.time = 0

    @staticmethod
    def fieldnames_matching(detector_hit, info):
        """Returns the variable in the correct format
        :param detector_hit: list of str values with information about the detector hit
        :param info: specify value that will be set

        :return:
            value in specified datatype
        """
        if info not in DetectorHit.fieldnames:
            return None
        if info in ['hit_ID', 'layer_ID', 'module_id', 'particle_info']:
            return detector_hit[DetectorHit.fieldnames.index(info)]
        elif info in ['particle_ID', 'cell_ID']:
            return int(detector_hit[DetectorHit.fieldnames.index(info)])
        else:
            return float(detector_hit[DetectorHit.fieldnames.index(info)])
