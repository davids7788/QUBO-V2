class DetectorHit:
    """Class for storing information about particle-detector interactions.
    """

    def __init__(self,
                 detector_hit: dict):
        """Set fields according to the information provided by the given dictionary
        :param detector_hit: dictionary with information about a detector hit
        """
        self.hit_dictionary = detector_hit

        self.hit_id = None
        self.x = None
        self.y = None
        self.z = None
        self.layer_id = None
        self.module_id = None
        self.cell_id = None
        self.particle_id = None
        self.is_signal = None
        self.particle_energy = None

        self.fill_fields()

    def fill_fields(self) -> None:
        """Fills the fields with the dictionary with information from the hit_dictionary.
        """

        for key, value in self.hit_dictionary.items():
            if key == 'hit_ID':
                self.hit_id = str(value)
            elif key == 'x':
                self.x = float(value)
            elif key == 'y':
                self.y = float(value)
            elif key == 'z':
                self.z = float(value)
            elif key == 'layer_ID':
                self.layer_id = int(value)
            elif key == 'module_ID':
                self.module_id = int(value)
            elif key == 'cell_ID':
                self.cell_id = int(value)
            elif key == 'particle_ID':
                if value is not None:
                    self.particle_id = int(value)
                else:
                    self.particle_id = None
            elif key == 'is_signal':
                self.is_signal = bool(value)
            elif key == 'particle_energy':
                self.particle_energy = float(value)
            else:
                print(f'key: {key} not found!')
                print('Exiting...')
                exit()
