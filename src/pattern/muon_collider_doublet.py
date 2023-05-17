class MuCoDoublet:
    def __init__(self,
                 hit_1_pdg: str = '',
                 hit_2_pdg: str = '',
                 hit_1_position: tuple[float, float, float] = (0.0, 0.0, 0.0),
                 hit_2_position: tuple[float, float, float] = (0.0, 0.0, 0.0),
                 hit_1_momentum: tuple[float, float, float] = (0.0, 0.0, 0.0),
                 hit_2_momentum: tuple[float, float, float] = (0.0, 0.0, 0.0),
                 hit_1_id: str = '',
                 hit_2_id: str = '',
                 time_1: float = -1e10,
                 time_2: float = -1e10):
        """Class for doublet objects, consisting of two hits on consecutive detector layers.
        :param hit_1_pdg: particle pdg (simulation only), default = ''
        :param hit_2_pdg: particle pdg number (simulation only), default = ''
        :param hit_1_position: particle position (x, y, z) of hit 1, default: (0.0, 0.0, 0.0)
        :param hit_2_position: particle position (x, y, z) of hit 2, default: (0.0, 0.0, 0.0)
        :param hit_1_position: particle position (x, y, z) of hit 1, default: (0.0, 0.0, 0.0)
        :param hit_2_position: particle position (x, y, z) of hit 2, default: (0.0, 0.0, 0.0)
        :param hit_1_id: unique integer hit id on detector of hit 1, default: ''
        :param hit_2_id: unique integer hit id on detector of hit 2, default: ''
        :param time_1: absolute time of the hit 1 [s], default: -1e10
        :param time_2: absolute time of the hit 2 [s], default: -1e10
        """
        self.hit_1_pdg = hit_1_pdg
        self.hit_2_pdg = hit_2_pdg
        self.hit_1_position = hit_1_position
        self.hit_2_position = hit_2_position
        self.hit_1_momentum = hit_1_momentum
        self.hit_2_momentum = hit_2_momentum
        self.hit_1_id = hit_1_id
        self.hit_2_id = hit_2_id
        self.time_1 = time_1
        self.time_2 = time_2

    def is_correct_match(self) -> bool:
        """Checks if doublet hits stem from the same particle (Note: momentum at vertex used for that).
        :return
            True if created from same particle, else False.
        """
        if self.hit_1_momentum == (-1000.0, -1000.0, -1000.0) or \
                self.hit_2_momentum == (-1000.0, -1000.0, -1000.0):
            return False

        if self.hit_1_momentum == self.hit_2_momentum:
            return True
        return False
