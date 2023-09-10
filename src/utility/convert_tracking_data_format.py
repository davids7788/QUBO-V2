import numpy as np


def get_key4hep_csv_format() -> list[str]:
    """Returns list of strings matching the key4hep_csv header
    :return
        header of key4hep_csv as list of strings
    """
    return ['particle_id',
            'geometry_id',
            'system_id',
            'layer_id',
            'side_id',
            'module_id',
            'sensor_id',
            'tx',
            'ty',
            'tz',
            'tt',
            'tpx',
            'tpy',
            'tpz',
            'te',
            'deltapx',
            'deltapy',
            'deltapz',
            'deltae',
            'index']


def from_key4hep_csv(key4hep_csv_entry: list[str],
                     key4hep_csv_format: list[str]) -> list[str]:
    """Converts row entry of a key4hep tracking data file
    :param key4hep_csv_entry: one row from a key4hep_csv tracking file
    :param key4hep_csv_format: list of strings representing the header of the key4hep_csv format

    :return
        list of strings rearranged to match the DetectorHit class format
    """

    def get_cell_id(csv_entry: list[str]) -> str:
        """Returns the cell id, used for ACTS Kalman Filter for LUXE inside key4hep environment for.
        :param csv_entry: list of str with information about particle hit in detector

        :return
            concatenated string of module and stave to identify region where hit stems from
        """
        return f'{bin(int(csv_entry[key4hep_csv_format.index("module_id")])).zfill(3)[2:]}' \
               f'{bin(int(csv_entry[key4hep_csv_format.index("layer_id")])).zfill(3)[2:]}10'

    row_trimmed = [key4hep_csv_entry[key4hep_csv_format.index('index')],                          # hit_id
                   np.round(1e-3 * float(key4hep_csv_entry[key4hep_csv_format.index('tx')]), 3),  # x-coordinate
                   np.round(1e-3 * float(key4hep_csv_entry[key4hep_csv_format.index('ty')]), 3),  # y-coordinate
                   np.round(1e-3 * float(key4hep_csv_entry[key4hep_csv_format.index('tz')]), 3),  # z-coordinate
                   key4hep_csv_entry[key4hep_csv_format.index('layer_id')],                       # layer_id
                   key4hep_csv_entry[key4hep_csv_format.index('module_id')],                      # module_id
                   int(get_cell_id(key4hep_csv_entry), base=2),                        # cell_id for ACTS tracking
                   key4hep_csv_entry[key4hep_csv_format.index('particle_id')],         # particle_id
                   True,                                                               # not provided in key4hep csv
                   key4hep_csv_entry[key4hep_csv_format.index('te')],                  # energy of particle
                   key4hep_csv_entry[key4hep_csv_format.index('tt')]]                  # time of particle
    return row_trimmed
