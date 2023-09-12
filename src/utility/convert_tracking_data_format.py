import numpy as np
from pattern.detector_hit import DetectorHit


def get_simplified_simulation_csv_format() -> list[str]:
    """Returns list of strings matching the key4hep_csv header
    :return
        header of key4hep_csv as list of strings
    """
    return ['hit_ID',
            'x',
            'y',
            'z',
            'layer_ID',
            'particle_ID',
            'particle_energy']


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


def from_simplified_simulation(simplified_sim_entry: list[str]) -> DetectorHit:
    """Converts row entry of a key4hep tracking data file
    :param simplified_sim_entry one row from a simplified simulation .csv tracking file

    :return
        DetectorHit object
    """
    # get structure of key4hep_csv
    fieldnames = get_simplified_simulation_csv_format()

    hit_dictionary = {}
    hit_dictionary.update({'hit_ID': simplified_sim_entry[fieldnames.index('hit_ID')]})
    hit_dictionary.update({'x': simplified_sim_entry[fieldnames.index('x')]})
    hit_dictionary.update({'y': simplified_sim_entry[fieldnames.index('y')]})
    hit_dictionary.update({'z': simplified_sim_entry[fieldnames.index('z')]})
    hit_dictionary.update({'layer_ID': -1})
    hit_dictionary.update({'module_ID': -1})
    hit_dictionary.update({'cell_ID': -1})
    hit_dictionary.update({'particle_ID': simplified_sim_entry[fieldnames.index('particle_ID')]})
    hit_dictionary.update({'is_signal': True})
    # TODO: provide more detailed info to match key4hep tracking files

    # create detector hit object
    hit = DetectorHit(hit_dictionary)

    return hit


def from_key4hep_csv(key4hep_csv_entry: list[str]) -> DetectorHit:
    """Converts row entry of a key4hep tracking data file
    :param key4hep_csv_entry: one row from a key4hep_csv tracking file

    :return
        DetectorHit object
    """
    # get structure of key4hep_csv
    fieldnames = get_key4hep_csv_format()

    hit_dictionary = {}
    hit_dictionary.update({'hit_ID': key4hep_csv_entry[fieldnames.index('index')]})
    hit_dictionary.update({'x': np.round(1e-3 * float(key4hep_csv_entry[fieldnames.index('tx')]), 3)})
    hit_dictionary.update({'y': np.round(1e-3 * float(key4hep_csv_entry[fieldnames.index('ty')]), 3)})
    hit_dictionary.update({'z': np.round(1e-3 * float(key4hep_csv_entry[fieldnames.index('tz')]), 3)})
    hit_dictionary.update({'layer_ID': key4hep_csv_entry[fieldnames.index('layer_id')]})
    hit_dictionary.update({'module_ID': key4hep_csv_entry[fieldnames.index('module_id')]})
    hit_dictionary.update({'cell_ID': int(get_cell_id_key4hep(key4hep_csv_entry), base=2)})
    hit_dictionary.update({'particle_ID': key4hep_csv_entry[fieldnames.index('particle_id')]})
    hit_dictionary.update({'is_signal': True})  # TODO: provide info of key4hep.csv file

    # create detector hit object
    hit = DetectorHit(hit_dictionary)

    return hit


def get_cell_id_key4hep(csv_entry: list[str]) -> str:
    """Returns the cell id, used for ACTS Kalman Filter for LUXE inside key4hep environment for.
    :param csv_entry: list of str with information about particle hit in detector

    :return
        concatenated string of module and stave to identify region where hit stems from
    """
    return f'{bin(int(csv_entry[get_key4hep_csv_format().index("module_id")])).zfill(3)[2:]}' \
           f'{bin(int(csv_entry[get_key4hep_csv_format().index("layer_id")])).zfill(3)[2:]}10'
