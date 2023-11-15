import numpy as np
import csv
from pattern.detector_hit import DetectorHit

try:
    from pyLCIO import IOIMPL, UTIL, EVENT
except ImportError:
    print('pyLCIO not available in the environment!')


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
            'module_ID',
            'is_signal',
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


def load_data(tracking_data_file: str,
              tracking_data_format: str) -> list[DetectorHit]:
    """Handles the data loading and picks the correct function according to the data format
    
    :param tracking_data_file:    tracking dat file
    :param tracking data format:  format of the tracking data 
                                  --> 'key4hep slcio', 'key4hep csv', 'simplified simulation csv'
        
    """

    if tracking_data_format == 'simplified simulation csv':
        return load_tracking_data_from_simplified_simulation_csv(tracking_data_file)
    elif tracking_data_format == 'key4hep csv':
        return load_tracking_data_from_key4hep_csv(tracking_data_file)
    elif tracking_data_format == 'key4hep slcio':
        return load_tracking_data_from_slcio(tracking_data_file)


def load_tracking_data_from_simplified_simulation_csv(tracking_data_file: str) -> list[DetectorHit]:
    """Loads tracking data from .csv file and returns a list of DetectorHit objects created from the file entries.
    :param tracking_data_file: .csv tracking data file

    :return:
        list of DetectorHit objects
    """
    detector_hits = []
    with open(tracking_data_file, 'r') as file:
        csv_reader = csv.reader(file)
        _ = next(csv_reader)  # access header, csv files should consist of one line of header

        for row in csv_reader:
            hit = from_simplified_simulation(row)
            detector_hits.append(hit)

    return detector_hits


def load_tracking_data_from_key4hep_csv(tracking_data_file: str) -> list[DetectorHit]:
    """Loads tracking data from .csv file and returns a list of DetectorHit objects created from the file entries.
    :param tracking_data_file: .csv tracking data file

    :return:
        list of DetectorHit objects
    """
    detector_hits = []
    with open(tracking_data_file, 'r') as file:
        csv_reader = csv.reader(file)
        _ = next(csv_reader)  # access header, csv files should consist of one line of header

        for row in csv_reader:
            hit = from_key4hep_csv(row)
            detector_hits.append(hit)

    return detector_hits


def load_tracking_data_from_slcio(tracking_data_file: str) -> list[DetectorHit]:
    """Loads tracking data from .csv file and returns a list of DetectorHit objects created from the file entries.
    :param tracking_data_file: .csv tracking data file

    :return:
        list of DetectorHit objects
    """
    detector_hits = []
    
    # file reader
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(tracking_data_file)
    
    # there should only be one event in the file
    for i, event in enumerate(reader):
        
        tracker_hits = event.getCollection('SiTrackerHits')
        relation_collection = event.getCollection('SiTrackerHitRelations')
        relation = UTIL.LCRelationNavigator(relation_collection)
        
        encoding = tracker_hits.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder = UTIL.BitField64(encoding)
        
        for i, tracker_hit in enumerate(tracker_hits):
 
            if relation.getRelatedToObjects(tracker_hit)[0].isProducedBySecondary():
                rel_object = None
            else:
                rel_object = relation.getRelatedToObjects(tracker_hit)[0].getMCParticle()
            cell_id = int(tracker_hit.getCellID0())
            decoder.setValue(cell_id)
            layer = int(decoder['layer'].value())
            
            position = [float(1e-3 * tracker_hit.getPosition()[0]), 
                        float(1e-3 * tracker_hit.getPosition()[1]),
                        np.around(1e-3 * float(tracker_hit.getPosition()[2]), 4)]

 
            if rel_object:
                particle_id = id(rel_object)
                is_signal = True
            else:
                particle_id = None
                is_signal = False
                                           

            hit = from_slcio(i,
                             position,
                             layer,
                             particle_id,
                             cell_id,
                             is_signal)
            detector_hits.append(hit)

    return detector_hits


def from_slcio(hit_id: str,
               position: list[float, float, float],
               layer_id: int,
               particle_id: int | None,
               cell_id: int,
               is_signal: bool) -> DetectorHit:
    """Converts row entry of a key4hep slcio tracking data file
    
    :param simplified_sim_entry one row from a key4hep .slcio tracking file
    

    :return
        DetectorHit object
    """
    hit_dictionary = {}
    hit_dictionary.update({'hit_ID': hit_id})
    hit_dictionary.update({'x': position[0]})
    hit_dictionary.update({'y': position[1]})
    hit_dictionary.update({'z': position[2]})

    hit_dictionary.update({'layer_ID': layer_id})
    
    # not provided, set to -1
    hit_dictionary.update({'module_ID': -1})
    hit_dictionary.update({'is_signal': particle_id})

    hit_dictionary.update({'cell_ID': cell_id})
    hit_dictionary.update({'particle_ID': particle_id})
    hit_dictionary.update({'is_signal': True})

    # create detector hit object
    hit = DetectorHit(hit_dictionary)
    
    return hit


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

    # old versions don't have layer_ID and module_ID matching key4hep sample --> rework
    try:
        hit_dictionary.update({'layer_ID': simplified_sim_entry[fieldnames.index('layer_ID')]})
    except IndexError:
        hit_dictionary.update({'layer_ID': -1})
    try:
        hit_dictionary.update({'module_ID': simplified_sim_entry[fieldnames.index('module_ID')]})
    except IndexError:
        hit_dictionary.update({'module_ID': -1})
    hit_dictionary.update({'is_signal': simplified_sim_entry[fieldnames.index('is_signal')]})
    hit_dictionary.update({'cell_ID': get_cell_id_simplified_simulation(simplified_sim_entry)})
    hit_dictionary.update({'particle_ID': simplified_sim_entry[fieldnames.index('particle_ID')]})
    hit_dictionary.update({'is_signal': True})

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
    hit_dictionary.update({'x': np.around(1e-3 * float(key4hep_csv_entry[fieldnames.index('tx')]), 8)})
    hit_dictionary.update({'y': np.around(1e-3 * float(key4hep_csv_entry[fieldnames.index('ty')]), 8)})
    hit_dictionary.update({'z': np.around(1e-3 * float(key4hep_csv_entry[fieldnames.index('tz')]), 8)})
    hit_dictionary.update({'layer_ID': key4hep_csv_entry[fieldnames.index('layer_id')]})
    hit_dictionary.update({'module_ID': key4hep_csv_entry[fieldnames.index('module_id')]})
    hit_dictionary.update({'cell_ID': int(get_cell_id_key4hep(key4hep_csv_entry), base=2)})
    hit_dictionary.update({'particle_ID': key4hep_csv_entry[fieldnames.index('particle_id')]})
    hit_dictionary.update({'is_signal': True})

    # create detector hit object
    hit = DetectorHit(hit_dictionary)

    return hit


def get_cell_id_key4hep(csv_entry: list[str]) -> str:
    """Returns the cell id, used for ACTS Kalman Filter for LUXE inside key4hep environment for.
    :param csv_entry: list of str with information about particle hit in detector

    :return
        concatenated string of module and stave to identify region where hit stems from
    """
    return f'{bin(int(csv_entry[get_key4hep_csv_format().index("module_id")]))[2:].zfill(3)}' \
           f'{bin(int(csv_entry[get_key4hep_csv_format().index("layer_id")]))[2:].zfill(3)}10'


def get_cell_id_simplified_simulation(csv_entry: list[str]) -> str:
    """Returns the cell id, used for ACTS Kalman Filter for LUXE inside key4hep environment for.
    :param csv_entry: list of str with information about particle hit in detector

    :return
        concatenated string of module and stave to identify region where hit stems from
    """
    return f'{bin(int(csv_entry[get_simplified_simulation_csv_format().index("module_ID")]))[2:].zfill(3)}' \
           f'{bin(int(csv_entry[get_simplified_simulation_csv_format().index("layer_ID")]))[2:].zfill(3)}10'
