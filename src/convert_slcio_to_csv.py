import csv
import argparse
import time
import os

from pathlib import Path
from pyLCIO import IOIMPL, UTIL, EVENT

parser = argparse.ArgumentParser(description='convert .slcio to .csv information',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--slcio_file',
                    action='store',
                    type=str,
                    default=None,
                    help='Marlin created reco .slcio file')

parser.add_argument('--events',
                    action='store',
                    type=int,
                    default=0,
                    help='Number of events to store')

parser.add_argument('--endcaps',
                    action='store_true',
                    help='Hits in endcap regions')

parser.add_argument('--VXDTrackerHits',
                    action='store_true',
                    help='Hits in the VXD')

parser.add_argument('--InnerTrackerBarrelHits',
                    action='store_true',
                    help='Hits in the inner tracker barrel')

parser.add_argument('--OuterTrackerBarrelHits',
                    action='store_true',
                    help='Hits in the outer tracker barrel')

parser.add_argument('--DLFilteredHits',
                    action='store_true',
                    help='DoubleLayerFilteredHits')

parser.add_argument('--all',
                    action='store_true',
                    help='Hits in the VXD, inner and outer barrel')




args = parser.parse_args()

slcio_file = args.slcio_file
events = args.events
vxd_hits = args.VXDTrackerHits
endcaps = args.endcaps
inner_tracker_hits = args.InnerTrackerBarrelHits
outer_tracker_hits = args.OuterTrackerBarrelHits
dl_filtered_hits = args.DLFilteredHits
all_hits = args.all
if all_hits:
    vxd_hits = True
    inner_tracker_hits = True
    outer_tracker_hit = True
    dl_filtered_hits = True

fieldnames = ['hit_ID', 'x', 'y', 'z', 'layer', 'px', 'py', 'pz', 'time', 'PDG', 'event']

dl_filtered = []
vxd_tracker_barrel = []
vxd_tracker_endcap = []
inner_tracker_barrel = []
inner_tracker_endcap = []
outer_tracker_barrel = []
outer_tracker_endcap = []

reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(slcio_file)


def write_to_array(tracker_hits,
                   relation, 
                   decoder,
                   event_number,
                   endcap,
                   array):
    """Writes information about a hit into the currently opened .csv file.
    """
    hit_id = 0
    for tracker_hit in tracker_hits:
        rel_object = relation.getRelatedToObjects(tracker_hit)[0].getMCParticle()
        cellID = int(tracker_hit.getCellID0())
        decoder.setValue(cellID)
        if endcap:
            if tracker_hit.getPosition()[2] < 0:
                layer = f'{decoder["layer"].value()}_left'
            else:
                layer = f'{decoder["layer"].value()}_right'
        else:
            layer = decoder['layer'].value()
        try: 
            px = rel_object.getMomentum()[0]
            py = rel_object.getMomentum()[1]
            pz = rel_object.getMomentum()[2]
            pdg = rel_object.getPDG()
        except ReferenceError:
            px = -1000
            py = -1000
            pz = -1000
            pdg = -1000

        array.append({'hit_ID': hit_id,
                      'x': tracker_hit.getPosition()[0],
                      'y': tracker_hit.getPosition()[1],
                      'z': tracker_hit.getPosition()[2],
                      'layer': layer,
                      'px': px,
                      'py': py,
                      'pz': pz,
                      'time': tracker_hit.getTime(),
                      'PDG': pdg,
                      'event': i})
        hit_id += 1
                    

def get_hit_and_relation_information(event, 
                                     collection_name, 
                                     relation_collection_name):
    """Returns information about the hit and the relation to the particle causing the hit.
    """
    relation_collection = event.getCollection(relation_collection_name)
    relation = UTIL.LCRelationNavigator(relation_collection)
    tracker_hits = event.getCollection(collection_name)
    encoding = tracker_hits.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
    decoder = UTIL.BitField64(encoding)
    return tracker_hits, relation, decoder

for i, event in enumerate(reader):
    if Path(f'event_{i}').is_dir():
        pass
    else:
        os.mkdir(f'event_{i}')
    if i == events:
        break
    
    # DLFiltered
    if dl_filtered_hits:
        print(f'Processing event: {event.getEventNumber()}', end='\r')

        tracker_hits, relation, decoder = get_hit_and_relation_information(event, 
                                                                           "VXDTrackerHits_DLFiltered",      
                                                                           "VXDTrackerHitRelations")
        write_to_array(tracker_hits, relation, decoder, i, False, dl_filtered)
        with open(f"event_{i}/{'.'.join(slcio_file.split('.')[0:-1])}_DLFiltered_{event.getEventNumber()}.csv", 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for line in dl_filtered:
                writer.writerow(line)
        dl_filtered = []
        
    # VXDHits
    if vxd_hits:
        tracker_hits, relation, decoder = get_hit_and_relation_information(event, 
                                                                           "VXDTrackerHits",      
                                                                           "VXDTrackerHitRelations")
        write_to_array(tracker_hits, relation, decoder, i, False, vxd_tracker_barrel)
        with open(f"event_{i}/{'.'.join(slcio_file.split('.')[0:-1])}_VXDTracker_{event.getEventNumber()}.csv", 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for line in vxd_tracker_barrel:
                writer.writerow(line) 
        vxd_tracker_barrel = []
        
        if endcaps:
            tracker_hits, relation, decoder = get_hit_and_relation_information(event, 
                                                                               "VXDEndcapTrackerHits",      
                                                                               "VXDEndcapTrackerHitRelations")
            write_to_array(tracker_hits, relation, decoder, i, True, vxd_tracker_endcap)
            with open(f"event_{i}/{'.'.join(slcio_file.split('.')[0:-1])}_VXDTrackerEndcap_{event.getEventNumber()}.csv", 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                for line in vxd_tracker_endcap:
                    writer.writerow(line)  
            vxd_tracker_endcap = []
            
    # InnerTrackerHits
    if inner_tracker_hits:
        tracker_hits, relation, decoder = get_hit_and_relation_information(event, 
                                                                           "ITrackerHits", 
                                                                           "InnerTrackerBarrelHitsRelations")
        write_to_array(tracker_hits, relation, decoder, i, False, inner_tracker_barrel)
        with open(f"event_{i}/{'.'.join(slcio_file.split('.')[0:-1])}_ITracker_{event.getEventNumber()}.csv", 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for line in inner_tracker_barrel:
                writer.writerow(line)  
        inner_tracker_barrel = []
        
        if endcaps:
            tracker_hits, relation, decoder = get_hit_and_relation_information(event, 
                                                                               "ITrackerEndcapHits", 
                                                                               "InnerTrackerEndcapHitsRelations")
            write_to_array(tracker_hits, relation, decoder, i, True, inner_tracker_endcap)
            with open(f"event_{i}/{'.'.join(slcio_file.split('.')[0:-1])}_ITrackerEndcap_{event.getEventNumber()}.csv", 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                for line in inner_tracker_endcap:
                    writer.writerow(line)  
            inner_tracker_endcap = []
            
    # OuterTrackerHits
    if inner_tracker_hits:
        tracker_hits, relation, decoder = get_hit_and_relation_information(event, 
                                                                           "OTrackerHits", 
                                                                           "OuterTrackerBarrelHitsRelations")
        write_to_array(tracker_hits, relation, decoder, i, False, outer_tracker_barrel)
        with open(f"event_{i}/{'.'.join(slcio_file.split('.')[0:-1])}_OTracker_{event.getEventNumber()}.csv", 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for line in outer_tracker_barrel:
                writer.writerow(line)  
        outer_tracker_barrel = []
        
        if endcaps:
            tracker_hits, relation, decoder = get_hit_and_relation_information(event, 
                                                                               "OTrackerEndcapHits", 
                                                                               "OuterTrackerEndcapHitsRelations")
            write_to_array(tracker_hits, relation, decoder, i, True, outer_tracker_endcap)  
            with open(f"event_{i}/{'.'.join(slcio_file.split('.')[0:-1])}_OTrackerEndcap_{event.getEventNumber()}.csv", 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                for line in outer_tracker_endcap:
                    writer.writerow(line)
            outer_tracker_endcap = []

reader.close()
   
print('\n\nConversion done successfully!')