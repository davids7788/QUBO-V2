import argparse
import yaml
import os

from pathlib import Path
from pattern_building.MuCo_triplet_creator import MuCoTripletCreator
from pattern_building.MuCo_segment_manager import MuCoSegmentManager

parser = argparse.ArgumentParser(description='Muon Collider Pattern Builder',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--event_folder',
                    action='store',
                    type=str,
                    default=None,
                    help='Muon Collider hits')

parser.add_argument('--geometry_folder',
                    action='store',
                    type=str,
                    default=None,
                    help='Folder with Muon Collider geometry files')

parser.add_argument('--configuration',
                    action='store',
                    type=str,
                    default=None,
                    help='File with Muon Collider pattern building configuration')

parser.add_argument('--no_endcaps',
                    action='store_true',
                    help='Only barrel region used for tracking')


args = parser.parse_args()
event_folder = args.event_folder
geometry_folder = args.geometry_folder
with open(args.configuration, 'r') as f:
    configuration = yaml.safe_load(f)

no_endcaps = args.no_endcaps

new_folder = '/'.join([event_folder, args.configuration.split('/')[-1]])
if Path(new_folder).is_dir():
    pass
else:
    os.mkdir(new_folder)


s_manager = MuCoSegmentManager(configuration, geometry_folder)
s_manager.create_MuCo_segments()
s_manager.layer_mapping(no_endcaps=no_endcaps)

mu_co_creator = MuCoTripletCreator(event_folder)
mu_co_creator.load_tracking_data(dlfilter=False,
                                 s_manager=s_manager,
                                 no_endcaps=no_endcaps)

s_manager.create_segment_mapping(no_endcaps=no_endcaps)

mu_co_creator.create_xplet_list(s_manager,
                                configuration)

mu_co_creator.set_qubo_coefficients(s_manager,
                                    configuration,
                                    save_to_folder=new_folder)
mu_co_creator.write_info_file(save_to_folder=new_folder,
                              configuration=configuration)
