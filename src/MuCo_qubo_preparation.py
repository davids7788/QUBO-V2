import argparse
import yaml

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

parser.add_argument('--target_folder',
                    action='store',
                    type=str,
                    default=None,
                    help='Folder in which the triplet list will be stored')


args = parser.parse_args()
event_folder = args.event_folder
geometry_folder = args.geometry_folder
with open(args.configuration, 'r') as f:
    configuration = yaml.safe_load(f)
target_folder = args.target_folder

s_manager = MuCoSegmentManager(configuration, geometry_folder)
s_manager.create_MuCo_segments()
s_manager.layer_mapping()

mu_co_creator = MuCoTripletCreator()
mu_co_creator.load_tracking_data(event_folder, dlfilter=False, s_manager=s_manager)

s_manager.create_segment_mapping()

mu_co_creator.create_xplet_list(s_manager)
exit()
mu_co_creator.set_qubo_coefficients(target_folder)

