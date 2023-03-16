import argparse
from pattern_building.MuCo_triplet_creator import MuCoTripletCreator

parser = argparse.ArgumentParser(description='Muon Collider Pattern Builder',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--Muon_Collider_Tracker_File',
                    action='store',
                    type=str,
                    default=None,
                    help='Muon Collider hits')

parser.add_argument('--target_folder',
                    action='store',
                    type=str,
                    default=None,
                    help='Folder in which the triplet list will be stored')


args = parser.parse_args()
tracking_file = args.Muon_Collider_Tracker_File
target_folder = args.target_folder

mu_co_creator = MuCoTripletCreator()
mu_co_creator.load_tracking_data(tracking_file)
mu_co_creator.create_xplet_list()
mu_co_creator.set_qubo_coefficients(target_folder)

