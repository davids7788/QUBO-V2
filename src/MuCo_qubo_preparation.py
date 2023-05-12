import argparse
import yaml

from pattern_building.MuCo_triplet_creator import MuCoTripletCreator
from pattern_building.MuCo_segment_manager import MuCoSegmentManager

parser = argparse.ArgumentParser(description='Muon Collider Pattern Builder',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--Muon_Collider_Tracker_File',
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
tracking_file = args.Muon_Collider_Tracker_File
geometry_folder = args.geometry_folder
with open(args.configuration, 'r') as f:
    configuration = yaml.safe_load(f)
target_folder = args.target_folder
fieldnames = ['hit_ID', 'x', 'y', 'z', 'layer', 'MC_particle_ID', 'px', 'py', 'pz', 'time', 'PDG', 'event']

x = [-1.6092182473447714, -6.8953764319931565, -462.8257255902779]
y = [-30.165, -129.4685, 115.298]
z = [-9.247761474434203, -5.73165722903271, 1659.8515]
layer = ["0", "0", "4_right"]
file_name = ["Output_REC.000_VXD_0.csv", "Output_REC.000_ITracker_0.csv","Output_REC.000_ITrackerEndcap_0.csv"]


s_manager = MuCoSegmentManager(configuration, geometry_folder)
s_manager.create_MuCo_segments()
for i in range(3):
    segment = s_manager.get_segment_at_known_xyz_value(x[i],
                                                       y[i],
                                                       z[i],
                                                       file_name[i],
                                                       layer[i])
    print(f'phi_start: {segment.phi_start}')
    print(f'phi_end: {segment.phi_end}')
    print(f'r_start: {segment.r_start}')
    print(f'r_end: {segment.r_end}')
    print(f'z_start: {segment.z_start}')
    print(f'z_end: {segment.z_end}')
    print()
    print()
exit()
# s_manager.segment_mapping_MuCo

mu_co_creator = MuCoTripletCreator()
mu_co_creator.load_tracking_data(tracking_file)
mu_co_creator.create_xplet_list()
mu_co_creator.set_qubo_coefficients(target_folder)

