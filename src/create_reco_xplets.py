import sys
import numpy as np
import csv

from pattern.x_plet import Xplet

# sys argv[1]: qubo folder
# sys argv[2]: geometry file

with open(sys.argv[2], 'r') as file:
    csv_reader_geometry_file = csv.reader(file)
    csv_header_geometry_file = next(csv_reader_geometry_file)
    z = csv_header_geometry_file.index("z")
    z_values = []
    for row in csv_reader_geometry_file:
        z_values.append(float(row[z]))
    z_values_unique = list(set(z_values))
    z_values_unique.sort()
    if len(z_values_unique) == 4:
        last_layer = [z_values_unique[-1], z_values_unique[-1]]
        first_layer = [z_values_unique[0], z_values_unique[0]]
    else:
        last_layer = [z_values_unique[-2], z_values_unique[-1]]
        first_layer = [z_values_unique[0], z_values_unique[1]]

folder = sys.argv[1]
triplet_list_folder = "/".join(folder.split("/")[0:-1])

log_file = np.load(f"{folder}/qubo_log.npy", allow_pickle=True)
triplet_list = np.load(f"{triplet_list_folder}/triplet_list.npy",
                       allow_pickle=True)

computed_solution = log_file[()]["computed solution vector"]


def triplet_start(t):
    """Returns z value of first triplet hit
    """
    return t.doublet_1.hit_1_position[2]


reco_x_plets = []
for triplet_1_index, triplet_1 in enumerate(triplet_list):
    if computed_solution[triplet_1_index] == 0 or \
            triplet_start(triplet_1) not in first_layer:
        continue
    for triplet_2_index, connectivity in triplet_1.interactions.items():
        if (connectivity < 0 and computed_solution[triplet_2_index] == 1 and
                triplet_2_index > triplet_1_index):
            if triplet_start(triplet_list[triplet_2_index]) in first_layer:
                continue
            x_plet = Xplet()
            x_plet.add_triplet(triplet_1)
            x_plet.add_triplet(triplet_list[triplet_2_index])
            reco_x_plets.append(x_plet)
np.save(f"{folder}/reco_xplet_list", reco_x_plets)
