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
triplet_list = np.load(f"{triplet_list_folder}/triplet_list.npy", allow_pickle=True)

computed_solution = log_file[()]["computed solution vector"]


def triplet_start(t):
    """Returns z value of first triplet hit
    """
    return t.doublet_1.hit_1_position[2]


kept_triplets = []
for solution, triplet in zip(computed_solution, triplet_list):
    if solution == 1:
        kept_triplets.append(triplet)
kept_triplets.sort(key=triplet_start)


reco_x_plets = []
for i, t1 in enumerate(kept_triplets):
    if t1.doublet_1.hit_1_position[2] in first_layer:
        t_list = [t1]
        for t2 in kept_triplets[i + 1:]:
            if t_list[-1].doublet_2 == t2.doublet_1:
                t_list.append(t2)

        reco_pattern = Xplet()
        for triplet in t_list:
            reco_pattern.add_triplet(triplet)
        reco_x_plets.append(reco_pattern)

np.save(f"{folder}/reco_xplet_list", reco_x_plets)
