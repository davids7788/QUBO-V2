import sys
import csv
import numpy as np

from pattern.x_plet import Xplet
from pattern.doublet import Doublet
from pattern.triplet import Triplet

# sys argv[1]: tracking data file
# sys.argv[2]: geometry file
# sys.argv[3]: save folder


tracking_data_file = sys.argv[1]
geometry_file = sys.argv[2]
save_to_folder = sys.argv[3]
output_name = ".".join(tracking_data_file.split("/")[-1].split(".")[0:-1])

# extract information about last and first detector layer from .csv file
with open(geometry_file, 'r') as file:
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


generated_x_plets = []
rows = []

with open(tracking_data_file, 'r') as file:
    csv_reader_tracking_file = csv.reader(file)
    csv_header_tracking_file = next(csv_reader_tracking_file)

    # preparing indices
    x_index = csv_header_tracking_file.index("x")
    y_index = csv_header_tracking_file.index("y")
    z_index = csv_header_tracking_file.index("z")
    hit_id_index = csv_header_tracking_file.index("hit_ID")
    particle_id_index = csv_header_tracking_file.index("particle_ID")
    layer_id_index = csv_header_tracking_file.index("layer_ID")
    particle_energy_index = csv_header_tracking_file.index("particle_energy")
    
    for row in csv_reader_tracking_file:
        rows.append(row)

    def sort_row_by_particle_id(row_entry):
        return int(row_entry[particle_id_index])

    def sort_row_by_z_value(row_entry):
        return float(row_entry[z_index])

    rows.sort(key=sort_row_by_particle_id)
    rows.append([-1, -1, -1, -1, -1, -1, -1])
    

    x_plet_pieces = []
    current_particle_id = 0
    for entry in rows:
        if int(entry[particle_id_index]) == current_particle_id:
                x_plet_pieces.append(entry)

        else:
            x_plet_pieces.sort(key=sort_row_by_z_value)
            doublet_list = []
            for j in range(len(x_plet_pieces[0:-1])):
                doublet_list.append(Doublet(int(x_plet_pieces[j][particle_id_index]),
                                            int(x_plet_pieces[j + 1][particle_id_index]),
                                            (float(x_plet_pieces[j][x_index]),
                                             float(x_plet_pieces[j][y_index]),
                                             float(x_plet_pieces[j][z_index])),
                                            (float(x_plet_pieces[j + 1][x_index]),
                                             float(x_plet_pieces[j + 1][y_index]),
                                             float(x_plet_pieces[j + 1][z_index])),
                                            int(x_plet_pieces[j][hit_id_index]),
                                            int(x_plet_pieces[j + 1][hit_id_index]),
                                            float(x_plet_pieces[j][particle_energy_index]),
                                            float(x_plet_pieces[j + 1][particle_energy_index])))
            triplet_list = []
            for k in range(len(doublet_list) - 1):
                triplet_list.append(Triplet(doublet_list[k], doublet_list[k + 1], -1))
            truth_pattern = Xplet()
            for triplet in triplet_list:
                truth_pattern.add_triplet(triplet)
            if not truth_pattern.is_empty:
                generated_x_plets.append(truth_pattern)
            x_plet_pieces = [entry]
            current_particle_id = int(entry[particle_id_index])
     
    
    np.save(f"{save_to_folder}/{output_name}_gen_xplet_list", generated_x_plets)
    print(f"Number of generated xplets: {len(generated_x_plets)}\n")
