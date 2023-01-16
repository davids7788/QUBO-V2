import csv
import numpy as np
import os

from pattern.doublet import Doublet
from pattern.triplet import Triplet
from pattern.x_plet import Xplet


def gen_xplets(tracking_data_file, save_to_folder):
    """Creates truth generated xplet list.
        :param tracking_data_file: .csv tracking data file
        :param save_to_folder: location to save the gen_xplet list
    """

    print("Creating Xplets on generator level with truth information...")
    output_name = ".".join(tracking_data_file.split("/")[-1].split(".")[0:-1])
    if os.path.isfile(output_name):
        return
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
        # layer_id_index = csv_header_tracking_file.index("layer_ID")
        particle_energy_index = csv_header_tracking_file.index("particle_energy")
        # particle_time_index = csv_header_tracking_file.index("time")

        for row in csv_reader_tracking_file:
            rows.append(row)

        def sort_row_by_particle_id(row_entry):
            return int(row_entry[particle_id_index])

        def sort_row_by_z_value(row_entry):
            return float(row_entry[z_index])

        rows.sort(key=sort_row_by_particle_id)

        # bug fix, needs to be remade in a more sophisticated way...
        rows.append([-1, -1, -1, -1, -1, -1, -1, -1])

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
                                                x_plet_pieces[j][hit_id_index],
                                                x_plet_pieces[j + 1][hit_id_index],
                                                float(x_plet_pieces[j][particle_energy_index]),
                                                float(x_plet_pieces[j + 1][particle_energy_index])))
                triplet_list = []
                for k in range(len(doublet_list) - 1):
                    triplet_list.append(Triplet(doublet_list[k], doublet_list[k + 1]))
                truth_pattern = Xplet()
                for triplet in triplet_list:
                    truth_pattern.add_triplet(triplet)
                if not truth_pattern.is_empty:
                    generated_x_plets.append(truth_pattern)
                x_plet_pieces = [entry]
                current_particle_id = int(entry[particle_id_index])

        np.save(f"{save_to_folder}/{output_name}_gen_xplet_list", generated_x_plets)
        print(f"Number of generated Xplets: {len(generated_x_plets)}\n")
        print("-----------------------------------\n")
