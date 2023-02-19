import csv
import numpy as np
import os

from pattern.doublet import Doublet
from pattern.triplet import Triplet
from pattern.x_plet import Xplet


def gen_xplets_simplified_LUXE(tracking_data_file: str,
                               save_to_folder: str) -> None:
    """Creates truth generated xplet list.
        :param tracking_data_file: .csv tracking data file
        :param save_to_folder: location to save the gen_xplet list
    """

    print("Creating Xplets on generator level with truth information...")
    output_name = ".".join(tracking_data_file.split("/")[-1].split(".")[0:-1])
    if os.path.isfile(output_name):
        return
    generated_x_plets = []

    with open(tracking_data_file, 'r') as file:
        csv_reader_tracking_file = csv.reader(file)
        csv_header_tracking_file = next(csv_reader_tracking_file)

        # preparing indices
        x_index = csv_header_tracking_file.index("x")
        y_index = csv_header_tracking_file.index("y")
        z_index = csv_header_tracking_file.index("z")
        hit_id_index = csv_header_tracking_file.index("hit_ID")
        particle_id_index = csv_header_tracking_file.index("particle_ID")
        particle_energy_index = csv_header_tracking_file.index("particle_energy")
        try:
            time = csv_header_tracking_file.index("time")
        except ValueError:
            time = -1

        def sort_row_by_z_value(row_entry: list[str]) -> float:
            return float(row_entry[z_index])

        def create_doublet(first_hit: list[str],
                           second_hit: list[str]):
            if not time:
                first_hit.append(time)
                second_hit.append(time)
            doublet = Doublet(int(first_hit[particle_id_index]),
                              int(second_hit[particle_id_index]),
                              (float(first_hit[x_index]),
                               float(first_hit[y_index]),
                               float(first_hit[z_index])),
                              (float(second_hit[x_index]),
                               float(second_hit[y_index]),
                               float(second_hit[z_index])),
                              first_hit[hit_id_index],
                              second_hit[hit_id_index],
                              float(first_hit[particle_energy_index]),
                              float(second_hit[particle_energy_index]),
                              float(first_hit[-1]),
                              float(second_hit[-1]))
            return doublet

        x_plet_dict = {}
        for row in csv_reader_tracking_file:
            if row[particle_id_index] not in x_plet_dict.keys():
                x_plet_dict.update({row[particle_id_index]: [row]})
            else:
                x_plet_dict[row[particle_id_index]].append(row)

        for x_plet_pieces in x_plet_dict.values():
            x_plet_pieces.sort(key=sort_row_by_z_value)
            doublet_list = []
            for j in range(len(x_plet_pieces[0:-1])):
                doublet_list.append(create_doublet(x_plet_pieces[j], x_plet_pieces[j + 1]))
            triplet_list = []
            for k in range(len(doublet_list) - 1):
                triplet_list.append(Triplet(doublet_list[k], doublet_list[k + 1]))
            truth_pattern = Xplet()
            for triplet in triplet_list:
                truth_pattern.add_triplet(triplet)
            if not truth_pattern.is_empty:
                generated_x_plets.append(truth_pattern)

        np.save(f"{save_to_folder}/{output_name}_gen_xplet_list", generated_x_plets)
        print(f"Number of generated Xplets: {len(generated_x_plets)}\n")
        print("-----------------------------------\n")
