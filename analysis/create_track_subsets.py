import csv
from scipy.stats import binned_statistic_2d
import numpy as np

for slicing in range(10):
    tracking_data_file = \
        f"/nfs/dust/luxe/user/spatarod/towards_paper/e-laser/phase-0/gpc/7.0/smeared/e0gpc_7.0_000{slicing}_sl.csv"

    file_prefix = "/nfs/dust/luxe/user/spatarod/towards_paper/e-laser/phase-0/gpc"
    outfile_folder = "/".join([file_prefix,
                               tracking_data_file.split("/")[-3],
                               tracking_data_file.split("/")[-2]]) + "_sliced"
    file_split = tracking_data_file.split("/")[-1].split(".")

    t_x = []
    t_y = []
    p_numbers = []
    rows = []

    x_min = 0.05275912
    x_max = 0.55430088
    y_min = - 0.00688128
    y_max = 0.00688128
    z1 = 3.9560125
    z2 = 4.0560125
    z3 = 4.1560125
    z4 = 4.2560125

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
            if float(row[z_index]) == z1:
                t_x.append(float(row[x_index]))
                t_y.append(float(row[y_index]))
                p_numbers.append(row[particle_energy_index])

    x_bins = 18 * 32
    y_bins = 16
    detector_stats = binned_statistic_2d(t_x, t_y, None, "count",
                                         bins=[x_bins, y_bins],
                                         range=[[x_min, x_max], [y_min, y_max]], expand_binnumbers=True)
    print(max([max(s) for s in detector_stats.statistic]))
    # exit()
    max_occupancy = max([max(s) for s in detector_stats.statistic])
    file_name = file_split[0] + "." + file_split[1] + "_" + str(int(max_occupancy)) + "_tracks." + file_split[2]

    bin_number = np.where(detector_stats.statistic == max_occupancy)

    if len(bin_number[0]) > 1:
        a = bin_number[0]
        b = bin_number[1]
        bin_number = [a[0], b[0]]

    number_list = []

    for i, b_number_x, b_number_y in zip([0 + i for i in range(len(detector_stats.binnumber[0]))],
                                         detector_stats.binnumber[0],
                                         detector_stats.binnumber[1]):
        if b_number_x == bin_number[0] + 1 and b_number_y == bin_number[1] + 1:
            number_list.append(i)

    z1_x_range = []
    z1_y_range = []
    z2_x_range = []
    z2_y_range = []
    z3_x_range = []
    z3_y_range = []
    z4_x_range = []
    z4_y_range = []

    with open(outfile_folder + "/" + file_name, 'w', newline='') as csvfile:
        fieldnames = ['hit_ID', 'x', 'y', 'z', 'layer_ID', 'particle_ID', 'particle_energy']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in rows:
            row[particle_id_index]
            if int(row[particle_id_index]) in number_list:
                writer.writerow({'hit_ID': row[hit_id_index],
                                 'x': row[x_index],
                                 'y': row[y_index],
                                 'z': row[z_index],
                                 'layer_ID': row[layer_id_index],
                                 'particle_ID': row[particle_id_index],
                                 'particle_energy': row[particle_energy_index]})
            if row[z_index] == z1:
                if row[y_index] < z1_y_range[0]:
                    z1_y_range[0] = row[y_index]
                if row[y_index] > z1_y_range[1]:
                    z1_y_range[1] = row[y_index]
                if row[x_index] < z1_x_range[0]:
                    z1_x_range[0] = row[x_index]
                if row[x_index] > z1_x_range[1]:
                    z1_x_range[1] = row[x_index]

            if row[z_index] == z2:
                if row[y_index] < z2_y_range[0]:
                    z2_y_range[0] = row[y_index]
                if row[y_index] > z2_y_range[1]:
                    z2_y_range[1] = row[y_index]
                if row[x_index] < z2_x_range[0]:
                    z2_x_range[0] = row[x_index]
                if row[x_index] > z2_x_range[1]:
                    z2_x_range[1] = row[x_index]

            if row[z_index] == z3:
                if row[y_index] < z3_y_range[0]:
                    z3_y_range[0] = row[y_index]
                if row[y_index] > z3_y_range[1]:
                    z3_y_range[1] = row[y_index]
                if row[x_index] < z3_x_range[0]:
                    z3_x_range[0] = row[x_index]
                if row[x_index] > z3_x_range[1]:
                    z3_x_range[1] = row[x_index]

            if row[z_index] == z4:
                if row[y_index] < z4_y_range[0]:
                    z4_y_range[0] = row[y_index]
                if row[y_index] > z4_y_range[1]:
                    z4_y_range[1] = row[y_index]
                if row[x_index] < z4_x_range[0]:
                    z4_x_range[0] = row[x_index]
                if row[x_index] > z4_x_range[1]:
                    z4_x_range[1] = row[x_index]

        for row in rows:
            if int(row[particle_id_index]) not in number_list:
                if row[z_index] == z1 and z1_x_range[0] < row[x_index] < z1_x_range[1]\
                        and z1_y_range[0] < row[y_index] < z1_y_range[1]:
                    writer.writerow({'hit_ID': row[hit_id_index],
                                     'x': row[x_index],
                                     'y': row[y_index],
                                     'z': row[z_index],
                                     'layer_ID': row[layer_id_index],
                                     'particle_ID': row[particle_id_index],
                                     'particle_energy': row[particle_energy_index]})
                if row[z_index] == z2 and z2_x_range[0] < row[x_index] < z2_x_range[1]\
                        and z2_y_range[0] < row[y_index] < z2_y_range[1]:
                    writer.writerow({'hit_ID': row[hit_id_index],
                                     'x': row[x_index],
                                     'y': row[y_index],
                                     'z': row[z_index],
                                     'layer_ID': row[layer_id_index],
                                     'particle_ID': row[particle_id_index],
                                     'particle_energy': row[particle_energy_index]})
                if row[z_index] == z3 and z3_x_range[0] < row[x_index] < z3_x_range[1]\
                        and z3_y_range[0] < row[y_index] < z3_y_range[1]:
                    writer.writerow({'hit_ID': row[hit_id_index],
                                     'x': row[x_index],
                                     'y': row[y_index],
                                     'z': row[z_index],
                                     'layer_ID': row[layer_id_index],
                                     'particle_ID': row[particle_id_index],
                                     'particle_energy': row[particle_energy_index]})
                if row[z_index] == z4 and z4_x_range[0] < row[x_index] < z4_x_range[1]\
                        and z4_y_range[0] < row[y_index] < z4_y_range[1]:
                    writer.writerow({'hit_ID': row[hit_id_index],
                                     'x': row[x_index],
                                     'y': row[y_index],
                                     'z': row[z_index],
                                     'layer_ID': row[layer_id_index],
                                     'particle_ID': row[particle_id_index],
                                     'particle_energy': row[particle_energy_index]})
