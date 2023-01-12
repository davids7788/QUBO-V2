import csv
import numpy as np

with open("../../simplified_simulation_files/e0gpc_7.0/smeared/e0gpc_7.0_0000_sl.csv", 'r') as file:
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


    def fill_statistics(track_candidate):
        for entry1, entry2 in zip(track_candidate[0:-1], track_candidate[1:]):
            x0 = x0_at_z_ref(float(entry2[x_index]),
                             float(entry1[x_index]),
                             float(entry2[z_index]),
                             float(entry1[z_index]),
                             3.9560125)
            dx_x0_values.append((float(entry2[x_index]) - float(entry1[x_index])) / x0)
            dy.append(abs(float(entry2[y_index]) - float(entry1[y_index])) / x0)
        for entry1, entry2, entry3 in zip(track_candidate[0:-2], track_candidate[1:-1], track_candidate[2:]):
            xz_angle_d1 = np.arctan2((float(entry2[x_index]) - float(entry1[x_index])),
                                     (float(entry2[z_index]) - float(entry1[z_index])))
            xz_angle_d2 = np.arctan2((float(entry3[x_index]) - float(entry2[x_index])),
                                     (float(entry3[z_index]) - float(entry2[z_index])))
            yz_angle_d1 = np.arctan2((float(entry2[y_index]) - float(entry1[y_index])),
                                     (float(entry2[z_index]) - float(entry1[z_index])))
            yz_angle_d2 = np.arctan2((float(entry3[y_index]) - float(entry2[y_index])),
                                     (float(entry3[z_index]) - float(entry2[z_index])))

            t_angle.append(np.sqrt((xz_angle_d2 - xz_angle_d1)**2 + (yz_angle_d2 - yz_angle_d1)**2))

    def x0_at_z_ref(x_end: float,
                    x_start: float,
                    z_end: float,
                    z_start: float,
                    z_ref: float):
        dx = x_end - x_start
        dz = z_end - z_start
        return x_end - dx * abs(z_end - z_ref) / dz


    def sort_key(entry):
        return int(entry[particle_id_index])

    rows = []
    for row in csv_reader_tracking_file:
        rows.append(row)

    rows.sort(key=sort_key)
    dy = []
    dx_x0_values = []
    t_angle = []

    track = []
    track_id = None
    for i in range(len(rows)):
        if not track:
            track_id = rows[i][particle_id_index]
            track.append(rows[i])
        else:
            if rows[i][particle_id_index] != track_id:
                if len(track) < 4:
                    track = [rows[i]]
                    track_id = rows[i][particle_id_index]
                    continue
                elif len(track) == 4:
                    fill_statistics(track)

                    track = [rows[i]]
                    track_id = rows[i][particle_id_index]
                else:
                    print(len(track_id))
                    for t in track:
                        print(t)
                    exit()
            else:
                track.append(rows[i])
    print(np.std(dx_x0_values))
    print(np.around(np.mean(dx_x0_values), 6), np.around(np.std(dx_x0_values), 6))
    t_survived = 0
    for t in t_angle:
        if t < 0.001:
            t_survived += 1
    print(100 * np.around(t_survived / len(t_angle), 4))
    dy_survived = 0
    for j in dy:
        if j < 0.00085:
            dy_survived += 1
    print(100 * np.around(dy_survived / len(dy), 4))
    print(np.around(np.mean(dy), 6), np.around(np.std(dy), 6))
