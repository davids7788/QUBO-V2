import numpy as np
import csv
import os
import matplotlib.pyplot as plt
from numpy.linalg import norm


def convert_to_csv_true(path_to_file):
    """Converts the .npy file into a .csv tracking file containing truth information
    :param path_to_file: path to .npy file
    """
    folder = "/".join(path_to_file.split("/")[0:-1])
    file = path_to_file.split("/")[-1]
    if not os.path.isdir(f"{folder}/true"):
        os.mkdir(f"{folder}/true")

    file_name = None
    if "sl" in file:
        file_name = '_'.join(file.split("_")[0:-2]) + "_sl"
    elif "fl" in file:
        file_name = '_'.join(file.split("_")[0:-2]) + "_fl"
    else:
        print('Check LUXE geometry file!')
        print('Exiting...')
        exit()

    data = np.load(f"{path_to_file}.npy", allow_pickle=True)[()]

    hit_id = 0
    with open(f"{folder}/true/{file_name}.csv", 'w', newline='') as csvfile:
        fieldnames = ['hit_ID', 'x', 'y', 'z', 'layer_ID', 'module_ID', 'is_signal', 'particle_ID', 'particle_energy']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for key, plane in zip(data['Plane list'].keys(), data['Plane list'].values()):
            for particle, hit in zip(plane.true_hits_dictionary.keys(), plane.true_hits_dictionary.values()):
                layer, module = get_layer_and_module(hit[0], hit[1], hit[2])
                writer.writerow({'hit_ID': hit_id,
                                 'x': hit[0],
                                 'y': hit[1],
                                 'z': hit[2],
                                 'layer_ID': layer,
                                 'module_ID': module,
                                 'is_signal': True,
                                 'particle_ID': particle,
                                 'particle_energy': norm(data["Particle momentum history log"][particle][key])})
                hit_id += 1


def convert_to_csv_smeared(path_to_file):
    """Converts the .npy file into a .csv tracking file containing smeared hits
    :param path_to_file: path to .npy file
    """
    folder = "/".join(path_to_file.split("/")[0:-1])
    file = path_to_file.split("/")[-1]
    if not os.path.isdir(f"{folder}/smeared"):
        os.mkdir(f"{folder}/smeared")

    file_name = None
    if "sl" in file:
        file_name = '_'.join(file.split("_")[0:-2]) + "_sl"
    elif "fl" in file:
        file_name = '_'.join(file.split("_")[0:-2]) + "_fl"
    else:
        print('Check LUXE geometry file!')
        print('Exiting...')
        exit()

    data = np.load(f"{path_to_file}.npy", allow_pickle=True)[()]

    hit_id = 0
    with open(f"{folder}/smeared/{file_name}.csv", 'w', newline='') as csvfile:
        fieldnames = ['hit_ID', 'x', 'y', 'z', 'layer_ID', 'module_ID', 'is_signal', 'particle_ID', 'particle_energy']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for key, plane in zip(data['Plane list'].keys(), data['Plane list'].values()):
            # sigma_x_resolution = (plane.limits_x[1] - plane.limits_x[0]) / (plane.num_bins[0]) / np.sqrt(12)
            # sigma_y_resolution = (plane.limits_y[1] - plane.limits_y[0]) / (plane.num_bins[1]) / np.sqrt(12)
            sigma_x_resolution = 5e-6
            sigma_y_resolution = 5e-6
            for particle, hit in zip(plane.true_hits_dictionary.keys(), plane.true_hits_dictionary.values()):
                smeared_x = np.random.normal(hit[0], sigma_x_resolution)
                smeared_y = np.random.normal(hit[1], sigma_y_resolution)
                layer, module = get_layer_and_module(smeared_x, smeared_y, hit[2])
                writer.writerow({'hit_ID': hit_id,
                                 'x': smeared_x,
                                 'y': smeared_y,
                                 'z': hit[2],
                                 'layer_ID': layer,
                                 'module_ID': module,
                                 'is_signal': True,
                                 'particle_ID': particle,
                                 'particle_energy': norm(data["Particle momentum history log"][particle][key])})
                hit_id += 1


def convert_to_csv_pixel(path_to_file):
    """Converts the .npy file into a .csv tracking file containing pixelated hits, the mid of a pixel is provided
    :param path_to_file: path to .npy file
    """
    folder = "/".join(path_to_file.split("/")[0:-1])
    file = path_to_file.split("/")[-1]

    if not os.path.isdir(f"{folder}/occupancy"):
        os.mkdir(f"{folder}/occupancy")
    # if not os.path.isdir(f"{folder}/train_images_NN"):
    #     os.mkdir(f"{folder}/train_images_NN")
    if not os.path.isdir(f"{folder}/pixel"):
        os.mkdir(f"{folder}/pixel")

    file_name = None
    if "sl" in file:
        file_name = '_'.join(file.split("_")[0:-2]) + "_sl"
    elif "fl" in file:
        file_name = '_'.join(file.split("_")[0:-2]) + "_fl"
    else:
        print('Check LUXE geometry file!')
        print('Exiting...')
        exit()

    # single: one particle triggers exactly one pixel
    # multi: one particle may trigger more than one pixel

    # num hits per detector
    single_pixel_hits_occupancy = {}
    # num triggered pixels per detector
    net_single_pixel_hits_occupancy = {}

    multi_pixel_hits_occupancy = {}
    net_multi_pixel_hits_occupancy = {}

    # mid of bin as information for tracking
    single_hits_pixel_for_tracking = {}
    multi_hits_pixel_for_tracking = {}

    # tracking number of hits as histogram data for each pixel
    hist_data_single_hits = {}
    hist_data_multi_hits = {}

    data = np.load(f"{path_to_file}.npy", allow_pickle=True)[()]
    for plane_name, plane in zip(data['Plane list'].keys(), data['Plane list'].values()):

        # mid to edge of bin
        d_mid_edge_x = 0.5 * plane.bin_width_in_x()
        d_mid_edge_y = 0.5 * plane.bin_width_in_y()

        multi_pixel_hits = {}

        # every particle on each plan
        for particle, hit in zip(plane.true_hits_dictionary.keys(), plane.true_hits_dictionary.values()):
            multi_pixel_hits.update({particle: hit})

            momentum = data['Particle momentum history log'][particle][plane_name]

            angle_xz = np.arctan2(momentum[0], momentum[2])
            angle_yz = np.arctan2(momentum[1], momentum[2])

            momentum_factor_x = None
            if angle_xz != np.pi:
                momentum_factor_x = 1 / abs(1 - angle_xz / np.pi)

            momentum_factor_y = None
            if angle_yz != np.pi:
                momentum_factor_y = 1 / abs(1 - angle_yz / np.pi)

            # bin_edges = [(x_min, x_max),
            #              (y_min, y_max)]

            bin_edges = plane.bin_edges_for_hit(hit)
            single_hits_pixel_for_tracking.update({plane_name + "_" + particle: [bin_edges[0][0] + d_mid_edge_x,
                                                                                 bin_edges[0][0] + d_mid_edge_y,
                                                                                 hit[2]]})
            multi_hits_pixel_for_tracking.update({plane_name + "_" + particle: [bin_edges[0][0] + d_mid_edge_x,
                                                                                bin_edges[0][0] + d_mid_edge_y,
                                                                                hit[2]]})

            dx_to_lower_bin_edge = hit[0] - bin_edges[0][0]
            dx_to_upper_bin_edge = bin_edges[0][1] - hit[0]

            dy_to_lower_bin_edge = hit[1] - bin_edges[1][0]
            dy_to_upper_bin_edge = bin_edges[1][1] - hit[1]

            # upper---------lower
            # |-----|-----|-----|upper
            # |  1  |  2  |  3  |  |
            # |-----|-----|-----|  |
            # |  4  | hit |  5  |  |
            # |-----|-----|-----|  |
            # |  6  |  7  |  8  |  |
            # |-----|-----|-----|lower

            check_2 = False
            check_4 = False
            check_5 = False
            check_7 = False

            # check 5
            if angle_xz > 0:
                momentum_factor_x = 1 / momentum_factor_x
            if abs(d_mid_edge_x - dx_to_lower_bin_edge) / d_mid_edge_x * momentum_factor_x > 0.95:
                if bin_edges[0][0] > plane.limits_x[0]:
                    additional_pixel = [bin_edges[0][0] - d_mid_edge_x,
                                        bin_edges[0][0] + d_mid_edge_y,
                                        hit[2]]
                    multi_pixel_hits.update({particle + '_5': additional_pixel})
                    multi_hits_pixel_for_tracking.update({plane_name + "_" + particle + '_5': additional_pixel})
                    check_5 = True

            elif angle_xz < 0:
                momentum_factor_x = 1 / momentum_factor_x
                # check 4
            if abs((d_mid_edge_x - dx_to_upper_bin_edge) / d_mid_edge_x) * momentum_factor_x > 0.95:
                if bin_edges[0][1] < plane.limits_x[1]:
                    additional_pixel = [bin_edges[0][1] + d_mid_edge_x,
                                        bin_edges[0][0] + d_mid_edge_y,
                                        hit[2]]
                    multi_pixel_hits.update({particle + '_4': additional_pixel})
                    multi_hits_pixel_for_tracking.update({plane_name + "_" + particle + '_4': additional_pixel})
                    check_4 = True

            if angle_yz > 0:
                momentum_factor_y = 1 / momentum_factor_y
                # check 7
            if abs(d_mid_edge_y - dy_to_lower_bin_edge) / d_mid_edge_y * momentum_factor_y > 0.95:
                if bin_edges[1][0] > plane.limits_y[0]:
                    additional_pixel = [bin_edges[0][0] + d_mid_edge_x,
                                        bin_edges[1][0] - d_mid_edge_y,
                                        hit[2]]
                    multi_pixel_hits.update({particle + '_7': additional_pixel})
                    multi_hits_pixel_for_tracking.update({plane_name + "_" + particle + '_7': additional_pixel})
                    check_7 = True

            if angle_yz < 0:
                momentum_factor_y = 1 / momentum_factor_y
                # check 2
            if abs((d_mid_edge_y - dy_to_upper_bin_edge) / d_mid_edge_y) * momentum_factor_y > 0.95:
                if bin_edges[1][1] < plane.limits_y[1]:
                    additional_pixel = [bin_edges[0][0] + d_mid_edge_x,
                                        bin_edges[1][1] + d_mid_edge_y,
                                        hit[2]]
                    multi_pixel_hits.update({particle + '_2': additional_pixel})
                    multi_hits_pixel_for_tracking.update({plane_name + "_" + particle + '_2': additional_pixel})
                    check_2 = True

            # check 1
            if check_2 and check_4:
                additional_pixel = [bin_edges[0][1] + d_mid_edge_x,
                                    bin_edges[1][1] + d_mid_edge_y,
                                    hit[2]]
                multi_pixel_hits.update({particle + '_1': additional_pixel})
                multi_hits_pixel_for_tracking.update({plane_name + "_" + particle + '_1': additional_pixel})

            # check 6
            if check_4 and check_7:
                additional_pixel = [bin_edges[0][1] + d_mid_edge_x,
                                    bin_edges[1][0] - d_mid_edge_y,
                                    hit[2]]
                multi_pixel_hits.update({particle + '_6': additional_pixel})
                multi_hits_pixel_for_tracking.update({plane_name + "_" + particle + '_6': additional_pixel})

            # check 3
            if check_2 and check_5:
                additional_pixel = [bin_edges[0][0] - d_mid_edge_x,
                                    bin_edges[1][1] + d_mid_edge_y,
                                    hit[2]]
                multi_pixel_hits.update({particle + '_3': additional_pixel})
                multi_hits_pixel_for_tracking.update({plane_name + "_" + particle + '_3': additional_pixel})

            # check 8
            if check_5 and check_7:
                additional_pixel = [bin_edges[0][0] - d_mid_edge_x,
                                    bin_edges[1][0] - d_mid_edge_y,
                                    hit[2]]
                multi_pixel_hits.update({particle + '_8': additional_pixel})
                multi_hits_pixel_for_tracking.update({plane_name + "_" + particle + '_8': additional_pixel})

        # histogram part
        correction_factor = None
        if "sl" in file:
            correction_factor = 2 ** 3
            #     label = r"occupancy / $2^6$ pixel"
        elif "fl" in file:
            correction_factor = 1
            #     label = r"occupancy / pixel"

        # single pixel hit case
        plt.figure(figsize=(16, 9), dpi=200)
        h_single, x_single, y_single, im_single = plt.hist2d(
            [value[0] for value in plane.true_hits_dictionary.values()],
            [value[1] for value in plane.true_hits_dictionary.values()],
            bins=[int(1 / correction_factor * plane.num_bins[0]),
                  int(1 / correction_factor * plane.num_bins[1])],
            range=[[plane.limits_x[0], plane.limits_x[1]],
                   [plane.limits_y[0], plane.limits_y[1]]])
        plt.close()
        # plt.plot([plane.limits_x[0], plane.limits_x[1]], [plane.limits_y[0], plane.limits_y[0]], "k", linewidth=15)
        # plt.vlines(plane.limits_x[0], plane.limits_y[0], plane.limits_y[1], "k", linewidth=15)
        # plt.plot([plane.limits_x[0], plane.limits_x[1]], [plane.limits_y[1], plane.limits_y[1]], "k", linewidth=15)
        # plt.vlines(plane.limits_x[1], plane.limits_y[0], plane.limits_y[1], "k", linewidth=15,
        #            label="detector edges")
        # plt.xlabel("x [m]", fontsize=18)
        # plt.ylabel("y [m]", fontsize=18)
        # plt.xticks(fontsize=18)
        # plt.yticks(fontsize=18)
        # plt.xlim([plane.limits_x[0], plane.limits_x[1]])
        # plt.ylim([plane.limits_y[0], plane.limits_y[1]])
        # plt.legend(loc="best", fontsize=18)
        # cbar = plt.colorbar()
        # cbar.set_label(label, size=18)
        # cbar.ax.tick_params(labelsize=18)
        # plt.title(r"Pixel occupancy", fontsize=18, loc="left")

        # plt.savefig(f"{folder}/occupancy/" + file_name + "/" + plane_name.replace(" ", "_") +
        #             "_single_pixel_hit.pdf")
        # plt.close()

        # flatten to make a histogram out of the data
        flattened_2d_single = h_single.reshape(1, int(1 / correction_factor * plane.num_bins[0]) *
                                               int(1 / correction_factor * plane.num_bins[1]))[0].tolist()

        # hits on plane
        single_pixel_hits_occupancy.update({f"Single hits on {plane_name}": sum(flattened_2d_single)})

        # histogram for looking at hit density
        n_pixel_single, bins_pixel_single, patches_single = plt.hist(flattened_2d_single,
                                                                     bins=int(max(flattened_2d_single)) + 1,
                                                                     range=(0, max(flattened_2d_single)))

        hist_data_single_hits.update({f"Histogram data of {plane_name}": n_pixel_single})

        # getting rid of double hits
        for j in range(np.shape(h_single)[0]):
            if max(h_single[j]) > 1:
                for k in range(np.shape(h_single)[1]):
                    if h_single[j, k] > 1:
                        h_single[j, k] = 1

        # save as train image for NN
        # np.save(f"{folder}/train_images_NN/{file}_single_pixel_hit_{plane_name.replace(" ", "_")}", h_single)

        # net triggered pixels
        net_single_pixel_hits_occupancy.update({f"Net single hits on {plane_name}":
                                                sum(h_single.reshape(1, int(1 / correction_factor * plane.num_bins[0]) *
                                                                     int(1 / correction_factor * plane.num_bins[1]))
                                                    [0].tolist())})

        # multi pixel hit case
        plt.figure(figsize=(16, 9), dpi=200)
        h_multi, x_multi, y_multi, im_multi = plt.hist2d([value[0] for value in multi_pixel_hits.values()],
                                                         [value[1] for value in multi_pixel_hits.values()],
                                                         bins=[int(1 / correction_factor * plane.num_bins[0]),
                                                               int(1 / correction_factor * plane.num_bins[1])],
                                                         range=[[plane.limits_x[0], plane.limits_x[1]],
                                                                [plane.limits_y[0], plane.limits_y[1]]])

        # plt.plot([plane.limits_x[0], plane.limits_x[1]], [plane.limits_y[0], plane.limits_y[0]], "k", linewidth=15)
        # plt.vlines(plane.limits_x[0], plane.limits_y[0], plane.limits_y[1], "k", linewidth=15)
        # plt.plot([plane.limits_x[0], plane.limits_x[1]], [plane.limits_y[1], plane.limits_y[1]], "k", linewidth=15)
        # plt.vlines(plane.limits_x[1], plane.limits_y[0], plane.limits_y[1], "k", linewidth=15,
        #            label="detector edges")
        # plt.xlabel("x [m]", fontsize=18)
        # plt.ylabel("y [m]", fontsize=18)
        # plt.xticks(fontsize=18)
        # plt.yticks(fontsize=18)
        # plt.xlim([plane.limits_x[0], plane.limits_x[1]])
        # plt.ylim([plane.limits_y[0], plane.limits_y[1]])
        # plt.legend(loc="best", fontsize=18)
        # cbar = plt.colorbar()
        # cbar.set_label(label, size=18)
        # cbar.ax.tick_params(labelsize=18)
        # plt.title(r"Pixel occupancy", fontsize=18, loc="left")
        # plt.savefig(f"{folder}/occupancy/" + file_name + "/" + plane_name.replace(" ", "_") +
        #             "_multi_pixel_hit.pdf")
        plt.close()

        # flatten to make a histogram out of the data
        flattened_2d_multi = h_multi.reshape(1, int(1 / correction_factor * plane.num_bins[0]) *
                                             int(1 / correction_factor * plane.num_bins[1]))[0].tolist()
        # hits on plane
        multi_pixel_hits_occupancy.update({f"Multi hits on {plane_name}": sum(flattened_2d_multi)})

        # histogram for looking at hit density
        n_pixel_multi, bins_pixel_multi, patches_multi = plt.hist(flattened_2d_multi,
                                                                  range=(0, max(flattened_2d_multi)),
                                                                  bins=int(max(flattened_2d_multi)) + 1)

        hist_data_multi_hits.update({f"Histogram data of {plane_name}": n_pixel_multi})

        # getting rid of double hits
        for j in range(np.shape(h_multi)[0]):
            if max(h_multi[j]) > 1:
                for k in range(np.shape(h_multi)[1]):
                    if h_multi[j, k] > 1:
                        h_multi[j, k] = 1

        # save as train image for NN
        # np.save(f"{folder}/train_images_NN/" + file_name + "/multi_pixel_hit_" + plane_name.replace(" ", "_"),
        #         h_multi)

        # net triggered pixels
        net_multi_pixel_hits_occupancy.update({f"Net multi hits on {plane_name}":
                                               sum(h_multi.reshape(1, int(1 / correction_factor * plane.num_bins[0]) *
                                                                   int(1 / correction_factor * plane.num_bins[
                                                                               1]))[0].tolist())})

        # writing csv files to make the data set compatible with the ML track challenge
    hit_id = 0
    with open(f'{folder}/pixel/single_pixel_hits_{file_name}.csv', 'w', newline='') as csv_file:
        fieldnames = ['hit_ID', 'x', 'y', 'z', 'layer_ID', 'module_ID', 'is_signal', 'particle_ID', 'particle_energy']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for key, hit in zip(single_hits_pixel_for_tracking.keys(), single_hits_pixel_for_tracking.values()):
            layer, module = get_layer_and_module(hit[0], hit[1], hit[2])
            writer.writerow({'hit_ID': hit_id,
                             'x': hit[0],
                             'y': hit[1],
                             'z': hit[2],
                             'layer_ID': layer,
                             'module_ID': module,
                             'is_signal': True,
                             'particle_ID': key.split("_")[1],
                             'particle_energy': norm(data["Particle momentum history log"][key.split("_")[1]]
                                                     [key.split("_")[0]])})
            hit_id += 1

    hit_id = 0
    with open(f'{folder}/pixel/multi_pixel_hits_{file_name}.csv', 'w', newline='') as csv_file:
        fieldnames = ['hit_ID', 'x', 'y', 'z', 'layer_ID', 'module_ID', 'is_signal', 'particle_ID', 'particle_energy']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for key, hit in zip(multi_hits_pixel_for_tracking.keys(), multi_hits_pixel_for_tracking.values()):
            layer, module = get_layer_and_module(hit[0], hit[1], hit[2])
            writer.writerow({'hit_ID': hit_id,
                             'x': hit[0],
                             'y': hit[1],
                             'z': hit[2],
                             'layer_ID': layer,
                             'module_ID': module,
                             'is_signal': True,
                             'particle_ID': key.split("_")[1],
                             'particle_energy': norm(data["Particle momentum history log"][key.split("_")[1]]
                                                     [key.split("_")[0]])})
            hit_id += 1

    np.save(f"{folder}/occupancy/occupancy_info_" + file_name,
            np.array([{"single hit": single_pixel_hits_occupancy,
                       "single hit net": net_single_pixel_hits_occupancy,
                       "multi hit": multi_pixel_hits_occupancy,
                       "multi hit net": net_multi_pixel_hits_occupancy,
                       "histogram single": hist_data_single_hits,
                       "histogram multi": hist_data_single_hits}]))


def get_layer_and_module(x, y, z):
    """Returns module and layer for given x,y,z value of hits.
    :param x: x value of hit
    :param y: y value of hit
    :param z: z value of hit

    :return
        module und layer of specified hit x,y,z value
    """
    if z == 3.962:
        layer = 0
    elif z == 3.950:
        layer = 1
    elif z == 4.062:
        layer = 2
    elif z == 4.050:
        layer = 3
    elif z == 4.162:
        layer = 4
    elif z == 4.150:
        layer = 5
    elif z == 4.262:
        layer = 6
    elif z == 4.250:
        layer = 7
    else:
        return -1, -1  # fake value if data is corrupted

    module = 0
    if layer % 2 == 0:
        for key, value in get_close_to_beamline_stave().items():
            if value[0] <= x <= value[1] and value[2] <= y <= value[3]:
                module = int(key)
                break
    else:
        for key, value in get_not_close_to_beamline_stave().items():
            if value[0] <= x <= value[1] and value[2] <= y <= value[3]:
                module = int(key)
                break

    return layer, module


def get_close_to_beamline_stave():
    """Returns sensors on the stave closest to beamline.
    :return
        dictionary of detector module dimensions
    """
    return {"0": [0.05275912, 0.08270088, -0.00688128, 0.00688128],
            "1": [0.08285912, 0.11280088, -0.00688128, 0.00688128],
            "2": [0.11295912, 0.14290088, -0.00688128, 0.00688128],
            "3": [0.14305912, 0.17300088, -0.00688128, 0.00688128],
            "4": [0.17315912, 0.20310088, -0.00688128, 0.00688128],
            "5": [0.20325912, 0.23320088, -0.00688128, 0.00688128],
            "6": [0.23335912, 0.26330088, -0.00688128, 0.00688128],
            "7": [0.26345912, 0.29340088, -0.00688128, 0.00688128],
            "8": [0.29355912, 0.32350088, -0.00688128, 0.00688128]}


def get_not_close_to_beamline_stave():
    """Returns sensors on the stave which is not closest to beamline.
    :return
        dictionary of detector module dimensions
    """
    return {"0": [0.28355912, 0.31350088, -0.00688128, 0.00688128],
            "1": [0.31365912, 0.34360088, -0.00688128, 0.00688128],
            "2": [0.34375912, 0.37370088, -0.00688128, 0.00688128],
            "3": [0.37385912, 0.40380088, -0.00688128, 0.00688128],
            "4": [0.40395912, 0.43390088, -0.00688128, 0.00688128],
            "5": [0.43405912, 0.46400088, -0.00688128, 0.00688128],
            "6": [0.46415912, 0.49410088, -0.00688128, 0.00688128],
            "7": [0.49425912, 0.52420088, -0.00688128, 0.00688128],
            "8": [0.52435912, 0.55430088, -0.00688128, 0.00688128]}
