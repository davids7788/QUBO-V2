import os
import sys
from src.detector_plane import DetectorPlane
from src.experimental_results_MC_toy import ExperimentalResults
import matplotlib.pyplot as plt
import numpy as np
import csv

folder = sys.argv[1] 
os.chdir(folder)

if not os.path.isdir("occupancy"):
    os.mkdir("occupancy")

if not os.path.isdir("train_images_NN"):
    os.mkdir("train_images_NN")
    
if not os.path.isdir("pixel"):
    os.mkdir("pixel")
    
    
for file in os.listdir():
    # make sure only raw SimplifiedSimulation files as npy files are in this directory!
    if ".npy" not in file:
        continue

    # every file gets a folder for train images data -> up to 72 per simulation
    if not os.path.isdir("train_images_NN/" + ".".join(file.split(".")[0:-1])):
        os.mkdir("train_images_NN/" + ".".join(file.split(".")[0:-1]))

    if not os.path.isdir("occupancy/" + ".".join(file.split(".")[0:-1])):
        os.mkdir("occupancy/" + ".".join(file.split(".")[0:-1]))

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

    data = np.load(file, allow_pickle=True)[()]
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

            if angle_xz != np.pi:
                momentum_factor_x = 1 / abs(1 - angle_xz / np.pi)

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

            if angle_xz < 0:
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
        if "sl" in file:
            label = r"occupancy / $2^6$ pixel"
            correction_factor = 2**3
        elif "fl" in file:
            label = r"occupancy / pixel"
            correction_factor = 1

        # single pixel hit case
        plt.figure(figsize=(16, 9), dpi=200)
        h_single, x_single, y_single, im_single = plt.hist2d([value[0] for value in
                                                              plane.true_hits_dictionary.values()],
                                                             [value[1] for value in
                                                              plane.true_hits_dictionary.values()],
                                                             bins=[int(1 / correction_factor * plane.num_bins[0]),
                                                                   int(1 / correction_factor * plane.num_bins[1])],
                                                             range=[[plane.limits_x[0], plane.limits_x[1]],
                                                                    [plane.limits_y[0], plane.limits_y[1]]])

        plt.plot([plane.limits_x[0], plane.limits_x[1]], [plane.limits_y[0], plane.limits_y[0]], "k", linewidth=15)
        plt.vlines(plane.limits_x[0], plane.limits_y[0], plane.limits_y[1], "k", linewidth=15)
        plt.plot([plane.limits_x[0], plane.limits_x[1]], [plane.limits_y[1], plane.limits_y[1]], "k", linewidth=15)
        plt.vlines(plane.limits_x[1], plane.limits_y[0], plane.limits_y[1], "k", linewidth=15, label="detector edges")
        plt.xlabel("x [m]", fontsize=18)
        plt.ylabel("y [m]", fontsize=18)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlim([plane.limits_x[0], plane.limits_x[1]])
        plt.ylim([plane.limits_y[0], plane.limits_y[1]])
        plt.legend(loc="best", fontsize=18)
        cbar = plt.colorbar()
        cbar.set_label(label, size=18)
        cbar.ax.tick_params(labelsize=18)
        plt.title(r"Pixel occupancy", fontsize=18, loc="left")
        plt.savefig("occupancy/"
                    + ".".join(file.split(".")[0:-1]) +
                    "/" + plane_name.replace(" ", "_") +
                    "_single_pixel_hit.pdf")
        plt.close()

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
        np.save("train_images_NN/" +
                ".".join(file.split(".")[0:-1]) +
                "/single_pixel_hit_" +
                plane_name.replace(" ", "_"), h_single)

        # net triggered pixels
        net_single_pixel_hits_occupancy.update({f"Net single hits on {plane_name}":
                                                sum(h_single.reshape(1, int(1 / correction_factor *
                                                                            plane.num_bins[0]) *
                                                                     int(1 / correction_factor *
                                                                         plane.num_bins[1]))[0].tolist())})

        # multi pixel hit case
        plt.figure(figsize=(16, 9), dpi=200)
        h_multi, x_multi, y_multi, im_multi = plt.hist2d([value[0] for value in multi_pixel_hits.values()],
                                                         [value[1] for value in multi_pixel_hits.values()],
                                                         bins=[int(1 / correction_factor * plane.num_bins[0]),
                                                               int(1 / correction_factor * plane.num_bins[1])],
                                                         range=[[plane.limits_x[0], plane.limits_x[1]],
                                                                [plane.limits_y[0], plane.limits_y[1]]])

        plt.plot([plane.limits_x[0], plane.limits_x[1]], [plane.limits_y[0], plane.limits_y[0]], "k", linewidth=15)
        plt.vlines(plane.limits_x[0], plane.limits_y[0], plane.limits_y[1], "k", linewidth=15)
        plt.plot([plane.limits_x[0], plane.limits_x[1]], [plane.limits_y[1], plane.limits_y[1]], "k", linewidth=15)
        plt.vlines(plane.limits_x[1], plane.limits_y[0], plane.limits_y[1], "k", linewidth=15, label="detector edges")
        plt.xlabel("x [m]", fontsize=18)
        plt.ylabel("y [m]", fontsize=18)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlim([plane.limits_x[0], plane.limits_x[1]])
        plt.ylim([plane.limits_y[0], plane.limits_y[1]])
        plt.legend(loc="best", fontsize=18)
        cbar = plt.colorbar()
        cbar.set_label(label, size=18)
        cbar.ax.tick_params(labelsize=18)
        plt.title(r"Pixel occupancy", fontsize=18, loc="left")
        plt.savefig("occupancy/" +
                    ".".join(file.split(".")[0:-1]) + "/" +
                    plane_name.replace(" ", "_") +
                    "_multi_pixel_hit.pdf")
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
        np.save("train_images_NN/" +
                ".".join(file.split(".")[0:-1]) +
                "/multi_pixel_hit_" +
                plane_name.replace(" ", "_"), h_multi)

        # net triggered pixels
        net_multi_pixel_hits_occupancy.update({f"Net multi hits on {plane_name}":
                                               sum(h_multi.reshape(1, int(1 / correction_factor *
                                                                          plane.num_bins[0]) *
                                                                   int(1 / correction_factor *
                                                                       plane.num_bins[1]))[0].tolist())})

    # writing csv files to make the data set compatible with the ML track challenge
    hit_id = 0
    with open('pixel/single_pixel_hits_' + '.'.join(file.split('.')[0:-1]) + '.csv', 'w', newline='') as csvfile:
        fieldnames = ['hit_ID', 'x', 'y', 'z', 'layer_ID', 'particle_ID']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for key, hit in zip(single_hits_pixel_for_tracking.keys(), single_hits_pixel_for_tracking.values()):
            writer.writerow({'hit_ID': hit_id,
                             'x': hit[0],
                             'y': hit[1],
                             'z': hit[2],
                             'layer_ID': key.split("_")[0],
                             'particle_ID': key.split("_")[1]})
            hit_id += 1

    hit_id = 0
    with open('pixel/multi_pixel_hits_' + '.'.join(file.split('.')[0:-1]) + '.csv', 'w', newline='') as csvfile:
        fieldnames = ['hit_ID', 'x', 'y', 'z', 'layer_ID', 'particle_ID']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for key, hit in zip(multi_hits_pixel_for_tracking.keys(), multi_hits_pixel_for_tracking.values()):
            writer.writerow({'hit_ID': hit_id,
                             'x': hit[0],
                             'y': hit[1],
                             'z': hit[2],
                             'layer_ID': key.split("_")[0],
                             'particle_ID': key.split("_")[1]})
            hit_id += 1

    np.save("occupancy/occupancy_info_" +
            '.'.join(file.split('.')[0:-1]),
            np.array([{"single hit": single_pixel_hits_occupancy,
                       "single hit net": net_single_pixel_hits_occupancy,
                       "multi hit": multi_pixel_hits_occupancy,
                       "multi hit net": net_multi_pixel_hits_occupancy,
                       "histogram single": hist_data_single_hits,
                       "histogram multi": hist_data_single_hits}]))
