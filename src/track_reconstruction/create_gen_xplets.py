import csv
import numpy as np
import os

from pattern.detector_hit import DetectorHit
from pattern.multiplet import Multiplet


def gen_xplets_simplified_LUXE(tracking_data_file: str,
                               save_to_folder: str) -> None:
    """Creates truth generated xplet list.
        :param tracking_data_file: .csv tracking data file
        :param save_to_folder: location to save the gen_xplet list
    """

    print("Create X-plets on generator level with truth information...")
    output_name = ".".join(tracking_data_file.split("/")[-1].split(".")[0:-1])
    if os.path.isfile(output_name):
        return

    multiplets_dict = {}
    gen_multiplets = []

    with open(tracking_data_file, 'r') as file:
        csv_reader_tracking_file = csv.reader(file)
        _ = next(csv_reader_tracking_file)

        for row in csv_reader_tracking_file:
            # Create new DetectorHit object
            hit = DetectorHit(row)

            # check if particle ID was already seen
            if hit.particle_id in multiplets_dict.keys():
                multiplets_dict[hit.particle_id].append(hit)
            else:
                multiplets_dict.update({hit.particle_id: [hit]})

        # Loop over dictionary entries, which are sorted as particle_number: [detector_hit_0, detector_hit_1, ...]
        for particle_entry in multiplets_dict.values():
            particle_entry.sort(key=lambda entry: entry.z)

            truth_multiplet = Multiplet()
            for value in particle_entry:
                truth_multiplet.add_hit(value)

            gen_multiplets.append(truth_multiplet)

        np.save(f"{save_to_folder}/{output_name}_gen_xplet_list", gen_multiplets)
        print(f"Number of generated Xplets: {len(gen_multiplets)}\n")
        print("-----------------------------------\n")
