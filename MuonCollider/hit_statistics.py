import numpy as np
import sys
import csv
import matplotlib.pyplot as plt

sys.path.append("../src")
t_list = np.load("triplet_list.npy", allow_pickle=True)

data = []

tracking_data_file = "/nfs/dust/luxe/user/spatarod/Muon-Collider/DPG2023/Muon1GeV/Output_REC.000_VXD.csv"


with open(tracking_data_file, 'r') as file:
    csv_reader = csv.reader(file)
    csv_header = next(csv_reader)  # access header, csv files should consist of one line of header
    x_index = csv_header.index("x")
    y_index = csv_header.index("y")
    z_index = csv_header.index("z")
    hit_id_index = csv_header.index("hit_ID")
    mc_particle_id_index = csv_header.index("MC_particle_ID")
    pdg_index = csv_header.index("PDG")
    time_index = csv_header.index("time")
    for row in csv_reader:
        data.append([float(r) if row.index(r) not in [time_index, hit_id_index] else r for r in row])


plt.hist([np.sqrt(float(entry[x_index])**2 + float(entry[x_index])**2) for entry in data], bins=100)
plt.savefig("x-y_distribution.pdf")
plt.show()
