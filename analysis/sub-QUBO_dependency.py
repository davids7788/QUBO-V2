import matplotlib.pyplot as plt
import os
import sys
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chisquare, chi2
sys.path.append("../src")
from pattern.doublet import Doublet

xi = 6.0
prefix = f"/nfs/dust/luxe/user/spatarod/ConnectionListStudy/e0gpc/{xi}/smeared"

BX = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

sub_qubo_size = 7

num_BX_impact_list_completed = 0
num_BX_connection_list_completed = 0

for folder in os.listdir(prefix):
    if "example" in folder:
        for subfolder in os.listdir(f"{prefix}/{folder}"):
            if f"eigensolver_{sub_qubo_size}q_impact_list_reverse" in subfolder:
                is_complete = False
                for file in os.listdir(f"{prefix}/{folder}/{subfolder}"):
                    if "reco_xplet_list_ambiguity_solved" in file:
                        is_complete = True
                        num_BX_impact_list_completed +=1
                        
            if f"eigensolver_{sub_qubo_size}q_connection_list" in subfolder:
                is_complete = False
                for file in os.listdir(f"{prefix}/{folder}/{subfolder}"):
                    if "reco_xplet_list_ambiguity_solved" in file:
                        is_complete = True
                        num_BX_connection_list_completed +=1                        



def get_eff_and_frate(reco_xplets, gen_xplets):

    matched_tracks = 0
    fake_tracks = 0

    matched_tracks_bookkeeping = set()
    
    fake_chi_squared = []
    matched_chi_squared = []
    
    for track in reco_xplets:
        matched = False
        p_id = None
        ids = set(track.particle_ids.values())
        for test_id in ids:
            count = 0
            for particle_id in track.particle_ids.values():
                if test_id == particle_id:
                    count += 1
            if count >= 3:
                p_id = test_id
                matched = True
        if matched:
            if p_id not in matched_tracks_bookkeeping:
                matched_tracks += 1
                matched_tracks_bookkeeping.add(p_id)
                matched_chi_squared.append(fit_lin_track(track))
            else:
                pass
        else:
            fake_tracks += 1
            fake_chi_squared.append(fit_lin_track(track))
   
    return matched_tracks / len(gen_xplets), fake_tracks / len(reco_xplets), matched_chi_squared, fake_chi_squared


def lin_func(x, a, b):
    """Linear function with slope a and bias b
    :param x: data points
    :param a: slope
    :param b: bias
    """
    return a * x + b




def fit_lin_track(track):
    """Linear fit of the track in xz and yz direction. Chi squared values are averaged. Chi squared and p-value
    attributes are set for the xplet.
    """
    x = [value[0] for value in track.coordinates.values()]
    y = [value[1] for value in track.coordinates.values()]
    z = [value[2] for value in track.coordinates.values()]

    popt_xz, _ = curve_fit(lin_func, z, x)
    popt_yz, _ = curve_fit(lin_func, z, y)

    f_exp=[popt_xz[0] * z_i + popt_xz[1] for z_i in z]
    chi_xz = sum([((x_i - e) / 5e-6)**2  for x_i, e in zip(x, f_exp)]) / (len(x) - 2)

    f_exp=[popt_yz[0] * z_i + popt_yz[1] for z_i in z]
    chi_yz = sum([((y_i - e) / 5e-6)**2  for y_i, e in zip(y, f_exp)]) / (len(y) - 2)
 
    return 0.5 * (chi_xz + chi_yz), chi2.sf(0.5 * (chi_xz + chi_yz), df=len(x) - 1 - (len(x) - 2))
       
xi_7_eff_connection_list = []
xi_7_connection_list_chi_squared_matched = []
xi_7_connection_list_chi_squared_fake = []

xi_7_eff_impact_list = []
xi_7_impact_list_chi_squared_matched = []
xi_7_impact_list_chi_squared_fake = []


xi_7_frate_connection_list = []
xi_7_frate_impact_list = []

for folder in os.listdir(prefix):
    if "example" in folder:
        for subfolder in os.listdir(f"{prefix}/{folder}"):
            if f"eigensolver_{sub_qubo_size}q_impact_list_reverse" in subfolder:
                print(f"{prefix}/{folder}/{subfolder}")
                reco_xplets = np.load(f"{prefix}/{folder}/{subfolder}/reco_xplet_list_ambiguity_solved.npy", allow_pickle=True)
                gen_xplets = np.load(f"{prefix}/" + "_".join(["_".join(folder.split("_")[0:3]),"sl", "gen_xplet_list.npy"]),
                                     allow_pickle=True)
                statistics = get_eff_and_frate(reco_xplets, gen_xplets)
                xi_7_eff_impact_list.append(statistics[0])
                xi_7_frate_impact_list.append(statistics[1])
                
                for chi in statistics[2]:
                    xi_7_impact_list_chi_squared_matched.append(chi[0])
                for chi in statistics[3]:
                    xi_7_impact_list_chi_squared_fake.append(chi[0])

           
            if f"eigensolver_{sub_qubo_size}q_connection_list" in subfolder:
                print(f"{prefix}/{folder}/{subfolder}")
                reco_xplets = np.load(f"{prefix}/{folder}/{subfolder}/reco_xplet_list_ambiguity_solved.npy",
                                    allow_pickle=True)
                gen_xplets = np.load(f"{prefix}/" + "_".join(["_".join(folder.split("_")[0:3]),"sl", "gen_xplet_list.npy"]),
                                     allow_pickle=True)
                statistics = get_eff_and_frate(reco_xplets, gen_xplets)
                xi_7_eff_connection_list.append(statistics[0])
                xi_7_frate_connection_list.append(statistics[1])
                
                for chi in statistics[2]:
                    xi_7_connection_list_chi_squared_matched.append(chi[0])
                for chi in statistics[3]:
                    xi_7_connection_list_chi_squared_fake.append(chi[0])

#####################################################
plt.figure(figsize=(8,6), dpi=400)
plt.hist(xi_7_connection_list_chi_squared_matched, 
         bins=50,
         range=(0, 50),
         color="black",
         histtype="step",
         linewidth=2,
         label="matched")
plt.hist(xi_7_connection_list_chi_squared_fake, 
         bins=50, 
         range=(0, 50),
         color="red", 
         histtype="step",
         linewidth=2,
         label="fake")
plt.title(r"$\xi$= 6, 10BX, connection list", fontsize=18, loc="left")
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlabel("Reduced $\chi^2$", loc="right", fontsize=16)
plt.legend(loc="best", fontsize=16)
plt.ylabel("Number of tracks", loc="top", fontsize=16)
plt.savefig("chi_squared_connection_list_only.pdf", bbox_inches="tight")
plt.savefig("chi_squared_connection_list_only.jpg", bbox_inches="tight")
           
#####################################################
plt.figure(figsize=(8,6), dpi=400)
plt.hist(xi_7_impact_list_chi_squared_matched, 
         bins=50,
         range=(0, 50),
         color="royalblue",
         histtype="step",
         linewidth=2,
         label="matched")
plt.hist(xi_7_impact_list_chi_squared_fake, 
         bins=50, 
         range=(0, 50),
         color="darkgreen", 
         histtype="step",
         linewidth=2,
         label="fake")
plt.title(r"$\xi$= 6, 10BX, impact list algorithm", fontsize=18, loc="left")
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(loc="best", fontsize=16)
plt.xlabel("Reduced $\chi^2$", loc="right", fontsize=16)
plt.ylabel("Number of tracks", loc="top", fontsize=16)
plt.savefig("chi_squared_impact_list_only.pdf", bbox_inches="tight")
plt.savefig("chi_squared_impact_list_only.jpg", bbox_inches="tight")
           
#####################################################
           
           
           

