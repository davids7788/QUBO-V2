import matplotlib.pyplot as plt
import os
import sys
import numpy as np
sys.path.append("../src")


xi = 5.0


sub_qubo_size_array = [5, 7, 10, 12, 16]


sub_qubo_size_5_impact_list_truth_emergy = []
sub_qubo_size_7_impact_list_truth_emergy = []
sub_qubo_size_10_impact_list_truth_emergy = []
sub_qubo_size_12_impact_list_truth_emergy = []
sub_qubo_size_16_impact_list_truth_emergy = []

sub_qubo_size_5_impact_list_found_emergy = []
sub_qubo_size_7_impact_list_found_emergy = []
sub_qubo_size_10_impact_list_found_emergy = []
sub_qubo_size_12_impact_list_found_emergy = []
sub_qubo_size_16_impact_list_found_emergy = []


sub_qubo_size_5_connection_list_truth_emergy = []
sub_qubo_size_7_connection_list_truth_emergy = []
sub_qubo_size_10_connection_list_truth_emergy = []
sub_qubo_size_12_connection_list_truth_emergy = []
sub_qubo_size_16_connection_list_truth_emergy = []

sub_qubo_size_5_connection_list_found_emergy = []
sub_qubo_size_7_connection_list_found_emergy = []
sub_qubo_size_10_connection_list_found_emergy = []
sub_qubo_size_12_connection_list_found_emergy = []
sub_qubo_size_16_connection_list_found_emergy = []






prefix = f"/nfs/dust/luxe/user/spatarod/ConnectionListStudy/e0gpc/{xi}/smeared"
for sub_qubo_size in sub_qubo_size_array:
    for folder in os.listdir(prefix):
        if "example" in folder:
            for subfolder in os.listdir(f"{prefix}/{folder}"):
                if f"eigensolver_{sub_qubo_size}q_impact_list_reverse" in subfolder:
                    print(f"{prefix}/{folder}/{subfolder}")
                        qubo_log = np.load(f"{prefix}/{folder}/{subfolder}/qubo_log.npy", allow_pickle=True)[()]

                        min_energy = qubo_log["truth minimum energy]
                        found_minimum_energy = qubo_log["computed minimum energy]  

                        if subqubo_size == 5:
                            xi_3_matched_impact += statistics[0]
                            xi_3_fake_impact += statistics[1]
                            xi_3_gen_impact += statistics[2]
                            xi_3_reco_impact += statistics[3]
                        if xi == 4.0:
                            xi_4_matched_impact += statistics[0]
                            xi_4_fake_impact += statistics[1]
                            xi_4_gen_impact += statistics[2]
                            xi_4_reco_impact += statistics[3]                        
                        if xi == 5.0:
                            xi_5_matched_impact += statistics[0]
                            xi_5_fake_impact += statistics[1]
                            xi_5_gen_impact += statistics[2]
                            xi_5_reco_impact += statistics[3]  
                        if xi == 6.0:
                            xi_6_matched_impact += statistics[0]
                            xi_6_fake_impact += statistics[1]
                            xi_6_gen_impact += statistics[2]
                            xi_6_reco_impact += statistics[3]  
                        if xi == 7.0:
                            xi_7_matched_impact += statistics[0]
                            xi_7_fake_impact += statistics[1]
                            xi_7_gen_impact += statistics[2]
                            xi_7_reco_impact += statistics[3]                          
                        
                        

                if f"eigensolver_{sub_qubo_size}q_connection_list" in subfolder:
                    print(f"{prefix}/{folder}/{subfolder}")
                    reco_xplets = np.load(f"{prefix}/{folder}/{subfolder}/reco_xplet_list_ambiguity_solved.npy",
                                        allow_pickle=True)
                    gen_xplets = np.load(f"{prefix}/" + "_".join(["_".join(folder.split("_")[0:3]),"sl", "gen_xplet_list.npy"]),
                                         allow_pickle=True)
                    statistics = get_eff_and_frate(reco_xplets, gen_xplets)
                    if xi == 3.0:
                        xi_3_matched_connection += statistics[0]
                        xi_3_fake_connection += statistics[1]
                        xi_3_gen_connection += statistics[2]
                        xi_3_reco_connection += statistics[3]
                    if xi == 4.0:
                        xi_4_matched_connection += statistics[0]
                        xi_4_fake_connection += statistics[1]
                        xi_4_gen_connection += statistics[2]
                        xi_4_reco_connection += statistics[3]                        
                    if xi == 5.0:
                        xi_5_matched_connection += statistics[0]
                        xi_5_fake_connection += statistics[1]
                        xi_5_gen_connection += statistics[2]
                        xi_5_reco_connection += statistics[3]  
                    if xi == 6.0:
                        xi_6_matched_connection += statistics[0]
                        xi_6_fake_connection += statistics[1]
                        xi_6_gen_connection += statistics[2]
                        xi_6_reco_connection += statistics[3]  
                    if xi == 7.0:
                        xi_7_matched_connection += statistics[0]
                        xi_7_fake_connection += statistics[1]
                        xi_7_gen_connection += statistics[2]
                        xi_7_reco_connection += statistics[3] 

                if f"VQE_IdealQasmSim_{sub_qubo_size}q_TwoLocal_linear_NFT_impact_list_reverse" in subfolder and xi < 6:
                    print(f"{prefix}/{folder}/{subfolder}")
                    reco_xplets = np.load(f"{prefix}/{folder}/{subfolder}/reco_xplet_list_ambiguity_solved.npy", allow_pickle=True)
                    gen_xplets = np.load(f"{prefix}/" + "_".join(["_".join(folder.split("_")[0:3]),"sl", "gen_xplet_list.npy"]),
                                         allow_pickle=True)
                    statistics = get_eff_and_frate(reco_xplets, gen_xplets)
                    if xi == 3.0:
                        xi_3_matched_impact_VQE += statistics[0]
                        xi_3_fake_impact_VQE += statistics[1]
                        xi_3_gen_impact_VQE += statistics[2]
                        xi_3_reco_impact_VQE += statistics[3]
                    if xi == 4.0:
                        xi_4_matched_impact_VQE += statistics[0]
                        xi_4_fake_impact_VQE += statistics[1]
                        xi_4_gen_impact_VQE += statistics[2]
                        xi_4_reco_impact_VQE += statistics[3]                        
                    if xi == 5.0:
                        xi_5_matched_impact_VQE += statistics[0]
                        xi_5_fake_impact_VQE += statistics[1]
                        xi_5_gen_impact_VQE += statistics[2]
                        xi_5_reco_impact_VQE += statistics[3]  
                                    

                        
                if f"VQE_IdealQasmSim_{sub_qubo_size}q_TwoLocal_linear_NFT_connection_list" in subfolder and xi < 6:
                    print(f"{prefix}/{folder}/{subfolder}")
                    reco_xplets = np.load(f"{prefix}/{folder}/{subfolder}/reco_xplet_list_ambiguity_solved.npy",
                                        allow_pickle=True)
                    gen_xplets = np.load(f"{prefix}/" + "_".join(["_".join(folder.split("_")[0:3]),"sl", "gen_xplet_list.npy"]),
                                         allow_pickle=True)
                    statistics = get_eff_and_frate(reco_xplets, gen_xplets)
                    if xi == 3.0:
                        xi_3_matched_connection_VQE += statistics[0]
                        xi_3_fake_connection_VQE += statistics[1]
                        xi_3_gen_connection_VQE += statistics[2]
                        xi_3_reco_connection_VQE += statistics[3]
                    if xi == 4.0:
                        xi_4_matched_connection_VQE += statistics[0]
                        xi_4_fake_connection_VQE += statistics[1]
                        xi_4_gen_connection_VQE += statistics[2]
                        xi_4_reco_connection_VQE += statistics[3]                        
                    if xi == 5.0:
                        xi_5_matched_connection_VQE += statistics[0]
                        xi_5_fake_connection_VQE += statistics[1]
                        xi_5_gen_connection_VQE += statistics[2]
                        xi_5_reco_connection_VQE += statistics[3]  
                      
                        
                        
def get_standard_deviation(k, n):
    if k == 0:
        return 0
    if k > n:
        k = n
    return np.sqrt(((k + 1) * (k + 2)) / ((n + 2) * (n + 3)) - ((k + 1)**2 / (n + 2)**2))                        

                        
plt.figure(figsize=(8,6), dpi=400)
plt.errorbar(xi_values, [xi_3_matched_impact / xi_3_gen_impact, 
                         xi_4_matched_impact / xi_4_gen_impact,
                         xi_5_matched_impact / xi_5_gen_impact,
                         xi_6_matched_impact / xi_6_gen_impact,
                         xi_7_matched_impact / xi_7_gen_impact],
             yerr=[get_standard_deviation(xi_3_matched_impact, xi_3_gen_impact),
                   get_standard_deviation(xi_4_matched_impact, xi_4_gen_impact),
                   get_standard_deviation(xi_5_matched_impact, xi_5_gen_impact),
                   get_standard_deviation(xi_6_matched_impact, xi_6_gen_impact),
                   get_standard_deviation(xi_7_matched_impact, xi_7_gen_impact)],
             color="royalblue",
             alpha=1,
             linewidth=2,
             marker="o",
             markersize=6,
             label="impact list eigensolver")

plt.errorbar(xi_values, [xi_3_matched_connection / xi_3_gen_connection, 
                         xi_4_matched_connection / xi_4_gen_connection,
                         xi_5_matched_connection / xi_5_gen_connection,
                         xi_6_matched_connection / xi_6_gen_connection,
                         xi_7_matched_connection / xi_7_gen_connection],
             yerr=[get_standard_deviation(xi_3_matched_connection, xi_3_gen_connection),
                   get_standard_deviation(xi_5_matched_connection, xi_4_gen_connection),
                   get_standard_deviation(xi_5_matched_connection, xi_5_gen_connection),
                   get_standard_deviation(xi_6_matched_connection, xi_6_gen_connection),
                   get_standard_deviation(xi_7_matched_connection, xi_7_gen_connection)],
             color="black",
             alpha=1,
             linewidth=2,
             marker="v",
             markersize=6,
             label="connection list eigensolver")

plt.errorbar(xi_values[0:3], [xi_3_matched_impact_VQE / xi_3_gen_impact_VQE, 
                         xi_4_matched_impact_VQE / xi_4_gen_impact_VQE,
                         xi_5_matched_impact_VQE / xi_5_gen_impact_VQE],
             yerr=[get_standard_deviation(xi_3_matched_impact_VQE, xi_3_gen_impact_VQE),
                   get_standard_deviation(xi_4_matched_impact_VQE, xi_4_gen_impact_VQE),
                   get_standard_deviation(xi_5_matched_impact_VQE, xi_5_gen_impact_VQE)],
             color="royalblue",
             alpha=0.8,
             linewidth=2,
             marker="d",
             markersize=6,
             label="impact list VQE")


plt.errorbar(xi_values[0:3], [xi_3_matched_connection_VQE / xi_3_gen_connection_VQE, 
                         xi_4_matched_connection_VQE / xi_4_gen_connection_VQE,
                         xi_5_matched_connection_VQE / xi_5_gen_connection_VQE],
             yerr=[get_standard_deviation(xi_3_matched_connection_VQE, xi_3_gen_connection_VQE),
                   get_standard_deviation(xi_5_matched_connection_VQE, xi_4_gen_connection_VQE),
                   get_standard_deviation(xi_5_matched_connection_VQE, xi_5_gen_connection_VQE)],
             color="black",
             alpha=0.8,
             linewidth=2,
             marker="s",
             markersize=6,
             label="connection list VQE")

plt.grid()
plt.xticks([3, 4, 5, 6, 7], fontsize=16)
plt.yticks(fontsize=16)
plt.legend(loc="best", fontsize=16)
plt.xlabel(r"$\xi$", loc="right", fontsize=16)
plt.ylabel("Efficiency", loc="top", fontsize=16)
plt.savefig("Efficiency_vs_xi.pdf", bbox_inches="tight")
plt.savefig("Efficiency_vs_xi.jpg", bbox_inches="tight")


plt.figure(figsize=(8,6), dpi=400)
plt.errorbar(xi_values, [xi_3_fake_impact / xi_3_reco_impact, 
                         xi_4_fake_impact / xi_4_reco_impact,
                         xi_5_fake_impact / xi_5_reco_impact,
                         xi_6_fake_impact / xi_6_reco_impact,
                         xi_7_fake_impact / xi_7_reco_impact],
             yerr=[get_standard_deviation(xi_3_fake_impact, xi_3_reco_impact),
                   get_standard_deviation(xi_4_fake_impact, xi_4_reco_impact),
                   get_standard_deviation(xi_5_fake_impact, xi_5_reco_impact),
                   get_standard_deviation(xi_6_fake_impact, xi_6_reco_impact),
                   get_standard_deviation(xi_7_fake_impact, xi_7_reco_impact)],
             color="royalblue",
             alpha=1,
             linewidth=2,
             marker="o",
             markersize=6,
             label="impact list eigensolver")


plt.errorbar(xi_values[0:3], [xi_3_fake_impact_VQE / xi_3_reco_impact_VQE, 
                         xi_4_fake_impact_VQE / xi_4_reco_impact_VQE,
                         xi_5_fake_impact_VQE / xi_5_reco_impact_VQE],
             yerr=[get_standard_deviation(xi_3_fake_impact_VQE, xi_3_reco_impact_VQE),
                   get_standard_deviation(xi_4_fake_impact_VQE, xi_4_reco_impact_VQE),
                   get_standard_deviation(xi_5_fake_impact_VQE, xi_5_reco_impact_VQE)],
             color="royalblue",
             alpha=0.8,
             linewidth=2,
             marker="d",
             markersize=6,
             label="impact list VQE")

plt.errorbar(xi_values, [xi_3_fake_connection / xi_3_reco_connection, 
                         xi_4_fake_connection / xi_4_reco_connection,
                         xi_5_fake_connection / xi_5_reco_connection,
                         xi_6_fake_connection / xi_6_reco_connection,
                         xi_7_fake_connection / xi_7_reco_connection],
             yerr=[get_standard_deviation(xi_3_fake_connection, xi_3_reco_connection),
                   get_standard_deviation(xi_5_fake_connection, xi_4_reco_connection),
                   get_standard_deviation(xi_5_fake_connection, xi_5_reco_connection),
                   get_standard_deviation(xi_6_fake_connection, xi_6_reco_connection),
                   get_standard_deviation(xi_7_fake_connection, xi_7_reco_connection)],
             color="black",
             alpha=1,
             linewidth=2,
             marker="v",
             markersize=6,
             label="connection list eigensolver")

plt.errorbar(xi_values[0:3], [xi_3_fake_connection_VQE / xi_3_reco_connection_VQE, 
                         xi_4_fake_connection_VQE / xi_4_reco_connection_VQE,
                         xi_5_fake_connection_VQE / xi_5_reco_connection_VQE],
             yerr=[get_standard_deviation(xi_3_fake_connection_VQE, xi_3_reco_connection_VQE),
                   get_standard_deviation(xi_5_fake_connection_VQE, xi_4_reco_connection_VQE),
                   get_standard_deviation(xi_5_fake_connection_VQE, xi_5_reco_connection_VQE)],
             color="black",
             alpha=0.8,
             linewidth=2,
             marker="s",
             markersize=6,
             label="connection list VQE")

plt.grid()
plt.xticks([3, 4, 5, 6, 7], fontsize=16)
plt.yticks(fontsize=16)
plt.legend(loc="best", fontsize=16)
plt.xlabel(r"$\xi$", loc="right", fontsize=16)
plt.ylabel("Fake rate", loc="top", fontsize=16)
plt.savefig("Fake_rate_vs_xi.pdf", bbox_inches="tight")
plt.savefig("Fake_rate_vs_xi.jpg", bbox_inches="tight")