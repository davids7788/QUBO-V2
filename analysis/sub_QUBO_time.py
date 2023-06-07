import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../src")

vqe_5q_0000 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0000_sl-towards_paper/"
                      "407733599_VQE_IdealQasmSim_5q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                      allow_pickle=True)[()]["time tracking subQUBOs"].values()

len_vqe_5q_0000 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0000_sl-towards_paper/"
                          "407733599_VQE_IdealQasmSim_5q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                          allow_pickle=True)[()]["truth solution vector"]

vqe_5q_0001 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0001_sl-towards_paper/"
                      "911290885_VQE_IdealQasmSim_5q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                      allow_pickle=True)[()]["time tracking subQUBOs"].values()
len_vqe_5q_0001 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0001_sl-towards_paper/"
                          "911290885_VQE_IdealQasmSim_5q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                          allow_pickle=True)[()]["truth solution vector"]

vqe_5q_0002 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0002_sl-towards_paper/"
                      "702678135_VQE_IdealQasmSim_5q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                      allow_pickle=True)[()]["time tracking subQUBOs"].values()
len_vqe_5q_0002 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0002_sl-towards_paper/"
                          "702678135_VQE_IdealQasmSim_5q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                          allow_pickle=True)[()]["truth solution vector"]

vqe_5q_0003 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0003_sl-towards_paper/"
                      "847723440_VQE_IdealQasmSim_5q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                      allow_pickle=True)[()]["time tracking subQUBOs"].values()
len_vqe_5q_0003 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0003_sl-towards_paper/"
                          "847723440_VQE_IdealQasmSim_5q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                          allow_pickle=True)[()]["truth solution vector"]

vqe_5q_0004 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0004_sl-towards_paper/"
                      "380007307_VQE_IdealQasmSim_5q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                      allow_pickle=True)[()]["time tracking subQUBOs"].values()
len_vqe_5q_0004 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0004_sl-towards_paper/"
                          "380007307_VQE_IdealQasmSim_5q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                          allow_pickle=True)[()]["truth solution vector"]

######

vqe_7q_0000 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0000_sl-towards_paper/"
                      "219886525_VQE_IdealQasmSim_7q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                      allow_pickle=True)[()]["time tracking subQUBOs"].values()
vqe_7q_0001 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0001_sl-towards_paper/"
                      "316702619_VQE_IdealQasmSim_7q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                      allow_pickle=True)[()]["time tracking subQUBOs"].values()
vqe_7q_0002 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0002_sl-towards_paper/"
                      "311781405_VQE_IdealQasmSim_7q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                      allow_pickle=True)[()]["time tracking subQUBOs"].values()
vqe_7q_0003 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0003_sl-towards_paper/"
                      "551052192_VQE_IdealQasmSim_7q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                      allow_pickle=True)[()]["time tracking subQUBOs"].values()
vqe_7q_0004 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0004_sl-towards_paper/"
                      "678230537_VQE_IdealQasmSim_7q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                      allow_pickle=True)[()]["time tracking subQUBOs"].values()

len_vqe_7q_0000 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0000_sl-towards_paper/"
                          "219886525_VQE_IdealQasmSim_7q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                          allow_pickle=True)[()]["truth solution vector"]
len_vqe_7q_0001 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0001_sl-towards_paper/"
                          "316702619_VQE_IdealQasmSim_7q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                          allow_pickle=True)[()]["truth solution vector"]
len_vqe_7q_0002 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0002_sl-towards_paper/"
                          "311781405_VQE_IdealQasmSim_7q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                          allow_pickle=True)[()]["truth solution vector"]
len_vqe_7q_0003 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0003_sl-towards_paper/"
                          "551052192_VQE_IdealQasmSim_7q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                          allow_pickle=True)[()]["truth solution vector"]
len_vqe_7q_0004 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0004_sl-towards_paper/"
                          "678230537_VQE_IdealQasmSim_7q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                          allow_pickle=True)[()]["truth solution vector"]

#########

vqe_10q_0000 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0000_sl-towards_paper/"
                       "959787474_VQE_IdealQasmSim_10q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                       allow_pickle=True)[()]["time tracking subQUBOs"].values()
vqe_10q_0001 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0001_sl-towards_paper/"
                       "722418773_VQE_IdealQasmSim_10q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                       allow_pickle=True)[()]["time tracking subQUBOs"].values()
vqe_10q_0002 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0002_sl-towards_paper/"
                       "812656486_VQE_IdealQasmSim_10q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                       allow_pickle=True)[()]["time tracking subQUBOs"].values()
vqe_10q_0003 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0003_sl-towards_paper/"
                       "508903724_VQE_IdealQasmSim_10q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                       allow_pickle=True)[()]["time tracking subQUBOs"].values()
vqe_10q_0004 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0004_sl-towards_paper/"
                       "666973413_VQE_IdealQasmSim_10q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                       allow_pickle=True)[()]["time tracking subQUBOs"].values()

len_vqe_10q_0000 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0000_sl-towards_paper/"
                           "959787474_VQE_IdealQasmSim_10q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                           allow_pickle=True)[()]["truth solution vector"]
len_vqe_10q_0001 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0001_sl-towards_paper/"
                           "722418773_VQE_IdealQasmSim_10q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                           allow_pickle=True)[()]["truth solution vector"]
len_vqe_10q_0002 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0002_sl-towards_paper/"
                           "812656486_VQE_IdealQasmSim_10q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                           allow_pickle=True)[()]["truth solution vector"]
len_vqe_10q_0003 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0003_sl-towards_paper/"
                           "508903724_VQE_IdealQasmSim_10q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                           allow_pickle=True)[()]["truth solution vector"]
len_vqe_10q_0004 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0004_sl-towards_paper/"
                           "666973413_VQE_IdealQasmSim_10q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                           allow_pickle=True)[()]["truth solution vector"]

####################

vqe_12q_0000 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0000_sl-towards_paper/"
                       "473785829_VQE_IdealQasmSim_12q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                       allow_pickle=True)[()]["time tracking subQUBOs"].values()
vqe_12q_0001 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0001_sl-towards_paper/"
                       "401240514_VQE_IdealQasmSim_12q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                       allow_pickle=True)[()]["time tracking subQUBOs"].values()
vqe_12q_0002 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0002_sl-towards_paper/"
                       "561247530_VQE_IdealQasmSim_12q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                       allow_pickle=True)[()]["time tracking subQUBOs"].values()
vqe_12q_0003 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0003_sl-towards_paper/"
                       "814069472_VQE_IdealQasmSim_12q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                       allow_pickle=True)[()]["time tracking subQUBOs"].values()
vqe_12q_0004 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0004_sl-towards_paper/"
                       "235838502_VQE_IdealQasmSim_12q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                       allow_pickle=True)[()]["time tracking subQUBOs"].values()

len_vqe_12q_0000 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0000_sl-towards_paper/"
                           "473785829_VQE_IdealQasmSim_12q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                           allow_pickle=True)[()]["truth solution vector"]
len_vqe_12q_0001 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0001_sl-towards_paper/"
                           "401240514_VQE_IdealQasmSim_12q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                           allow_pickle=True)[()]["truth solution vector"]
len_vqe_12q_0002 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0002_sl-towards_paper/"
                           "561247530_VQE_IdealQasmSim_12q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                           allow_pickle=True)[()]["truth solution vector"]
len_vqe_12q_0003 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0003_sl-towards_paper/"
                           "814069472_VQE_IdealQasmSim_12q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                           allow_pickle=True)[()]["truth solution vector"]
len_vqe_12q_0004 = np.load("../qubo_files/xi_5/smeared/e0gpc_5.0_0004_sl-towards_paper/"
                           "235838502_VQE_IdealQasmSim_12q_TwoLocal_linear_NFT_impact_list_reverse/qubo_log.npy",
                           allow_pickle=True)[()]["truth solution vector"]

timing_vqe_5 = []
timing_vqe_7 = []
timing_vqe_10 = []
timing_vqe_12 = []

len_vqe_5 = []
len_vqe_7 = []
len_vqe_10 = []
len_vqe_12 = []

for t1, t2, t3, t4, t5 in zip(vqe_5q_0000, vqe_5q_0001, vqe_5q_0002, vqe_5q_0003, vqe_5q_0004):
    timing_vqe_5.append(t1)
    timing_vqe_5.append(t2)
    timing_vqe_5.append(t3)
    timing_vqe_5.append(t4)
    timing_vqe_5.append(t5)

for t1, t2, t3, t4, t5 in zip(vqe_7q_0000, vqe_7q_0001, vqe_7q_0002, vqe_7q_0003, vqe_7q_0004):
    timing_vqe_7.append(t1)
    timing_vqe_7.append(t2)
    timing_vqe_7.append(t3)
    timing_vqe_7.append(t4)
    timing_vqe_7.append(t5)

for t1, t2, t3, t4, t5 in zip(vqe_10q_0000, vqe_10q_0001, vqe_10q_0002, vqe_10q_0003, vqe_10q_0004):
    timing_vqe_10.append(t1)
    timing_vqe_10.append(t2)
    timing_vqe_10.append(t3)
    timing_vqe_10.append(t4)
    timing_vqe_10.append(t5)

for t1, t2, t3, t4, t5 in zip(vqe_12q_0000, vqe_12q_0001, vqe_12q_0002, vqe_12q_0003, vqe_12q_0004):
    timing_vqe_12.append(t1)
    timing_vqe_12.append(t2)
    timing_vqe_12.append(t3)
    timing_vqe_12.append(t4)
    timing_vqe_12.append(t5)

####
len_vqe_5.append(len(len_vqe_5q_0000) / 5 )
len_vqe_5.append(len(len_vqe_5q_0001) / 5 )
len_vqe_5.append(len(len_vqe_5q_0002) / 5 )
len_vqe_5.append(len(len_vqe_5q_0003) / 5 )
len_vqe_5.append(len(len_vqe_5q_0004) / 5 )

len_vqe_7.append(len(len_vqe_7q_0000) / 7)
len_vqe_7.append(len(len_vqe_7q_0001) / 7)
len_vqe_7.append(len(len_vqe_7q_0002) / 7)
len_vqe_7.append(len(len_vqe_7q_0003) / 7)
len_vqe_7.append(len(len_vqe_7q_0004) / 7)

len_vqe_10.append(len(len_vqe_10q_0000) / 10)
len_vqe_10.append(len(len_vqe_10q_0001) / 10)
len_vqe_10.append(len(len_vqe_10q_0002) / 10)
len_vqe_10.append(len(len_vqe_10q_0003) / 10)
len_vqe_10.append(len(len_vqe_10q_0004) / 10)

len_vqe_12.append(len(len_vqe_12q_0000) / 12)
len_vqe_12.append(len(len_vqe_12q_0001) / 12)
len_vqe_12.append(len(len_vqe_12q_0002) / 12)
len_vqe_12.append(len(len_vqe_12q_0003) / 12)
len_vqe_12.append(len(len_vqe_12q_0004) / 12)

# print(np.mean(timing_vqe_5), np.std(timing_vqe_5))
# print(np.mean(timing_vqe_7), np.std(timing_vqe_7))
# print(np.mean(timing_vqe_10), np.std(timing_vqe_10))
# print(np.mean(timing_vqe_12), np.std(timing_vqe_12))

sub_QUBO_size = [5, 7, 10, 12]
fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
plt.title("xi=5, BX0-4", fontsize=18, loc="left")
# ax.errorbar(sub_QUBO_size,
#             [np.mean(timing_vqe_5), np.mean(timing_vqe_7), np.mean(timing_vqe_10), np.mean(timing_vqe_12)],
#             yerr=[np.std(timing_vqe_5) / np.sqrt(len(timing_vqe_5)),
#                   np.std(timing_vqe_7) / np.sqrt(len(timing_vqe_7)),
#                   np.std(timing_vqe_10) / np.sqrt(len(timing_vqe_10)),
#                   np.std(timing_vqe_12)  / np.sqrt(len(timing_vqe_12))],
#             color="royalblue",
#             linestyle=" ",
#             marker="o",
#             markersize=6,
#             label="sub-QUBO time")
# ax.legend(loc="lower center", fontsize=14)
# ax.set_xlabel("sub-QUBO size", fontsize=16, loc="right")
# ax.set_ylabel("time [s]", fontsize=16, loc="top")
# ax.set_xticks([5, 7, 10, 12])
# ax.tick_params(axis="both", which="major", labelsize=16)
ax.errorbar(sub_QUBO_size,
             [np.mean(len_vqe_5), np.mean(len_vqe_7), np.mean(len_vqe_10), np.mean(len_vqe_12)],
             yerr=[np.std(len_vqe_5) / np.sqrt(len(len_vqe_5)),
                   np.std(len_vqe_7) / np.sqrt(len(len_vqe_7)),
                   np.std(len_vqe_10) / np.sqrt(len(len_vqe_10)),
                   np.std(len_vqe_12) / np.sqrt(len(len_vqe_12))],
             color="black",
             linestyle=" ",
             marker="^",
             markersize=6,
             label="circuit executions")
ax.set_xlabel("sub-QUBO size", fontsize=16, loc="right")
ax.set_xticks([5, 7, 10, 12])
ax.legend(loc="upper center", fontsize=16)
ax.set_ylabel("circuit executions / iteration", fontsize=16, loc="top")
ax.tick_params(axis="both", which="major", labelsize=16)
plt.savefig("subqubo_timing.pdf", bbox_inches="tight")
