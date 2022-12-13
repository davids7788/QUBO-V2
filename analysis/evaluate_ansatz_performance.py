from ROOT import *
import numpy as np
import sys
import matplotlib.pyplot as plt

sys.path.append("../src")

vqe_NFT_7q_linear = np.load("../test_qubos/e0gpc_5.0_0000_sl-c_4/"
                            "332757277_VQE_IdealQasmSim_7q_TwoLocal_NFT_linear/qubo_log.npy", allow_pickle=True)

vqe_NFT_7q_linear_deep = np.load("../test_qubos/e0gpc_5.0_0000_sl-c_4/"
                                 "710478200_VQE_IdealQasmSim_7q_TwoLocal_NFT_linear_deep/qubo_log.npy", allow_pickle=True)

vqe_NFT_7q_circular = np.load("../test_qubos/e0gpc_5.0_0000_sl-c_4/"
                              "961945021_VQE_IdealQasmSim_7q_TwoLocal_NFT_circular/qubo_log.npy", allow_pickle=True)

vqe_NFT_7q_circular_deep = np.load("../test_qubos/e0gpc_5.0_0000_sl-c_4/"
                                   "241986228_VQE_IdealQasmSim_7q_TwoLocal_NFT_circular_deep/qubo_log.npy", allow_pickle=True)

vqe_NFT_7q_full = np.load("../test_qubos/e0gpc_5.0_0000_sl-c_4/"
                          "853664635_VQE_IdealQasmSim_7q_TwoLocal_NFT_full/qubo_log.npy", allow_pickle=True)

vqe_NFT_7q_full_deep = np.load("../test_qubos/e0gpc_5.0_0000_sl-c_4/"
                               "250965046_VQE_IdealQasmSim_7q_TwoLocal_NFT_full_deep/qubo_log.npy", allow_pickle=True)


def evaluate_sample(sample, description):
    num_subqubos = 0
    num_correctly_solved_subqubos = 0
    processing_time = []

    for subQUBO in sample[()]["compare to analytical solution"].values():
        if num_subqubos == 6188:
            break
        if subQUBO:
            num_correctly_solved_subqubos += 1
        num_subqubos += 1
    for subQUBO_time in sample[()]["time tracking subQUBOs"].values():

        processing_time.append(subQUBO_time)
    print(f"Results for {description}:")
    print(f"{np.around(np.mean(processing_time), 2)} +/- {np.around(np.std(processing_time), 2) }")
    print(f"{np.around(100 * num_correctly_solved_subqubos / num_subqubos, 2)} +/-  "
          f"{np.around(100 * np.sqrt(num_correctly_solved_subqubos) / num_subqubos, 2)} % "
          f"----  evaluated for {num_subqubos} sub-QUBOs")
    print("------------------------------------------------------------------")


evaluate_sample(vqe_NFT_7q_linear, "linear")
evaluate_sample(vqe_NFT_7q_circular, "circular")
evaluate_sample(vqe_NFT_7q_full, "full")
evaluate_sample(vqe_NFT_7q_linear_deep, "linear deep")
evaluate_sample(vqe_NFT_7q_circular_deep, "circular deep")
evaluate_sample(vqe_NFT_7q_full_deep, "full deep")