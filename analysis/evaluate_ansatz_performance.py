from ROOT import *
import numpy as np
import sys

sys.path.append("../src")

vqe_NFT_7q = np.load("../test_results_qubo/e0gpc_4.0_0000_sl-c_4/"
                     "443062997_VQE_IdealQasmSim_7q_TwoLocal_NFT/qubo_log.npy", allow_pickle=True)

vqe_NFT_7q_deep = np.load("../test_results_qubo/e0gpc_4.0_0000_sl-c_4/"
                          "967526882_VQE_IdealQasmSim_7q_TwoLocal_NFT_deep/qubo_log.npy", allow_pickle=True)

vqe_NFT_7q_full = np.load("../test_results_qubo/e0gpc_4.0_0000_sl-c_4/"
                          "537331657_VQE_IdealQasmSim_7q_TwoLocal_NFT_full/qubo_log.npy", allow_pickle=True)

vqe_NFT_7q_circular = np.load("../test_results_qubo/e0gpc_4.0_0000_sl-c_4/"
                              "103357245_VQE_IdealQasmSim_7q_TwoLocal_NFT_circular/qubo_log.npy", allow_pickle=True)


def evaluate_sample(sample, description):
    num_subqubos = 0
    num_correctly_solved_subqubos = 0

    for subQUBO in sample[()]["compare to analytical solution"].values():
        if num_subqubos == 5000:
            break
        if subQUBO:
            num_correctly_solved_subqubos += 1
        num_subqubos += 1

    print(f"Results for {description}:")
    print(f"{np.around(100 * num_correctly_solved_subqubos / num_subqubos, 2)} +/-  "
          f"{np.around(100 * np.sqrt(num_correctly_solved_subqubos) / num_subqubos, 2)} % "
          f"----  evaluated for {num_subqubos} sub-QUBOs")
    print("------------------------------------------------------------------")


evaluate_sample(vqe_NFT_7q, "baseline")
evaluate_sample(vqe_NFT_7q_circular, "circular")
evaluate_sample(vqe_NFT_7q_full, "full")
evaluate_sample(vqe_NFT_7q_deep, "deep")
