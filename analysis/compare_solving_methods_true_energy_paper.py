import matplotlib.pyplot as plt
import ROOT
from ROOT import *
import numpy as np
import sys
sys.path.append("../src")

ROOT.gROOT.LoadMacro("macros/LuxeStyle.C")
ROOT.gROOT.LoadMacro("macros/LuxeLabels.C")

SetLuxeStyle()

mismatched_def = int(sys.argv[1])
folder_eigensolver_7q = sys.argv[2]
folder_vqe_ideal_qasm_sim_TwoLocal = sys.argv[3]
folder_eigensolver_12q = sys.argv[4]
folder_eigensolver_14q = sys.argv[5]
folder_eigensolver_18q = sys.argv[6]
folder_bit_flip = sys.argv[7]

gen_folder = "-".join("/".join(folder_bit_flip.split("/")[0:-1]).split("-")[0:-1])
xi = folder_bit_flip.split("e0gpc_")[1].split("_")[0]


# Load data sets
gen = np.load(gen_folder + "_gen_xplet_list.npy", allow_pickle=True)
print(f"Generated tracks: {len(gen)}\n")
edges_left = [1.5]
edges_mid = [2.0 + 0.4 * i for i in range(15)]
edges_right = [8.0, 8.5, 9.5]
edges = np.array(edges_left + edges_mid + edges_right)

gen_xplets = TH1F('gen xplets', 'energy', len(edges) - 1, edges)

gStyle.SetPadTickY(0)
canv = TCanvas("example","xplet efficiency ", 800, 600)

h_frame = TH1F('frame', 'x', len(edges) - 1, edges)
h_frame.GetYaxis().SetTitle("fraction of counts")
h_frame.GetXaxis().SetTitle("energy [GeV]")
h_frame.GetXaxis().SetNdivisions(105)
h_frame.SetMaximum(1.18)
h_frame.GetXaxis().SetLabelOffset(0.02)
h_frame.GetXaxis().SetRangeUser(0.0, 11.0)
h_frame.Draw()


for g_x in gen:
    gen_xplets.Fill(g_x.energy[0] / 1000)


if folder_eigensolver_7q != "-1":
    reco_eigensolver_7q = np.load(f"{folder_eigensolver_7q}/reco_xplet_list.npy",
                                  allow_pickle=True)
    matched_eigensolver_7q = []
    fully_mismatched_eigensolver_7q = []
    for xplet in reco_eigensolver_7q:
        if len(set(xplet.particle_ids.values())) == 1:
            matched_eigensolver_7q.append(xplet)
        elif len(set(xplet.particle_ids.values())) > mismatched_def:
            fully_mismatched_eigensolver_7q.append(xplet)
    print(f"Efficiency on xplets for Numpy Eigensolver 7 qubits: "
          f"{np.around(len(matched_eigensolver_7q) / len(gen), 3)}")
    print(f"Mismatched (> {mismatched_def} particle IDs / xplet) for Numpy Eigensolver 7 qubits: "
          f"{np.around(len(fully_mismatched_eigensolver_7q) / len(reco_eigensolver_7q), 2)}")
    matched_xplets_eigensolver_7q = TH1F('matched xplets Numpy Eigensolver 7q', 'x', len(edges) - 1, edges)
    fully_mismatched_xplets_eigensolver_7q = TH1F('mismatched Numpy Eigensolver 7q', 'x', len(edges) - 1, edges)
    reco_xplets_eigensolver_7q = TH1F('reco xplets Numpy Eigensolver 7q', 'x', len(edges) - 1, edges)
    for m_x in matched_eigensolver_7q:
        matched_xplets_eigensolver_7q.Fill(m_x.energy[0] / 1000)
    for m_x in fully_mismatched_eigensolver_7q:
        fully_mismatched_xplets_eigensolver_7q.Fill(m_x.energy[0] / 1000)
    for r_x in reco_eigensolver_7q:
        reco_xplets_eigensolver_7q.Fill(r_x.energy[0] / 1000)
    Teff_eigensolver_7q = TEfficiency(matched_xplets_eigensolver_7q, gen_xplets)
    Teff_fully_mis_eigensolver_7q = TEfficiency(fully_mismatched_xplets_eigensolver_7q, reco_xplets_eigensolver_7q)
    Teff_eigensolver_7q.SetMarkerStyle(20)
    Teff_eigensolver_7q.SetMarkerColorAlpha(2, 0.85)
    Teff_eigensolver_7q.SetMarkerSize(1.5)
    Teff_eigensolver_7q.Draw("PSAME")

    Teff_fully_mis_eigensolver_7q.SetMarkerStyle(24)
    Teff_fully_mis_eigensolver_7q.SetMarkerColorAlpha(2, 0.85)
    Teff_fully_mis_eigensolver_7q.SetMarkerSize(1.5)
    Teff_fully_mis_eigensolver_7q.Draw("PSAME")




if folder_eigensolver_12q != "-1":
    reco_eigensolver_12q = np.load(f"{folder_eigensolver_12q}/reco_xplet_list.npy",
                                   allow_pickle=True)
    matched_eigensolver_12q = []
    fully_mismatched_eigensolver_12q = []
    for xplet in reco_eigensolver_12q:
        if len(set(xplet.particle_ids.values())) == 1:
            matched_eigensolver_12q.append(xplet)
        elif len(set(xplet.particle_ids.values())) > mismatched_def:
            fully_mismatched_eigensolver_12q.append(xplet)
    print(f"Efficiency on xplets for Numpy Eigensolver 12 qubits: "
          f"{np.around(len(matched_eigensolver_12q) / len(gen), 3)}")
    print(f"Mismatched (> {mismatched_def} particle IDs / xplet) for Numpy Eigensolver 12 qubits: "
          f"{np.around(len(fully_mismatched_eigensolver_12q) / len(reco_eigensolver_12q), 3)}")
    matched_xplets_eigensolver_12q = TH1F('matched xplets Numpy Eigensolver 12q', 'x', len(edges) - 1, edges)
    fully_mismatched_xplets_eigensolver_12q = TH1F('mismatched Numpy Eigensolver 12q', 'x', len(edges) - 1, edges)
    reco_xplets_eigensolver_12q = TH1F('reco xplets Numpy Eigensolver 12q', 'x', len(edges) - 1, edges)
    for m_x in matched_eigensolver_12q:
        matched_xplets_eigensolver_12q.Fill(m_x.energy[0] / 1000)
    for m_x in fully_mismatched_eigensolver_12q:
        fully_mismatched_xplets_eigensolver_12q.Fill(m_x.energy[0] / 1000)
    for r_x in reco_eigensolver_12q:
        reco_xplets_eigensolver_12q.Fill(r_x.energy[0] / 1000)
    Teff_eigensolver_12q = TEfficiency(matched_xplets_eigensolver_12q, gen_xplets)
    Teff_fully_mis_eigensolver_12q = TEfficiency(fully_mismatched_xplets_eigensolver_12q, reco_xplets_eigensolver_12q)

    Teff_eigensolver_12q.SetMarkerStyle(21)
    Teff_eigensolver_12q.SetMarkerColorAlpha(3, 0.85)
    Teff_eigensolver_12q.SetMarkerSize(1.5)
    Teff_eigensolver_12q.Draw("PSAME")

    Teff_fully_mis_eigensolver_12q.SetMarkerStyle(25)
    Teff_fully_mis_eigensolver_12q.SetMarkerColorAlpha(3, 0.85)
    Teff_fully_mis_eigensolver_12q.SetMarkerSize(1.5)
    Teff_fully_mis_eigensolver_12q.Draw("PSAME")


if folder_eigensolver_14q != "-1":
    reco_eigensolver_14q = np.load(f"{folder_eigensolver_14q}/reco_xplet_list.npy",
                                   allow_pickle=True)
    matched_eigensolver_14q = []
    fully_mismatched_eigensolver_14q = []
    for xplet in reco_eigensolver_14q:
        if len(set(xplet.particle_ids.values())) == 1:
            matched_eigensolver_14q.append(xplet)
        elif len(set(xplet.particle_ids.values())) > mismatched_def:
            fully_mismatched_eigensolver_14q.append(xplet)
    print(f"Efficiency on xplets for Numpy Eigensolver 14 qubits: "
          f"{np.around(len(matched_eigensolver_14q) / len(gen), 3)}")
    print(f"Mismatched (> {mismatched_def} particle IDs / xplet) for Numpy Eigensolver 14 qubits: "
          f"{np.around(len(fully_mismatched_eigensolver_14q) / len(reco_eigensolver_14q), 3)}")
    matched_xplets_eigensolver_14q = TH1F('matched xplets Numpy Eigensolver 14q', 'x', len(edges) - 1, edges)
    fully_mismatched_xplets_eigensolver_14q = TH1F('mismatched Numpy Eigensolver 14q', 'x', len(edges) - 1, edges)
    reco_xplets_eigensolver_14q = TH1F('reco xplets Numpy Eigensolver 14q', 'x', len(edges) - 1, edges)
    for m_x in matched_eigensolver_14q:
        matched_xplets_eigensolver_14q.Fill(m_x.energy[0] / 1000)
    for m_x in fully_mismatched_eigensolver_14q:
        fully_mismatched_xplets_eigensolver_14q.Fill(m_x.energy[0] / 1000)
    for r_x in reco_eigensolver_14q:
        reco_xplets_eigensolver_14q.Fill(r_x.energy[0] / 1000)
    Teff_eigensolver_14q = TEfficiency(matched_xplets_eigensolver_14q, gen_xplets)
    Teff_fully_mis_eigensolver_14q = TEfficiency(fully_mismatched_xplets_eigensolver_14q, reco_xplets_eigensolver_14q)

    Teff_eigensolver_14q.SetMarkerStyle(22)
    Teff_eigensolver_14q.SetMarkerColorAlpha(4, 0.85)
    Teff_eigensolver_14q.SetMarkerSize(1.5)
    Teff_eigensolver_14q.Draw("PSAME")

    Teff_fully_mis_eigensolver_14q.SetMarkerStyle(26)
    Teff_fully_mis_eigensolver_14q.SetMarkerColorAlpha(4, 0.85)
    Teff_fully_mis_eigensolver_14q.SetMarkerSize(1.5)
    Teff_fully_mis_eigensolver_14q.Draw("PSAME")
if folder_eigensolver_18q != "-1":
    reco_eigensolver_18q = np.load(f"{folder_eigensolver_18q}/reco_xplet_list.npy",
                                   allow_pickle=True)
    matched_eigensolver_18q = []
    fully_mismatched_eigensolver_18q = []
    for xplet in reco_eigensolver_18q:
        if len(set(xplet.particle_ids.values())) == 1:
            matched_eigensolver_18q.append(xplet)
        elif len(set(xplet.particle_ids.values())) > mismatched_def:
            fully_mismatched_eigensolver_18q.append(xplet)
    print(f"Efficiency on xplets for Numpy Eigensolver 18 qubits: "
          f"{np.around(len(matched_eigensolver_18q) / len(gen), 3)}")
    print(f"Mismatched (> {mismatched_def} particle IDs / xplet) for Numpy Eigensolver 18 qubits: "
          f"{np.around(len(fully_mismatched_eigensolver_18q) / len(reco_eigensolver_18q), 3)}")
    matched_xplets_eigensolver_18q = TH1F('matched xplets Numpy Eigensolver 18q', 'x', len(edges) - 1, edges)
    fully_mismatched_xplets_eigensolver_18q = TH1F('mismatched Numpy Eigensolver 18q', 'x', len(edges) - 1, edges)
    reco_xplets_eigensolver_18q = TH1F('reco xplets Numpy Eigensolver 18q', 'x', len(edges) - 1, edges)
    for m_x in matched_eigensolver_18q:
        matched_xplets_eigensolver_18q.Fill(m_x.energy[0] / 1000)
    for m_x in fully_mismatched_eigensolver_18q:
        fully_mismatched_xplets_eigensolver_18q.Fill(m_x.energy[0] / 1000)
    for r_x in reco_eigensolver_18q:
        reco_xplets_eigensolver_18q.Fill(r_x.energy[0] / 1000)
    Teff_eigensolver_18q = TEfficiency(matched_xplets_eigensolver_18q, gen_xplets)
    Teff_fully_mis_eigensolver_18q = TEfficiency(fully_mismatched_xplets_eigensolver_18q, reco_xplets_eigensolver_18q)

    Teff_eigensolver_18q.SetMarkerStyle(23)
    Teff_eigensolver_18q.SetMarkerColorAlpha(6, 0.85)
    Teff_eigensolver_18q.SetMarkerSize(1.5)
    Teff_eigensolver_18q.Draw("PSAME")

    Teff_fully_mis_eigensolver_18q.SetMarkerStyle(32)
    Teff_fully_mis_eigensolver_18q.SetMarkerColorAlpha(6, 0.85)
    Teff_fully_mis_eigensolver_18q.SetMarkerSize(1.5)
    Teff_fully_mis_eigensolver_18q.Draw("PSAME")


if folder_vqe_ideal_qasm_sim_TwoLocal != "-1":
    reco_vqe_ideal_qasm_sim_TwoLocal = np.load(f"{folder_vqe_ideal_qasm_sim_TwoLocal}/reco_xplet_list.npy",
                                               allow_pickle=True)
    matched_vqe_ideal_qasm_sim_TwoLocal = []
    fully_mismatched_vqe_ideal_qasm_sim_TwoLocal = []
    for xplet in reco_vqe_ideal_qasm_sim_TwoLocal:
        if len(set(xplet.particle_ids.values())) == 1:
            matched_vqe_ideal_qasm_sim_TwoLocal.append(xplet)
        elif len(set(xplet.particle_ids.values())) > mismatched_def:
            fully_mismatched_vqe_ideal_qasm_sim_TwoLocal.append(xplet)
    print(f"Efficiency on xplets for VQE Ideal Sim TwoLocal: "
          f"{np.around(len(matched_vqe_ideal_qasm_sim_TwoLocal) / len(gen), 3)}")
    print(f"Mismatched (> {mismatched_def} particle IDs / xplet) for VQE Ideal Sim TwoLocal "
          f"{np.around(len(fully_mismatched_vqe_ideal_qasm_sim_TwoLocal) / len(reco_vqe_ideal_qasm_sim_TwoLocal), 3)}")
    matched_xplets_vqe_ideal_qasm_sim_TwoLocal = TH1F('matched xplets VQE Ideal Sim TwoLocal', 'x', len(edges) - 1,
                                                      edges)
    fully_mismatched_xplets_vqe_ideal_qasm_sim_TwoLocal = TH1F('mismatched VQE Ideal Sim TwoLocal', 'x', len(edges) - 1,
                                                               edges)
    reco_xplets_vqe_ideal_qasm_sim_TwoLocal = TH1F('reco xplets VQE Ideal Sim TwoLocal', 'x', len(edges) - 1, edges)
    for m_x in matched_vqe_ideal_qasm_sim_TwoLocal:
        matched_xplets_vqe_ideal_qasm_sim_TwoLocal.Fill(m_x.energy[0] / 1000)
    for m_x in fully_mismatched_vqe_ideal_qasm_sim_TwoLocal:
        fully_mismatched_xplets_vqe_ideal_qasm_sim_TwoLocal.Fill(m_x.energy[0] / 1000)
    for r_x in reco_vqe_ideal_qasm_sim_TwoLocal:
        reco_xplets_vqe_ideal_qasm_sim_TwoLocal.Fill(r_x.energy[0] / 1000)
    Teff_vqe_ideal_qasm_sim_TwoLocal = TEfficiency(matched_xplets_vqe_ideal_qasm_sim_TwoLocal, gen_xplets)
    Teff_fully_mis_vqe_ideal_qasm_sim_TwoLocal = TEfficiency(fully_mismatched_xplets_vqe_ideal_qasm_sim_TwoLocal,
                                                             reco_xplets_vqe_ideal_qasm_sim_TwoLocal)

    Teff_vqe_ideal_qasm_sim_TwoLocal.SetMarkerStyle(33)
    Teff_vqe_ideal_qasm_sim_TwoLocal.SetMarkerColorAlpha(7, 0.85)
    Teff_vqe_ideal_qasm_sim_TwoLocal.SetMarkerSize(1.5)
    Teff_vqe_ideal_qasm_sim_TwoLocal.Draw("PSAME")

    Teff_fully_mis_vqe_ideal_qasm_sim_TwoLocal.SetMarkerStyle(27)
    Teff_fully_mis_vqe_ideal_qasm_sim_TwoLocal.SetMarkerColorAlpha(7, 0.85)
    Teff_fully_mis_vqe_ideal_qasm_sim_TwoLocal.SetMarkerSize(1.5)
    Teff_fully_mis_vqe_ideal_qasm_sim_TwoLocal.Draw("PSAME")

if folder_bit_flip != "-1":
    reco_bit_flip = np.load(f"{folder_bit_flip}/reco_xplet_list.npy", allow_pickle=True)
    matched_bit_flip = []
    fully_mismatched_bit_flip = []
    for xplet in reco_bit_flip:
        if len(set(xplet.particle_ids.values())) == 1:
            matched_bit_flip.append(xplet)
        elif len(set(xplet.particle_ids.values())) > mismatched_def:
            fully_mismatched_bit_flip.append(xplet)
    print(f"Efficiency on xplets for bit flip search: "
          f"{np.around(len(matched_bit_flip) / len(gen), 3)}")
    print(f"Mismatched (> {mismatched_def} particle IDs / xplet) for bit flip search "
          f"{np.around(len(fully_mismatched_bit_flip) / len(reco_bit_flip), 3)}")
    matched_xplets_bit_flip = TH1F('matched xplets bit flip search', 'x', len(edges) - 1, edges)
    fully_mismatched_xplets_bit_flip = TH1F('mismatched bit flip search', 'x', len(edges) - 1, edges)
    reco_xplets_bit_flip = TH1F('reco xplets bit flip search', 'x', len(edges) - 1, edges)
    for m_x in matched_bit_flip:
        matched_xplets_bit_flip.Fill(m_x.energy[0] / 1000)
    for m_x in fully_mismatched_bit_flip:
        fully_mismatched_xplets_bit_flip.Fill(m_x.energy[0] / 1000)
    for r_x in reco_bit_flip:
        reco_xplets_bit_flip.Fill(r_x.energy[0] / 1000)
    Teff_bit_flip = TEfficiency(matched_xplets_bit_flip, gen_xplets)
    Teff_fully_mis_bit_flip = TEfficiency(fully_mismatched_xplets_bit_flip, reco_xplets_bit_flip)

    Teff_bit_flip.SetMarkerStyle(34)
    Teff_bit_flip.SetMarkerColorAlpha(1, 0.85)
    Teff_bit_flip.SetMarkerSize(1.5)
    Teff_bit_flip.Draw("PSAME")

    Teff_fully_mis_bit_flip.SetMarkerStyle(28)
    Teff_fully_mis_bit_flip.SetMarkerColorAlpha(1, 0.85)
    Teff_fully_mis_bit_flip.SetMarkerSize(1.5)
    Teff_fully_mis_bit_flip.Draw("PSAME")



gPad.SetRightMargin(0.13)
gPad.Draw()

rightmax = gen_xplets.GetMaximum() * 1.18
scale = gPad.GetUymax() / rightmax

gen_xplets.SetLineColorAlpha(33, 0.99)
gen_xplets.SetLineWidth(2)
gen_xplets.Scale(scale)
gen_xplets.Draw("HISTSAME")

axis = TGaxis(gPad.GetUxmax(),
              gPad.GetUymin(),
              gPad.GetUxmax(),
              gPad.GetUymax(),
              0,
              rightmax,
              510,
              "+L")
axis.SetTitle("counts")
axis.SetTitleSize(0.05)
axis.SetLabelFont(42)
axis.SetTextFont(42)
axis.SetLabelSize(0.05)
axis.SetTitleOffset(1.4)
axis.Draw()


LUXELabel(0.18, 0.87)

leg = TLegend(0.63, 0.23, 0.78, 0.65)

latex = TLatex()
latex.SetTextSize(0.05)
latex.SetTextFont(42)
latex.DrawLatex(3.7, 1.07, f"e-laser, #xi = {xi}, N_{{xplets,gen}}={len(gen)}")

if folder_eigensolver_7q != "-1":
    leg.AddEntry(Teff_eigensolver_7q, "Eigensolver (Q7)", "p")
    leg.AddEntry(Teff_fully_mis_eigensolver_7q, "Eigensolver (Q7)", "p")

if folder_eigensolver_12q != "-1":
    leg.AddEntry(Teff_eigensolver_12q, "Eigensolver (Q12)", "p")
    leg.AddEntry(Teff_fully_mis_eigensolver_12q, "Eigensolver (Q12)", "p")

if folder_eigensolver_14q != "-1":
    leg.AddEntry(Teff_eigensolver_14q, "Eigensolver (Q14)", "p")
    leg.AddEntry(Teff_fully_mis_eigensolver_14q, "Eigensolver (Q14)", "p")

if folder_eigensolver_18q != "-1":
    leg.AddEntry(Teff_eigensolver_18q, "Eigensolver (Q18)", "p")
    leg.AddEntry(Teff_fully_mis_eigensolver_18q, "Eigensolver (Q18)", "p")

if folder_vqe_ideal_qasm_sim_TwoLocal != "-1":
    leg.AddEntry(Teff_vqe_ideal_qasm_sim_TwoLocal, "VQE (Q7)", "p")
    leg.AddEntry(Teff_fully_mis_vqe_ideal_qasm_sim_TwoLocal, "VQE (Q7)", "p")

if folder_bit_flip != "-1":
    leg.AddEntry(Teff_bit_flip, "bit flip search", "p")
    leg.AddEntry(Teff_fully_mis_bit_flip, "bit flip search", "p")

leg.AddEntry(gen_xplets, "gen xplets", "l")
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetTextSize(0.03)
leg.Draw()



gen_out = "/".join(gen_folder.split("/")[0:-1])
canv.SaveAs(f"{gen_out}/1Dxplet_efficiency_comparison_xi_{xi}_energy.pdf")