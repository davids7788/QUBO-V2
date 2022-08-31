import matplotlib.pyplot as plt
import ROOT
from ROOT import *
import numpy as np
import sys
sys.path.append("../src")

ROOT.gROOT.LoadMacro("macros/LuxeStyle.C")
ROOT.gROOT.LoadMacro("macros/LuxeLabels.C")

SetLuxeStyle()

folder_bit_flip = sys.argv[1]
folder_eigensolver = sys.argv[2]
folder_qaoa_ideal_qasm_sim = sys.argv[3]
folder_vqe_ideal_qasm_sim_TwoLocal = sys.argv[4]
# folder_vqe_ideal_qasm_sim_HamiltonianDriven = sys.argv[5]

gen_folder = "-".join("/".join(folder_bit_flip.split("/")[0:-1]).split("-")[0:-1])
xi = folder_bit_flip.split("e0gpc_")[1].split("_")[0]


# Load data sets
gen = np.load(gen_folder + "_gen_xplet_list.npy", allow_pickle=True)

reco_bit_flip = np.load(f"{folder_bit_flip}/reco_xplet_list.npy", allow_pickle=True)
reco_eigensolver = np.load(f"{folder_eigensolver}/reco_xplet_list.npy", allow_pickle=True)
reco_qaoa_ideal_qasm_sim = np.load(f"{folder_qaoa_ideal_qasm_sim}/reco_xplet_list.npy", allow_pickle=True)
reco_vqe_ideal_qasm_sim_TwoLocal = np.load(f"{folder_vqe_ideal_qasm_sim_TwoLocal}/reco_xplet_list.npy", allow_pickle=True)
# reco_vqe_ideal_qasm_sim_HamiltonianDriven = np.load(f"{folder_vqe_ideal_qasm_sim_HamiltonianDriven}/reco_xplet_list.npy", allow_pickle=True)

# Prepare matched and fully mismatched
matched_bit_flip = []
matched_eigensolver = []
matched_qaoa_ideal_qasm_sim = []
matched_vqe_ideal_qasm_sim_TwoLocal = []
# matched_vqe_ideal_qasm_sim_HamiltonianDriven = []

fully_mismatched_bit_flip = []
fully_mismatched_eigensolver = []
fully_mismatched_qaoa_ideal_qasm_sim = []
fully_mismatched_vqe_ideal_qasm_sim_TwoLocal = []
# fully_mismatched_vqe_ideal_qasm_sim_HamiltonianDriven = []

for xplet in reco_bit_flip:
    if len(set(xplet.particle_ids.values())) == 1:
        matched_bit_flip.append(xplet)
    elif len(set(xplet.particle_ids.values())) >= len(xplet.hit_ids) - 1:
        fully_mismatched_bit_flip.append(xplet)

for xplet in reco_eigensolver:
    if len(set(xplet.particle_ids.values())) == 1:
        matched_eigensolver.append(xplet)
    elif len(set(xplet.particle_ids.values())) >= len(xplet.hit_ids) - 1:
        fully_mismatched_eigensolver.append(xplet)

for xplet in reco_qaoa_ideal_qasm_sim:
    if len(set(xplet.particle_ids.values())) == 1:
        matched_qaoa_ideal_qasm_sim.append(xplet)
    elif len(set(xplet.particle_ids.values())) >= len(xplet.hit_ids):
        fully_mismatched_qaoa_ideal_qasm_sim.append(xplet)        
        
for xplet in reco_vqe_ideal_qasm_sim_TwoLocal:
    if len(set(xplet.particle_ids.values())) == 1:
        matched_vqe_ideal_qasm_sim_TwoLocal.append(xplet)
    elif len(set(xplet.particle_ids.values())) >= len(xplet.hit_ids) - 1:
        fully_mismatched_vqe_ideal_qasm_sim_TwoLocal.append(xplet)
        
# for xplet in reco_vqe_ideal_qasm_sim_HamiltonianDriven:
#     if len(set(xplet.particle_ids.values())) == 1:
#         matched_vqe_ideal_qasm_sim_HamiltonianDriven.append(xplet)
#     elif len(set(xplet.particle_ids.values())) >= len(xplet.hit_ids) - 1:
#         fully_mismatched_vqe_ideal_qasm_sim_HamiltonianDriven.append(xplet)        


print(f"Generated tracks: {len(gen)}\n")
print(f"Matched xplets bit flip: {len(matched_bit_flip)}")
print(f"Matched xplets Numpy Eigensolver: {len(matched_eigensolver)}")
print(f"Matched xplets QAOA Ideal Sim: {len(matched_qaoa_ideal_qasm_sim)}")
print(f"Matched xplets VQE Ideal Sim TwoLocal: {len(matched_vqe_ideal_qasm_sim_TwoLocal)}")
# print(f"Matched xplets VQE Ideal Sim HamiltonianDriven: {len(matched_vqe_ideal_qasm_sim_HamiltonianDriven)}\n")

print(f"Fully mismatched xplets bit flip: {len(fully_mismatched_bit_flip)}")
print(f"Fully mismatched xplets Numpy Eigensolver: {len(fully_mismatched_eigensolver)}")
print(f"Fully mismatched xplets QAOA Ideal Sim: {len(fully_mismatched_qaoa_ideal_qasm_sim)}")
print(f"Fully mismatched xplets VQE Ideal Sim TwoLocal: {len(fully_mismatched_vqe_ideal_qasm_sim_TwoLocal)}")
# print(f"Fully mismatched xplets VQE Ideal Sim HamiltonianDriven: {len(fully_mismatched_vqe_ideal_qasm_sim_HamiltonianDriven)}\n")


gen_xplets = TH1F('gen xplets', 
                  'x', 50, 0.05275912, 0.55430088)

matched_xplets_bit_flip = TH1F('matched xplets bit flip', 'x', 50, 0.05275912, 0.55430088)
matched_xplets_eigensolver = TH1F('matched xplets Numpy Eigensolver', 'x', 50, 0.05275912, 0.55430088)
matched_xplets_qaoa_ideal_qasm_sim = TH1F('matched xplets QAOA Ideal Sim', 'x', 50, 0.05275912, 0.55430088)
matched_xplets_vqe_ideal_qasm_sim_TwoLocal = TH1F('matched xplets VQE Ideal Sim TwoLocal', 'x', 50, 0.05275912, 0.55430088)
# matched_xplets_vqe_ideal_qasm_sim_HamiltonianDriven = TH1F('matched xplets VQE Ideal Sim HamiltonianDriven', 
#                                                            'x', 50, 0.05275912, 0.55430088)

fully_mismatched_xplets_bit_flip = TH1F('fully mismatched bit flip', 'x', 50, 0.05275912, 0.55430088)
fully_mismatched_xplets_eigensolver = TH1F('fully mismatched Numpy Eigensolver', 'x', 50, 0.05275912, 0.55430088)
fully_mismatched_xplets_qaoa_ideal_qasm_sim = TH1F('fully mismatched QAOA Ideal Sim', 'x', 50, 0.05275912, 0.55430088)
fully_mismatched_xplets_vqe_ideal_qasm_sim_TwoLocal = TH1F('fully mismatched VQE Ideal Sim TwoLocal', 'x', 50, 0.05275912, 0.55430088)
# fully_mismatched_xplets_vqe_ideal_qasm_sim_HamiltonianDriven = TH1F('fully mismatched VQE Ideal Sim HamiltonianDriven', 
#                                                                     'x', 50, 0.05275912, 0.55430088)

reco_xplets_bit_flip = TH1F('reco xplets bit flip', 'x', 50, 0.05275912, 0.55430088)
reco_xplets_eigensolver = TH1F('reco xplets Numpy Eigensolver', 'x', 50, 0.05275912, 0.55430088)
reco_xplets_qaoa_ideal_qasm_sim = TH1F('reco xplets QAOA Ideal Sim', 'x', 50, 0.05275912, 0.55430088)
reco_xplets_vqe_ideal_qasm_sim_TwoLocal = TH1F('reco xplets VQE Ideal Sim TwoLocal', 'x', 50, 0.05275912, 0.55430088)
# reco_xplets_vqe_ideal_qasm_sim_HamiltonianDriven = TH1F('reco xplets VQE Ideal Sim HamiltonianDriven', 'x', 50, 0.05275912, 0.55430088)

h_frame = TH1F('frame', 
               'x', 50, 0.05275912, 0.55430088)


# gen xplets
for g_x in gen:
    gen_xplets.Fill(g_x.coordinates[0][0])


# matched xplets
for m_x in matched_bit_flip:
    matched_xplets_bit_flip.Fill(m_x.coordinates[0][0])
for m_x in matched_eigensolver:
    matched_xplets_eigensolver.Fill(m_x.coordinates[0][0])
for m_x in matched_qaoa_ideal_qasm_sim:
    matched_xplets_qaoa_ideal_qasm_sim.Fill(m_x.coordinates[0][0])
for m_x in matched_vqe_ideal_qasm_sim_TwoLocal:
    matched_xplets_vqe_ideal_qasm_sim_TwoLocal.Fill(m_x.coordinates[0][0])
# for m_x in matched_vqe_ideal_qasm_sim_HamiltonianDriven:
#     matched_xplets_vqe_ideal_qasm_sim_HamiltonianDriven.Fill(m_x.coordinates[0][0])
    

# fully mismatched xplets
for m_x in fully_mismatched_bit_flip:
    fully_mismatched_xplets_bit_flip.Fill(m_x.coordinates[0][0])
for m_x in fully_mismatched_eigensolver:
    fully_mismatched_xplets_eigensolver.Fill(m_x.coordinates[0][0])
for m_x in fully_mismatched_qaoa_ideal_qasm_sim:
    fully_mismatched_xplets_qaoa_ideal_qasm_sim.Fill(m_x.coordinates[0][0])
for m_x in fully_mismatched_vqe_ideal_qasm_sim_TwoLocal:
    fully_mismatched_xplets_vqe_ideal_qasm_sim_TwoLocal.Fill(m_x.coordinates[0][0])
# for m_x in fully_mismatched_vqe_ideal_qasm_sim_HamiltonianDriven:
#     fully_mismatched_xplets_vqe_ideal_qasm_sim_HamiltonianDriven.Fill(m_x.coordinates[0][0])
    
# reco_xplets
for r_x in reco_bit_flip:
    reco_xplets_bit_flip.Fill(r_x.coordinates[0][0])
for r_x in reco_eigensolver:
    reco_xplets_eigensolver.Fill(r_x.coordinates[0][0])
for r_x in reco_qaoa_ideal_qasm_sim:
    reco_xplets_qaoa_ideal_qasm_sim.Fill(r_x.coordinates[0][0])
for r_x in reco_vqe_ideal_qasm_sim_TwoLocal:
    reco_xplets_vqe_ideal_qasm_sim_TwoLocal.Fill(r_x.coordinates[0][0])
# for r_x in reco_vqe_ideal_qasm_sim_HamiltonianDriven:
#     reco_xplets_vqe_ideal_qasm_sim_HamiltonianDriven.Fill(r_x.coordinates[0][0])
    
Teff_bit_flip = TEfficiency(matched_xplets_bit_flip, gen_xplets)
Teff_eigensolver = TEfficiency(matched_xplets_eigensolver, gen_xplets)
Teff_qaoa_ideal_qasm_sim = TEfficiency(matched_xplets_qaoa_ideal_qasm_sim, gen_xplets)
Teff_vqe_ideal_qasm_sim_TwoLocal = TEfficiency(matched_xplets_vqe_ideal_qasm_sim_TwoLocal, gen_xplets)
# Teff_vqe_ideal_qasm_sim_HamiltonianDriven = TEfficiency(matched_xplets_vqe_ideal_qasm_sim_HamiltonianDriven, gen_xplets)

Teff_fully_mis_bit_flip = TEfficiency(fully_mismatched_xplets_bit_flip, reco_xplets_bit_flip)
Teff_fully_mis_eigensolver = TEfficiency(fully_mismatched_xplets_eigensolver, reco_xplets_eigensolver)
Teff_fully_mis_qaoa_ideal_qasm_sim = TEfficiency(fully_mismatched_xplets_qaoa_ideal_qasm_sim, reco_xplets_qaoa_ideal_qasm_sim)
Teff_fully_mis_vqe_ideal_qasm_sim_TwoLocal = TEfficiency(fully_mismatched_xplets_vqe_ideal_qasm_sim_TwoLocal,
                                                         reco_xplets_vqe_ideal_qasm_sim_TwoLocal)
# Teff_fully_mis_vqe_ideal_qasm_sim_HamiltonianDrive = TEfficiency(fully_mismatched_xplets_vqe_ideal_qasm_sim_HamiltonianDriven, 
#                                                                  reco_xplets_vqe_ideal_qasm_sim_HamiltonianDriven)

gStyle.SetPadTickY(0)   
canv = TCanvas("example","xplet efficiency ", 800, 600)

h_frame.GetYaxis().SetTitle("fraction of counts")
h_frame.GetXaxis().SetTitle("x [m]")
h_frame.GetXaxis().SetNdivisions(105)
h_frame.SetMaximum(1.05)
h_frame.GetXaxis().SetLabelOffset(0.02)
h_frame.GetXaxis().SetRangeUser(0.05, 0.51)
h_frame.Draw()

# bit flip
Teff_bit_flip.SetMarkerStyle(20)
Teff_bit_flip.SetMarkerColorAlpha(1, 0.85)
Teff_bit_flip.SetMarkerSize(1.5)
Teff_bit_flip.Draw("PSAME")

Teff_fully_mis_bit_flip.SetMarkerStyle(24)
Teff_fully_mis_bit_flip.SetMarkerColorAlpha(1, 0.85)
Teff_fully_mis_bit_flip.SetMarkerSize(1.5)
Teff_fully_mis_bit_flip.Draw("PSAME")

# Eigensolver
Teff_eigensolver.SetMarkerStyle(21)
Teff_eigensolver.SetMarkerColorAlpha(2, 0.85)
Teff_eigensolver.SetMarkerSize(1.5)
Teff_eigensolver.Draw("PSAME")

Teff_fully_mis_eigensolver.SetMarkerStyle(25)
Teff_fully_mis_eigensolver.SetMarkerColorAlpha(2, 0.85)
Teff_fully_mis_eigensolver.SetMarkerSize(1.5)
Teff_fully_mis_eigensolver.Draw("PSAME")

# QAOA
Teff_qaoa_ideal_qasm_sim.SetMarkerStyle(22)
Teff_qaoa_ideal_qasm_sim.SetMarkerColorAlpha(3, 0.85)
Teff_qaoa_ideal_qasm_sim.SetMarkerSize(1.95)
Teff_qaoa_ideal_qasm_sim.Draw("PSAME")

Teff_fully_mis_qaoa_ideal_qasm_sim.SetMarkerStyle(26)
Teff_fully_mis_qaoa_ideal_qasm_sim.SetMarkerColorAlpha(3, 0.85)
Teff_fully_mis_qaoa_ideal_qasm_sim.SetMarkerSize(1.95)
Teff_fully_mis_qaoa_ideal_qasm_sim.Draw("PSAME")

# VQE TwoLocal
Teff_vqe_ideal_qasm_sim_TwoLocal.SetMarkerStyle(23)
Teff_vqe_ideal_qasm_sim_TwoLocal.SetMarkerColorAlpha(4, 0.85)
Teff_vqe_ideal_qasm_sim_TwoLocal.SetMarkerSize(1.5)
Teff_vqe_ideal_qasm_sim_TwoLocal.Draw("PSAME")

Teff_fully_mis_vqe_ideal_qasm_sim_TwoLocal.SetMarkerStyle(32)
Teff_fully_mis_vqe_ideal_qasm_sim_TwoLocal.SetMarkerColorAlpha(4, 0.85)
Teff_fully_mis_vqe_ideal_qasm_sim_TwoLocal.SetMarkerSize(1.5)
Teff_fully_mis_vqe_ideal_qasm_sim_TwoLocal.Draw("PSAME")

gPad.SetRightMargin(0.13)
gPad.Draw()

rightmax = gen_xplets.GetMaximum()
scale = gPad.GetUymax() / rightmax

gen_xplets.SetLineColorAlpha(33, 0.85)
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


LUXELabel(0.2, 0.65);

leg = TLegend(0.2, 0.3, 0.6, 0.5)
leg.SetNColumns(2)

leg.AddEntry(Teff_bit_flip, "Bit Flip", "p")
leg.AddEntry(Teff_fully_mis_bit_flip, "Bit Flip", "p")

leg.AddEntry(Teff_eigensolver, "Eigensolver", "p")
leg.AddEntry(Teff_fully_mis_eigensolver, "Eigensolver", "p")

leg.AddEntry(Teff_qaoa_ideal_qasm_sim, "QAOA", "p")
leg.AddEntry(Teff_fully_mis_qaoa_ideal_qasm_sim, "QAOA", "p")

leg.AddEntry(Teff_vqe_ideal_qasm_sim_TwoLocal,"VQE", "p")
leg.AddEntry(Teff_fully_mis_vqe_ideal_qasm_sim_TwoLocal,"VQE", "p")

# leg.AddEntry(Teff_vqe_ideal_qasm_sim_HamiltonianDriven, "VQE HamiltonianDriven", "p")
# leg.AddEntry(Teff_fully_mis_vqe_ideal_qasm_sim_HamiltonianDriven, "VQE HamiltonianDriven", "p")

leg.AddEntry(gen_xplets, "gen xplets", "l")
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetTextSize(0.04)
leg.SetHeader(f"#xi = {xi} eff vs. mismatched")
leg.Draw()


gen_out = "/".join(gen_folder.split("/")[0:-1])
canv.SaveAs(f"{gen_out}/1Dxplet_efficiency_comparison_xi_{xi}_detector_position.pdf")