import ROOT
from ROOT import *
import numpy as np
import sys
sys.path.append("../src")

ROOT.gROOT.LoadMacro("macros/LuxeStyle.C")
ROOT.gROOT.LoadMacro("macros/LuxeLabels.C")

SetLuxeStyle()

folder = sys.argv[1]
if "eigensolver" in folder:
    label = "Numpy Eigensolver"
if "bit_flip_only" in folder:
    label = "bit flip search"
if "VQE" in folder:
    label = "VQE"
if "QAOA" in folder:
    label = "QAOA"

gen = "-".join("/".join(folder.split("/")[0:-1]).split("-")[0:-1])
xi = folder.split("e0gpc_")[1].split("_")[0]

gen_x = np.load(gen + "_gen_xplet_list.npy", allow_pickle=True)
reco_x = np.load(f"{folder}/reco_xplet_list.npy", allow_pickle=True)

matched = []
fake = []

for t in reco_x:
    if len(set(t.particle_ids.values())) == 1:
        matched.append(t)
    else:
        fake.append(t)

print(f"Generated tracks: {len(gen_x)}")
print(f"Matched tracks: {len(matched)}")
print(f"Fake tracks: {len(fake)}")

num_pixels_x = 18 * 1024
num_pixels_y = 512
bins_x = 18 * 2 * 4
bins_y = 4

bin_size_x = int(num_pixels_x / bins_x)
bin_size_y = int(num_pixels_y / bins_y)

gen_xplets = ROOT.TH2F('gen',
                       'x vs y', bins_x, 0.05275912, 0.55430088, bins_y, -0.00688128, 0.00688128)

matched_xplets = ROOT.TH2F('matched',
                           'x vs y', bins_x, 0.05275912, 0.55430088, bins_y, -0.00688128, 0.00688128)

fake_xplets = ROOT.TH2F('fake_xplets', 
                        'x vs y', bins_x, 0.05275912, 0.55430088, bins_y, -0.00688128, 0.00688128)

reco_xplets = ROOT.TH2F('reco_xplets', 
                        'x vs y', bins_x, 0.05275912, 0.55430088, bins_y, -0.00688128, 0.00688128)

gen_xplets_energy = ROOT.TH2F('gen energy',
                              'gen energy', bins_x, 0.05275912, 0.55430088, bins_y, -0.00688128, 0.00688128)

gen_xplets_energy2 = ROOT.TH2F('gen energy2',
                               'gen energy2', bins_x, 0.05275912, 0.55430088, bins_y, -0.00688128, 0.00688128)

matched_xplets_energy = ROOT.TH2F('matched energy',
                                  'matched energy', bins_x, 0.05275912, 0.55430088, bins_y, -0.00688128, 0.00688128)

matched_xplets_energy2 = ROOT.TH2F('matched energy2',
                                   'matched energy2', bins_x, 0.05275912, 0.55430088, bins_y, -0.00688128, 0.00688128)


e_lim_1 = 2000
e_lim_2 = 2500

counter = 0
counter_energy_1 = 0
counter_energy_2 = 0
for g_x in gen_x:
    gen_xplets.Fill(g_x.coordinates[0][0],
                    g_x.coordinates[0][1])
    counter += 1
    if g_x.energy[0] > e_lim_1:
        gen_xplets_energy.Fill(g_x.coordinates[0][0],
                               g_x.coordinates[0][1])
        counter_energy_1 += 1
    if g_x.energy[0] > e_lim_2:
        gen_xplets_energy2.Fill(g_x.coordinates[0][0],
                                g_x.coordinates[0][1])
        counter_energy_2 += 1

print(counter)
print(counter_energy_1)
print(counter_energy_2)
exit()

        
for f_x in fake:
    fake_xplets.Fill(f_x.coordinates[0][0],
                     f_x.coordinates[0][1])
    
for r_x in reco_x:
    reco_xplets.Fill(r_x.coordinates[0][0],
                     r_x.coordinates[0][1])

for m_x in matched:
    matched_xplets.Fill(m_x.coordinates[0][0],
                        m_x.coordinates[0][1])
    if m_x.energy[0] > e_lim_1:
        matched_xplets_energy.Fill(m_x.coordinates[0][0],
                                   m_x.coordinates[0][1])
    if m_x.energy[0] > e_lim_2:
        matched_xplets_energy2.Fill(m_x.coordinates[0][0],
                                    m_x.coordinates[0][1])


def get_standard_deviation(k, n):
    if k == 0:
        return 0
    if k > n:
        k = n
    return np.sqrt(((k + 1) * (k + 2)) / ((n + 2) * (n + 3)) - ((k + 1)**2 / (n + 2)**2))


eff_standard_deviation_e_lower = []
eff_standard_deviation_e_upper = []

eff_standard_deviation_e_lower_energy = []
eff_standard_deviation_e_upper_energy = []

eff_standard_deviation_e_lower_energy2 = []
eff_standard_deviation_e_upper_energy2 = []


matched_xplets.Divide(gen_xplets)
matched_xplets_energy.Divide(gen_xplets_energy)
matched_xplets_energy2.Divide(gen_xplets_energy2)


occupancy = []
efficiency = []
for i in range(matched_xplets.GetNbinsX() + 1):
    for j in range(matched_xplets.GetNbinsY() + 1):
        if gen_xplets.GetBinContent(i, j) > 1:
            occupancy.append(gen_xplets.GetBinContent(i, j))
            efficiency.append(matched_xplets.GetBinContent(i, j))
            eff_standard_deviation_e_lower.append(
                get_standard_deviation(matched_xplets.GetBinContent(i, j),
                                       gen_xplets.GetBinContent(i, j))
                if matched_xplets.GetBinContent(i, j) / gen_xplets.GetBinContent(i, j) -
                get_standard_deviation(matched_xplets.GetBinContent(i, j),
                                       gen_xplets.GetBinContent(i, j)) >= 0
                else matched_xplets.GetBinContent(i, j) / gen_xplets.GetBinContent(i, j))
            eff_standard_deviation_e_upper.append(
                get_standard_deviation(matched_xplets.GetBinContent(i, j),
                                       gen_xplets.GetBinContent(i, j))
                if min(matched_xplets.GetBinContent(i, j) / gen_xplets.GetBinContent(i, j), 1) +
                get_standard_deviation(matched_xplets.GetBinContent(i, j),
                                       gen_xplets.GetBinContent(i, j)) <= 1
                else 1 - min(matched_xplets.GetBinContent(i, j) / gen_xplets.GetBinContent(i, j), 1))

occupancy_energy = []
efficiency_energy = []
for i in range(matched_xplets_energy.GetNbinsX() + 1):
    for j in range(matched_xplets_energy.GetNbinsY() + 1):
        if gen_xplets_energy.GetBinContent(i, j) > 1:
            occupancy_energy.append(gen_xplets_energy.GetBinContent(i, j))
            efficiency_energy.append(matched_xplets_energy.GetBinContent(i, j))
            eff_standard_deviation_e_lower_energy.append(
                get_standard_deviation(matched_xplets_energy.GetBinContent(i, j),
                                       gen_xplets_energy.GetBinContent(i, j))
                if matched_xplets_energy.GetBinContent(i, j) / gen_xplets_energy.GetBinContent(i, j) -
                get_standard_deviation(matched_xplets_energy.GetBinContent(i, j),
                                       gen_xplets_energy.GetBinContent(i, j)) >= 0
                else matched_xplets_energy.GetBinContent(i, j) / gen_xplets_energy.GetBinContent(i, j))
            eff_standard_deviation_e_upper_energy.append(
                get_standard_deviation(matched_xplets_energy.GetBinContent(i, j),
                                       gen_xplets_energy.GetBinContent(i, j))
                if min(matched_xplets_energy.GetBinContent(i, j) / gen_xplets_energy.GetBinContent(i, j), 1) +
                get_standard_deviation(matched_xplets_energy.GetBinContent(i, j),
                                       gen_xplets_energy.GetBinContent(i, j)) <= 1
                else 1 - min(matched_xplets_energy.GetBinContent(i, j) / gen_xplets_energy.GetBinContent(i, j), 1))

occupancy_energy2 = []
efficiency_energy2 = []
for i in range(matched_xplets_energy2.GetNbinsX() + 1):
    for j in range(matched_xplets_energy2.GetNbinsY() + 1):
        if gen_xplets_energy2.GetBinContent(i, j) > 1:
            occupancy_energy2.append(gen_xplets_energy2.GetBinContent(i, j))
            efficiency_energy2.append(matched_xplets_energy2.GetBinContent(i, j))
            eff_standard_deviation_e_lower_energy2.append(
                get_standard_deviation(matched_xplets_energy2.GetBinContent(i, j),
                                       gen_xplets_energy2.GetBinContent(i, j))
                if matched_xplets_energy2.GetBinContent(i, j) / gen_xplets_energy2.GetBinContent(i, j) -
                get_standard_deviation(matched_xplets_energy2.GetBinContent(i, j),
                                       gen_xplets_energy2.GetBinContent(i, j)) >= 0
                else matched_xplets_energy2.GetBinContent(i, j) / gen_xplets_energy2.GetBinContent(i, j))
            eff_standard_deviation_e_upper_energy2.append(
                get_standard_deviation(matched_xplets_energy2.GetBinContent(i, j),
                                       gen_xplets_energy2.GetBinContent(i, j))
                if min(matched_xplets_energy2.GetBinContent(i, j) / gen_xplets_energy2.GetBinContent(i, j), 1) +
                get_standard_deviation(matched_xplets_energy2.GetBinContent(i, j),
                                       gen_xplets_energy2.GetBinContent(i, j)) <= 1
                else 1 - min(matched_xplets_energy2.GetBinContent(i, j) / gen_xplets_energy2.GetBinContent(i, j), 1))


canv = TCanvas("c1", "eff vs. occ plot", 800, 600)

gr0 = TGraphMultiErrors("e1",
                        "e1",
                        len(occupancy),
                        np.array(occupancy),
                        np.array(efficiency),
                        np.zeros(len(occupancy)),
                        np.zeros(len(occupancy)),
                        np.array(eff_standard_deviation_e_lower),
                        np.array(eff_standard_deviation_e_upper))
gr0.SetMarkerStyle(8)
gr0.SetMarkerSize(1.2)
gr0.SetMarkerColor(1)
gr0.SetLineColor(2)
gr0.SetLineWidth(0)
gr0.GetAttLine(0).SetLineColor(kBlack)
gr0.GetAttLine(0).SetLineWidth(1)
gr0.SetMaximum(1.05)
gr0.SetMinimum(0.0)
gr0.SetTitle(f"#xi = {xi}, label; occupancy / {bin_size_x} x {bin_size_y} pixel; efficiency")
gr0.Draw("APS; Z; s=0.5")


gr1 = TGraph(len(occupancy_energy),
             np.array(occupancy_energy),
             np.array(efficiency_energy))

gr1.SetMarkerStyle(8)
gr1.SetMarkerSize(1.2)
gr1.SetLineColor(2)
gr1.SetMarkerColor(2)
gr1.SetLineWidth(0)
gr1.Draw("PSAME2")

gr2 = TGraph(len(occupancy_energy2),
             np.array(occupancy_energy2),
             np.array(efficiency_energy2))

gr2.SetMarkerStyle(8)
gr2.SetLineColor(2)
gr2.SetMarkerSize(1.2)
gr2.SetMarkerColor(4)
gr2.SetLineWidth(0)
gr2.Draw("PSAME2")


leg=TLegend(0.3, 0.2, 0.6, 0.4)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetTextSize(0.05)
leg.AddEntry(gr0, f"#xi = {xi}", "P")
leg.AddEntry(gr1, f"#xi = {xi}, > {np.around(e_lim_1 / 1000, 2)} GeV", "P")
leg.AddEntry(gr2, f"#xi = {xi}, > {np.around(e_lim_2 / 1000, 2)} GeV", "P")

leg.Draw()

LUXELabel(0.8, 0.25)

canv.SaveAs(f"{folder}/eff_vs_occupancy_{xi}.pdf")
