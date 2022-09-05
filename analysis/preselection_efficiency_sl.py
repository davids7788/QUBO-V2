import ROOT
from ROOT import *
import numpy as np
import sys
sys.path.append("../src")

ROOT.gROOT.LoadMacro("macros/LuxeStyle.C")
ROOT.gROOT.LoadMacro("macros/LuxeLabels.C")

SetLuxeStyle()

preselection_folder = sys.argv[1]
gen_x = preselection_folder.split("-")[0] + "_gen_xplet_list.npy"
triplet_list = preselection_folder + "triplet_list.npy"
xi = preselection_folder.split("e0gpc_")[1].split("_")[0]

# only full gen particle tracks are considered --> all layers have to be hit
gen_particles = set()

for xplet in gen_x:
    for triplet_id in xplet.triplet_ids:
        gen_particles.add(triplet_list[triplet_id].doublet_1.hit_1_particle_key)
        gen_particles.add(triplet_list[triplet_id].doublet_1.hit_2_particle_key)
        gen_particles.add(triplet_list[triplet_id].doublet_2.hit_1_particle_key)
        gen_particles.add(triplet_list[triplet_id].doublet_2.hit_2_particle_key)


t_list_skimmed = []
for triplet in t_list_skimmed:
    if triplet.doublet_1.hit_1_particle_key not in gen_particles_triplets:
        continue
    if triplet.doublet_1.hit_2_particle_key not in gen_particles_triplets:
        continue
    if triplet.doublet_2.hit_1_particle_key not in gen_particles_triplets:
        continue
    if triplet.doublet_2.hit_2_particle_key not in gen_particles_triplets:
        continue
    else:
        t_list_skimmed.append(triplet)


edges_left = [0.05, 0.07, 0.085]
edges_mid = [0.1 + 0.01 * i for i in range(24)]
edges_right = [0.35, 0.38, 0.42]
edges = np.array(edges_left + edges_mid + edges_right)

gen_doublets = TH1F('gen doublets', 'energy', len(edges) - 1, edges)
gen_triplets = TH1F('gen triplets', 'energy', len(edges) - 1, edges)

selected_doublets = TH1F('gen doublets', 'energy', len(edges) - 1, edges)
selected_triplets = TH1F('gen triplets', 'energy', len(edges) - 1, edges)

for xplet in gen_xplets:
    gen_doublets.Fill(xplet.energy[0])
    gen_doublets.Fill(xplet.energy[1])
    gen_doublets.Fill(xplet.energy[2])
    gen_triplets.Fill(xplet.energy[0])
    gen_triplets.Fill(xplet.energy[1])

used_doublets = set()
for triplet in triplet_list:
    if triplet.doublet_1.hit_1_particle_key not in gen_particles_triplets:
        continue
    if triplet.doublet_1.hit_2_particle_key not in gen_particles_triplets:
        continue
    if triplet.doublet_2.hit_1_particle_key not in gen_particles_triplets:
        continue
    if triplet.doublet_2.hit_2_particle_key not in gen_particles_triplets:
        continue
    else:
        if triplet.doublet_1.hit_1_particle_key == triplet.doublet_1.hit_2_particle_key == \
                triplet.doublet_2.hit_1_particle_key == triplet.doublet_2.hit_2_particle_key:
            if triplet.doublet_1 in used_doublets:
                pass
            else:
                selected_doublets.Fill(triplet.doublet_1.energy)
            if triplet.doublet_2 in used_doublets:
                pass
            else:
                selected_doublets.Fill(triplet.doublet_2.energy)
            selected_triplets.Fill(triplet.doublet_1.energy_1)


gStyle.SetPadTickY(0)
canv = TCanvas("example", "doublet efficiency ", 800, 600)

h_frame = TH1F('frame', 'energy', len(edges) - 1, edges)
h_frame.GetYaxis().SetTitle("fraction of counts")
h_frame.GetXaxis().SetTitle("energy [GeV]")
h_frame.GetXaxis().SetNdivisions(105)
h_frame.SetMaximum(1.18)
h_frame.GetXaxis().SetLabelOffset(0.02)
h_frame.GetXaxis().SetRangeUser(0.05, 0.51)
h_frame.Draw()

eff_doublets = TEfficiency(selected_doublets, gen_doublets)

eff_doublets.SetMarkerStyle(8)
eff_doublets.SetMarkerColorAlpha(kBlack, 0.85)
eff_doublets.SetMarkerSize(0.95)
eff_doublets.Draw("PSAME")

LUXELabel(0.75, 0.85)

leg = TLegend(0.55, 0.7, 0.75, 0.9)
leg.AddEntry(eff_doublets, "doublet efficiency", "p")
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetTextSize(0.04)
leg.SetHeader(f"#xi = {xi}")
leg.Draw()
canv.SaveAs(f"{folder}/preselection_efficiency_doublets_{xi}.pdf")

gStyle.SetPadTickY(0)
canv2 = TCanvas("example", "doublet efficiency ", 800, 600)

h_frame_2 = TH1F('frame', 'energy', len(edges) - 1, edges)
h_frame_2.GetYaxis().SetTitle("fraction of counts")
h_frame_2.GetXaxis().SetTitle("energy [GeV]")
h_frame_2.GetXaxis().SetNdivisions(105)
h_frame_2.SetMaximum(1.18)
h_frame_2.GetXaxis().SetLabelOffset(0.02)
h_frame_2.GetXaxis().SetRangeUser(0.05, 0.51)
h_frame_2.Draw()

eff_triplets = TEfficiency(selected_triplets, gen_triplets)
eff_triplets.SetMarkerStyle(8)
eff_triplets.SetMarkerColorAlpha(kBlack, 0.85)
eff_triplets.SetMarkerSize(0.95)
eff_triplets.Draw("PSAME")

LUXELabel(0.75, 0.85)

leg = TLegend(0.55, 0.7, 0.75, 0.9)
leg.AddEntry(eff_doublets, "triplet efficiency", "p")
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetTextSize(0.04)
leg.SetHeader(f"#xi = {xi}")
leg.Draw()
canv2.SaveAs(f"{folder}/preselection_efficiency_triplets_{xi}.pdf")

