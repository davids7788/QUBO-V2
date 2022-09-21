import ROOT
from ROOT import *
import numpy as np
import sys
sys.path.append("../src")

ROOT.gROOT.LoadMacro("macros/LuxeStyle.C")
ROOT.gROOT.LoadMacro("macros/LuxeLabels.C")

SetLuxeStyle()

preselection_folder = sys.argv[1]
gen_x = np.load(preselection_folder.split("-")[0] + "_gen_xplet_list.npy", allow_pickle=True)
triplet_list = np.load(preselection_folder + "/triplet_list.npy", allow_pickle=True)
xi = preselection_folder.split("e0gpc_")[1].split("_")[0]

# only full gen particle tracks are considered --> all layers have to be hit
gen_particles = set()

edges_left = [1.5]
edges_mid = [2.0 + 0.4 * i for i in range(15)]
edges_right = [8.0, 8.5, 9.5]
edges = np.array(edges_left + edges_mid + edges_right)

gen_doublets = TH1F('gen doublets', 'energy', len(edges) - 1, edges)
gen_triplets = TH1F('gen triplets', 'energy', len(edges) - 1, edges)

selected_doublets = TH1F('selected doublets', 'energy', len(edges) - 1, edges)
matched_doublets = TH1F('matched doublets', 'energy', len(edges) - 1, edges)

selected_triplets = TH1F('selected triplets', 'energy', len(edges) - 1, edges)
matched_triplets = TH1F('matched triplets', 'energy', len(edges) - 1, edges)

for xplet in gen_x:
    gen_doublets.Fill(xplet.energy[0] / 1000)
    gen_doublets.Fill(xplet.energy[1] / 1000)
    gen_doublets.Fill(xplet.energy[2] / 1000)
    gen_triplets.Fill(xplet.energy[0] / 1000)
    gen_triplets.Fill(xplet.energy[1] / 1000)

counter = 0
matched_doublet_set = set()
doublet_set = set()

for triplet in triplet_list:
    selected_triplets.Fill(triplet.doublet_1.energy_1 / 1000)
    if triplet.doublet_1 not in doublet_set: 
        selected_doublets.Fill(triplet.doublet_1.energy_1 / 1000)
        doublet_set.add(triplet.doublet_1)
    if triplet.doublet_2 not in doublet_set: 
        selected_doublets.Fill(triplet.doublet_2.energy_1 / 1000)
        doublet_set.add(triplet.doublet_2)

    if triplet.doublet_1.hit_1_particle_key == triplet.doublet_1.hit_2_particle_key == \
                triplet.doublet_2.hit_1_particle_key == triplet.doublet_2.hit_2_particle_key:
        if triplet.doublet_1 not in matched_doublet_set: 
            matched_doublets.Fill(triplet.doublet_1.energy_1 / 1000)
            matched_doublet_set.add(triplet.doublet_1)
        if triplet.doublet_2 not in matched_doublet_set: 
            matched_doublets.Fill(triplet.doublet_2.energy_1 / 1000)
            matched_doublet_set.add(triplet.doublet_2)
            
        matched_triplets.Fill(triplet.doublet_1.energy_1 / 1000)
        

for i in range(len(edges)):
    if matched_doublets.GetBinContent(i) > gen_doublets.GetBinContent(i):
        matched_doublets.SetBinContent(i, gen_doublets.GetBinContent(i))
    if matched_triplets.GetBinContent(i) > gen_triplets.GetBinContent(i):
        matched_triplets.SetBinContent(i, gen_triplets.GetBinContent(i))
        
        
gStyle.SetPadTickY(0)
canv = TCanvas("example", "doublet efficiency ", 800, 600)
pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
pad1.SetBottomMargin(0.02)
pad1.Draw()
pad1.cd()


h_frame = TH1F('frame', 'energy', len(edges) - 1, edges)
h_frame.GetYaxis().SetTitle("fraction of counts")
h_frame.GetXaxis().SetTitle("energy [GeV]")
h_frame.GetYaxis().SetTitleSize(0.06)
h_frame.GetXaxis().SetNdivisions(105)
h_frame.SetMaximum(1.05)
h_frame.GetYaxis().SetTitleOffset(0.0)
h_frame.GetXaxis().SetLabelOffset(0.02)
h_frame.GetXaxis().SetRangeUser(0.0, 11.0)
h_frame.Draw()

# Teff
eff_doublets = TEfficiency(matched_doublets, gen_doublets)
eff_doublets.SetMarkerStyle(8)
eff_doublets.SetMarkerColor(kBlack)
eff_doublets.SetMarkerSize(0.95)
eff_doublets.Draw("PSAME")

pad1.SetRightMargin(0.13)
pad1.Draw()

rightmax = selected_doublets.GetMaximum() * 1.1
scale = pad1.GetUymax() / rightmax

selected_doublets.SetLineColorAlpha(kRed, 0.5)
selected_doublets.SetLineWidth(2)
selected_doublets.Scale(scale)
selected_doublets.Draw("HISTSAME")

matched_doublets.SetLineColorAlpha(kBlue, 0.5)
matched_doublets.SetLineWidth(2)
matched_doublets.Scale(scale)
matched_doublets.Draw("HISTSAME")

axis = TGaxis(gPad.GetUxmax(),
              gPad.GetUymin(),
              gPad.GetUxmax(),
              gPad.GetUymax(),
              0,
              rightmax,
              510,
              "+L")
axis.SetTitle("counts")
axis.SetTitleSize(0.06)
axis.SetLabelFont(42)
axis.SetTextFont(42)
axis.SetLabelSize(0.06)
axis.SetTitleOffset(1.1)
axis.Draw()
canv.Update()
canv.cd()



LUXELabel(0.7, 0.85)

leg = TLegend(0.55, 0.55, 0.75, 0.75)
leg.AddEntry(eff_doublets, "doublet efficiency", "p")
leg.AddEntry(selected_doublets, "selected doublets", "l")
leg.AddEntry(matched_doublets, "matched doublets", "l")
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetTextSize(0.04)
leg.SetHeader(f"#xi = {xi}")
leg.Draw()


pad2 = TPad("pad2", "pad2", 0.0, 0.05, 1., 0.3)
pad2.SetTopMargin(0.)
pad2.SetBottomMargin(0.35)
pad2.SetRightMargin(0.13)
pad2.SetGridy()
pad2.Draw()
pad2.cd()
pad2.Update()


matched_clone = matched_doublets.Clone()
matched_clone.Divide(selected_doublets)
matched_clone.SetMarkerStyle(8)
matched_clone.SetLineColor(kBlack)
matched_clone.SetMarkerSize(0.95)
matched_clone.GetYaxis().SetTitle("#frac{matched}{selected}")
matched_clone.GetYaxis().SetTitleSize(0.16)
matched_clone.GetYaxis().SetLabelSize(0.16)
matched_clone.GetYaxis().SetTitleOffset(0.4)
matched_clone.SetMinimum(0.0)
matched_clone.SetMaximum(1.1)
matched_clone.GetXaxis().SetTitle("energy [GeV]")
matched_clone.GetXaxis().SetTitleSize(0.16)
matched_clone.GetXaxis().SetLabelSize(0.16)
matched_clone.GetXaxis().SetTitleOffset(1)
matched_clone.GetYaxis().SetNdivisions(105)
matched_clone.Draw("HIST")

canv.SaveAs(f"{preselection_folder}/preselection_efficiency_doublets_{xi}.pdf")


##############################################
## Triplet plots
canv2 = TCanvas("example2", "triplets efficiency ", 800, 600)
pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
pad1.SetBottomMargin(0.02)
pad1.Draw()
pad1.cd()


h_frame = TH1F('frame', 'energy', len(edges) - 1, edges)
h_frame.GetYaxis().SetTitle("fraction of counts")
h_frame.GetXaxis().SetTitle("energy [GeV]")
h_frame.GetYaxis().SetTitleSize(0.06)
h_frame.GetXaxis().SetNdivisions(105)
h_frame.SetMaximum(1.05)
h_frame.GetYaxis().SetTitleOffset(0.0)
h_frame.GetXaxis().SetLabelOffset(0.02)
h_frame.GetXaxis().SetRangeUser(0.0, 11.0)
h_frame.Draw()

# Teff
eff_triplets = TEfficiency(matched_triplets, gen_triplets)
eff_triplets.SetMarkerStyle(8)
eff_triplets.SetMarkerColor(kBlack)
eff_triplets.SetMarkerSize(0.95)
eff_triplets.Draw("PSAME")

pad1.SetRightMargin(0.13)
pad1.Draw()

rightmax = selected_triplets.GetMaximum() * 1.1
scale = pad1.GetUymax() / rightmax

selected_triplets.SetLineColorAlpha(kRed, 0.5)
selected_triplets.SetLineWidth(2)
selected_triplets.Scale(scale)
selected_triplets.Draw("HISTSAME")

matched_triplets.SetLineColorAlpha(kBlue, 0.5)
matched_triplets.SetLineWidth(2)
matched_triplets.Scale(scale)
matched_triplets.Draw("HISTSAME")

axis = TGaxis(gPad.GetUxmax(),
              gPad.GetUymin(),
              gPad.GetUxmax(),
              gPad.GetUymax(),
              0,
              rightmax,
              510,
              "+L")
axis.SetTitle("counts")
axis.SetTitleSize(0.06)
axis.SetLabelFont(42)
axis.SetTextFont(42)
axis.SetLabelSize(0.06)
axis.SetTitleOffset(1.1)
axis.Draw()
canv2.Update()
canv2.cd()



LUXELabel(0.7, 0.85)

leg = TLegend(0.55, 0.55, 0.75, 0.75)
leg.AddEntry(eff_triplets, "triplets efficiency", "p")
leg.AddEntry(selected_triplets, "selected triplets", "l")
leg.AddEntry(matched_triplets, "matched triplets", "l")
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetTextSize(0.04)
leg.SetHeader(f"#xi = {xi}")
leg.Draw()


pad2 = TPad("pad2", "pad2", 0.0, 0.05, 1., 0.3)
pad2.SetTopMargin(0.)
pad2.SetBottomMargin(0.35)
pad2.SetRightMargin(0.13)
pad2.SetGridy()
pad2.Draw()
pad2.cd()
pad2.Update()


matched_clone = matched_triplets.Clone()
matched_clone.Divide(selected_doublets)
matched_clone.SetMarkerStyle(8)
matched_clone.SetLineColor(kBlack)
matched_clone.SetMarkerSize(0.95)
matched_clone.GetYaxis().SetTitle("#frac{matched}{selected}")
matched_clone.GetYaxis().SetTitleSize(0.16)
matched_clone.GetYaxis().SetLabelSize(0.16)
matched_clone.GetYaxis().SetTitleOffset(0.4)
matched_clone.SetMinimum(0.0)
matched_clone.SetMaximum(1.1)
matched_clone.GetXaxis().SetTitle("energy [GeV]")
matched_clone.GetXaxis().SetTitleSize(0.16)
matched_clone.GetXaxis().SetLabelSize(0.16)
matched_clone.GetXaxis().SetTitleOffset(1)
matched_clone.GetYaxis().SetNdivisions(105)
matched_clone.Draw("HIST")

canv2.SaveAs(f"{preselection_folder}/preselection_efficiency_triplets_{xi}.pdf")