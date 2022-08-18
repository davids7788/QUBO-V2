import matplotlib.pyplot as plt
import ROOT
from ROOT import *
import numpy as np
import sys
sys.path.append("../src")

ROOT.gROOT.LoadMacro("macros/LuxeStyle.C")
ROOT.gROOT.LoadMacro("macros/LuxeLabels.C")

SetLuxeStyle()

folder = sys.argv[1]
gen = "-".join("/".join(folder.split("/")[0:-1]).split("-")[0:-1])
xi = folder.split("e0gpc_")[1].split("_")[0]

region = int(sys.argv[2])
#  0: all xplets
#  1: y > 0
# -1: y < 0

if region == 0:
    save_extension = "_all"
elif region == 1:
    save_extension = "_upper"
elif region == -1:
    save_extension = "_lower"
else:
    print("No valid region selected!")
    exit()

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

gen_xplets = TH1F('gen xplets', 
                  'x', 50, 0.05275912, 0.55430088)

matched_xplets = TH1F('matched xplets', 
                      'x', 50, 0.05275912, 0.55430088)

fake_xplets = TH1F('fake xplets', 
                   'x', 50, 0.05275912, 0.55430088)

reco_xplets = TH1F('reco xplets', 
                   'x', 50, 0.05275912, 0.55430088)

h_frame = TH1F('frame', 
               'x', 50, 0.05275912, 0.55430088)



for g_x in gen_x:
    if region == -1:
        if g_x.coordinates[0][1] > 0:
            continue
    elif region == 1:
        if g_x.coordinates[0][1] < 0:
            continue  
    gen_xplets.Fill(g_x.coordinates[0][0])


for m_x in matched:
    if region == -1:
        if m_x.coordinates[0][1] > 0:
            continue
    if region == 1:
        if m_x.coordinates[0][1] < 0:
            continue
    matched_xplets.Fill(m_x.coordinates[0][0])


for f_x in fake:
    if region == -1:
        if f_x.coordinates[0][1] > 0:
            continue
    if region == 1:
        if f_x.coordinates[0][1] < 0:
            continue
    fake_xplets.Fill(f_x.coordinates[0][0])
    

for r_x in reco_x:
    if region == -1:
        if r_x.coordinates[0][1] > 0:
            continue
    if region == 1:
        if r_x.coordinates[0][1] < 0:
            continue
    reco_xplets.Fill(r_x.coordinates[0][0])
    
Teff_error = TEfficiency(matched_xplets, gen_xplets)
Tfake_error = TEfficiency(fake_xplets, reco_xplets)


gStyle.SetPadTickY(0)   
canv = TCanvas("example","xplet efficiency ", 800, 600)

h_frame.GetYaxis().SetTitle("fraction of counts")
h_frame.GetXaxis().SetTitle("x [m]")
h_frame.Draw()

Teff_error.SetMarkerStyle(8)
Teff_error.SetMarkerColorAlpha(kBlue, 0.85)
Teff_error.SetMarkerSize(0.95)
Teff_error.Draw("PSAME")

Tfake_error.SetMarkerStyle(8)
Tfake_error.SetMarkerColorAlpha(kRed, 0.85)
Tfake_error.SetMarkerSize(0.95)
Tfake_error.Draw("PSAME")

gPad.SetRightMargin(0.13)
gPad.Draw()

rightmax = gen_xplets.GetMaximum()
scale = gPad.GetUymax() / rightmax

gen_xplets.SetLineColorAlpha(33, 0.8)
gen_xplets.SetLineWidth(3)
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
axis.SetTextSize(0.05)
axis.SetTitleOffset(1.4)
axis.Draw()


LUXELabel(0.75, 0.85);

leg = TLegend(0.55, 0.7, 0.75, 0.9)
leg.AddEntry(Teff_error, "efficiency", "p")
leg.AddEntry(Tfake_error, "fake rate", "p")
leg.AddEntry(gen_xplets, "gen xplets", "l")
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetTextSize(0.04)
leg.SetHeader(f"#xi = {xi}")
leg.Draw()

canv.SaveAs(f"{folder}/1Dxplet_efficiency_xi_{xi}{save_extension}.pdf")