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

gen_x = np.load(gen + "_gen_xplet_list.npy", allow_pickle=True)
reco_x = np.load(f"{folder}/reco_xplet_list.npy", allow_pickle=True)

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

gen_xplets = ROOT.TH2F('gen_xplets', 
                       'x vs y', 100, 0.05275912, 0.55430088, 100, -0.00688128, 0.00688128)

matched_xplets = ROOT.TH2F('gen_xplets', 
                           'x vs y', 100, 0.05275912, 0.55430088, 100, -0.00688128, 0.00688128)

fake_xplets = ROOT.TH2F('fake_xplets', 
                        'x vs y', 100, 0.05275912, 0.55430088, 100, -0.00688128, 0.00688128)

reco_xplets = ROOT.TH2F('reco_xplets', 
                        'x vs y', 100, 0.05275912, 0.55430088, 100, -0.00688128, 0.00688128)



for g_x in gen_x:
    if region == -1:
        if g_x.coordinates[0][1] > 0:
            continue
    elif region == 1:
        if g_x.coordinates[0][1] < 0:
            continue  
    gen_xplets.Fill(g_x.coordinates[0][0],
                    g_x.coordinates[0][1])


for m_x in matched:
    if region == -1:
        if m_x.coordinates[0][1] > 0:
            continue
    if region == 1:
        if m_x.coordinates[0][1] < 0:
            continue
    matched_xplets.Fill(m_x.coordinates[0][0],
                        m_x.coordinates[0][1])


for f_x in fake:
    if region == -1:
        if f_x.coordinates[0][1] > 0:
            continue
    if region == 1:
        if f_x.coordinates[0][1] < 0:
            continue
    fake_xplets.Fill(f_x.coordinates[0][0],
                     f_x.coordinates[0][1])
    

for r_x in reco_x:
    if region == -1:
        if r_x.coordinates[0][1] > 0:
            continue
    if region == 1:
        if r_x.coordinates[0][1] < 0:
            continue
    reco_xplets.Fill(r_x.coordinates[0][0],
                     r_x.coordinates[0][1])


canv = TCanvas("example1", "xplet efficiency ", 600, 400)
matched_xplets.Divide(gen_xplets)

matched_xplets.GetYaxis().SetTitle("y [m]")
matched_xplets.GetXaxis().SetTitle("x [m]")

matched_xplets.Draw("COLZ")
LUXELabel(0.2, 0.85);
gPad.SetRightMargin(0.1)
gPad.Draw()

leg=TLegend(0.2,0.14,0.4,0.4)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetTextSize(0.04)
leg.SetHeader("#splitline{#xi = " + xi + "}{xplet efficiency}")
leg.Draw()

canv.Draw()
canv.SaveAs(f"{folder}/2Dxplet_efficiency_xi_{xi}.pdf")
try:
    canv.Close(); 
    gSystem.ProcessEvents()
except NameError:
    pass

#
canv2 = TCanvas("example2","fake x-plets ",600,400)

fake_xplets.GetYaxis().SetTitle("y [m]")
fake_xplets.GetXaxis().SetTitle("x [m]")

fake_xplets.Draw("COLZ")
LUXELabel(0.2,0.85);
gPad.SetRightMargin(0.1)
gPad.Draw()

leg=TLegend(0.2,0.14,0.4,0.4)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetTextSize(0.04)
leg.SetHeader("#splitline{#xi = " + xi + "}{fake xplets}")
leg.Draw()

canv2.Update()
canv2.Draw()
canv2.SaveAs(f"{folder}/2Dfake_xplets_xi_{xi}{save_extension}.pdf")

try:
    canv2.Close(); 
    gSystem.ProcessEvents()
except NameError:
    pass
##

fake_xplets.Divide(reco_xplets)
canv3 = TCanvas("example3","xplet fake rate ",600,400)

fake_xplets.GetYaxis().SetTitle("y [m]")
fake_xplets.GetXaxis().SetTitle("x [m]")

fake_xplets.Draw("COLZ")
LUXELabel(0.2,0.85);
gPad.SetRightMargin(0.1)
gPad.Draw()

leg=TLegend(0.2,0.14,0.4,0.4)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetTextSize(0.04)
leg.SetHeader("#splitline{#xi = " + xi + "}{xplets fake rate}")
leg.Draw()

canv3.Draw()
canv3.SaveAs(f"{folder}/2Dxplet_fake_rate_xi_{xi}{save_extension}.pdf")