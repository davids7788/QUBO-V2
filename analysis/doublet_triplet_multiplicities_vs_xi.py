import matplotlib.pyplot as plt
import ROOT
from ROOT import *
import numpy as np
import sys
from array import array

ROOT.gROOT.LoadMacro("macros/LuxeStyle.C")
ROOT.gROOT.LoadMacro("macros/LuxeLabels.C")

SetLuxeStyle()

p_4_0_array = [2124, 2047, 2065, 2022, 2099, 2064, 2060, 2113, 1992, 2139]
d_4_0_array = [8123, 7651, 7865, 7642, 7853, 7770, 7830, 8007, 7402, 8207]
t_4_0_array = [5279, 4945, 5206, 5058, 5210, 5124, 5131, 5294, 4870, 5376]

p_5_0_array = [10849, 10576, 10168, 10600, 10616, 10579, 10831, 10587, 10639, 10358]
d_5_0_array = [66083, 63032, 59393, 62627, 63692, 62656, 65532, 63041, 63945, 61531]
t_5_0_array = [47675, 44028, 41679, 43982, 45097, 44505, 46124, 43740, 44996, 42688]

p_7_0_array = [67052, 67429, 65927, 67310, 67520, 67759, 66383, 68196, 68341, 67221]
d_7_0_array = [1069207, 1077725, 1035783, 1072821, 1077949, 1088879, 1046797, 1093273, 1102585, 1070506]
t_7_0_array = [1202492, 1228155, 1162599, 1205280, 1215384, 1255901, 1171173, 1237625, 1264322, 1202798]


x = [4, 5, 7]
p = [int(np.mean(p_4_0_array)), int(np.mean(p_5_0_array)), int(np.mean(p_7_0_array))]
d = [int(np.mean(d_4_0_array)), int(np.mean(d_5_0_array)), int(np.mean(d_7_0_array))]
t = [int(np.mean(t_4_0_array)), int(np.mean(t_5_0_array)), int(np.mean(t_7_0_array))]

num_doublets = array("d")
num_triplets = array("d")
num_particles = array("d")
xi = array("d")

for i in range(3):
    num_doublets.append(d[i])
    num_triplets.append(t[i])
    num_particles.append(p[i])
    xi.append(x[i])

canvas = TCanvas( 'c1', 'd/t multiplicities', 800, 600)
gPad.SetTopMargin(0.15)
gPad.SetLogy()
gr = TGraph(3, xi, num_doublets)
gr.SetLineColor(kBlack)
gr.GetXaxis().SetNdivisions(305)
gr.SetMarkerColor(kBlack)
gr.SetMarkerSize(2.0)
gr.SetMinimum(1e3)
gr.SetMaximum(5e6)
gr.SetMarkerStyle( 21 )
gr.GetXaxis().SetTitle('#xi')
gr.GetYaxis().SetTitle('counts')
# gr.Fit("expo")
# ff = gr.GetFunction("expo")
# ff.SetLineColor(38)
# ff.SetLineStyle(2)
gr.Draw('ALP')

canvas.Update()

gr2 = TGraph(3, xi, num_triplets)
gr2.SetLineColor(kRed)
gr2.SetMarkerColor(kRed)
gr2.SetMarkerStyle(22)
gr2.SetMarkerSize(2.0)
gr2.GetXaxis().SetTitle('#xi')
gr2.GetYaxis().SetTitle('counts')
# gr2.Fit("expo")
# ff2 = gr2.GetFunction("expo")
# ff2.SetLineColor(46)
# ff2.SetLineStyle(2)
gr2.Draw('LPSAME')

canvas.Update()

LUXELabel(0.2, 0.7, "e-laser, phase-0")

latex_l1 = TLatex()
latex_l1.SetTextSize(0.05)
latex_l1.SetTextFont(42)
latex_l1.DrawLatex(3.825, 6e6, str(p[0]))

latex_l1 = TLatex()
latex_l1.SetTextSize(0.05)
latex_l1.SetTextFont(42)
latex_l1.DrawLatex(4.8, 6e6, str(p[1]))

latex_l1 = TLatex()
latex_l1.SetTextSize(0.05)
latex_l1.SetTextFont(42)
latex_l1.DrawLatex(6.8, 6e6, str(p[2]))

latex_l1 = TLatex()
latex_l1.SetTextSize(0.05)
latex_l1.SetTextFont(42)
latex_l1.DrawLatex(5.35, 1.5e7, "average number of positrons")

leg = TLegend(0.7, 0.35, 0.9, 0.5);
leg.AddEntry(gr, "doublets", "p");
leg.AddEntry(gr2, "triplets", "p");
# leg.AddEntry(ff, "a #cdot exp(#xi)", "p");
# leg.AddEntry(ff2, "a #cdot exp(#xi)", "p");
leg.Draw();

canvas.SaveAs(f"doublet_triplet_multiplicities_vs_xi.pdf")                      