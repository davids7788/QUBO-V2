import matplotlib.pyplot as plt
import ROOT
from ROOT import *
import numpy as np
import sys
from array import array

ROOT.gROOT.LoadMacro("macros/LuxeStyle.C")
ROOT.gROOT.LoadMacro("macros/LuxeLabels.C")

SetLuxeStyle()

[9293, 83056, 1527225]
[5393, 48446, 1692158]

x = [4, 5, 7]
d = [9293, 83056, 1527225]
t = [5393, 48446, 1692158]

num_doublets = array("d")
num_triplets = array("d")
xi = array("d")

for i in range(3):
    num_doublets.append(d[i])
    num_triplets.append(t[i])
    xi.append(x[i])

canvas = TCanvas( 'c1', 'd/t multiplicities', 800, 600)
gPad.SetLogy()
gr = TGraph(3, xi, num_doublets)
gr.SetLineColor(38)
gr.GetXaxis().SetNdivisions(305)
gr.SetMarkerColor(38)
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
gr2.SetLineColor(46)
gr2.SetMarkerColor(46)
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

LUXELabel(0.2, 0.85, "e-laser, phase-0")

leg = TLegend(0.65, 0.3, 0.8, 0.4);
leg.AddEntry(gr, "doublets", "p");
leg.AddEntry(gr2, "triplets", "p");
# leg.AddEntry(ff, "a #cdot exp(#xi)", "p");
# leg.AddEntry(ff2, "a #cdot exp(#xi)", "p");
leg.Draw();

canvas.SaveAs(f"doublet_triplet_multiplicities_vs_xi.pdf")                      