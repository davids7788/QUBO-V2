import matplotlib.pyplot as plt
import ROOT
from ROOT import *
import numpy as np
import sys
from array import array

ROOT.gROOT.LoadMacro("macros/LuxeStyle.C")
ROOT.gROOT.LoadMacro("macros/LuxeLabels.C")

SetLuxeStyle()


x = [4, 5, 7]
c_1_frate_baseline = [0.027, 0.082, 0.198]



c_1_frate_baseline_d = array("d")

xi = array("d")

for i in range(3):
    c_1_frate_baseline_d.append(c_1_frate_baseline[i])
    xi.append(x[i])
    
canvas = TCanvas( 'c1', 'eff/f_rate vs xi', 800, 600)
gr = TGraph(3, xi, c_1_frate_baseline_d)
gr.SetLineColor(38)
gr.GetXaxis().SetNdivisions(505)
# gr.GetYaxis().SetNdivisions(501)
gr.SetMarkerColor(1)
gr.SetMarkerSize(2.0)
gr.SetLineColor(1)
gr.SetLineWidth(2)
gr.SetMinimum(0.0)
gr.SetMaximum(0.5)
gr.SetMarkerStyle(20)
gr.GetXaxis().SetTitle('#xi')
gr.GetYaxis().SetTitle('fake rate')
gr.Draw('ALP')

canvas.Update()

# gr2 = TGraph(3, xi, frate_GNN_d)
# gr2.SetLineColor(46)
# gr2.SetMarkerColor(2)
# gr2.SetMarkerStyle(21)
# gr2.SetLineColor(2)
# gr2.SetLineWidth(2)
# gr2.SetMarkerSize(2.0)
# gr2.Draw('LPSAME')
# canvas.Update()

# gr3 = TGraph(3, xi, frate_ACTS_d)
# gr3.SetLineColor(46)
# gr3.SetMarkerColor(4)
# gr3.SetLineColor(4)
# gr3.SetLineWidth(2)
# gr3.SetMarkerStyle(22)
# gr3.SetMarkerSize(2.0)
# gr3.Draw('LPSAME')

# canvas.Update()

LUXELabel(0.55, 0.85, "e-laser, phase-0")
leg = TLegend(0.2, 0.6, 0.5, 0.9);
# leg.AddEntry(gr3, "ACTS placeholder", "p")
# leg.AddEntry(gr2, "GNN placeholder", "p")
leg.AddEntry(gr, "bit flip", "p")

leg.Draw()

canvas.SaveAs(f"frate_vs_xi_multi_setups.pdf")  