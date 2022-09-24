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
eff_bit_flip = [0.999, 0.986, 0.889]
eff_GNN = [0.55, 0.55, 0.55]
eff_ACTS = [0.6, 0.6, 0.6]
# frate = [0.004, 0.039, 0.136]

eff_bit_flip_d = array("d")
eff_GNN_d = array("d")
eff_ACTS_d = array("d")

xi = array("d")

for i in range(3):
    eff_bit_flip_d.append(eff_bit_flip[i])
    eff_GNN_d.append(eff_GNN[i])
    eff_ACTS_d.append(eff_ACTS[i])
    xi.append(x[i])
    
canvas = TCanvas( 'c1', 'eff/f_rate vs xi', 800, 600)
gr = TGraph(3, xi, eff_bit_flip_d)
gr.SetLineColor(38)
gr.GetXaxis().SetNdivisions(505)
# gr.GetYaxis().SetNdivisions(501)
gr.SetMarkerColor(1)
gr.SetMarkerSize(2.0)
gr.SetLineColor(1)
gr.SetLineWidth(2)
gr.SetMinimum(0.45)
gr.SetMaximum(1.05)
gr.SetMarkerStyle(20)
gr.GetXaxis().SetTitle('#xi')
gr.GetYaxis().SetTitle('efficiency')
gr.Draw('ALP')

canvas.Update()

gr2 = TGraph(3, xi, eff_GNN_d)
gr2.SetLineColor(46)
gr2.SetMarkerColor(2)
gr2.SetMarkerStyle(21)
gr2.SetLineColor(2)
gr2.SetLineWidth(2)
gr2.SetMarkerSize(2.0)
gr2.Draw('LPSAME')

canvas.Update()

gr3 = TGraph(3, xi, eff_ACTS_d)
gr3.SetLineColor(46)
gr3.SetMarkerColor(4)
gr3.SetLineColor(4)
gr3.SetLineWidth(2)
gr3.SetMarkerStyle(22)
gr3.SetMarkerSize(2.0)
gr3.Draw('LPSAME')

canvas.Update()

LUXELabel(0.55, 0.85, "e-laser, phase-0")
leg = TLegend(0.2, 0.5, 0.5, 0.8);
leg.AddEntry(gr3, "ACTS placeholder", "p")
leg.AddEntry(gr2, "GNN placeholder", "p")
leg.AddEntry(gr, "bit flip", "p")

leg.Draw()

canvas.SaveAs(f"eff_vs_xi.pdf")                      