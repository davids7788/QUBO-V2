import ROOT
from ROOT import *
import numpy as np
import sys
from array import array

ROOT.gROOT.LoadMacro("macros/LuxeStyle.C")
ROOT.gROOT.LoadMacro("macros/LuxeLabels.C")

SetLuxeStyle()


x = [4, 5, 7]
c_1_CTD_eff = [0.701, 0.472, 0.094]

# c_1_CTD_frate = [0.027, 0.082, 0.198]

eff_c_1_CTD_d = array("d")
eff_GNN_d = array("d")
eff_ACTS_d = array("d")

xi = array("d")

for i in range(3):
    eff_c_1_CTD_d.append(c_1_CTD_eff[i])
    # eff_GNN_d.append(eff_GNN[i])
    # eff_ACTS_d.append(eff_ACTS[i])
    xi.append(x[i])
    
canvas = TCanvas( 'c1', 'eff/f_rate vs xi', 800, 600)
gr = TGraph(3, xi, eff_c_1_CTD_d)
gr.SetLineColor(38)
gr.GetXaxis().SetNdivisions(505)
# gr.GetYaxis().SetNdivisions(501)
gr.SetMarkerColor(1)
gr.SetMarkerSize(2.0)
gr.SetLineColor(1)
gr.SetLineWidth(2)
gr.SetMinimum(0.02)
gr.SetMaximum(1.05)
gr.SetMarkerStyle(20)
gr.GetXaxis().SetTitle('#xi')
gr.GetYaxis().SetTitle('efficiency')
gr.Draw('ALP')

canvas.Update()

# gr2 = TGraph(3, xi, eff_GNN_d)
# gr2.SetLineColor(46)
# gr2.SetMarkerColor(2)
# gr2.SetMarkerStyle(21)
 #gr2.SetLineColor(2)
 #gr2.SetLineWidth(2)
 #gr2.SetMarkerSize(2.0)
# gr2.Draw('LPSAME')

canvas.Update()

 #gr3 = TGraph(3, xi, eff_ACTS_d)
 #gr3.SetLineColor(46)
 #gr3.SetMarkerColor(4)
 #gr3.SetLineColor(4)
 #gr3.SetLineWidth(2)
 #gr3.SetMarkerStyle(22)
# gr3.SetMarkerSize(2.0)
# gr3.Draw('LPSAME')

canvas.Update()

LUXELabel(0.55, 0.85, "e-laser, phase-0")
leg = TLegend(0.2, 0.2, 0.5, 0.5);
# leg.AddEntry(gr3, "ACTS placeholder", "p")
# leg.AddEntry(gr2, "GNN placeholder", "p")
leg.AddEntry(gr, "Eigensolver 7Q", "p")
leg.Draw()

canvas.SaveAs(f"eff_vs_xi_multi_setups.pdf")  