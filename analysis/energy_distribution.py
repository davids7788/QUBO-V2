import ROOT
from ROOT import *
import numpy as np
import sys
import h5py
sys.path.append("../src")

ROOT.gROOT.LoadMacro("macros/LuxeStyle.C")
ROOT.gROOT.LoadMacro("macros/LuxeLabels.C")

SetLuxeStyle()

file_xi_4 = h5py.File("/nfs/dust/luxe/group/MCProduction/Signal/ptarmigan-v0.8.1/e-laser/phase0/gpc/4.0/e0gpc_4.0_0000_particles.h5", 'r')
file_xi_5 = h5py.File("/nfs/dust/luxe/group/MCProduction/Signal/ptarmigan-v0.8.1/e-laser/phase0/gpc/5.0/e0gpc_5.0_0000_particles.h5", 'r')
file_xi_7 = h5py.File("/nfs/dust/luxe/user/kropfann/7.0_phase0_dt0.5/00/e0gpc_7.0_0000_particles.h5", 'r')

h_frame = TH1F('xi', 'energy', 40, 1, 11)
xi_4 = TH1F('xi 4', 'energy xi 4', 40, 1, 11)
xi_5 = TH1F('xi 5', 'energy xi 5', 40, 1, 11)
xi_7 = TH1F('xi 7', 'energy xi 7', 40, 1, 11)

sum_w_xi_4 = sum(file_xi_4["/final-state/positron/weight"])
sum_w_xi_5 = sum(file_xi_5["/final-state/positron/weight"])
sum_w_xi_7 = sum(file_xi_7["/final-state/positron/weight"])


for entry, weight in zip(file_xi_4["/final-state/positron/momentum"], file_xi_4["/final-state/positron/weight"]):
    xi_4.Fill(entry[0], weight)
for entry, weight in zip(file_xi_5["/final-state/positron/momentum"], file_xi_5["/final-state/positron/weight"]):
    xi_5.Fill(entry[0], weight)
for entry, weight in zip(file_xi_7["/final-state/positron/momentum"], file_xi_7["/final-state/positron/weight"]):
    xi_7.Fill(entry[0], weight)

canv = TCanvas("example","xplet efficiency ", 800, 600)    
h_frame.GetYaxis().SetTitle("fraction of counts / 250 MeV")
h_frame.GetXaxis().SetTitle("energy [GeV]")
h_frame.SetMaximum(0.085)
h_frame.Draw()

xi_4.SetLineColor(kBlue)
xi_4.SetMarkerSize(0)
xi_4.DrawNormalized("EHISTSAME")

xi_5.SetLineColor(kRed)
xi_5.SetMarkerSize(0)
xi_5.DrawNormalized("EHISTSAME")

xi_7.SetLineColor(kBlack)
xi_7.SetMarkerSize(0)
xi_7.DrawNormalized("EHISTSAME")

LUXELabel(0.2, 0.85)

leg = TLegend(0.6, 0.7, 0.8, 0.9)
leg.AddEntry(xi_4, f"#xi = 4.0", "l")
leg.AddEntry(xi_5, f"#xi = 5.0", "l")
leg.AddEntry(xi_7, f"#xi = 7.0", "l")
leg.SetHeader("e-laser, phase-0")
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetTextSize(0.05)
leg.Draw()

canv.SaveAs(f"energy_distribution.pdf")









                         
