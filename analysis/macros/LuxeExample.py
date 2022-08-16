import ROOT
import sys
import os
from ROOT import *

ROOT.gROOT.LoadMacro("LuxeStyle.C")
ROOT.gROOT.LoadMacro("LuxeLabels.C")

SetLuxeStyle()

canv=TCanvas("canv","A LUXE Plot",800,600)

h1=TH1F("h1", "A LUXE Histo", 100, -3, 3);
h1.FillRandom("gaus", 10000);

h1.SetYTitle("dN_{e}/dE_{e} [1/GeV]");
h1.SetXTitle("E_{e}  [GeV]");
h1.GetYaxis().SetTitleOffset(1.4);
h1.GetXaxis().SetTitleOffset(1.4);
h1.SetFillColor(kBlue)
h1.Draw()

LUXELabel(0.2,0.85,"CDR")

leg=TLegend(0.7,0.8,0.9,0.9)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetTextSize(0.04)
leg.AddEntry(h1,"A LUXE Histo","f")
leg.Draw()

input("Press enter to continue")

canv.SaveAs("luxe.pdf")
