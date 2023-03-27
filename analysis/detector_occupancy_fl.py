import ROOT
from ROOT import *
import numpy as np
import sys
sys.path.append("../src/simplified_simulation")
sys.path.append("../src")

ROOT.gROOT.LoadMacro("macros/LuxeStyle.C")
ROOT.gROOT.LoadMacro("macros/LuxeLabels.C")

SetLuxeStyle()

num_files = 10

data = []
for n in range(num_files):
    print(n)
    file = f"/nfs/dust/luxe/user/spatarod/occupancy-study-1.5T/e0gpc_7.0_000{n}_particles_fl.npy"
    for i in range(0, 9):
        single_chip_data = np.load(file, allow_pickle=True)[()]['True detector hits'][f'Plane {i}']
        for entry in single_chip_data.values():
            data.append(entry)

xi = 7

occupancy = ROOT.TH2F('detector occupancy',               # per mm²
                      'x vs y', 
                      500, 
                      50, 
                      550, 
                      14, 
                      - 7, 
                      + 7)

for value in data:
    occupancy.Fill(value[0] * 1000, value[1] * 1000, 1 / num_files)
     
canv = TCanvas("example1", "detector occupancy", 800, 600)


occupancy.GetYaxis().SetTitle("y [mm]")
occupancy.GetXaxis().SetTitle("x [mm]")
occupancy.GetYaxis().SetNdivisions(211)
occupancy.GetXaxis().SetNdivisions(211)

occupancy.Draw("COLZ")

latex = TLatex()
latex.SetTextSize(0.045)
latex.SetTextFont(42)
latex.DrawLatex(55, 7.2, f"40TW laser, e-laser, #xi = {xi}, inner stave first detector layer")

latex = TLatex()
latex.SetTextSize(0.045)
latex.SetTextFont(42)
latex.SetTextAngle(90)
latex.DrawLatex(645, -3, "hits / mm^{2}")


gPad.SetRightMargin(0.15)
gPad.SetTopMargin(0.1)
gPad.Draw()

canv.Draw()

canv.SaveAs(f"detector_occupancy_xi_{xi}.pdf")
canv.SaveAs(f"detector_occupancy_xi_{xi}.C")
