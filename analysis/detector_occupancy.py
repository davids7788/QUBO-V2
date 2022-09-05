import ROOT
from ROOT import *
import numpy as np
import sys
sys.path.append("../src")

ROOT.gROOT.LoadMacro("macros/LuxeStyle.C")
ROOT.gROOT.LoadMacro("macros/LuxeLabels.C")

SetLuxeStyle()

file = sys.argv[1]
data = np.load(file, allow_pickle=True)[()]['True detector hits']['Plane 0']
xi = file.split("e0gpc_")[1].split("_")[0]

# Setup as of beginning of September 2021
CHIPS_SIZE_X = 0.02994176
CHIPS_SIZE_Y = 0.01376256

X_MID_POSITION_CHIPS_LIST = [0.06773, 0.09783, 0.12793, 0.15803, 0.18813, 0.21823, 0.24833, 0.27843, 0.30853,
                             0.29853, 0.32863, 0.35873, 0.38883, 0.41893, 0.44903, 0.47913, 0.50923, 0.53933]

Z_POSITIONS = [3.9620125, 3.9500125, 4.0620125, 4.0500125, 4.1620125, 4.1500125, 4.2620125, 4.2500125]

PIXEL_X = 1024
PIXEL_Y = 512

# Simplified setup
CHIPS_SIZE_X_SIMPLE = 18 * CHIPS_SIZE_X
CHIPS_SIZE_Y_SIMPLE = CHIPS_SIZE_Y

Z_POSITIONS_SIMPLE = [(Z_POSITIONS[0] + Z_POSITIONS[1]) / 2,
                      (Z_POSITIONS[2] + Z_POSITIONS[3]) / 2,
                      (Z_POSITIONS[4] + Z_POSITIONS[5]) / 2,
                      (Z_POSITIONS[6] + Z_POSITIONS[7]) / 2]

PIXEL_X_SIMPLE = 18 * PIXEL_X
PIXEL_Y_SIMPLE = PIXEL_Y

X_START_POSITION_SIMPLE = X_MID_POSITION_CHIPS_LIST[0] - 0.5 * CHIPS_SIZE_X
X_END_POSITION_SIMPLE = X_MID_POSITION_CHIPS_LIST[-1] + 0.5 * CHIPS_SIZE_X

Y_START_POSITION_SIMPLE = - 0.5 * CHIPS_SIZE_Y_SIMPLE
Y_END_POSITION_SIMPLE = 0.5 * CHIPS_SIZE_Y_SIMPLE


occupancy = ROOT.TH2F('detector occupancy', 
                      'x vs y', 
                      1816, 
                      X_START_POSITION_SIMPLE, 
                      X_END_POSITION_SIMPLE, 
                      52, 
                      Y_START_POSITION_SIMPLE, 
                      Y_END_POSITION_SIMPLE)

for value in data.values():
    occupancy.Fill(value[0], value[1])
    
canv = TCanvas("example1", "xplet efficiency ", 600, 400)

occupancy.GetYaxis().SetTitle("y [m]")
occupancy.GetXaxis().SetTitle("x [m]")
occupancy.Draw("COLZ")
LUXELabel(0.2, 0.82)

latex = TLatex()
latex.SetTextSize(0.045)
latex.SetTextFont(42)
latex.DrawLatex(0.055, 0.0072, f"phase-0, e-laser, #xi = {xi}, B = 0.95T, multiple scattering")

latex = TLatex()
latex.SetTextSize(0.045)
latex.SetTextFont(42)
latex.SetTextAngle(90)
latex.DrawLatex(0.645, -0.0008, "occupancy / 100 pixel")


gPad.SetRightMargin(0.15)
gPad.SetTopMargin(0.1)
gPad.Draw()

canv.Draw()

canv.SaveAs(f"occupancy_{xi}.pdf")






















