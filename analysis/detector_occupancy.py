import ROOT
from ROOT import *
import numpy as np
import sys
sys.path.append("../src/simplified_simulation")
sys.path.append("../src")

ROOT.gROOT.LoadMacro("macros/LuxeStyle.C")
ROOT.gROOT.LoadMacro("macros/LuxeLabels.C")

SetLuxeStyle()

file = "/nfs/dust/luxe/user/spatarod/towards_paper/e-laser/phase-0/gpc/e0gpc_7.0/e0gpc_7.0_0000_particles_sl.npy"
data = np.load(file, allow_pickle=True)[()]['True detector hits']['Plane 0']
xi = 7


# Setup as of beginning of September 2021, information not needed, jsut added for completeness
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




#occupancy = ROOT.TH2F('detector occupancy',        # actual pixels taken into account
#                      'x vs y', 
#                      2304, 
#                      X_START_POSITION_SIMPLE, 
#                      X_END_POSITION_SIMPLE, 
#                      64, 
#                      Y_START_POSITION_SIMPLE, 
#                      Y_END_POSITION_SIMPLE)

occupancy = ROOT.TH2F('detector occupancy',               # per mm²
                      'x vs y', 
                      500, 
                      50, 
                      550, 
                      14, 
                      - 7, 
                      + 7)

for value in data.values():
    occupancy.Fill(value[0] * 1000, value[1] * 1000)
     
canv = TCanvas("example1", "xplet efficiency ", 800, 600)


occupancy.GetYaxis().SetTitle("y [mm]")
occupancy.GetXaxis().SetTitle("x [mm]")
occupancy.GetYaxis().SetNdivisions(211)
occupancy.GetXaxis().SetNdivisions(211)
# for i in range(15):
#     occupancy.GetYaxis().ChangeLabel(i,-1,-1,-1,-1,-1,f"{7-i}")

occupancy.Draw("COLZ")

latex = TLatex()
latex.SetTextSize(0.045)
latex.SetTextFont(42)
latex.DrawLatex(55, 7.2, f"40TW laser, e-laser, #xi = {xi}, first detector layer")

latex = TLatex()
latex.SetTextSize(0.045)
latex.SetTextFont(42)
latex.SetTextAngle(90)
latex.DrawLatex(645, -3, "hits / mm^{2}")


gPad.SetRightMargin(0.15)
gPad.SetTopMargin(0.1)
gPad.Draw()

canv.Draw()

canv.SaveAs(f"occupancy_{xi}.pdf")
