import ROOT
from ROOT import *
import numpy as np
import sys
sys.path.append("../src")
sys.path.append("../src/preselection")
sys.path.append("../src/pattern")
from doublet import Doublet
from triplet import Triplet

ROOT.gROOT.LoadMacro("macros/LuxeStyle.C")
ROOT.gROOT.LoadMacro("macros/LuxeLabels.C")

SetLuxeStyle()

preselection_folder = sys.argv[1]
gen_x = np.load(preselection_folder.split("e0gpc_")[0] + "/e0gpc_5.0_0000_sl_gen_xplet_list.npy", allow_pickle=True)
print(len(gen_x))
xi = preselection_folder.split("e0gpc_")[1].split("_")[0]

# only full gen particle tracks are considered --> all layers have to be hit
gen_particles = set()

doublet_angles = TH1F('gen doublets', 'energy', 50, 0.0526 - 0.0012, 0.0526 + 0.0012)
triplet_angles = TH1F('gen triplets passed doublets', 'energy', 50, 0.0, 1.1e-3)





def x0_at_z_ref(x_end: float,
                x_start: float,
                z_end: float,
                z_start: float,
                z_ref=3.9560125):
    """Function for calculation x position of a doublet at a z-reference value.
    :param x_end: x-position of second hit
    :param x_start: x-position of first hit
    :param z_end: z-position of second hit
    :param z_start: z-position of first hit
    :param z_ref: z-value reference layer
    :return:
        x_position at the reference layer
    """
    dx = x_end - x_start
    dz = z_end - z_start
    return x_end - dx * abs(z_end - z_ref) / dz    
    
def doublet_criteria_check(x1: float,
                           x2: float,
                           z1: float,
                           z2: float,
                           z_ref=3.9560125):
    """Checks if hits may be combined to doublets, applying dx/x0 criterion
    :param x1: x value first hit
    :param x2: x value second hit
    :param z1: z value first hit
    :param z2: z value second hit
    :param z_ref: z-value reference
    :return:
        True if criteria applies, else False
    """
    if abs(((x2 - x1) / x0_at_z_ref(x1, x2, z1, z2, z_ref) -
            0.0526)) > 0.0009:
        return False
    return True

for xplet in gen_x:
    d1 = Doublet(-1, 
                 -1,
                 xplet.coordinates[0],
                 xplet.coordinates[1],
                 -1,
                 -1,
                 xplet.energy[0] / 1000,
                 xplet.energy[1] / 1000)
    d2 = Doublet(-1, 
                 -1,
                 xplet.coordinates[1],
                 xplet.coordinates[2],
                 -1,
                 -1,
                 xplet.energy[1] / 1000,
                 xplet.energy[2] / 1000)

    d3 = Doublet(-1, 
                 -1,
                 xplet.coordinates[2],
                 xplet.coordinates[3],
                 -1,
                 -1,
                 xplet.energy[2] / 1000,
                 xplet.energy[3] / 1000)
    doublet_angles.Fill((d1.hit_2_position[0] - d1.hit_1_position[0]) / x0_at_z_ref(d1.hit_1_position[0],
                                                                                    d1.hit_2_position[0],
                                                                                    d1.hit_1_position[2],
                                                                                    d1.hit_2_position[2]))
    doublet_angles.Fill((d2.hit_2_position[0] - d2.hit_1_position[0]) / x0_at_z_ref(d2.hit_1_position[0],
                                                                                    d2.hit_2_position[0],
                                                                                    d2.hit_1_position[2],
                                                                                    d2.hit_2_position[2]))   
    doublet_angles.Fill((d3.hit_2_position[0] - d3.hit_1_position[0]) / x0_at_z_ref(d3.hit_1_position[0],
                                                                                    d3.hit_2_position[0],
                                                                                    d3.hit_1_position[2],
                                                                                    d3.hit_2_position[2]))                       
                        
    t1 = Triplet(d1, d2, -1)
    t2 = Triplet(d2, d3, -1)

    angles = t1.angles_between_doublets()
    triplet_angles.Fill(np.sqrt(angles[0]**2 + angles[1]**2))

    angles = t2.angles_between_doublets()
    triplet_angles.Fill(np.sqrt(angles[0]**2 + angles[1]**2))           

lower_bound_doublets = TLine(0.0517, 0.0, 0.0517, 0.065)
lower_bound_doublets.SetLineWidth(2)
lower_bound_doublets.SetLineColor(kRed)
lower_bound_doublets.SetLineStyle(7)
upper_bound_doublets = TLine(0.0535, 0.0, 0.0535, 0.065)
upper_bound_doublets.SetLineWidth(2)
upper_bound_doublets.SetLineColor(kRed)
upper_bound_doublets.SetLineStyle(7)
                        
upper_bound_triplets = TLine(1e-3, 0.0, 1e-3, 0.065)
upper_bound_triplets.SetLineWidth(2)
upper_bound_triplets.SetLineColor(kRed)
upper_bound_triplets.SetLineStyle(7)
                        
                        
gStyle.SetPadTickY(0)

canv = TCanvas("example", "doublet angle", 800, 600)

h_frame = TH1F('frame', 'energy', 50, 0.0526 - 0.0012, 0.0526 + 0.0012)
h_frame.GetYaxis().SetTitle("fraction of counts")
h_frame.GetXaxis().SetTitle("[a.u]")
h_frame.GetYaxis().SetTitleSize(0.055)
h_frame.GetXaxis().SetTitleSize(0.055)
h_frame.GetXaxis().SetLabelOffset(0.02)
h_frame.SetMaximum(0.065)
h_frame.Draw("P")

# Teff
doublet_angles.Scale(1 / doublet_angles.GetEntries())
doublet_angles.SetMarkerStyle(8)
doublet_angles.SetMarkerColor(kBlack)
doublet_angles.SetMarkerSize(1.95)
doublet_angles.Draw("HISTSAME")

lower_bound_doublets.Draw("SAME")
upper_bound_doublets.Draw("SAME")


LUXELabel(0.7, 0.8)

leg = TLegend(0.45, 0.25, 0.6, 0.4)
leg.AddEntry(doublet_angles, "dx/x_{0}", "l")
leg.AddEntry(lower_bound_doublets, "selection: #mu #pm 3 #sigma", "l")
leg.SetBorderSize(0)
leg.SetHeader(f"#xi = {xi}")
leg.SetFillColor(0)
leg.SetTextSize(0.04)
leg.Draw()


canv.SaveAs(f"doublet_angles_{xi}.pdf")
if canv: 
    canv.Close()
    gSystem.ProcessEvents()

canv2 = TCanvas("example", "triplet angle", 800, 600)

h_frame = TH1F('frame', 'energy', 50, -0.0001, 1.1e-3)
h_frame.GetYaxis().SetTitle("fraction of counts")
h_frame.GetXaxis().SetTitle("angle [rad]")
h_frame.GetYaxis().SetTitleSize(0.055)
h_frame.GetXaxis().SetTitleSize(0.055)
h_frame.GetXaxis().SetLabelOffset(0.02)
h_frame.SetMaximum(0.065)
h_frame.Draw("P")

triplet_angles.Scale(1 / doublet_angles.GetEntries())
triplet_angles.SetMarkerStyle(8)
triplet_angles.SetMarkerColor(kBlack)
triplet_angles.SetMarkerSize(1.95)
triplet_angles.Draw("HISTSAME")

upper_bound_triplets.Draw("SAME")


LUXELabel(0.7, 0.8)

leg = TLegend(0.55, 0.4, 0.7, 0.6)
leg.AddEntry(triplet_angles, "#sqrt{xz_{angle}^{2} + yz_{angle}^{2}}", "l")
leg.AddEntry(lower_bound_doublets, "selection: #leq 0.001", "l")
leg.SetBorderSize(0)
leg.SetHeader(f"#xi = {xi}")
leg.SetFillColor(0)
leg.SetTextSize(0.04)
leg.Draw()

canv2.SaveAs(f"triplet_angles_{xi}.pdf")

