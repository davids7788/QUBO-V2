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

gen_x = np.load("/nfs/dust/luxe/user/spatarod/towards_paper/e-laser/phase-0/gpc/7.0/e0gpc_7.0_0000_sl_gen_xplet_list.npy", allow_pickle=True)[()]

xi = 7

# only full gen particle tracks are considered --> all layers have to be hit
gen_particles = set()

dx_x0 = 0.0526
three_sigma = 0.0011

doublet_angles_low_energy = TH1F('gen doublets low', 'energy', 50, 0.0526 - 0.0013, 0.0526 + 0.0013)
doublet_angles_high_energy = TH1F('gen doublets high', 'energy', 50, 0.0526 - 0.0013, 0.0526 + 0.0013)

triplet_angles_low_energy = TH1F('gen triplets passed doublets low eergy', 'energy', 50, 0.0, 1.2e-3)
triplet_angles_high_energy = TH1F('gen triplets passed doublets high energy', 'energy', 50, 0.0, 1.2e-3)

energy_split = 3

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
            0.0526)) > three_sigma:
        return False
    return True

for xplet in gen_x:
    try:
        d1 = Doublet(-1, 
                 -1,
                 xplet.coordinates[0],
                 xplet.coordinates[1],
                 -1,
                 -1,
                 xplet.energy[0] / 1000,
                 xplet.energy[1] / 1000)
    except:
        pass
    try:
        d2 = Doublet(-1, 
                 -1,
                 xplet.coordinates[1],
                 xplet.coordinates[2],
                 -1,
                 -1,
                 xplet.energy[1] / 1000,
                 xplet.energy[2] / 1000)
    except:
        pass

    try:
        d3 = Doublet(-1, 
                 -1,
                 xplet.coordinates[2],
                 xplet.coordinates[3],
                 -1,
                 -1,
                 xplet.energy[2] / 1000,
                 xplet.energy[3] / 1000)
    except:
        pass

    if d1.energy_1 > 3:
        doublet_angles_low_energy.Fill((d1.hit_2_position[0] - d1.hit_1_position[0]) / x0_at_z_ref(d1.hit_1_position[0],
                                                                                        d1.hit_2_position[0],
                                                                                        d1.hit_1_position[2],
                                                                                        d1.hit_2_position[2])) 
    else:    
        doublet_angles_high_energy.Fill((d1.hit_2_position[0] - d1.hit_1_position[0]) / x0_at_z_ref(d1.hit_1_position[0],
                                                                                        d1.hit_2_position[0],
                                                                                        d1.hit_1_position[2],
                                                                                        d1.hit_2_position[2])) 
    
    if d2.energy_1 > 3:
        doublet_angles_low_energy.Fill((d2.hit_2_position[0] - d2.hit_1_position[0]) / x0_at_z_ref(d2.hit_1_position[0],
                                                                                       d2.hit_2_position[0],
                                                                                       d2.hit_1_position[2],
                                                                                       d2.hit_2_position[2])) 
    else:    
        doublet_angles_high_energy.Fill((d2.hit_2_position[0] - d2.hit_1_position[0]) / x0_at_z_ref(d2.hit_1_position[0],
                                                                                        d2.hit_2_position[0],
                                                                                        d2.hit_1_position[2],
                                                                                        d2.hit_2_position[2]))  
        
    if d3.energy_1 > 3:
        doublet_angles_low_energy.Fill((d3.hit_2_position[0] - d3.hit_1_position[0]) / x0_at_z_ref(d3.hit_1_position[0],
                                                                                       d3.hit_2_position[0],
                                                                                       d3.hit_1_position[2],
                                                                                       d3.hit_2_position[2])) 
    else:    
        doublet_angles_high_energy.Fill((d3.hit_2_position[0] - d3.hit_1_position[0]) / x0_at_z_ref(d3.hit_1_position[0],
                                                                                        d3.hit_2_position[0],
                                                                                        d3.hit_1_position[2],
                                                                                        d3.hit_2_position[2]))                        
                        
    t1 = Triplet(d1, d2, -1)
    t2 = Triplet(d2, d3, -1)

    angles = t1.angles_between_doublets()
    if d1.energy_1 < 3:
        triplet_angles_low_energy.Fill(np.sqrt(angles[0]**2 + angles[1]**2))
    else:
        triplet_angles_high_energy.Fill(np.sqrt(angles[0]**2 + angles[1]**2))

    angles = t2.angles_between_doublets()
    if d1.energy_1 < 3:
        triplet_angles_low_energy.Fill(np.sqrt(angles[0]**2 + angles[1]**2))
    else:
        triplet_angles_high_energy.Fill(np.sqrt(angles[0]**2 + angles[1]**2))          

lower_bound_doublets = TLine(0.0515, 0.0, 0.0515, 0.065)
lower_bound_doublets.SetLineWidth(2)
lower_bound_doublets.SetLineColor(kRed)
lower_bound_doublets.SetLineStyle(7)
upper_bound_doublets = TLine(0.0537, 0.0, 0.0537, 0.065)
upper_bound_doublets.SetLineWidth(2)
upper_bound_doublets.SetLineColor(kRed)
upper_bound_doublets.SetLineStyle(7)
                        
upper_bound_triplets = TLine(1e-3, 0.0, 1e-3, 0.065)
upper_bound_triplets.SetLineWidth(2)
upper_bound_triplets.SetLineColor(kRed)
upper_bound_triplets.SetLineStyle(7)
                        
                        
gStyle.SetPadTickY(0)

canv = TCanvas("example", "doublet angle", 800, 600)

h_frame = TH1F('frame', 'energy', 50, 0.0526 - 0.0013, 0.0526 + 0.0013)
h_frame.GetYaxis().SetTitle("fraction of counts")
h_frame.GetXaxis().SetTitle("dx/x_{0}")
h_frame.GetYaxis().SetTitleSize(0.055)
h_frame.GetXaxis().SetTitleSize(0.055)
h_frame.GetXaxis().SetLabelOffset(0.02)
h_frame.SetMaximum(0.065)
h_frame.Draw("P")

# Teff
doublet_angles_low_energy.Scale(1 / (doublet_angles_low_energy.GetEntries())) # + doublet_angles_high_energy.GetEntries()))
doublet_angles_low_energy.SetMarkerStyle(8)
doublet_angles_low_energy.SetLineColor(kBlack)
doublet_angles_low_energy.SetMarkerSize(1.95)
doublet_angles_low_energy.Draw("HISTSAME")

doublet_angles_high_energy.Scale(1 / (doublet_angles_high_energy.GetEntries())) # + doublet_angles_low_energy.GetEntries()))
doublet_angles_high_energy.SetMarkerStyle(8)
doublet_angles_high_energy.SetLineColor(kBlue)
doublet_angles_high_energy.SetMarkerSize(1.95)
doublet_angles_high_energy.Draw("HISTSAME")

lower_bound_doublets.Draw("SAME")
upper_bound_doublets.Draw("SAME")


LUXELabel(0.25, 0.85, "e-laser, phase-0")

leg = TLegend(0.7, 0.65, 0.85, 0.85)

leg.SetBorderSize(0)
leg.SetHeader(f"#xi = {xi}")
leg.AddEntry(doublet_angles_low_energy, " < 3GeV", "l")
leg.AddEntry(doublet_angles_high_energy, " > 3GeV", "l")
leg.SetFillColor(0)
leg.SetTextSize(0.05)
leg.Draw()


canv.SaveAs(f"doublet_angles_{xi}.pdf")
if canv: 
    canv.Close()
    gSystem.ProcessEvents()

canv2 = TCanvas("example", "triplet angle", 800, 600)

h_frame = TH1F('frame', 'energy', 50, -0.0001, 1.1e-3)
h_frame.GetYaxis().SetTitle("fraction of counts")
h_frame.GetXaxis().SetTitle("#sqrt{#theta_{xz}^{2} + #theta_{yz}^{2}}")
h_frame.GetYaxis().SetTitleSize(0.055)
h_frame.GetXaxis().SetTitleSize(0.055)
h_frame.GetXaxis().SetLabelOffset(0.02)
h_frame.SetMaximum(0.065)
h_frame.Draw("P")


triplet_angles_low_energy.Scale(1 / (triplet_angles_low_energy.GetEntries())) # + triplet_angles_high_energy.GetEntries()))
triplet_angles_low_energy.SetMarkerStyle(8)
triplet_angles_low_energy.SetLineColor(kBlack)
triplet_angles_low_energy.SetMarkerSize(1.95)
triplet_angles_low_energy.Draw("HISTSAME")

triplet_angles_high_energy.Scale(1 / (triplet_angles_high_energy.GetEntries())) # + triplet_angles_low_energy.GetEntries()))
triplet_angles_high_energy.SetMarkerStyle(8)
triplet_angles_high_energy.SetLineColor(kBlue)
triplet_angles_high_energy.SetMarkerSize(1.95)
triplet_angles_high_energy.Draw("HISTSAME")

upper_bound_triplets.Draw("SAME")

LUXELabel(0.45, 0.85, "e-laser, phase-0")
leg = TLegend(0.65, 0.6, 0.8, 0.8)

leg.SetBorderSize(0)
leg.SetHeader(f"#xi = {xi}")
leg.AddEntry(triplet_angles_low_energy, " < 3GeV", "l")
leg.AddEntry(triplet_angles_high_energy, " > 3GeV", "l")
leg.SetFillColor(0)
leg.SetTextSize(0.05)
leg.Draw()

canv2.SaveAs(f"triplet_angles_{xi}.pdf")

