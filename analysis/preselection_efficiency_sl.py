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
triplet_list = np.load(preselection_folder + "/triplet_list.npy", allow_pickle=True)
xi = preselection_folder.split("e0gpc_")[1].split("_")[0]

# only full gen particle tracks are considered --> all layers have to be hit
gen_particles = set()

edges_left = [1.5]
edges_mid = [2.0 + 0.4 * i for i in range(15)]
edges_right = [8.0, 8.5, 9.5]
edges = np.array(edges_left + edges_mid + edges_right)

gen_doublets = TH1F('gen doublets', 'energy', len(edges) - 1, edges)
matched_doublets = TH1F('matched doublets', 'energy', len(edges) - 1, edges)
gen_triplets = TH1F('gen triplets', 'energy', len(edges) - 1, edges)
gen_triplets_passed_doublets = TH1F('gen triplets passed doublets', 'energy', len(edges) - 1, edges)
matched_triplets = TH1F('matched triplets', 'energy', len(edges) - 1, edges)




def x0_at_z_ref(x_end: float,
                x_start: float,
                z_end: float,
                z_start: float,
                z_ref: float):
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
    gen_doublets.Fill(xplet.energy[0] / 1000)
    gen_doublets.Fill(xplet.energy[1] / 1000)
    gen_doublets.Fill(xplet.energy[2] / 1000)
    gen_triplets.Fill(xplet.energy[0] / 1000)
    gen_triplets.Fill(xplet.energy[1] / 1000)

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


    t1_pass = True
    t2_pass = True
    if doublet_criteria_check(xplet.coordinates[0][0],
                              xplet.coordinates[1][0],
                              xplet.coordinates[0][2],
                              xplet.coordinates[1][2]):
        matched_doublets.Fill(xplet.energy[0] / 1000)
    else:
        t1_pass = False
    
    if doublet_criteria_check(xplet.coordinates[1][0],
                                  xplet.coordinates[2][0],
                                  xplet.coordinates[1][2],
                                  xplet.coordinates[2][2]):
        matched_doublets.Fill(xplet.energy[1] / 1000)

    else:
        t1_pass = False
        t2_pass = False
    if doublet_criteria_check(xplet.coordinates[2][0],
                              xplet.coordinates[3][0],
                              xplet.coordinates[2][2],
                              xplet.coordinates[3][2]):

        matched_doublets.Fill(xplet.energy[2] / 1000)   
    else:
        t2_pass = False
    t1 = Triplet(d1, d2, -1)
    t2 = Triplet(d2, d3, -1)
    
    if t1_pass:
        gen_triplets_passed_doublets.Fill(t1.doublet_1.energy_1)
        angles = t1.angles_between_doublets()
        if np.sqrt(angles[0]**2 + angles[1]**2) < 0.001:
            matched_triplets.Fill(t1.doublet_1.energy_1)
    if t2_pass:
        gen_triplets_passed_doublets.Fill(t2.doublet_1.energy_1)
        angles = t2.angles_between_doublets()
        if np.sqrt(angles[0]**2 + angles[1]**2) < 0.001:
            matched_triplets.Fill(t2.doublet_1.energy_1)            
    
gStyle.SetPadTickY(0)

canv = TCanvas("example", "preselection efficiency ", 800, 600)

h_frame = TH1F('frame', 'energy', len(edges) - 1, edges)
h_frame.GetYaxis().SetTitle("fraction of counts")
h_frame.GetXaxis().SetTitle("true energy [GeV]")
h_frame.GetYaxis().SetTitleSize(0.055)
h_frame.GetXaxis().SetTitleSize(0.055)
h_frame.GetXaxis().SetNdivisions(105)
h_frame.GetYaxis().SetNdivisions(105)
h_frame.SetMaximum(1.01)
h_frame.SetMinimum(0.8)
h_frame.GetXaxis().SetLabelOffset(0.02)
h_frame.GetXaxis().SetRangeUser(0.0, 12.0)
h_frame.Draw("PSAME")

# Teff
eff_doublets = TEfficiency(matched_doublets, gen_doublets)
eff_doublets.SetMarkerStyle(8)
eff_doublets.SetMarkerColor(kRed)
eff_doublets.SetMarkerSize(1.95)
eff_doublets.Draw("PSAME")

matched_triplets_clone = matched_triplets.Clone()

eff_triplets = TEfficiency(matched_triplets, gen_triplets_passed_doublets)
eff_triplets.SetMarkerStyle(8)
eff_triplets.SetMarkerColor(kBlue)
eff_triplets.SetMarkerSize(1.95)
eff_triplets.Draw("PSAME")

# Teff total

eff_total = TEfficiency(matched_triplets_clone, gen_triplets)
eff_total.SetMarkerStyle(8)
eff_total.SetMarkerColor(kBlack)
eff_total.SetMarkerSize(1.95)
eff_total.Draw("PSAME")



LUXELabel(0.8, 0.2)

leg = TLegend(0.6, 0.4, 0.85, 0.65)
leg.AddEntry(eff_doublets, "doublet efficiency", "p")
leg.AddEntry(eff_triplets, "triplet efficiency", "p")
leg.AddEntry(eff_total, "total efficiency", "p")
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetTextSize(0.04)
leg.SetHeader(f"#xi = {xi}")
leg.Draw()


canv.SaveAs(f"preselection_efficiency_{xi}.pdf")