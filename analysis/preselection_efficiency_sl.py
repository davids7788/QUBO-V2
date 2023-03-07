import ROOT
from ROOT import *
import numpy as np
import sys
sys.path.append("../src")
sys.path.append("../src/pattern_building")
sys.path.append("../src/pattern")
from doublet import Doublet
from triplet import Triplet

ROOT.gROOT.LoadMacro("macros/LuxeStyle.C")
ROOT.gROOT.LoadMacro("macros/LuxeLabels.C")

SetLuxeStyle()

gen_x = np.load("/nfs/dust/luxe/user/spatarod/towards_paper/e-laser/phase-0/gpc/e0gpc_7.0/e0gpc_7.0_0000_sl_gen_xplet_list.npy", allow_pickle=True)[()]
xi = 7

for g in gen_x:
    if g.coordinates == {}:
        print("Yes")

# only full gen particle tracks are considered --> all layers have to be hit
gen_particles = set()

edges_left = [1]
edges_mid = [2.0 + 0.5 * i for i in range(15)]
edges_right = [10, 12]
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
            0.0526)) > 0.0011:
        return False
    return True

for xplet in gen_x:
    gen_doublets.Fill(xplet.energy[0] / 1000)
    gen_doublets.Fill(xplet.energy[1] / 1000)
    gen_doublets.Fill(xplet.energy[2] / 1000)
    gen_triplets.Fill(xplet.energy[0] / 1000)
    gen_triplets.Fill(xplet.energy[1] / 1000)

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

    try:
        t1_pass = True
        t2_pass = True
        if doublet_criteria_check(xplet.coordinates[0][0],
                                  xplet.coordinates[1][0],
                                  xplet.coordinates[0][2],
                                  xplet.coordinates[1][2]):
            matched_doublets.Fill(xplet.energy[0] / 1000)
        else:
            t1_pass = False
    except:
        t1_pass = False
    
    try:
        if doublet_criteria_check(xplet.coordinates[1][0],
                                  xplet.coordinates[2][0],
                                  xplet.coordinates[1][2],
                                  xplet.coordinates[2][2]):
            matched_doublets.Fill(xplet.energy[1] / 1000)

        else:
            t1_pass = False
            t2_pass = False
    except:
        t1_pass = False
        t2_pass = False
    
    try:
        if doublet_criteria_check(xplet.coordinates[2][0],
                              xplet.coordinates[3][0],
                              xplet.coordinates[2][2],
                              xplet.coordinates[3][2]):

            matched_doublets.Fill(xplet.energy[2] / 1000)   
        else:
            t2_pass = False
    except:
        t2_pass = False
    
    
    if t1_pass:
        t1 = Triplet(d1, d2, -1)
        gen_triplets_passed_doublets.Fill(t1.doublet_1.energy_1)
        angles = t1.angles_between_doublets()
        if np.sqrt(angles[0]**2 + angles[1]**2) < 0.001:
            matched_triplets.Fill(t1.doublet_1.energy_1)
    if t2_pass:
        t2 = Triplet(d2, d3, -1)
        gen_triplets_passed_doublets.Fill(t2.doublet_1.energy_1)
        angles = t2.angles_between_doublets()
        if np.sqrt(angles[0]**2 + angles[1]**2) < 0.001:
            matched_triplets.Fill(t2.doublet_1.energy_1)            
    
gStyle.SetPadTickY(0)

canv = TCanvas("example", "pattern_building efficiency ", 800, 600)

gPad.SetTopMargin(0.15)
h_frame = TH1F('frame', 'energy', len(edges) - 1, edges)
h_frame.GetYaxis().SetTitle("Selection efficiency")
h_frame.GetXaxis().SetTitle("True energy [GeV]")
h_frame.GetYaxis().SetTitleSize(0.055)
h_frame.GetXaxis().SetTitleSize(0.055)
# h_frame.GetXaxis().SetNdivisions(115)
# h_frame.GetYaxis().SetNdivisions(115)
h_frame.SetMaximum(1.01)
h_frame.SetMinimum(0.8)
h_frame.GetXaxis().SetLabelOffset(0.02)
h_frame.GetXaxis().SetRangeUser(0.0, 12.0)
gPad.SetTicky()
h_frame.Draw()

# Teff
eff_doublets = TEfficiency(matched_doublets, gen_doublets)
eff_doublets.SetMarkerStyle(21)
eff_doublets.SetMarkerColor(kRed)
eff_doublets.SetMarkerSize(1.5)
eff_doublets.Draw("PSAME X0")
gPad.Update()

eff_d_gr = eff_doublets.GetPaintedGraph().Clone()
for i in range(eff_d_gr.GetN()):
    eff_d_gr.GetEXlow()[i] = 0
    eff_d_gr.GetEXhigh()[i] = 0 
eff_d_gr.SetMarkerStyle(21)
eff_d_gr.SetMarkerColor(kRed)
eff_d_gr.SetLineColor(kRed)
eff_d_gr.SetMarkerSize(1.5)
eff_d_gr.Draw("PSAME")

matched_triplets_clone = matched_triplets.Clone()

eff_triplets = TEfficiency(matched_triplets, gen_triplets_passed_doublets)
eff_triplets.SetMarkerStyle(22)
eff_triplets.SetMarkerColor(kBlue)
eff_triplets.SetMarkerSize(1.5)
eff_triplets.Draw("PSAME X0")
gPad.Update()

eff_tr_gr = eff_triplets.GetPaintedGraph().Clone()
for i in range(eff_tr_gr.GetN()):
    eff_tr_gr.GetEXlow()[i] = 0
    eff_tr_gr.GetEXhigh()[i] = 0 
eff_tr_gr.SetMarkerStyle(22)
eff_tr_gr.SetLineColor(kBlue)
eff_tr_gr.SetMarkerSize(1.5)
eff_tr_gr.Draw("PSAME")


# Teff total

eff_total = TEfficiency(matched_triplets_clone, gen_triplets)
eff_total.SetMarkerStyle(20)
eff_total.SetMarkerColor(kBlack)
eff_total.SetMarkerSize(1.5)
eff_total.Draw("PSAME X0")
gPad.Update()

eff_t_gr = eff_total.GetPaintedGraph().Clone()
for i in range(eff_t_gr.GetN()):
    eff_t_gr.GetEXlow()[i] = 0
    eff_t_gr.GetEXhigh()[i] = 0 
eff_t_gr.SetMarkerStyle(20)
eff_t_gr.SetLineColor(kBlack)
eff_t_gr.SetMarkerSize(1.5)
eff_t_gr.Draw("PSAME")

latex = TLatex()
latex.SetTextSize(0.045)
latex.SetTextFont(42)
latex.DrawLatex(5.5, 0.82, f"40TW laser, e-laser, #xi = {xi}")

leg = TLegend(0.65, 0.5, 0.85, 0.7)
leg.AddEntry(eff_doublets, "doublet", "p")
leg.AddEntry(eff_triplets, "triplet", "p")
leg.AddEntry(eff_total, "total", "p")
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetTextSize(0.05)
leg.Draw()



canv.SaveAs(f"preselection_efficiency_{xi}.pdf")