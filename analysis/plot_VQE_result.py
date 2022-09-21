import matplotlib.pyplot as plt
import numpy as np
import ROOT
from ROOT import*
import numpy as np

ROOT.gROOT.LoadMacro("macros/LuxeStyle.C")
ROOT.gROOT.LoadMacro("macros/LuxeLabels.C")

SetLuxeStyle()

ideal_sim = np.load("../../SummerStudent/result_ideal_simH.npy", allow_pickle=True)[()].eigenstate
fake_nairobi = np.load("../../SummerStudent/result_vqe_fake_nairobiH.npy", allow_pickle=True)[()].eigenstate
ibm_nairobi = np.load("../../SummerStudent/result_vqe_ibmH.npy", allow_pickle=True)[()].eigenstate

num_states = 7

list_ideal_sim = list(ideal_sim.items())
list_ideal_sim.sort(key=lambda x: x[1], reverse=True)

list_fake_nairobi = list(fake_nairobi.items())
list_fake_nairobi.sort(key=lambda x: x[1], reverse=True)

list_ibm_nairobi = list(ibm_nairobi.items())
list_ibm_nairobi.sort(key=lambda x: x[1], reverse=True)

try:
    sorted_ideal_sim = list_ideal_sim[0: num_states]
except:
    sorted_ideal_sim = list_ideal_sim[0: len(list_ideal_sim)]
    
try:
    sorted_list_fake_nairobi = list_fake_nairobi[0: num_states]
except:
    sorted_list_fake_nairobi = list_fake_nairobi[0: len(list_fake_nairobi)]

try:
    sorted_list_ibm_nairobi = list_ibm_nairobi[0: num_states]
except:
    sorted_list_ibm_nairobi = list_ibm_nairobi[0: len(list_ibm_nairobi)]
    
states = []

for s1 in sorted_ideal_sim:
    if s1[0] not in states:
        states.append(s1[0])
for s2 in sorted_list_fake_nairobi:
    if s2[0] not in states:
        states.append(s2[0])
for s3 in sorted_list_ibm_nairobi:
    if s3[0] not in states:
        states.append(s3[0])

num_states = 2**len(states[0]) 


values_list_ideal_sim = []
values_list_fake_nairobi = []
values_list_ibm_nairobi = []

ideal_sim_dump = 0
ideal_sim_keys = [key for key, value in list_ideal_sim]
for state in states:
    if state in ideal_sim_keys:   
        for key, value in list_ideal_sim:
            if key == state:
                values_list_ideal_sim.append(value**2)
    else:
        values_list_ideal_sim.append(0.0)
for key, value in list_ideal_sim: 
    if key not in states:
        if value**2 > ideal_sim_dump:
            ideal_sim_dump = value**2


fake_nairobi_dump = 0
fake_nairobi_keys = [key for key, value in list_fake_nairobi]
for state in states:
    if state in fake_nairobi_keys:   
        for key, value in list_fake_nairobi:
            if key == state:
                values_list_fake_nairobi.append(value**2)
    else:
        values_list_fake_nairobi.append(0.0)
for key, value in list_fake_nairobi:
    if key not in states:
        if value**2 > fake_nairobi_dump:
            fake_nairobi_dump = value**2
values_list_fake_nairobi.append(fake_nairobi_dump)
        
ibm_nairobi_dump = 0
ibm_nairobi_keys = [key for key, value in list_ibm_nairobi]
for state in states:
    if state in ibm_nairobi_keys:   
        for key, value in list_ibm_nairobi:
            if key == state:
                values_list_ibm_nairobi.append(value**2)
    else:
        values_list_ibm_nairobi.append(0.0)
for key, value in list_ibm_nairobi:
    if key not in states:
        if value**2 > ibm_nairobi_dump:
            ibm_nairobi_dump = value**2
values_list_ibm_nairobi.append(ibm_nairobi_dump)
states.append('max(rest)')

values_list_ideal_sim.append(0.0)
values_list_fake_nairobi.append(0.0)
values_list_ibm_nairobi.append(0.0)


ideal_sim_hist = TH1F("h1b", "ideal sim hist", len(states), 0, len(states))
fake_nairobi_hist = TH1F("h2b", "fake nairobi", len(states), 0, len(states))
ibm_nairobi_hist = TH1F("h3b", "ibm nairobi", len(states), 0, len(states))

for i in range(1, len(states) + 1):
    ideal_sim_hist.SetBinContent(i, values_list_ideal_sim[i - 1])
    fake_nairobi_hist.SetBinContent(i, values_list_fake_nairobi[i - 1])
    ibm_nairobi_hist.SetBinContent(i, values_list_ibm_nairobi[i - 1])
    ideal_sim_hist.GetXaxis().SetBinLabel(i, states[i - 1])
    ideal_sim_hist.GetXaxis().ChangeLabel(i, 290, 0.05)
for j in range(7):
    ideal_sim_hist.GetYaxis().ChangeLabel(j, -1, 0.05)
    

ideal_sim_hist.SetFillColor(1)
ideal_sim_hist.SetBarWidth(0.3)
ideal_sim_hist.SetBarOffset(0.05)
ideal_sim_hist.SetStats(0)
ideal_sim_hist.SetMinimum(0.0)
ideal_sim_hist.SetMaximum(1.05)

ideal_sim_hist.GetYaxis().SetTitle("Probabilities")


fake_nairobi_hist.SetFillColor(8)
fake_nairobi_hist.SetBarWidth(0.3)
fake_nairobi_hist.SetBarOffset(0.35)

fake_nairobi_hist.SetStats(0)

ibm_nairobi_hist.SetFillColor(9)
ibm_nairobi_hist.SetBarWidth(0.3)
ibm_nairobi_hist.SetBarOffset(0.65)
ibm_nairobi_hist.SetStats(0)


canv = TCanvas("hbar plot", "hbar plot", 800, 600)
ideal_sim_hist.Draw("b")
fake_nairobi_hist.Draw("same b")
ibm_nairobi_hist.Draw("same b")
gPad.SetBottomMargin(0.2)
gPad.Draw()

latex = TLatex()
latex.SetTextSize(0.05)
latex.SetTextFont(42)
latex.DrawLatex(4.65,.45,"N_{shots} = 512")
latex = TLatex()
latex.SetTextSize(0.05)
latex.SetTextFont(42)
latex.DrawLatex(4.65,.35,"error mitigation applied")


leg = TLegend(0.55, 0.6, 0.9, 0.9)

leg.AddEntry(ideal_sim_hist, "Ideal Simulation", "f")
leg.AddEntry(fake_nairobi_hist, "Fake Nairobi", "f")
leg.AddEntry(ibm_nairobi_hist, "IBM Nairobi", "f")
leg.SetHeader("Backends, 7 qubit systems")
leg.Draw()
canv.Draw()

canv.SaveAs("VQE_result.pdf")
