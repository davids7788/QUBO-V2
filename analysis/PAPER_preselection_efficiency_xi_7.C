#ifdef __CLING__
#pragma cling optimize(0)
#endif
void PAPER_preselection_efficiency_xi_7()
{
//=========Macro generated from canvas: example/preselection efficiency 
//=========  (Mon Mar 27 11:20:02 2023) by ROOT version 6.26/00
   TCanvas *example = new TCanvas("example", "preselection efficiency ",122,114,800,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   example->SetHighLightColor(2);
   example->Range(-2.005957,0.7513044,12.73716,1.055652);
   example->SetFillColor(0);
   example->SetBorderMode(0);
   example->SetBorderSize(2);
   example->SetTickx(1);
   example->SetTicky(1);
   example->SetLeftMargin(0.16);
   example->SetRightMargin(0.05);
   example->SetTopMargin(0.15);
   example->SetBottomMargin(0.16);
   example->SetFrameBorderMode(0);
   example->SetFrameBorderMode(0);
   Double_t xAxis1[18] = {1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 10, 12}; 
   
   TH1F *frame__1 = new TH1F("frame__1","energy",17, xAxis1);
   frame__1->SetMinimum(0.8);
   frame__1->SetMaximum(1.01);
   frame__1->SetStats(0);
   frame__1->SetLineWidth(2);
   frame__1->SetMarkerStyle(20);
   frame__1->SetMarkerSize(1.2);
   frame__1->GetXaxis()->SetTitle("True energy [GeV]");
   frame__1->GetXaxis()->SetRange(0,17);
   frame__1->GetXaxis()->SetLabelFont(42);
   frame__1->GetXaxis()->SetLabelOffset(0.02);
   frame__1->GetXaxis()->SetLabelSize(0.05);
   frame__1->GetXaxis()->SetTitleSize(0.055);
   frame__1->GetXaxis()->SetTitleOffset(1.4);
   frame__1->GetXaxis()->SetTitleFont(42);
   frame__1->GetYaxis()->SetTitle("Selection efficiency");
   frame__1->GetYaxis()->SetLabelFont(42);
   frame__1->GetYaxis()->SetLabelSize(0.05);
   frame__1->GetYaxis()->SetTitleSize(0.055);
   frame__1->GetYaxis()->SetTitleOffset(1.4);
   frame__1->GetYaxis()->SetTitleFont(42);
   frame__1->GetZaxis()->SetLabelFont(42);
   frame__1->GetZaxis()->SetLabelSize(0.05);
   frame__1->GetZaxis()->SetTitleSize(0.05);
   frame__1->GetZaxis()->SetTitleOffset(1);
   frame__1->GetZaxis()->SetTitleFont(42);
   frame__1->Draw("");
   Double_t xAxis1[18] = {1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 10, 12}; 
   
   TEfficiency * gen doublets_clone1 = new TEfficiency("gen doublets_clone","",17,xAxis1);
   
   gen doublets_clone1->SetConfidenceLevel(0.6826895);
   gen doublets_clone1->SetBetaAlpha(1);
   gen doublets_clone1->SetBetaBeta(1);
   gen doublets_clone1->SetWeight(1);
   gen doublets_clone1->SetStatisticOption(0);
   gen doublets_clone1->SetPosteriorMode(0);
   gen doublets_clone1->SetShortestInterval(0);
   gen doublets_clone1->SetTotalEvents(0,0);
   gen doublets_clone1->SetPassedEvents(0,0);
   gen doublets_clone1->SetTotalEvents(1,7059);
   gen doublets_clone1->SetPassedEvents(1,6945);
   gen doublets_clone1->SetTotalEvents(2,16419);
   gen doublets_clone1->SetPassedEvents(2,16401);
   gen doublets_clone1->SetTotalEvents(3,23491);
   gen doublets_clone1->SetPassedEvents(3,23484);
   gen doublets_clone1->SetTotalEvents(4,25463);
   gen doublets_clone1->SetPassedEvents(4,25459);
   gen doublets_clone1->SetTotalEvents(5,25611);
   gen doublets_clone1->SetPassedEvents(5,25599);
   gen doublets_clone1->SetTotalEvents(6,23082);
   gen doublets_clone1->SetPassedEvents(6,23070);
   gen doublets_clone1->SetTotalEvents(7,19620);
   gen doublets_clone1->SetPassedEvents(7,19604);
   gen doublets_clone1->SetTotalEvents(8,16557);
   gen doublets_clone1->SetPassedEvents(8,16547);
   gen doublets_clone1->SetTotalEvents(9,12930);
   gen doublets_clone1->SetPassedEvents(9,12917);
   gen doublets_clone1->SetTotalEvents(10,10010);
   gen doublets_clone1->SetPassedEvents(10,9998);
   gen doublets_clone1->SetTotalEvents(11,7270);
   gen doublets_clone1->SetPassedEvents(11,7262);
   gen doublets_clone1->SetTotalEvents(12,4995);
   gen doublets_clone1->SetPassedEvents(12,4982);
   gen doublets_clone1->SetTotalEvents(13,3552);
   gen doublets_clone1->SetPassedEvents(13,3547);
   gen doublets_clone1->SetTotalEvents(14,2136);
   gen doublets_clone1->SetPassedEvents(14,2132);
   gen doublets_clone1->SetTotalEvents(15,1461);
   gen doublets_clone1->SetPassedEvents(15,1457);
   gen doublets_clone1->SetTotalEvents(16,1233);
   gen doublets_clone1->SetPassedEvents(16,1227);
   gen doublets_clone1->SetTotalEvents(17,240);
   gen doublets_clone1->SetPassedEvents(17,239);
   gen doublets_clone1->SetTotalEvents(18,0);
   gen doublets_clone1->SetPassedEvents(18,0);
   gen doublets_clone1->SetFillColor(19);
   gen doublets_clone1->SetFillStyle(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ff0000");
   gen doublets_clone1->SetMarkerColor(ci);
   gen doublets_clone1->SetMarkerStyle(21);
   gen doublets_clone1->SetMarkerSize(1.5);
   gen doublets_clone1->Draw("psame x0");
   
   Double_t eff_graph_fx3001[17] = {
   1.5,
   2.25,
   2.75,
   3.25,
   3.75,
   4.25,
   4.75,
   5.25,
   5.75,
   6.25,
   6.75,
   7.25,
   7.75,
   8.25,
   8.75,
   9.5,
   11};
   Double_t eff_graph_fy3001[17] = {
   0.9838504,
   0.9989037,
   0.999702,
   0.9998429,
   0.9995315,
   0.9994801,
   0.9991845,
   0.999396,
   0.9989946,
   0.9988012,
   0.9988996,
   0.9973974,
   0.9985923,
   0.9981273,
   0.9972621,
   0.9951338,
   0.9958333};
   Double_t eff_graph_felx3001[17] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t eff_graph_fely3001[17] = {
   0.001643821,
   0.000323846,
   0.0001604622,
   0.0001241923,
   0.0001779839,
   0.0001974778,
   0.0002589365,
   0.0002576017,
   0.0003630594,
   0.0004551495,
   0.0005422146,
   0.0009387887,
   0.0009511491,
   0.001478207,
   0.002159481,
   0.002895145,
   0.00951566};
   Double_t eff_graph_fehx3001[17] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t eff_graph_fehy3001[17] = {
   0.001500261,
   0.0002558951,
   0.0001098856,
   7.517974e-05,
   0.0001333357,
   0.0001479428,
   0.0002016765,
   0.0001877306,
   0.000275153,
   0.0003410819,
   0.0003807083,
   0.0007119603,
   0.0006078905,
   0.0008960171,
   0.001309847,
   0.00192855,
   0.003447118};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(17,eff_graph_fx3001,eff_graph_fy3001,eff_graph_felx3001,eff_graph_fehx3001,eff_graph_fely3001,eff_graph_fehy3001);
   grae->SetName("eff_graph");
   grae->SetTitle("");
   grae->SetFillColor(19);
   grae->SetFillStyle(0);

   ci = TColor::GetColor("#ff0000");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#ff0000");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(21);
   grae->SetMarkerSize(1.5);
   
   TH1F *Graph_eff_graph3001 = new TH1F("Graph_eff_graph3001","",100,0,13.1);
   Graph_eff_graph3001->SetMinimum(0.9804354);
   Graph_eff_graph3001->SetMaximum(1.001689);
   Graph_eff_graph3001->SetDirectory(0);
   Graph_eff_graph3001->SetStats(0);
   Graph_eff_graph3001->SetLineWidth(2);
   Graph_eff_graph3001->SetMarkerStyle(20);
   Graph_eff_graph3001->SetMarkerSize(1.2);
   Graph_eff_graph3001->GetXaxis()->SetLabelFont(42);
   Graph_eff_graph3001->GetXaxis()->SetLabelSize(0.05);
   Graph_eff_graph3001->GetXaxis()->SetTitleSize(0.05);
   Graph_eff_graph3001->GetXaxis()->SetTitleOffset(1.4);
   Graph_eff_graph3001->GetXaxis()->SetTitleFont(42);
   Graph_eff_graph3001->GetYaxis()->SetLabelFont(42);
   Graph_eff_graph3001->GetYaxis()->SetLabelSize(0.05);
   Graph_eff_graph3001->GetYaxis()->SetTitleSize(0.05);
   Graph_eff_graph3001->GetYaxis()->SetTitleOffset(1.4);
   Graph_eff_graph3001->GetYaxis()->SetTitleFont(42);
   Graph_eff_graph3001->GetZaxis()->SetLabelFont(42);
   Graph_eff_graph3001->GetZaxis()->SetLabelSize(0.05);
   Graph_eff_graph3001->GetZaxis()->SetTitleSize(0.05);
   Graph_eff_graph3001->GetZaxis()->SetTitleOffset(1);
   Graph_eff_graph3001->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_eff_graph3001);
   
   grae->Draw("p");
   Double_t xAxis2[18] = {1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 10, 12}; 
   
   TEfficiency * gen triplets passed doublets_clone2 = new TEfficiency("gen triplets passed doublets_clone","",17,xAxis2);
   
   gen triplets passed doublets_clone2->SetConfidenceLevel(0.6826895);
   gen triplets passed doublets_clone2->SetBetaAlpha(1);
   gen triplets passed doublets_clone2->SetBetaBeta(1);
   gen triplets passed doublets_clone2->SetWeight(1);
   gen triplets passed doublets_clone2->SetStatisticOption(0);
   gen triplets passed doublets_clone2->SetPosteriorMode(0);
   gen triplets passed doublets_clone2->SetShortestInterval(0);
   gen triplets passed doublets_clone2->SetTotalEvents(0,0);
   gen triplets passed doublets_clone2->SetPassedEvents(0,0);
   gen triplets passed doublets_clone2->SetTotalEvents(1,4591);
   gen triplets passed doublets_clone2->SetPassedEvents(1,3972);
   gen triplets passed doublets_clone2->SetTotalEvents(2,10928);
   gen triplets passed doublets_clone2->SetPassedEvents(2,10465);
   gen triplets passed doublets_clone2->SetTotalEvents(3,15652);
   gen triplets passed doublets_clone2->SetPassedEvents(3,15498);
   gen triplets passed doublets_clone2->SetTotalEvents(4,16972);
   gen triplets passed doublets_clone2->SetPassedEvents(4,16930);
   gen triplets passed doublets_clone2->SetTotalEvents(5,17062);
   gen triplets passed doublets_clone2->SetPassedEvents(5,17059);
   gen triplets passed doublets_clone2->SetTotalEvents(6,15376);
   gen triplets passed doublets_clone2->SetPassedEvents(6,15376);
   gen triplets passed doublets_clone2->SetTotalEvents(7,13063);
   gen triplets passed doublets_clone2->SetPassedEvents(7,13063);
   gen triplets passed doublets_clone2->SetTotalEvents(8,11027);
   gen triplets passed doublets_clone2->SetPassedEvents(8,11027);
   gen triplets passed doublets_clone2->SetTotalEvents(9,8606);
   gen triplets passed doublets_clone2->SetPassedEvents(9,8606);
   gen triplets passed doublets_clone2->SetTotalEvents(10,6660);
   gen triplets passed doublets_clone2->SetPassedEvents(10,6660);
   gen triplets passed doublets_clone2->SetTotalEvents(11,4839);
   gen triplets passed doublets_clone2->SetPassedEvents(11,4839);
   gen triplets passed doublets_clone2->SetTotalEvents(12,3317);
   gen triplets passed doublets_clone2->SetPassedEvents(12,3317);
   gen triplets passed doublets_clone2->SetTotalEvents(13,2362);
   gen triplets passed doublets_clone2->SetPassedEvents(13,2362);
   gen triplets passed doublets_clone2->SetTotalEvents(14,1420);
   gen triplets passed doublets_clone2->SetPassedEvents(14,1420);
   gen triplets passed doublets_clone2->SetTotalEvents(15,969);
   gen triplets passed doublets_clone2->SetPassedEvents(15,969);
   gen triplets passed doublets_clone2->SetTotalEvents(16,816);
   gen triplets passed doublets_clone2->SetPassedEvents(16,816);
   gen triplets passed doublets_clone2->SetTotalEvents(17,159);
   gen triplets passed doublets_clone2->SetPassedEvents(17,159);
   gen triplets passed doublets_clone2->SetTotalEvents(18,0);
   gen triplets passed doublets_clone2->SetPassedEvents(18,0);
   gen triplets passed doublets_clone2->SetFillColor(19);
   gen triplets passed doublets_clone2->SetFillStyle(0);

   ci = TColor::GetColor("#0000ff");
   gen triplets passed doublets_clone2->SetMarkerColor(ci);
   gen triplets passed doublets_clone2->SetMarkerStyle(22);
   gen triplets passed doublets_clone2->SetMarkerSize(1.5);
   gen triplets passed doublets_clone2->Draw("psame x0");
   
   Double_t eff_graph_fx3002[17] = {
   1.5,
   2.25,
   2.75,
   3.25,
   3.75,
   4.25,
   4.75,
   5.25,
   5.75,
   6.25,
   6.75,
   7.25,
   7.75,
   8.25,
   8.75,
   9.5,
   11};
   Double_t eff_graph_fy3002[17] = {
   0.865171,
   0.9576318,
   0.990161,
   0.9975253,
   0.9998242,
   1,
   1,
   1,
   1,
   1,
   1,
   1,
   1,
   1,
   1,
   1,
   1};
   Double_t eff_graph_felx3002[17] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t eff_graph_fely3002[17] = {
   0.005230682,
   0.00201569,
   0.0008538552,
   0.0004431525,
   0.0001710046,
   0.0001197263,
   0.0001409241,
   0.0001669419,
   0.0002139002,
   0.0002763915,
   0.0003803826,
   0.0005548721,
   0.0007791297,
   0.001295654,
   0.001898115,
   0.002253611,
   0.01151198};
   Double_t eff_graph_fehx3002[17] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t eff_graph_fehy3002[17] = {
   0.005068183,
   0.00192995,
   0.0007886849,
   0.0003799855,
   9.569097e-05,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   grae = new TGraphAsymmErrors(17,eff_graph_fx3002,eff_graph_fy3002,eff_graph_felx3002,eff_graph_fehx3002,eff_graph_fely3002,eff_graph_fehy3002);
   grae->SetName("eff_graph");
   grae->SetTitle("");
   grae->SetFillColor(19);
   grae->SetFillStyle(0);

   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#0000ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(22);
   grae->SetMarkerSize(1.5);
   
   TH1F *Graph_eff_graph3002 = new TH1F("Graph_eff_graph3002","",100,0,13.1);
   Graph_eff_graph3002->SetMinimum(0.8459343);
   Graph_eff_graph3002->SetMaximum(1.014006);
   Graph_eff_graph3002->SetDirectory(0);
   Graph_eff_graph3002->SetStats(0);
   Graph_eff_graph3002->SetLineWidth(2);
   Graph_eff_graph3002->SetMarkerStyle(20);
   Graph_eff_graph3002->SetMarkerSize(1.2);
   Graph_eff_graph3002->GetXaxis()->SetLabelFont(42);
   Graph_eff_graph3002->GetXaxis()->SetLabelSize(0.05);
   Graph_eff_graph3002->GetXaxis()->SetTitleSize(0.05);
   Graph_eff_graph3002->GetXaxis()->SetTitleOffset(1.4);
   Graph_eff_graph3002->GetXaxis()->SetTitleFont(42);
   Graph_eff_graph3002->GetYaxis()->SetLabelFont(42);
   Graph_eff_graph3002->GetYaxis()->SetLabelSize(0.05);
   Graph_eff_graph3002->GetYaxis()->SetTitleSize(0.05);
   Graph_eff_graph3002->GetYaxis()->SetTitleOffset(1.4);
   Graph_eff_graph3002->GetYaxis()->SetTitleFont(42);
   Graph_eff_graph3002->GetZaxis()->SetLabelFont(42);
   Graph_eff_graph3002->GetZaxis()->SetLabelSize(0.05);
   Graph_eff_graph3002->GetZaxis()->SetTitleSize(0.05);
   Graph_eff_graph3002->GetZaxis()->SetTitleOffset(1);
   Graph_eff_graph3002->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_eff_graph3002);
   
   grae->Draw("p");
   Double_t xAxis3[18] = {1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 10, 12}; 
   
   TEfficiency * gen triplets_clone3 = new TEfficiency("gen triplets_clone","",17,xAxis3);
   
   gen triplets_clone3->SetConfidenceLevel(0.6826895);
   gen triplets_clone3->SetBetaAlpha(1);
   gen triplets_clone3->SetBetaBeta(1);
   gen triplets_clone3->SetWeight(1);
   gen triplets_clone3->SetStatisticOption(0);
   gen triplets_clone3->SetPosteriorMode(0);
   gen triplets_clone3->SetShortestInterval(0);
   gen triplets_clone3->SetTotalEvents(0,0);
   gen triplets_clone3->SetPassedEvents(0,0);
   gen triplets_clone3->SetTotalEvents(1,4706);
   gen triplets_clone3->SetPassedEvents(1,3972);
   gen triplets_clone3->SetTotalEvents(2,10946);
   gen triplets_clone3->SetPassedEvents(2,10465);
   gen triplets_clone3->SetTotalEvents(3,15660);
   gen triplets_clone3->SetPassedEvents(3,15498);
   gen triplets_clone3->SetTotalEvents(4,16976);
   gen triplets_clone3->SetPassedEvents(4,16930);
   gen triplets_clone3->SetTotalEvents(5,17074);
   gen triplets_clone3->SetPassedEvents(5,17059);
   gen triplets_clone3->SetTotalEvents(6,15388);
   gen triplets_clone3->SetPassedEvents(6,15376);
   gen triplets_clone3->SetTotalEvents(7,13080);
   gen triplets_clone3->SetPassedEvents(7,13063);
   gen triplets_clone3->SetTotalEvents(8,11038);
   gen triplets_clone3->SetPassedEvents(8,11027);
   gen triplets_clone3->SetTotalEvents(9,8620);
   gen triplets_clone3->SetPassedEvents(9,8606);
   gen triplets_clone3->SetTotalEvents(10,6673);
   gen triplets_clone3->SetPassedEvents(10,6660);
   gen triplets_clone3->SetTotalEvents(11,4847);
   gen triplets_clone3->SetPassedEvents(11,4839);
   gen triplets_clone3->SetTotalEvents(12,3330);
   gen triplets_clone3->SetPassedEvents(12,3317);
   gen triplets_clone3->SetTotalEvents(13,2368);
   gen triplets_clone3->SetPassedEvents(13,2362);
   gen triplets_clone3->SetTotalEvents(14,1424);
   gen triplets_clone3->SetPassedEvents(14,1420);
   gen triplets_clone3->SetTotalEvents(15,974);
   gen triplets_clone3->SetPassedEvents(15,969);
   gen triplets_clone3->SetTotalEvents(16,822);
   gen triplets_clone3->SetPassedEvents(16,816);
   gen triplets_clone3->SetTotalEvents(17,160);
   gen triplets_clone3->SetPassedEvents(17,159);
   gen triplets_clone3->SetTotalEvents(18,0);
   gen triplets_clone3->SetPassedEvents(18,0);
   gen triplets_clone3->SetFillColor(19);
   gen triplets_clone3->SetFillStyle(0);
   gen triplets_clone3->SetMarkerStyle(20);
   gen triplets_clone3->SetMarkerSize(1.5);
   gen triplets_clone3->Draw("psame x0");
   
   Double_t eff_graph_fx3003[17] = {
   1.5,
   2.25,
   2.75,
   3.25,
   3.75,
   4.25,
   4.75,
   5.25,
   5.75,
   6.25,
   6.75,
   7.25,
   7.75,
   8.25,
   8.75,
   9.5,
   11};
   Double_t eff_graph_fy3003[17] = {
   0.8440289,
   0.956057,
   0.9896552,
   0.9972903,
   0.9991215,
   0.9992202,
   0.9987003,
   0.9990034,
   0.9983759,
   0.9980519,
   0.9983495,
   0.9960961,
   0.9974662,
   0.997191,
   0.9948665,
   0.9927007,
   0.99375};
   Double_t eff_graph_felx3003[17] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t eff_graph_fely3003[17] = {
   0.005469578,
   0.002047625,
   0.0008733579,
   0.0004605953,
   0.0002902565,
   0.0002961636,
   0.0003975002,
   0.0003998402,
   0.0005597568,
   0.0007030336,
   0.0008129317,
   0.001406933,
   0.0015103,
   0.00221545,
   0.003457823,
   0.004334246,
   0.01422439};
   Double_t eff_graph_fehx3003[17] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t eff_graph_fehy3003[17] = {
   0.00532045,
   0.001962357,
   0.0008083509,
   0.0003976744,
   0.0002242271,
   0.0002218997,
   0.0003119772,
   0.0002957769,
   0.0004286413,
   0.0005330213,
   0.0005709575,
   0.001067577,
   0.001004593,
   0.001343869,
   0.002215602,
   0.002891587,
   0.005170872};
   grae = new TGraphAsymmErrors(17,eff_graph_fx3003,eff_graph_fy3003,eff_graph_felx3003,eff_graph_fehx3003,eff_graph_fely3003,eff_graph_fehy3003);
   grae->SetName("eff_graph");
   grae->SetTitle("");
   grae->SetFillColor(19);
   grae->SetFillStyle(0);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.5);
   
   TH1F *Graph_eff_graph3003 = new TH1F("Graph_eff_graph3003","",100,0,13.1);
   Graph_eff_graph3003->SetMinimum(0.822471);
   Graph_eff_graph3003->SetMaximum(1.01553);
   Graph_eff_graph3003->SetDirectory(0);
   Graph_eff_graph3003->SetStats(0);
   Graph_eff_graph3003->SetLineWidth(2);
   Graph_eff_graph3003->SetMarkerStyle(20);
   Graph_eff_graph3003->SetMarkerSize(1.2);
   Graph_eff_graph3003->GetXaxis()->SetLabelFont(42);
   Graph_eff_graph3003->GetXaxis()->SetLabelSize(0.05);
   Graph_eff_graph3003->GetXaxis()->SetTitleSize(0.05);
   Graph_eff_graph3003->GetXaxis()->SetTitleOffset(1.4);
   Graph_eff_graph3003->GetXaxis()->SetTitleFont(42);
   Graph_eff_graph3003->GetYaxis()->SetLabelFont(42);
   Graph_eff_graph3003->GetYaxis()->SetLabelSize(0.05);
   Graph_eff_graph3003->GetYaxis()->SetTitleSize(0.05);
   Graph_eff_graph3003->GetYaxis()->SetTitleOffset(1.4);
   Graph_eff_graph3003->GetYaxis()->SetTitleFont(42);
   Graph_eff_graph3003->GetZaxis()->SetLabelFont(42);
   Graph_eff_graph3003->GetZaxis()->SetLabelSize(0.05);
   Graph_eff_graph3003->GetZaxis()->SetTitleSize(0.05);
   Graph_eff_graph3003->GetZaxis()->SetTitleOffset(1);
   Graph_eff_graph3003->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_eff_graph3003);
   
   grae->Draw("p");
   TLatex *   tex = new TLatex(5.5,0.82,"40TW laser, e-laser, #xi = 7");
   tex->SetTextFont(42);
   tex->SetTextSize(0.045);
   tex->SetLineWidth(2);
   tex->Draw();
   
   TLegend *leg = new TLegend(0.65,0.5,0.85,0.7,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("gen doublets_clone","doublet","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("gen triplets passed doublets_clone","triplet","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("gen triplets_clone","total","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   leg->Draw();
   example->Modified();
   example->cd();
   example->SetSelected(example);
}
