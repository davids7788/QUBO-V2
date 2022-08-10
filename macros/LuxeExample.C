
#include "TH1F.h"
#include "TPad.h"
#include "TRandom.h"

#include "LuxeLabels.h"
#include "LuxeStyle.h"

#ifdef __CLING__
// these are not headers - do not treat them as such - needed for ROOT6
#include "LuxeLabels.C"
#endif

void LuxeExample()
{
  SetLuxeStyle();

  // generate some random data
  TH1F* hpx  = new TH1F("hpx","This is the px distribution",100,-3,3);
  hpx->FillRandom("gaus",10000);
  hpx->GetXaxis()->SetTitle("E_{e}  [GeV]");
  hpx->GetYaxis()->SetTitle("dN_{e}/dE_{e} [1/GeV]");
  hpx->GetXaxis()->SetTitleOffset(1.4);
  hpx->GetYaxis()->SetTitleOffset(1.4);
  hpx->SetFillColor(4);
    
  // and plot it
  hpx->Draw();


  //make a legend
  TLegend* leg=new TLegend(0.7,0.8,0.9,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hpx,"A LUXE Histo","f");
  leg->Draw();


#ifdef __CINT__
  gROOT->LoadMacro("LuxeLabels.C");
#endif

  LUXELabel(0.2,0.85,"CDR");

  return;
}

#ifndef __CINT__
int main() { 
  LuxeExample();
  gPad->Print("luxe.png");
  return 0;
}
#endif
