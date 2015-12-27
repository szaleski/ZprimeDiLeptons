#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "ZZStyle.C"
#include "TFile.h"
#include "TColor.h"
#include "TPaveText.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPad.h"
#include "TMath.h"
#include "TSystem.h"
#include <libgen.h>

using namespace std;

class FinalPlotAndTable_SH{

public:
  FinalPlotAndTable_SH();
  TH1F *DrawOverflow(TH1F *h);

  int Nbins;
};


FinalPlotAndTable_SH::FinalPlotAndTable_SH(){
 cout << "Producing final plots and table" << endl;	

 TFile *_file1 = TFile::Open("Results2673pb/Data-Dibosons-TTbarandTTbarLike-MC-OS-2673pb.root");
 TH1F *htotaldata = (TH1F*)_file1->Get("DataHisto");
 TH1F *htotaldataOF=DrawOverflow(htotaldata);

 //TH1F *htotalHisto = (TH1F*)_file0->Get("htotalHisto");
 //TH1F *htotalHistoOF=DrawOverflow(htotalHisto);

 //TFile *_file3 = TFile::Open("Results2673pb/Data-Dibosons-TTbarandTTbarLike-MC-OS-2673pb.root");
 //TH1F *htotalHistoRatio = (TH1F*)_file7->Get("hDivideHisto");
 //TH1F *htotalHistoRatioOF=DrawOverflow(htotalHistoRatio);
 
 
 TFile *_file4 = TFile::Open("Results2673pb/DY-MuMu-MC-OS-allbins-MC-2673pb.root");
 TH1F *htotalDY = (TH1F*)_file4->Get("hMassDYAll6");
 TH1F *htotalDYOF=DrawOverflow(htotalDY);
 
 TFile *_file5 = TFile::Open("Results2673pb/Data-Dibosons-TTbarandTTbarLike-MC-OS-2673pb.root");
 TH1F *htotaldiBoson = (TH1F*)_file5->Get("DiBosonBG");
 TH1F *htotaldiBosonOF=DrawOverflow(htotaldiBoson);

 TFile *_file6 = TFile::Open("Results2673pb/Data-Dibosons-TTbarandTTbarLike-MC-OS-2673pb.root");
 TH1F *htotalTlike = (TH1F*)_file6->Get("TTbarAndTTbarlikeHisto");
 TH1F *htotalTlikeOF=DrawOverflow(htotalTlike);

 TFile *_file7 = TFile::Open("Results2673pb/DiJets-Data-OS-2673pb-FR.root");
 TH1F *htotalDiJets = (TH1F*)_file7->Get("DataSub");
 TH1F *htotalDiJetsOF=DrawOverflow(htotalDiJets);

 TFile *_file8 = TFile::Open("Results2673pb/Wjets-25nsMC-OS-allbins-2673pb.root");
 TH1F *htotalWJets = (TH1F*)_file8->Get("WjetsHisto");
 TH1F *htotalWJetsOF=DrawOverflow(htotalWJets);
 

 char htotal_root_OF[300]; 
 sprintf(htotal_root_OF,"plots/htotal_root_ZprimeRecomass_OF_SH.root");
 TFile *file1 = new TFile(htotal_root_OF, "RECREATE");
 file1->cd();
 htotaldataOF->Write();
 //htotalHistoOF->Write();
 //htotalHistoRatioOF->Write();
 htotalDYOF->Write();
 htotaldiBosonOF->Write();
 htotalTlikeOF->Write();
 htotalDiJetsOF->Write();
 htotalWJetsOF->Write();
 
 file1->Write();
 file1->Close();
 

}


TH1F* FinalPlotAndTable_SH::DrawOverflow(TH1F *h)
{
  // This function paint the histogram h with an extra bin for overflows                                                                                                           
  UInt_t nx    = h->GetNbinsX()+1;
  Double_t *xbins= new Double_t[nx+1];
  for (UInt_t i=0;i<nx;i++)
    xbins[i]=h->GetBinLowEdge(i+1);
  xbins[nx]=xbins[nx-1]+h->GetBinWidth(nx);
  char *tempName= new char[strlen(h->GetName())+10];
  sprintf(tempName,"%swtOverFlow",h->GetName());
  // Book a temporary histogram having ab extra bin for overflows                                                                                                                  
  TH1F *htmp = new TH1F(tempName, h->GetTitle(), nx, xbins);
  // Reset the axis labels                                                                                                                                                         
  htmp->SetXTitle(h->GetXaxis()->GetTitle());
  htmp->SetYTitle(h->GetYaxis()->GetTitle());
  // Fill the new hitogram including the extra bin for overflows                                                                                                                   
  for (UInt_t i=1; i<=nx; i++)
    htmp->Fill(htmp->GetBinCenter(i), h->GetBinContent(i));
  // Fill the underflows                                                                                                                                                           
  htmp->Fill(h->GetBinLowEdge(1)-1, h->GetBinContent(0));
  // Restore the number of entries                                                                                                                                                 
  htmp->SetEntries(h->GetEntries());
  // FillStyle and color                                                                                                                                                           
  htmp->SetFillStyle(h->GetFillStyle());
  htmp->SetFillColor(h->GetFillColor());
  return htmp;
}
