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

class FinalPlotAndTable{

public:
  FinalPlotAndTable();
  TH1F *DrawOverflow(TH1F *h);
  void doOverflow();
  void doOF();
  int Nbins;
};


FinalPlotAndTable::FinalPlotAndTable(){

  /*
  TFile* inFile = new TFile("plots/htotal_root_ZprimeRecomass.root","READ");
  TFile* outFile = new TFile("plots/htotal_root_ZprimeRecomass_binned.root","RECREATE");

  std::vector<TH1F*> hists,hists_OF;

  hists.push_back((TH1F*) inFile->Get("htotaldata"));
  hists.push_back((TH1F*) inFile->Get("htotalHisto"));
  hists.push_back((TH1F*) inFile->Get("htotalHistoRatio"));
  hists.push_back((TH1F*) inFile->Get("hfourlepbestmass_4l_afterSel_new_DY"));
  hists.push_back((TH1F*) inFile->Get("hfourlepbestmass_4l_afterSel_new_diBoson"));
  hists.push_back((TH1F*) inFile->Get("hfourlepbestmass_4l_afterSel_new_Tlike"));
  hists.push_back((TH1F*) inFile->Get("hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData"));
  hists.push_back((TH1F*) inFile->Get("hfourlepbestmass_4l_afterSel_WJetsFromFakeRateFromMC"));

  const int NMBINS = 100;
  const double MMIN = 60., MMAX = 3000.;
  double logMbins[NMBINS+1];
  float binNormNr=0.;

  for (int ibin = 0; ibin <= NMBINS; ibin++) {
    logMbins[ibin] = exp(log(MMIN) + (log(MMAX)-log(MMIN))*ibin/NMBINS);
    cout << logMbins[ibin] << endl;
  }

  TH1F* h = new TH1F("n","t", NMBINS, logMbins);
  
  for(size_t histNr=0;histNr<hists.size();histNr++){
    hists[histNr]->SetBins(100,60.,3000.);
    //for(int binNr=1;binNr<=hists[histNr]->GetNbinsX();binNr++){
    for(int binNr=1;binNr<=100;binNr++){
      logMbins[binNr] = exp(log(MMIN) + (log(MMAX)-log(MMIN))*binNr/NMBINS);
      cout << logMbins[binNr] << endl;
      //if (binNr>1) binNormNr=logMbins[binNr]-logMbins[binNr-1];
      //else 
      binNormNr=logMbins[binNr];
      hists[histNr]->SetBinContent(binNr,hists[histNr]->GetBinContent(binNr)*binNormNr/hists[histNr]->GetBinWidth(binNr));
      hists[histNr]->SetBinError(binNr,hists[histNr]->GetBinError(binNr)*binNormNr/hists[histNr]->GetBinWidth(binNr));
    }
    //TH1F *tmphist=DrawOverflow(hists[histNr]);
    //tmphist->SetDirectory(outFile);
    hists[histNr]->SetDirectory(outFile);
    //delete tmphist;
  }
 

  outFile->Write();

  delete inFile;
  delete outFile;
  */
  //doOverflow();
  doOF();

}

void FinalPlotAndTable::doOF(){

  TFile *_file000 = TFile::Open("/lustre/cms/store/user/defilip/ZprimeAnalysis/histos/histosZprimeMuMu_combined/output_ZprimeToMuMu_M-5000_TuneCUETP8M1_13TeV-pythia8.root");
  TH1F *htotalSignal = (TH1F*)_file000->Get("ZprimeRecomassBinWidth");
  TH1F *htotalSignalOF=DrawOverflow(htotalSignal);

  char htotal_root_OF_0000[300];
  sprintf(htotal_root_OF_0000,"plots/h_ZprimeRecomass_Signal_OF.root");
  TFile *file0000 = new TFile(htotal_root_OF_0000, "RECREATE");
  file0000->cd();
  htotalSignalOF->Write();
  file0000->Write();
  file0000->Close();

  //

  TFile *_file0 = TFile::Open("plots/h_ZprimeRecomass_data.root");
  TH1F *htotaldata = (TH1F*)_file0->Get("ZprimeRecomassBinWidth");
  TH1F *htotaldataOF=DrawOverflow(htotaldata);

  char htotal_root_OF[300];
  sprintf(htotal_root_OF,"plots/h_ZprimeRecomass_data_OF.root");
  TFile *file00 = new TFile(htotal_root_OF, "RECREATE");
  file00->cd();
  htotaldataOF->Write();
  file00->Write();
  file00->Close();
  
  // 

  TFile *_file1 = TFile::Open("plots/h_ZprimeRecomass_DY.root");
  TH1F *htotalDY = (TH1F*)_file1->Get("ZprimeRecomassBinWidth");
  TH1F *htotalDYOF=DrawOverflow(htotalDY);

  char htotal_root_OF_1[300];
  sprintf(htotal_root_OF_1,"plots/h_ZprimeRecomass_DY_OF.root");
  TFile *file11 = new TFile(htotal_root_OF_1, "RECREATE");
  file11->cd();
  htotalDYOF->Write();
  file11->Write();
  file11->Close();

  TFile *_file2 = TFile::Open("plots/h_ZprimeRecomass_DiBoson.root");
  TH1F *htotaldiBoson = (TH1F*)_file2->Get("ZprimeRecomassBinWidth");
  TH1F *htotaldiBosonOF=DrawOverflow(htotaldiBoson);

  char htotal_root_OF_2[300];
  sprintf(htotal_root_OF_2,"plots/h_ZprimeRecomass_DiBoson_OF.root");
  TFile *file22 = new TFile(htotal_root_OF_2, "RECREATE");
  file22->cd();
  htotaldiBosonOF->Write();
  file22->Write();
  file22->Close();

  TFile *_file3 = TFile::Open("plots/h_ZprimeRecomass_Tlike.root");
  TH1F *htotalTlike = (TH1F*)_file3->Get("ZprimeRecomassBinWidth");
  TH1F *htotalTlikeOF=DrawOverflow(htotalTlike);

  char htotal_root_OF_3[300];
  sprintf(htotal_root_OF_3,"plots/h_ZprimeRecomass_Tlike_OF.root");
  TFile *file33 = new TFile(htotal_root_OF_3, "RECREATE");
  file33->cd();
  htotalTlikeOF->Write();
  file33->Write();
  file33->Close();
  
  
  TFile *_file4 = TFile::Open("plots/FR-DiJets-Data-OS-2800pb-BinWidth.root");
  TH1F *htotalDiJets = (TH1F*)_file4->Get("DataSub");
  TH1F *htotalDiJetsOF=DrawOverflow(htotalDiJets);

  char htotal_root_OF_4[300];
  sprintf(htotal_root_OF_4,"plots/FR-DiJets-Data-OS-2800pb-BinWidth_OF.root");
  TFile *file44 = new TFile(htotal_root_OF_4, "RECREATE");
  file44->cd();
  htotalDiJetsOF->Write();
  file44->Write();
  file44->Close();

  TFile *_file5 = TFile::Open("plots/FR-Wjets-25nsMC-OS-BinWidth-2800pb.root");
  TH1F *htotalWJets = (TH1F*)_file5->Get("WjetsHisto");
  TH1F *htotalWJetsOF=DrawOverflow(htotalWJets);
  
  char htotal_root_OF_5[300];
  sprintf(htotal_root_OF_5,"plots/FR-Wjets-25nsMC-OS-BinWidth-2800pb_OF.root");
  TFile *file55 = new TFile(htotal_root_OF_5, "RECREATE");
  file55->cd();
  htotalWJetsOF->Write();
  file55->Write();
  file55->Close();
  
  /*
  TFile *_file6 = TFile::Open("plots/FakeRate-DiJets-Wjets-OS-2800pb-BinWidth.root");
  TH1F *htotalFRJets = (TH1F*)_file6->Get("WjetsHisto");
  TH1F *htotalFRJetsOF=DrawOverflow(htotalWJets);
  char htotal_root_OF_6[300];
  sprintf(htotal_root_OF_6,"plots/FR-Wjets-25nsMC-OS-BinWidth-2673pb_OF.root");
  TFile *file66 = new TFile(htotal_root_OF_6, "RECREATE");
  file66->cd();
  htotalWJetsOF->Write();
  file66->Write();
  file66->Close();
  */

  /*
  char htotal_root_OF[300];
  sprintf(htotal_root_OF,"plots/htotal_root_ZprimeRecomass_OF_BinWidth.root");
  TFile *file1 = new TFile(htotal_root_OF, "RECREATE");
  file1->cd();
  htotaldataOF->Write();
  htotalDYOF->Write();
  htotaldiBosonOF->Write();
  htotalTlikeOF->Write();
  htotalDiJetsOF->Write();
  htotalWJetsOF->Write();
  file1->Write();
  file1->Close();
  */

}

void FinalPlotAndTable::doOverflow(){

 
 cout << "Producing final plots and table" << endl;	

 TFile *_file0 = TFile::Open("plots/htotal_root_ZprimeRecomass_binned.root");

 TH1F *htotaldata = (TH1F*)_file0->Get("htotaldata");
 htotaldata->Rebin(1);
 TH1F *htotaldataOF=DrawOverflow(htotaldata);

 TH1F *htotalHisto = (TH1F*)_file0->Get("htotalHisto");
 htotalHisto->Rebin(1);
 TH1F *htotalHistoOF=DrawOverflow(htotalHisto);

 TH1F *htotalHistoRatio = (TH1F*)_file0->Get("htotalHistoRatio");
 htotalHistoRatio->Rebin(1);
 TH1F *htotalHistoRatioOF=DrawOverflow(htotalHistoRatio);
 
 TH1F *htotalDY = (TH1F*)_file0->Get("hfourlepbestmass_4l_afterSel_new_DY");
 htotalDY->Rebin(1);
 TH1F *htotalDYOF=DrawOverflow(htotalDY);
 
 TH1F *htotaldiBoson = (TH1F*)_file0->Get("hfourlepbestmass_4l_afterSel_new_diBoson");
 htotaldiBoson->Rebin(1);
 TH1F *htotaldiBosonOF=DrawOverflow(htotaldiBoson);

 TH1F *htotalTlike = (TH1F*)_file0->Get("hfourlepbestmass_4l_afterSel_new_Tlike");
 htotalTlike->Rebin(1);
 TH1F *htotalTlikeOF=DrawOverflow(htotalTlike);
 
 TH1F *htotalDiJets = (TH1F*)_file0->Get("hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData");
 htotalDiJets->Rebin(1);
 TH1F *htotalDiJetsOF=DrawOverflow(htotalDiJets);

 TH1F *htotalWJets = (TH1F*)_file0->Get("hfourlepbestmass_4l_afterSel_WJetsFromFakeRateFromMC");
 htotalWJets->Rebin(1);
 TH1F *htotalWJetsOF=DrawOverflow(htotalWJets);
 

 char htotal_root_OF[300]; 
 sprintf(htotal_root_OF,"plots/htotal_root_ZprimeRecomass_OF.root");
 TFile *file1 = new TFile(htotal_root_OF, "RECREATE");
 file1->cd();
 htotaldataOF->Write();
 htotalHistoOF->Write();
 htotalHistoRatioOF->Write();
 htotalDYOF->Write();
 htotaldiBosonOF->Write();
 htotalTlikeOF->Write();
 htotalDiJetsOF->Write();
 htotalWJetsOF->Write();
 file1->Write();
 file1->Close();
 

}


TH1F* FinalPlotAndTable::DrawOverflow(TH1F *h)
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
