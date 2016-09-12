#include "TProfile.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TF2.h"
#include <string>
using std::string;
void MassResolution(){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetStatFormat("5.3f");
  gStyle->SetFitFormat("5.3f");
  //gStyle->SetStatFormat("5.3f");
  //gStyle->SetFitFormat("5.3f");
  //int BoxValue = 6611; //4680;  
  int BoxValue = 11111111; //4680;  
  gStyle->SetOptFit(11);
  gStyle->SetOptDate(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(BoxValue);
  //gStyle->SetOptStat(kFALSE);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(0); //(10);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPalette(0);
  TPaveLabel pl;
  TLatex lt;
  //lt.SetTextFont(62);
  lt.SetTextFont(70);
  lt.SetTextAlign(12);
  lt.SetTextSize(0.07);
  lt.SetTextColor(1);
  int NbBins   = 10;
  float MinBin = -1.0;
  float MaxBin =  1.0;
  int NbBinsTheta   = 100;
  float MinBinTheta =  -1000.0;  //0.0;
  float MaxBinTheta =  1000.0;  //180.0;
  float FitMin      = -5.;
  float FitMax      = 5.;
  //float FitMin      = -0.2;
  //float FitMax      = 0.2;
  float norm = 1.0;
  //==================================================================  
  //==================================================================  
  //================================================================== 
  //==================================================================
  TCanvas *c111 = new TCanvas("c111","MassResolution",600,600);
  c111->SetFillColor(0);
  TH1F *MassResolution1 = new TH1F("MassResolution1","MassResolution1",100,-100.,100.);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx();
  gPad->SetGridy();
  //========================================================== 
  //                                                              
  //               MC Z -> mumu                                    
  //==========================================================
  TFile *f1 = new TFile("/lustre/cms/store/user/defilip/ZprimeAnalysis/80X/histos/histosZprimeMuMu/output_ZprimeToMuMu_M-5000_TuneCUETP8M1_13TeV-pythia8_tree.root","READ");
  TH1F *massres; f1->GetObject("MassRecoGenDif",massres);
  MassResolution1->Add(massres);

  MassResolution1->Print();
  MassResolution1->SetLineStyle(0);
  MassResolution1->SetLineColor(4);
  MassResolution1->SetLineWidth(2);
  MassResolution1->SetMarkerColor(4); 
  //MassResolution->SetMarkerStyle(20);
  MassResolution1->SetMarkerStyle(2);
  //MassResolution->SetMarkerSize(0.70);
  MassResolution1->SetTitle("");
  MassResolution1->GetXaxis()->SetTitle("(reco dil. mass - gen dil. mass)");
  MassResolution1->GetYaxis()->SetTitle("#Events");
  MassResolution1->GetXaxis()->SetTitleOffset(1.3);
  MassResolution1->GetYaxis()->SetTitleOffset(2.0);
  MassResolution1->GetXaxis()->SetTitleSize(0.05);
  MassResolution1->Draw();
  //MassResolution->Sumw2();
  //MassResolution->GetXaxis()->SetRangeUser(0.95,1.05);
  TF1* fn1 = new TF1("fn1","gaus",FitMin,FitMax);
  fn1->SetLineColor(2);
  MassResolution1->Fit("fn1","R");
  MassResolution1->Draw("sames");
  //MassResolution1->Draw("Dot,sames");
  //MassResolution->Sumw2();
  //=================================================================
  
  gStyle->SetOptFit(1111);
  TH1F *h=(TH1F*)gROOT->FindObject("MassResolution1");
  cout<<"histogramme "<<h<<endl;
  cout<<"Liste des fonctions histogramme "<<h->GetListOfFunctions()<<endl;
  cout<<"Objet histogramme "<<h->GetListOfFunctions()->FindObject("stats")<<endl;
  TPaveStats *st2 =(TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
  st2->SetX1NDC(0.60); //(0.16);
  st2->SetY1NDC(0.50); //(0.40);
  st2->SetX2NDC(0.88); //(0.60);
  st2->SetY2NDC(0.85); //(0.80);
  st2->SetLineColor(kBlack); st2->SetFillColor(0);
  st2->Draw();
   
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText2 = new TPaveText(0.20, 0.70, 0.50, 0.87, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("#sqrt{s} = 13 TeV, CMS, MC");
  TText *t2 = tText2->AddText("Z' #rightarrow #mu #mu [5 TeV]");
  TText *t3 = tText2->AddText("Spring 16");
  tText2->Draw(); 
  //==========================================================
  TLegend *leg = new TLegend(0.60, 0.60, 0.87, 0.80);
  //leg->AddEntry(MassResolution1,"CMSSW745","f");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.04;
  leg->SetTextSize(tsize2); 
  //leg->Draw();  
  //======================================================================= 
  c111->Print("MC-Zprime5000-CMSSW80X-MassResolution-fit.png","png");
  //c111->Print("Pl
  
  //==================================================================  
  //==================================================================  
  //================================================================== 
  //==================================================================
  TCanvas *c1111 = new TCanvas("c1111","MassResolution",600,600);
  c1111->SetFillColor(0);
  TH1F *MassResolutionPull = new TH1F("MassResolutionPull","MassResolutionPull",100,-0.5,0.5);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx();
  gPad->SetGridy();
  //========================================================== 
  //                                                              
  //               MC Z -> mumu                                    
  //==========================================================
  TFile *f2 = new TFile("/lustre/cms/store/user/defilip/ZprimeAnalysis/80X/histos/histosZprimeMuMu/output_ZprimeToMuMu_M-5000_TuneCUETP8M1_13TeV-pythia8_tree.root","READ");
  TH1F *massrespull; f2->GetObject("MassRecoGenDifPull",massrespull);
  MassResolutionPull->Add(massrespull);

  MassResolutionPull->Print();
  MassResolutionPull->SetLineStyle(0);
  MassResolutionPull->SetLineColor(4);
  MassResolutionPull->SetLineWidth(2);
  MassResolutionPull->SetMarkerColor(4); 
  //MassResolutionPull->SetMarkerStyle(20);
  MassResolutionPull->SetMarkerStyle(2);
  //MassResolutionPull->SetMarkerSize(0.70);
  MassResolutionPull->SetTitle("");
  MassResolutionPull->GetXaxis()->SetTitle("(reco dil. mass - gen dil. mass)/(gen dil. mass)");
  MassResolutionPull->GetYaxis()->SetTitle("#Events");
  MassResolutionPull->GetXaxis()->SetTitleOffset(1.3);
  MassResolutionPull->GetYaxis()->SetTitleOffset(2.0);
  MassResolutionPull->GetXaxis()->SetTitleSize(0.05);
  MassResolutionPull->Draw();
  //MassResolutionPull->Sumw2();
  //MassResolutionPull->GetXaxis()->SetRangeUser(0.95,1.05);
  float FitMinPull      = -0.03;
  float FitMaxPull      = 0.03;
  TF1* fn11 = new TF1("fn11","gaus",FitMinPull,FitMaxPull);
  fn11->SetLineColor(2);
  MassResolutionPull->Fit("fn11","R");
  MassResolutionPull->Draw("sames");
  //MassResolutionPull->Draw("Dot,sames");
  //MassResolutionPull->Sumw2();
  //=================================================================
  
  gStyle->SetOptFit(1111);
  TH1F *hh=(TH1F*)gROOT->FindObject("MassResolutionPull");
  cout<<"histogramme "<<hh<<endl;
  cout<<"List of functions histogram "<<hh->GetListOfFunctions()<<endl;
  cout<<"Object histogramme "<<hh->GetListOfFunctions()->FindObject("stats")<<endl;
  TPaveStats *st22 =(TPaveStats*)hh->GetListOfFunctions()->FindObject("stats");
  st22->SetX1NDC(0.60); //(0.16);
  st22->SetY1NDC(0.50); //(0.40);
  st22->SetX2NDC(0.88); //(0.60);
  st22->SetY2NDC(0.85); //(0.80);
  st22->SetLineColor(kBlack); st22->SetFillColor(0);
  st22->Draw();
   
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText22 = new TPaveText(0.20, 0.70, 0.50, 0.87, "brNDC");
  tText22->SetBorderSize(0);
  tText22->SetFillColor(0);
  tText22->SetFillStyle(0);
  TText *t11 = tText22->AddText("#sqrt{s} = 13 TeV, CMS, MC");
  TText *t22 = tText22->AddText("Z' #rightarrow #mu #mu [5 TeV]");
  TText *t33 = tText22->AddText("Spring 16");
  tText22->Draw(); 
  //==========================================================
  TLegend *legg = new TLegend(0.60, 0.60, 0.87, 0.80);
  //legg->AddEntry(MassResolution1,"CMSSW745","f");
  legg->SetBorderSize(0.0);
  legg->SetMargin(0.3);
  legg->SetFillColor(0);
  legg->SetFillStyle(10);
  legg->SetLineColor(0);
  Float_t tsize22 = 0.04;
  legg->SetTextSize(tsize22); 
  //leg->Draw();  
  //======================================================================= 
  c1111->Print("MC-Zprime5000-CMSSW80X-MassResolutionPull-fit.png","png");
  //c1111->Print("Pl


}
