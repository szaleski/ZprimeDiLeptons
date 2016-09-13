//=========================================================================
//          Analysis code for Z' boson to Mu Mu analysis                  =  
//                  Tag and Proble code for Muon ID                       = 
//                        To run over MC                                  =  
//                    Author: Sherif Elgammal                             =
//                         1/1/2015                                       =
//=========================================================================
#include "TProfile.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TF2.h"
#include <string>
using std::string;
void TagProbeMuonPlots(){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetStatFormat("5.3f");
  gStyle->SetFitFormat("5.3f");
  int BoxValue = 11111111; //4680;  
  gStyle->SetOptFit(11);
  gStyle->SetOptDate(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(BoxValue);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(0); //(10);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPalette(0);
  TPaveLabel pl;
  TLatex lt;
  lt.SetTextFont(70);
  lt.SetTextAlign(12);
  lt.SetTextSize(0.07);
  lt.SetTextColor(1);
  //-----------------------------------
  float Eff1    = 0.992541;
  float Eff2    = 0.997441;
  float Eff3    = 0.998276;
  float Eff4    = 0.993334;
  float Eff5    = 0.994771;
  float Eff6    = 0.975608;
  float Eff7    = 0.984109;
  float Eff8    = 0.964307;
  float Eff9    = 0.99461;
  float Eff10   = 0.992715;
  float Eff11   = 0.997227;
  float Eff12   = 0.997424;
  float Eff13   = 0.99038;

  float error1     = 0.000632637;
  float error2     = 0.000245095;
  float error3     = 0.000196413;
  float error4     = 0.000414215;
  float error5     = 0.000248457;
  float error6     = 0.00135095;
  float error7     = 0.000522484;
  float error8     = 0.00163942;
  float error9     = 0.000251847;
  float error10    = 0.00043307;
  float error11    = 0.000251077;
  float error12    = 0.000247884;
  float error13    = 0.000725104;

  /*
    [0.0,0.9]   = 98.97
    [0.9,1.2]   = 99.30
    [1.2,2.1]   = 99.76
    [2.1,2.4]   = 99.15
    error[0.0,0.9]    = 0.02
    error[0.9,1.2]    = 0.03
    error[1.2,2.1]    = 0.01
    error[2.1,2.4]    = 0.05
  */
 /*
   [0.0,0.9]   = 98.90 (+/-) 0.01
   [0.9,1.2]   = 99.30 (+/-) 0.02
   [1.2,2.1]   = 99.80 (+/-) 0.0
   [2.1,2.4]   = 0.991507
   error[0.0,0.9]    = 0.01
   error[0.9,1.2]    = 0.02
   error[1.2,2.1]    = 0.0
   error[2.1,2.4]    = 0.000296783
 */



  //==============================================
  //  Data
  //==============================================
 float EffD1    = 0.974227;
 float EffD2    = 0.990483;
 float EffD3    = 0.995014;
 float EffD4    = 0.983568;
 float EffD5    = 0.99011;
 float EffD6    = 0.970588;
 float EffD7    = 0.981162;
 float EffD8    = 0.952618;
 float EffD9    = 0.987393;
 float EffD10   = 0.979916;
 float EffD11   = 0.994776;
 float EffD12   = 0.992913;
 float EffD13   = 0.981324;

 float errorD1     = 0.00665461;
 float errorD2     = 0.0026395;
 float errorD3     = 0.00188444;
 float errorD4     = 0.00358574;
 float errorD5     = 0.00190335;
 float errorD6     = 0.00815736;
 float errorD7     = 0.00313972;
 float errorD8     = 0.0108701;
 float errorD9     = 0.00216201;
 float errorD10    = 0.00409956;
 float errorD11    = 0.00197444;
 float errorD12    = 0.0023622;
 float errorD13    = 0.00563094;
 /*
   [0.0,0.9]   = 98.42 (+/-) 0.14
   [0.9,1.2]   = 98.18 (+/-) 0.27
   [1.2,2.1]   = 99.33 (+/-) 0.11
   [2.1,2.4]   = 0.977797
   error[0.0,0.9]    = 0.14
   error[0.9,1.2]    = 0.27
   error[1.2,2.1]    = 0.11
   error[2.1,2.4]    = 0.00435441
 */
 
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  TCanvas *c1 = new TCanvas("c1","Eff. pt",500,500);
  char textpro1[100],textNDF1[100],textRatio1[100];           
  c1->cd(1);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  gPad->SetGridx();
  TH1F* MassProbeTTandProbeTF1    = new TH1F("MassProbeTTandProbeTF1","",10,0.0,100.0);
  TH1F* MassProbeTTandProbeTF2    = new TH1F("MassProbeTTandProbeTF2","",10,0.0,100.0);
  TH1F* MassProbeTTandProbe720TF1 = new TH1F("MassProbeTTandProbe720TF1","",10,0.0,100.0);
  TH1F* MassProbeTTandProbe720TF2 = new TH1F("MassProbeTTandProbe720TF2","",10,0.0,100.0);
  TGraphAsymmErrors *MassEffForTagProbe539  = new TGraphAsymmErrors;
  TGraphAsymmErrors *MassEffForTagProbe720  = new TGraphAsymmErrors;
  TMultiGraph *mg1 = new TMultiGraph;
  //========================================================== 
  //                      Mass 539                         
  //==========================================================
  TFile *f1 = new TFile("Analyse_TagandProbe_CMSSW8012_MC_reHLT_DYtoMuMu120to200_13TeV_tree.root","READ");
  //---------------------------------------------------
  TH1F *PtProbeTT; f1->GetObject("PtProbeTT",PtProbeTT);
  TH1F *PtProbeTP; f1->GetObject("PtProbeTP",PtProbeTP);
  MassProbeTTandProbeTF1->Add(PtProbeTT);
  MassProbeTTandProbeTF2->Add(PtProbeTP);
  cout<<"====================== T&P1 =============================="<<endl;
  int nbEleTopMass1                             = MassProbeTTandProbeTF1->Integral();
  int nbEleBottomMass1                          = MassProbeTTandProbeTF2->Integral();
  float EffChargedHadMass1                      = (float)nbEleTopMass1/nbEleBottomMass1;
  cout<<"nbEleTopMass1                          = "<<nbEleTopMass1<<endl;
  cout<<"nbEleBottomMass1                       = "<<nbEleBottomMass1<<endl;
  cout<<"EffChargedHadMass1(electron)           = "<<EffChargedHadMass1<<endl;
  float StatErrorMass1                          = sqrt(nbEleBottomMass1 - nbEleTopMass1)/nbEleBottomMass1;
  cout<<"StatErrorMass1(electron)               = "<<StatErrorMass1<<endl;
  //---------------------------------------------------
  MassEffForTagProbe539->BayesDivide(MassProbeTTandProbeTF1,MassProbeTTandProbeTF2);
  MassEffForTagProbe539->SetLineStyle(0);
  MassEffForTagProbe539->SetLineColor(1);
  MassEffForTagProbe539->SetLineWidth(2);
  MassEffForTagProbe539->SetMarkerColor(1); 
  MassEffForTagProbe539->SetMarkerStyle(20);
  MassEffForTagProbe539->SetMarkerSize(0.875);
  //-----------------------------------------------------------------
  //========================================================== 
  //                      Mass 720                                                  
  //==========================================================
  TFile *f100 = new TFile("Analyse_TagandProbe_CMSSW8012_MC_reHLT_DYtoMuMu120to200_13TeV_tree.root","READ");  
  //---------------------------------------------------
  TH1F *PtProbeTTb; f1->GetObject("PtProbeTT",PtProbeTTb);
  TH1F *PtProbeTPb; f1->GetObject("PtProbeTP",PtProbeTPb);
  MassProbeTTandProbe720TF1->Add(PtProbeTTb);
  MassProbeTTandProbe720TF2->Add(PtProbeTPb);
  cout<<"====================== T&P1 =============================="<<endl;
  int nbEleTopMass2                             = MassProbeTTandProbe720TF1->Integral();
  int nbEleBottomMass2                          = MassProbeTTandProbe720TF2->Integral();
  float EffChargedHadMass2                      = (float)nbEleTopMass2/nbEleBottomMass2;
  cout<<"nbEleTopMass2                          = "<<nbEleTopMass2<<endl;
  cout<<"nbEleBottomMass2                       = "<<nbEleBottomMass2<<endl;
  cout<<"EffChargedHadMass2(electron)           = "<<EffChargedHadMass2<<endl;
  float StatErrorMass2                          = sqrt(nbEleBottomMass2 - nbEleTopMass2)/nbEleBottomMass2;
  cout<<"StatErrorMass2(electron)               = "<<StatErrorMass2<<endl;
  //---------------------------------------------------
  MassEffForTagProbe720->BayesDivide(MassProbeTTandProbe720TF1,MassProbeTTandProbe720TF2);
  MassEffForTagProbe720->SetLineStyle(0);
  MassEffForTagProbe720->SetLineColor(2);
  MassEffForTagProbe720->SetLineWidth(2);
  MassEffForTagProbe720->SetMarkerColor(2); 
  MassEffForTagProbe720->SetMarkerStyle(24);
  MassEffForTagProbe720->SetMarkerSize(0.875);
  //----------------------------------------------------------------- 
  mg1->Add(MassEffForTagProbe539);
  mg1->Add(MassEffForTagProbe720);
  mg1->Draw("AP");
  mg1->SetTitle("title;muon pt [GeV/c]; Eff.(High pt muons ID)");
  mg1->GetXaxis()->SetTitleOffset(1.7);
  mg1->GetYaxis()->SetTitleOffset(1.7);
  mg1->GetXaxis()->SetLabelSize(0.03);
  mg1->GetYaxis()->SetLabelSize(0.03);
  mg1->GetYaxis()->SetRangeUser(0.7,1.05);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText1 = new TPaveText(0.70, 0.90, 0.90, 0.95, "brNDC");
  tText1->SetBorderSize(0);
  tText1->SetFillColor(0);
  tText1->SetFillStyle(0);
  TText *t1 = tText1->AddText("22.2 pb^{-1} (13 TeV)");
  tText1->SetTextSize(0.035);
  tText1->Draw(); 
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TPaveText* tText2 = new TPaveText(0.2, 0.90, 0.4, 0.95, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t11 = tText2->AddText("CMS Spring16 Simulation");
  tText2->SetTextSize(0.035);
  tText2->Draw(); 
  //==========================================================
  TPaveText* tText11 = new TPaveText(0.50, 0.25, 0.80, 0.45, "brNDC");
  tText11->SetBorderSize(0);
  tText11->SetTextColor(1);
  tText11->SetFillStyle(0);
  Float_t tsize = 0.035;
  tText11->SetTextSize(tsize); 
  TText *t111 = tText1->AddText("Pt > 20 GeV/c"); 
  TText *t2 = tText1->AddText("|#eta|<2.4"); 
  tText11->Draw();
  //========================================================== 
  TLegend *leg = new TLegend(0.20, 0.20, 0.40, 0.35);
  leg->AddEntry(MassEffForTagProbe539,"MC","p");
  leg->AddEntry(MassEffForTagProbe720,"Data (PromptReco)","p");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize22 = 0.035;
  leg->SetTextSize(tsize22); 
  leg->Draw();
  //======================================================================= 
  c1->Print("Eff_TagProbe_ID_MC_13TeV_PtProbe_cmssw80X.png","png");
  c1->Print("Eff_TagProbe_ID_MC_13TeV_PtProbe_cmssw80X.pdf","pdf");
  //c1->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.eps","eps");
  //=======================================================================
  
  //==================================================================   
  //================================================================== 
  //==================================================================  
  //================================================================== 
  TCanvas *c4 = new TCanvas("c4","Eff. PVs",500,500);
  char textpro11[100],textNDF11[100],textRatio11[100];           
  c4->cd(1);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  gPad->SetGridx();
  TH1F* NbPVProbeTTandProbeTF1            = new TH1F("NbPVProbeTTandProbeTF1","",15,0.0,50.0);
  TH1F* NbPVProbeTTandProbeTF2            = new TH1F("NbPVProbeTTandProbeTF2","",15,0.0,50.0);
  TGraphAsymmErrors *TagProbePhi539       = new TGraphAsymmErrors;
  TH1F* NbPVProbeTTandProbe720TF1         = new TH1F("NbPVProbeTTandProbe720TF1","",15,0.0,50.0);
  TH1F* NbPVProbeTTandProbe720TF2         = new TH1F("NbPVProbeTTandProbe720TF2","",15,0.0,50.0);
  TGraphAsymmErrors *TagProbePhi720       = new TGraphAsymmErrors;
  TMultiGraph *mg4 = new TMultiGraph;
  //========================================================== 
  //                     Phi                                                  
  //==========================================================
  TFile *f5 = new TFile("Analyse_TagandProbe_CMSSW8012_MC_reHLT_DYtoMuMu120to200_13TeV_tree.root","READ");
  //---------------------------------------------------
  TH1F *PVProbeTT; f1->GetObject("PVProbeTT",PVProbeTT);
  TH1F *PVProbeTP; f1->GetObject("PVProbeTP",PVProbeTP);
  NbPVProbeTTandProbeTF1->Add(PVProbeTT);
  NbPVProbeTTandProbeTF2->Add(PVProbeTP);
  cout<<"====================== dBeta =============================="<<endl;
  int nbEleTopPV1                             = NbPVProbeTTandProbeTF1->Integral();
  int nbEleBottomPV1                          = NbPVProbeTTandProbeTF2->Integral();
  float EffChargedHadPV1                      = (float)nbEleTopPV1/nbEleBottomPV1;
  cout<<"nbEleTopPV1                          = "<<nbEleTopPV1<<endl;
  cout<<"nbEleBottomPV1                       = "<<nbEleBottomPV1<<endl;
  cout<<"EffChargedHadPV1(electron)           = "<<EffChargedHadPV1<<endl;
  float StatErrorPV1                          = sqrt(nbEleBottomPV1 - nbEleTopPV1)/nbEleBottomPV1;
  cout<<"StatErrorPV1(electron)               = "<<StatErrorPV1<<endl;
  //---------------------------------------------------
  TagProbePhi539->BayesDivide(NbPVProbeTTandProbeTF1,NbPVProbeTTandProbeTF2);
  TagProbePhi539->SetLineStyle(0);
  TagProbePhi539->SetLineColor(1);
  TagProbePhi539->SetLineWidth(2);
  TagProbePhi539->SetMarkerColor(1); 
  TagProbePhi539->SetMarkerStyle(20);
  TagProbePhi539->SetMarkerSize(0.875);
  //========================================================== 
  //                     Phi                                                  
  //==========================================================
  TFile *f50 = new TFile("Analyse_TagandProbe_CMSSW8012_MC_reHLT_DYtoMuMu120to200_13TeV_tree.root","READ");
  //---------------------------------------------------
  TH1F *PVProbeTTb; f1->GetObject("PVProbeTT",PVProbeTTb);
  TH1F *PVProbeTPb; f1->GetObject("PVProbeTP",PVProbeTPb);
  NbPVProbeTTandProbe720TF1->Add(PVProbeTTb);
  NbPVProbeTTandProbe720TF2->Add(PVProbeTPb);
  cout<<"====================== dBeta =============================="<<endl;
  int nbEleTopPV2                             = NbPVProbeTTandProbe720TF1->Integral();
  int nbEleBottomPV2                          = NbPVProbeTTandProbe720TF2->Integral();
  float EffChargedHadPV2                      = (float)nbEleTopPV2/nbEleBottomPV2;
  cout<<"nbEleTopPV2                          = "<<nbEleTopPV2<<endl;
  cout<<"nbEleBottomPV2                       = "<<nbEleBottomPV2<<endl;
  cout<<"EffChargedHadPV2(electron)           = "<<EffChargedHadPV2<<endl;
  float StatErrorPV2                          = sqrt(nbEleBottomPV2 - nbEleTopPV2)/nbEleBottomPV2;
  cout<<"StatErrorPV2(electron)               = "<<StatErrorPV2<<endl;
  //---------------------------------------------------
  TagProbePhi720->BayesDivide(NbPVProbeTTandProbe720TF1,NbPVProbeTTandProbe720TF2);
  TagProbePhi720->SetLineStyle(0);
  TagProbePhi720->SetLineColor(2);
  TagProbePhi720->SetLineWidth(2);
  TagProbePhi720->SetMarkerColor(2); 
  TagProbePhi720->SetMarkerStyle(24);
  TagProbePhi720->SetMarkerSize(0.875);
  //============================================================================
  mg4->Add(TagProbePhi539);
  mg4->Add(TagProbePhi720);
  mg4->Draw("AP");
  mg4->SetTitle("title;Number of vertices; Eff.(High pt muons ID)");
  mg4->GetXaxis()->SetTitleOffset(1.7);
  mg4->GetYaxis()->SetTitleOffset(1.7);
  mg4->GetXaxis()->SetLabelSize(0.03);
  mg4->GetYaxis()->SetLabelSize(0.03);
  mg4->GetXaxis()->SetRangeUser(0.0,37.0);
  mg4->GetYaxis()->SetRangeUser(0.7,1.05);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText1111 = new TPaveText(0.70, 0.90, 0.90, 0.95, "brNDC");
  tText1111->SetBorderSize(0);
  tText1111->SetFillColor(0);
  tText1111->SetFillStyle(0);
  TText *t1a = tText1111->AddText("22.2 pb^{-1} (13 TeV)");
  tText1111->SetTextSize(0.035);
  tText1111->Draw(); 
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TPaveText* tText22 = new TPaveText(0.2, 0.90, 0.4, 0.95, "brNDC");
  tText22->SetBorderSize(0);
  tText22->SetFillColor(0);
  tText22->SetFillStyle(0);
  TText *t1b = tText22->AddText("CMS Spring16 Simulation");
  tText22->SetTextSize(0.035);
  tText22->Draw(); 
  //========================================================== 
  TPaveText* tText1b = new TPaveText(0.50, 0.25, 0.80, 0.45, "brNDC");
  tText1b->SetBorderSize(0);
  tText1b->SetTextColor(1);
  tText1b->SetFillStyle(0);
  Float_t tsizee = 0.035;
  tText1b->SetTextSize(tsizee); 
  TText *t11b = tText1b->AddText("Pt > 20 GeV/c"); 
  TText *t22b = tText1b->AddText("|#eta|<2.4"); 
  tText1b->Draw();
  //========================================================== 
  TLegend *leg0 = new TLegend(0.20, 0.20, 0.40, 0.35);
  leg0->AddEntry(TagProbePhi539,"MC","p");
  leg0->AddEntry(TagProbePhi720,"Data (PromptReco)","p");
  leg0->SetBorderSize(0.0);
  leg0->SetMargin(0.3);
  leg0->SetFillColor(0);
  leg0->SetFillStyle(10);
  leg0->SetLineColor(0);
  Float_t tsize2a = 0.035;
  leg0->SetTextSize(tsize2a); 
  leg0->Draw();
  //======================================================================= 
  c4->Print("Eff_TagProbe_ID_MC_13TeV_vtxProbe_cmssw80X.png","png");
  c4->Print("Eff_TagProbe_ID_MC_13TeV_vtxProbe_cmssw80X.pdf","pdf");
  //c4->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Phi_EBEE.eps","eps");
  //=======================================================================

  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================
  //TCanvas *c5 = new TCanvas("c5","A Simple Graph with error bars",200,10,700,500);
  TCanvas *c5 = new TCanvas("c5","Eff. Eta",500,500);
  c5->SetFillColor(0);
  c5->SetGrid();
  c5->GetFrame()->SetFillColor(21);
  c5->GetFrame()->SetBorderSize(12);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  gPad->SetGridx();
  const Int_t n = 13;
  Float_t x[n]  = {-2.25,-1.85,-1.4,-1.05,-0.6,-0.25,0.0,0.25,0.6,1.05,1.4,1.85,2.25};
  Float_t y[n]  = {Eff1,
		   Eff2,
		   Eff3,
		   Eff4,
		   Eff5,
		   Eff6,
		   Eff7,
		   Eff8,
		   Eff9,
		   Eff10,
		   Eff11,
		   Eff12,
		   Eff13};
  Float_t ex[n] = {0.15,0.25,0.2,0.15,0.3,0.05,0.2,0.05,0.3,0.15,0.2,0.25,0.15};
  Float_t ey[n] = {error1,error2,error3,error4,error5,error6,error7,error8,error9,error10,error11,error12,error13};
  TGraphErrors *gr = new TGraphErrors(n,x,y,ex,ey);
  gr->SetMarkerColor(1);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.9);
  gr->Draw("AP");
  gr->GetYaxis()->SetNdivisions(10,kTRUE);
  gr->GetYaxis()->SetDecimals();
  //gr->GetYaxis()->SetRangeUser(0.0,0.4);
  //gr->GetXaxis()->SetNdivisions(510,kFALSE);
  gr->GetXaxis()->SetRangeUser(-2.4,2.4);
  gr->GetYaxis()->SetRangeUser(0.70,1.05);
  gr->SetTitle("title;muon #eta; Eff.(High pt muons ID)");
  gr->GetXaxis()->SetTitleOffset(1.7);
  gr->GetYaxis()->SetTitleOffset(1.7);
  gr->GetXaxis()->SetLabelSize(0.03);
  gr->GetYaxis()->SetLabelSize(0.03);
  gr->GetXaxis()->SetRangeUser(-2.1,2.1);
  Float_t x2[n]  = {-2.25,-1.85,-1.4,-1.05,-0.6,-0.25,0.0,0.25,0.6,1.05,1.4,1.85,2.25};
  Float_t y2[n]  = {EffD1,
		   EffD2,
		   EffD3,
		   EffD4,
		   EffD5,
		   EffD6,
		   EffD7,
		   EffD8,
		   EffD9,
		   EffD10,
		   EffD11,
		   EffD12,
		   EffD13};
  Float_t ex2[n] = {0.15,0.25,0.2,0.15,0.3,0.05,0.2,0.05,0.3,0.15,0.2,0.25,0.15};
  Float_t ey2[n] = {errorD1,errorD2,errorD3,errorD4,errorD5,errorD6,errorD7,errorD8,errorD9,errorD10,errorD11,errorD12,errorD13};
  TGraphErrors *gr2 = new TGraphErrors(n,x2,y2,ex2,ey2);
  gr2->SetMarkerColor(2);
  gr2->SetLineColor(2);
  gr2->SetMarkerStyle(24);
  gr2->SetMarkerSize(0.9);
  gr2->GetXaxis()->SetRangeUser(-2.1,2.1);
  gr2->Draw("p,sames");
  //-----------------------------------------------------------------
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  double LikelihoodCut1    = 1.2 ;
  TLine *line1 = new TLine( LikelihoodCut1 , 0.7,  LikelihoodCut1 , 1.05);
  line1->SetLineColor(2);
  line1->SetLineWidth(3);
  line1->SetLineStyle(2);
  line1->Draw();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  double LikelihoodCut2    = -1.2 ;
  TLine *line2 = new TLine( LikelihoodCut2 , 0.7,  LikelihoodCut2 , 1.05);
  line2->SetLineColor(2);
  line2->SetLineWidth(3);
  line2->SetLineStyle(2);
  line2->Draw();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText1c = new TPaveText(0.70, 0.90, 0.90, 0.95, "brNDC");
  tText1c->SetBorderSize(0);
  tText1c->SetFillColor(0);
  tText1c->SetFillStyle(0);
  TText *t11c = tText1c->AddText("22.2 pb^{-1} (13 TeV)");
  tText1c->SetTextSize(0.035);
  tText1c->Draw(); 
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TPaveText* tText2c = new TPaveText(0.2, 0.90, 0.4, 0.95, "brNDC");
  tText2c->SetBorderSize(0);
  tText2c->SetFillColor(0);
  tText2c->SetFillStyle(0);
  TText *t11d = tText2c->AddText("CMS Spring16 Simulation");
  tText2c->SetTextSize(0.035);
  tText2c->Draw(); 
  //========================================================== 
  TPaveText* tText11b = new TPaveText(0.50, 0.25, 0.80, 0.45, "brNDC");
  tText11b->SetBorderSize(0);
  tText11b->SetTextColor(1);
  tText11b->SetFillStyle(0);
  Float_t tsizeee = 0.035;
  tText11b->SetTextSize(tsizeee); 
  TText *t1bc = tText11b->AddText("Pt > 20 GeV/c"); 
  TText *t2bc = tText11b->AddText("|#eta|<2.4"); 
  tText11b->Draw();
  //==========================================================
  TLegend *leg1 = new TLegend(0.20, 0.20, 0.40, 0.35);
  leg1->AddEntry(gr,"MC","p");
  leg1->AddEntry(gr2,"Data (PromptReco)","p");
  leg1->SetBorderSize(0.0);
  leg1->SetMargin(0.3);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(10);
  leg1->SetLineColor(0);
  Float_t tsize2 = 0.03;
  leg1->SetTextSize(tsize2); 
  leg1->Draw();  
  //======================================================================= 
  c5->Print("Eff_TagProbe_ID_MC_13TeV_EtaProbe_cmssw80X.png","png");
  c5->Print("Eff_TagProbe_ID_MC_13TeV_EtaProbe_cmssw80X.pdf","pdf");
  c5->Update();
  //================================================================== 
  

}









