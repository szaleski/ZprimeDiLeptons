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
void TagProbeMuonPlotsTrackIso(){
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
  //------------------------------------------------------------------
  float Eff1    = 0.978225;
  float Eff2    = 0.978065;
  float Eff3    = 0.976409;
  float Eff4    = 0.979615;
  float Eff5    = 0.980686;
  float Eff6    = 0.977283;
  float Eff7    = 0.980582;
  float Eff8    = 0.978432;
  float Eff9    = 0.98169;
  float Eff10   = 0.980613;
  float Eff11   = 0.97379;
  float Eff12   = 0.975439;
  float Eff13   = 0.978401;

  float error1     = 0.00210502;
  float error2     = 0.0012664;
  float error3     = 0.00127943;
  float error4     = 0.00113777;
  float error5     = 0.000752352;
  float error6     = 0.00221699;
  float error7     = 0.000954362;
  float error8     = 0.0021149;
  float error9     = 0.000725463;
  float error10    = 0.00111376;
  float error11    = 0.00134457;
  float error12    = 0.00132619;
  float error13    = 0.00208803;
  /*
  [0.0,0.9]   = 98.1 (+/-) 0.04
  [0.9,1.2]   = 98.0 (+/-) 0.08
  [1.2,2.1]   = 97.6 (+/-) 0.07
  [2.1,2.4]   = 97.8 (+/-) 0.15
    
  error[0.0,0.9]    = 0.000439492
  error[0.9,1.2]    = 0.00079615
  error[1.2,2.1]    = 0.00065283
  error[2.1,2.4]    = 0.00148244
  */  
    
  //==============================================
  //  Data
  //==============================================
  float EffD1    = 0.966102;
  float EffD2    = 0.969136;
  float EffD3    = 0.972376;
  float EffD4    = 0.980392;
  float EffD5    = 0.983516;
  float EffD6    = 0.978261;
  float EffD7    = 0.979339;
  float EffD8    = 1;
  float EffD9    = 0.980769;
  float EffD10   = 0.974843;
  float EffD11   = 1;
  float EffD12   = 0.969925;
  float EffD13   = 1;
  
  float errorD1     = 0.0239697;
  float errorD2     = 0.0138029;
  float errorD3     = 0.012354;
  float errorD4     = 0.0113206;
  float errorD5     = 0.00672937;
  float errorD6     = 0.0217391;
  float errorD7     = 0.00923995;
  float errorD8     = 0;
  float errorD9     = 0.00726855;
  float errorD10    = 0.0125786;
  float errorD11    = 0;
  float errorD12    = 0.0150376;
  float errorD13    = 0;
  /*
  [0.0,0.9]   = 98.2 (+/-) 0.4
  [0.9,1.2]   = 97.7 (+/-) 0.8
  [1.2,2.1]   = 97.8 (+/-) 0.6
  [2.1,2.4]   = 98.1 (+/-) 1.4
  error[0.0,0.9]    = 0.00412775
  error[0.9,1.2]    = 0.00847997
  error[1.2,2.1]    = 0.00587387
  error[2.1,2.4]    = 0.0135982
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
  TFile *f1 = new TFile("Analyse_TagandProbe_MC_DY_60_Mass_120_RelTrackIso_ID_pt45.root","READ");
  //---------------------------------------------------
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
  TFile *f100 = new TFile("Analyse_TagandProbe_Data_DY_60_Mass_120_RelTrackIso_ID_pt45.root","READ");
  //---------------------------------------------------
  MassProbeTTandProbe720TF1->Add(PtProbeTT);
  MassProbeTTandProbe720TF2->Add(PtProbeTP);
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
  mg1->SetTitle("title;muon pt [GeV/c]; Eff.(Rel. Tracker iso)");
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
  TText *t1 = tText1->AddText("50.421 pb^{-1} (13 TeV)");
  tText1->SetTextSize(0.035);
  tText1->Draw(); 
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TPaveText* tText2 = new TPaveText(0.2, 0.90, 0.4, 0.95, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("CMS Spring15 Simulation");
  tText2->SetTextSize(0.035);
  tText2->Draw(); 
  //==========================================================
  TPaveText* tText1 = new TPaveText(0.50, 0.25, 0.80, 0.45, "brNDC");
  tText1->SetBorderSize(0);
  tText1->SetTextColor(1);
  tText1->SetFillStyle(0);
  Float_t tsize = 0.035;
  tText1->SetTextSize(tsize); 
  TText *t1 = tText1->AddText("Pt > 45 GeV/c"); 
  TText *t2 = tText1->AddText("|#eta|<2.4"); 
  tText1->Draw();
  //========================================================== 
  TLegend *leg = new TLegend(0.20, 0.20, 0.40, 0.35);
  leg->AddEntry(MassEffForTagProbe539,"MC","p");
  leg->AddEntry(MassEffForTagProbe720,"Data","p");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.035;
  leg->SetTextSize(tsize2); 
  leg->Draw();
  //======================================================================= 
  c1->Print("Eff_TagProbe_ID_MC_13TeV_PtProbe_cmssw720-trackiso.png","png");
  c1->Print("Eff_TagProbe_ID_MC_13TeV_PtProbe_cmssw720-trackiso.pdf","pdf");
  //c1->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.eps","eps");
  //=======================================================================
  
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  TCanvas *c4 = new TCanvas("c4","Eff. PVs",500,500);
  char textpro1[100],textNDF1[100],textRatio1[100];           
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
  TFile *f5 = new TFile("Analyse_TagandProbe_MC_DY_60_Mass_120_RelTrackIso_ID_pt45.root","READ");
  //---------------------------------------------------
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
  TFile *f50 = new TFile("Analyse_TagandProbe_Data_DY_60_Mass_120_RelTrackIso_ID_pt45.root","READ");
  //---------------------------------------------------
  NbPVProbeTTandProbe720TF1->Add(PVProbeTT);
  NbPVProbeTTandProbe720TF2->Add(PVProbeTP);
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
  mg4->SetTitle("title;Number of vertices; Eff.(Rel. Tracker iso)");
  mg4->GetXaxis()->SetTitleOffset(1.7);
  mg4->GetYaxis()->SetTitleOffset(1.7);
  mg4->GetXaxis()->SetLabelSize(0.03);
  mg4->GetYaxis()->SetLabelSize(0.03);
  mg4->GetXaxis()->SetRangeUser(0.0,37.0);
  mg4->GetYaxis()->SetRangeUser(0.7,1.05);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText1 = new TPaveText(0.70, 0.90, 0.90, 0.95, "brNDC");
  tText1->SetBorderSize(0);
  tText1->SetFillColor(0);
  tText1->SetFillStyle(0);
  TText *t1 = tText1->AddText("50.421 pb^{-1} (13 TeV)");
  tText1->SetTextSize(0.035);
  tText1->Draw(); 
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TPaveText* tText2 = new TPaveText(0.2, 0.90, 0.4, 0.95, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("CMS Spring15 Simulation");
  tText2->SetTextSize(0.035);
  tText2->Draw(); 
  //========================================================== 
  TPaveText* tText1 = new TPaveText(0.50, 0.25, 0.80, 0.45, "brNDC");
  tText1->SetBorderSize(0);
  tText1->SetTextColor(1);
  tText1->SetFillStyle(0);
  Float_t tsize = 0.035;
  tText1->SetTextSize(tsize); 
  TText *t1 = tText1->AddText("Pt > 45 GeV/c"); 
  TText *t2 = tText1->AddText("|#eta|<2.4"); 
  tText1->Draw();
  //========================================================== 
  TLegend *leg = new TLegend(0.20, 0.20, 0.40, 0.35);
  leg->AddEntry(TagProbePhi539,"MC","p");
  leg->AddEntry(TagProbePhi720,"Data","p");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.035;
  leg->SetTextSize(tsize2); 
  leg->Draw();
  //======================================================================= 
  c4->Print("Eff_TagProbe_ID_MC_13TeV_vtxProbe_cmssw720-trackiso.png","png");
  c4->Print("Eff_TagProbe_ID_MC_13TeV_vtxProbe_cmssw720-trackiso.pdf","pdf");
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
  gr->SetTitle("title;muon #eta; Eff.(Rel. Tracker iso)");
  gr->GetXaxis()->SetTitleOffset(1.7);
  gr->GetYaxis()->SetTitleOffset(1.7);
  gr->GetXaxis()->SetLabelSize(0.03);
  gr->GetYaxis()->SetLabelSize(0.03);
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
   TPaveText* tText1 = new TPaveText(0.70, 0.90, 0.90, 0.95, "brNDC");
  tText1->SetBorderSize(0);
  tText1->SetFillColor(0);
  tText1->SetFillStyle(0);
  TText *t1 = tText1->AddText("50.421 pb^{-1} (13 TeV)");
  tText1->SetTextSize(0.035);
  tText1->Draw(); 
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TPaveText* tText2 = new TPaveText(0.2, 0.90, 0.4, 0.95, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("CMS Spring15 Simulation");
  tText2->SetTextSize(0.035);
  tText2->Draw(); 
  //========================================================== 
  TPaveText* tText1 = new TPaveText(0.50, 0.25, 0.80, 0.45, "brNDC");
  tText1->SetBorderSize(0);
  tText1->SetTextColor(1);
  tText1->SetFillStyle(0);
  Float_t tsize = 0.035;
  tText1->SetTextSize(tsize); 
  TText *t1 = tText1->AddText("Pt > 45 GeV/c"); 
  TText *t2 = tText1->AddText("|#eta|<2.4"); 
  tText1->Draw();
  //==========================================================
  TLegend *leg = new TLegend(0.20, 0.20, 0.40, 0.35);
  leg->AddEntry(gr,"MC","p");
  leg->AddEntry(gr2,"Data","p");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.03;
  leg->SetTextSize(tsize2); 
  leg->Draw();  
  //======================================================================= 
  c5->Print("Eff_TagProbe_ID_MC_13TeV_EtaProbe_cmssw720-trackiso.png","png");
  c5->Print("Eff_TagProbe_ID_MC_13TeV_EtaProbe_cmssw720-trackiso.pdf","pdf");
  c5->Update();
  //================================================================== 
  

}









