#include "TProfile.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TF2.h"
#include <string>
using std::string;
void MuonsEffEta(){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  //gStyle->SetStatFormat("4.2f");
  //gStyle->SetFitFormat("4.2f");
  int BoxValue = 1110; //6610; //4680;  
  gStyle->SetOptFit(11);
  gStyle->SetOptDate(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(kFALSE);
  //gStyle->SetOptStat(111110);
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
  double norm = 1.0;
  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  TCanvas *c1 = new TCanvas("c1","Et in EB",600,600);
  char textpro1[100],textNDF1[100],textRatio1[100];           
  c1->cd(1);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  gPad->SetGridx();
  TGraphAsymmErrors *EleEtDeltaPhiEB     = new TGraphAsymmErrors;
  TGraphAsymmErrors *PhotonEtDeltaPhiEB  = new TGraphAsymmErrors;
  TGraphAsymmErrors *QCDEtDeltaPhiEB     = new TGraphAsymmErrors;
  TMultiGraph *mg1 = new TMultiGraph;
  //========================================================== 
  //                      DeltaPhi in EB                                                  
  //==========================================================
  TFile *f1 = new TFile("/lustre/cms/store/user/defilip/ZprimeAnalysis/80X/histos/histosZprimeMuMu/output_ZprimeToMuMu_M-5000_TuneCUETP8M1_13TeV-pythia8_tree.root","READ");
  //---------------------------------------------------
  cout<<"====================== error on pt =============================="<<endl;
  TH1F *EtaID; f1->GetObject("EtaID",EtaID);
  TH1F *EtaEffpterror; f1->GetObject("EtaEffpterror",EtaEffpterror);
  int nbEleTopEB1                             = EtaID->Integral();
  int nbEleBottomEB1                          = EtaEffpterror->Integral();
  float EffChargedHadEB1                      = (float)nbEleTopEB1/nbEleBottomEB1;
  cout<<"nbEleTopEB1                          = "<<nbEleTopEB1<<endl;
  cout<<"nbEleBottomEB1                       = "<<nbEleBottomEB1<<endl;
  cout<<"EffChargedHadEB1(electron)           = "<<EffChargedHadEB1<<endl;
  float StatErrorEB1                          = sqrt(nbEleBottomEB1 - nbEleTopEB1)/nbEleBottomEB1;
  cout<<"StatErrorEB1(electron)               = "<<StatErrorEB1<<endl;
  //---------------------------------------------------
  EleEtDeltaPhiEB->BayesDivide(EtaID,EtaEffpterror);
  EleEtDeltaPhiEB->SetLineStyle(0);
  EleEtDeltaPhiEB->SetLineColor(1);
  EleEtDeltaPhiEB->SetLineWidth(2);
  EleEtDeltaPhiEB->SetMarkerColor(1); 
  EleEtDeltaPhiEB->SetMarkerStyle(20);
  EleEtDeltaPhiEB->SetMarkerSize(0.875);
  /*
  //========================================================== 
  //                 DeltaPhi in EB                                                  
  //==========================================================
  TFile *f2 = new TFile("CMSSW701-Analyse_MC_GJet_pt_3000_V5_withPF.root","READ");
  PhotonEtDeltaPhiEB->BayesDivide(EtDeltaPhiEB2,EtDeltaPhiEB1);
  //---------------------------------------------------
  cout<<"====================== Gamma Jet EB =============================="<<endl;
  int nbEleTopEB2                             = EtDeltaPhiEB2->Integral();
  int nbEleBottomEB2                          = EtDeltaPhiEB1->Integral();
  float EffChargedHadEB3                      = (float)nbEleTopEB2/nbEleBottomEB2;
  cout<<"nbEleTopEB2                          = "<<nbEleTopEB2<<endl;
  cout<<"nbEleBottomEB2                       = "<<nbEleBottomEB2<<endl;
  cout<<"EffChargedHadEB3(photon)             = "<<EffChargedHadEB3<<endl;
  float StatErrorEB2                          = sqrt(nbEleBottomEB2 - nbEleTopEB2)/nbEleBottomEB2;
  cout<<"StatErrorEB2(photon)                 = "<<StatErrorEB2<<endl;
  //---------------------------------------------------
  PhotonEtDeltaPhiEB->SetLineStyle(0);
  PhotonEtDeltaPhiEB->SetLineColor(2);
  PhotonEtDeltaPhiEB->SetLineWidth(2);
  PhotonEtDeltaPhiEB->SetMarkerColor(2); 
  PhotonEtDeltaPhiEB->SetMarkerStyle(24);
  PhotonEtDeltaPhiEB->SetMarkerSize(0.875);
  //========================================================== 
  //                 QCD pfChargedHadron in EB                                                  
  //==========================================================
  TFile *f3 = new TFile("CMSSW701-Analyse_MC_QCD_15_pt_3000_V5_withPF.root","READ");
  QCDEtDeltaPhiEB->BayesDivide(EtDeltaPhiEB2,EtDeltaPhiEB1);
  //---------------------------------------------------
  cout<<"====================== QCD EB =============================="<<endl;
  int nbEleTopEB3                             = EtDeltaPhiEB2->Integral();
  int nbEleBottomEB3                          = EtDeltaPhiEB1->Integral();
  float EffChargedHadEB3                      = (float)nbEleTopEB3/nbEleBottomEB3;
  cout<<"nbEleTopEB3                          = "<<nbEleTopEB3<<endl;
  cout<<"nbEleBottomEB3                       = "<<nbEleBottomEB3<<endl;
  cout<<"EffChargedHadEB3(QCD)                = "<<EffChargedHadEB3<<endl;
  float StatErrorEB3                          = sqrt(nbEleBottomEB3 - nbEleTopEB3)/nbEleBottomEB3;
  cout<<"StatErrorEB3(QCD)                    = "<<StatErrorEB3<<endl;
  //---------------------------------------------------
  QCDEtDeltaPhiEB->SetLineStyle(0);
  QCDEtDeltaPhiEB->SetLineColor(4);
  QCDEtDeltaPhiEB->SetLineWidth(2);
  QCDEtDeltaPhiEB->SetMarkerColor(4); 
  QCDEtDeltaPhiEB->SetMarkerStyle(22);
  QCDEtDeltaPhiEB->SetMarkerSize(0.875);
  //-----------------------------------------------------------------
  */                                    
  mg1->Add(EleEtDeltaPhiEB);
  //mg1->Add(PhotonEtDeltaPhiEB);
  //mg1->Add(QCDEtDeltaPhiEB);
  mg1->Draw("AP");
  mg1->SetTitle("title;|#eta_{muons}|; Efficiency(N-1)");
  mg1->GetXaxis()->SetTitleOffset(1.7);
  mg1->GetYaxis()->SetTitleOffset(1.7);
  //mg1->GetXaxis()->SetRangeUser(0.0,800.0);
  //mg1->GetYaxis()->SetRangeUser(0.0,1.05);
  mg1->GetYaxis()->SetRangeUser(0.9,1.01);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText2 = new TPaveText(0.40, 0.90, 0.90, 0.92, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t2 = tText2->AddText("DeltaPhi");
  TText *t3 = tText2->AddText("#delta pt/pt<0.3; CMSSW701 [1/fb]"); 
  tText2->Draw();
  //========================================================== 
  //TPaveText* tText2 = new TPaveText(0.20, 0.30, 0.60, 0.45, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  //Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t1 = tText2->AddText("Eff(N-1) = 0.9969 #pm 0.0002");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9684 #pm 0.0005");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8921 #pm 0.0057");

  TText *t1 = tText2->AddText("Eff(N-1) = 0.9999 #pm 3.09149e-05");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9198 #pm 0.0008");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8316 #pm 0.0071");
  tText2->Draw();
  //========================================================== 
  TLegend *leg = new TLegend(0.60, 0.15, 0.80, 0.30);
  leg->AddEntry(EleEtDeltaPhiEB,"Z' [M = 1 TeV/c^{2}]","p");
  //leg->AddEntry(PhotonEtDeltaPhiEB,"#gamma + jets (EB)","p");
  //leg->AddEntry(QCDEtDeltaPhiEB,"QCD (EB)","p");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.03;
  leg->SetTextSize(tsize2); 
  leg->Draw();
  //======================================================================= 
  c1->Print("Eff_DeltaPt_HEEP_ID_MC_CMSSW701_Eta.png","png");
  //c1->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.pdf","pdf");
  //c1->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.eps","eps");
  //=======================================================================
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  
  TCanvas *c2 = new TCanvas("c2","Et in EB",600,600);
  char textpro11[100],textNDF11[100],textRatio11[100];           
  c2->cd(1);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  gPad->SetGridx();

  EleEtDeltaPhiEB     = new TGraphAsymmErrors;
  PhotonEtDeltaPhiEB  = new TGraphAsymmErrors;
  QCDEtDeltaPhiEB     = new TGraphAsymmErrors;
  TMultiGraph *mg2 = new TMultiGraph;
  //========================================================== 
  //                      DeltaPhi in EB                                                  
  //==========================================================
  TFile *f0 = new TFile("/lustre/cms/store/user/defilip/ZprimeAnalysis/80X/histos/histosZprimeMuMu/output_ZprimeToMuMu_M-5000_TuneCUETP8M1_13TeV-pythia8_tree.root","READ");
  //---------------------------------------------------
  cout<<"====================== NumberOftrackerLayers =============================="<<endl;
  //TH1F *EtaID; 
  f0->GetObject("EtaID",EtaID);
  TH1F *EtaEffptnumberOftrackerLayers; f0->GetObject("EtaEffptnumberOftrackerLayers",EtaEffptnumberOftrackerLayers);
  int nbEleTopEB2                             = EtaID->Integral();
  int nbEleBottomEB2                          = EtaEffptnumberOftrackerLayers->Integral();
  float EffChargedHadEB2                      = (float)nbEleTopEB2/nbEleBottomEB2;
  cout<<"nbEleTopEB2                          = "<<nbEleTopEB2<<endl;
  cout<<"nbEleBottomEB2                       = "<<nbEleBottomEB2<<endl;
  cout<<"EffChargedHadEB2(electron)           = "<<EffChargedHadEB2<<endl;
  float StatErrorEB2                          = sqrt(nbEleBottomEB2 - nbEleTopEB2)/nbEleBottomEB2;
  cout<<"StatErrorEB2(electron)               = "<<StatErrorEB2<<endl;
  //---------------------------------------------------
  EleEtDeltaPhiEB->BayesDivide(EtaID,EtaEffptnumberOftrackerLayers);
  EleEtDeltaPhiEB->SetLineStyle(0);
  EleEtDeltaPhiEB->SetLineColor(1);
  EleEtDeltaPhiEB->SetLineWidth(2);
  EleEtDeltaPhiEB->SetMarkerColor(1); 
  EleEtDeltaPhiEB->SetMarkerStyle(20);
  EleEtDeltaPhiEB->SetMarkerSize(0.875);
  /*
  //========================================================== 
  //                 DeltaPhi in EB                                                  
  //==========================================================
  TFile *f2 = new TFile("CMSSW701-Analyse_MC_GJet_pt_3000_V5_withPF.root","READ");
  PhotonEtDeltaPhiEB->BayesDivide(EtDeltaPhiEB2,EtDeltaPhiEB1);
  //---------------------------------------------------
  cout<<"====================== Gamma Jet EB =============================="<<endl;
  int nbEleTopEB2                             = EtDeltaPhiEB2->Integral();
  int nbEleBottomEB2                          = EtDeltaPhiEB1->Integral();
  float EffChargedHadEB3                      = (float)nbEleTopEB2/nbEleBottomEB2;
  cout<<"nbEleTopEB2                          = "<<nbEleTopEB2<<endl;
  cout<<"nbEleBottomEB2                       = "<<nbEleBottomEB2<<endl;
  cout<<"EffChargedHadEB3(photon)             = "<<EffChargedHadEB3<<endl;
  float StatErrorEB2                          = sqrt(nbEleBottomEB2 - nbEleTopEB2)/nbEleBottomEB2;
  cout<<"StatErrorEB2(photon)                 = "<<StatErrorEB2<<endl;
  //---------------------------------------------------
  PhotonEtDeltaPhiEB->SetLineStyle(0);
  PhotonEtDeltaPhiEB->SetLineColor(2);
  PhotonEtDeltaPhiEB->SetLineWidth(2);
  PhotonEtDeltaPhiEB->SetMarkerColor(2); 
  PhotonEtDeltaPhiEB->SetMarkerStyle(24);
  PhotonEtDeltaPhiEB->SetMarkerSize(0.875);
  //========================================================== 
  //                 QCD pfChargedHadron in EB                                                  
  //==========================================================
  TFile *f3 = new TFile("CMSSW701-Analyse_MC_QCD_15_pt_3000_V5_withPF.root","READ");
  QCDEtDeltaPhiEB->BayesDivide(EtDeltaPhiEB2,EtDeltaPhiEB1);
  //---------------------------------------------------
  cout<<"====================== QCD EB =============================="<<endl;
  int nbEleTopEB3                             = EtDeltaPhiEB2->Integral();
  int nbEleBottomEB3                          = EtDeltaPhiEB1->Integral();
  float EffChargedHadEB3                      = (float)nbEleTopEB3/nbEleBottomEB3;
  cout<<"nbEleTopEB3                          = "<<nbEleTopEB3<<endl;
  cout<<"nbEleBottomEB3                       = "<<nbEleBottomEB3<<endl;
  cout<<"EffChargedHadEB3(QCD)                = "<<EffChargedHadEB3<<endl;
  float StatErrorEB3                          = sqrt(nbEleBottomEB3 - nbEleTopEB3)/nbEleBottomEB3;
  cout<<"StatErrorEB3(QCD)                    = "<<StatErrorEB3<<endl;
  //---------------------------------------------------
  QCDEtDeltaPhiEB->SetLineStyle(0);
  QCDEtDeltaPhiEB->SetLineColor(4);
  QCDEtDeltaPhiEB->SetLineWidth(2);
  QCDEtDeltaPhiEB->SetMarkerColor(4); 
  QCDEtDeltaPhiEB->SetMarkerStyle(22);
  QCDEtDeltaPhiEB->SetMarkerSize(0.875);
  //-----------------------------------------------------------------
  */
  mg2->Add(EleEtDeltaPhiEB);
  //mg2->Add(PhotonEtDeltaPhiEB);
  //mg2->Add(QCDEtDeltaPhiEB);
  mg2->Draw("AP");
  mg2->SetTitle("title;|#eta_{muons}|; Efficiency(N-1)");
  mg2->GetXaxis()->SetTitleOffset(1.7);
  mg2->GetYaxis()->SetTitleOffset(1.7);
  //mg2->GetXaxis()->SetRangeUser(0.0,800.0);
  //mg2->GetYaxis()->SetRangeUser(0.0,1.05);
  mg2->GetYaxis()->SetRangeUser(0.9,1.01);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  //TPaveText* tText2 = new TPaveText(0.40, 0.90, 0.90, 0.92, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  //Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t2 = tText2->AddText("DeltaPhi");
  //TText *t3 = tText2->AddText("NbTrackLayers>5.0; CMSSW701 [1/fb]"); 
  tText2->Draw();
  //========================================================== 
  //TPaveText* tText2 = new TPaveText(0.20, 0.30, 0.60, 0.45, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  //Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t1 = tText2->AddText("Eff(N-1) = 0.9969 #pm 0.0002");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9684 #pm 0.0005");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8921 #pm 0.0057");

  //TText *t1 = tText2->AddText("Eff(N-1) = 0.9998 #pm 3.82788e-05");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9198 #pm 0.0008");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8316 #pm 0.0071");
  tText2->Draw();
  //========================================================== 
  TLegend *leg1 = new TLegend(0.60, 0.15, 0.80, 0.30);
  leg1->AddEntry(EleEtDeltaPhiEB,"Z' [M = 1 TeV/c^{2}]","p");
  //leg1->AddEntry(PhotonEtDeltaPhiEB,"#gamma + jets (EB)","p");
  //leg1->AddEntry(QCDEtDeltaPhiEB,"QCD (EB)","p");
  leg1->SetBorderSize(0.0);
  leg1->SetMargin(0.3);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(10);
  leg1->SetLineColor(0);
  //Float_t tsize2 = 0.03;
  leg1->SetTextSize(tsize2); 
  leg1->Draw();
  //======================================================================= 
  c2->Print("Eff_numberOftrackerLayers_HEEP_ID_MC_CMSSW701_Eta.png","png");
  //c2->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.pdf","pdf");
  //c2->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.eps","eps");
  //=======================================================================
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  TCanvas *c3 = new TCanvas("c3","Et in EB",600,600);
  char textpro111[100],textNDF111[100],textRatio111[100];           
  c3->cd(1);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  gPad->SetGridx();
  EleEtDeltaPhiEB     = new TGraphAsymmErrors;
  PhotonEtDeltaPhiEB  = new TGraphAsymmErrors;
  QCDEtDeltaPhiEB     = new TGraphAsymmErrors;
  TMultiGraph *mg3 = new TMultiGraph;
  //========================================================== 
  //                      DeltaPhi in EB                                                  
  //==========================================================
  TFile *f11 = new TFile("/lustre/cms/store/user/defilip/ZprimeAnalysis/80X/histos/histosZprimeMuMu/output_ZprimeToMuMu_M-5000_TuneCUETP8M1_13TeV-pythia8_tree.root","READ");
  //---------------------------------------------------
  //TH1F *EtaID; 
  f11->GetObject("EtaID",EtaID);
  TH1F *EtaEffptnumberOfPixelHits; f11->GetObject("EtaEffptnumberOfPixelHits",EtaEffptnumberOfPixelHits);
 
  cout<<"====================== NumberOfPixelHits =============================="<<endl;
  int nbEleTopEB3                             = EtaID->Integral();
  int nbEleBottomEB3                          = EtaEffptnumberOfPixelHits->Integral();
  float EffChargedHadEB3                      = (float)nbEleTopEB3/nbEleBottomEB3;
  cout<<"nbEleTopEB3                          = "<<nbEleTopEB3<<endl;
  cout<<"nbEleBottomEB3                       = "<<nbEleBottomEB3<<endl;
  cout<<"EffChargedHadEB3(electron)           = "<<EffChargedHadEB3<<endl;
  float StatErrorEB3                          = sqrt(nbEleBottomEB3 - nbEleTopEB3)/nbEleBottomEB3;
  cout<<"StatErrorEB3(electron)               = "<<StatErrorEB3<<endl;
  //---------------------------------------------------
  EleEtDeltaPhiEB->BayesDivide(EtaID,EtaEffptnumberOfPixelHits);
  EleEtDeltaPhiEB->SetLineStyle(0);
  EleEtDeltaPhiEB->SetLineColor(1);
  EleEtDeltaPhiEB->SetLineWidth(2);
  EleEtDeltaPhiEB->SetMarkerColor(1); 
  EleEtDeltaPhiEB->SetMarkerStyle(20);
  EleEtDeltaPhiEB->SetMarkerSize(0.875);
  /*
  //========================================================== 
  //                 DeltaPhi in EB                                                  
  //==========================================================
  TFile *f2 = new TFile("CMSSW701-Analyse_MC_GJet_pt_3000_V5_withPF.root","READ");
  PhotonEtDeltaPhiEB->BayesDivide(EtDeltaPhiEB2,EtDeltaPhiEB1);
  //---------------------------------------------------
  cout<<"====================== Gamma Jet EB =============================="<<endl;
  int nbEleTopEB2                             = EtDeltaPhiEB2->Integral();
  int nbEleBottomEB2                          = EtDeltaPhiEB1->Integral();
  float EffChargedHadEB3                      = (float)nbEleTopEB2/nbEleBottomEB2;
  cout<<"nbEleTopEB2                          = "<<nbEleTopEB2<<endl;
  cout<<"nbEleBottomEB2                       = "<<nbEleBottomEB2<<endl;
  cout<<"EffChargedHadEB3(photon)             = "<<EffChargedHadEB3<<endl;
  float StatErrorEB2                          = sqrt(nbEleBottomEB2 - nbEleTopEB2)/nbEleBottomEB2;
  cout<<"StatErrorEB2(photon)                 = "<<StatErrorEB2<<endl;
  //---------------------------------------------------
  PhotonEtDeltaPhiEB->SetLineStyle(0);
  PhotonEtDeltaPhiEB->SetLineColor(2);
  PhotonEtDeltaPhiEB->SetLineWidth(2);
  PhotonEtDeltaPhiEB->SetMarkerColor(2); 
  PhotonEtDeltaPhiEB->SetMarkerStyle(24);
  PhotonEtDeltaPhiEB->SetMarkerSize(0.875);
  //========================================================== 
  //                 QCD pfChargedHadron in EB                                                  
  //==========================================================
  TFile *f3 = new TFile("CMSSW701-Analyse_MC_QCD_15_pt_3000_V5_withPF.root","READ");
  QCDEtDeltaPhiEB->BayesDivide(EtDeltaPhiEB2,EtDeltaPhiEB1);
  //---------------------------------------------------
  cout<<"====================== QCD EB =============================="<<endl;
  int nbEleTopEB3                             = EtDeltaPhiEB2->Integral();
  int nbEleBottomEB3                          = EtDeltaPhiEB1->Integral();
  float EffChargedHadEB3                      = (float)nbEleTopEB3/nbEleBottomEB3;
  cout<<"nbEleTopEB3                          = "<<nbEleTopEB3<<endl;
  cout<<"nbEleBottomEB3                       = "<<nbEleBottomEB3<<endl;
  cout<<"EffChargedHadEB3(QCD)                = "<<EffChargedHadEB3<<endl;
  float StatErrorEB3                          = sqrt(nbEleBottomEB3 - nbEleTopEB3)/nbEleBottomEB3;
  cout<<"StatErrorEB3(QCD)                    = "<<StatErrorEB3<<endl;
  //---------------------------------------------------
  QCDEtDeltaPhiEB->SetLineStyle(0);
  QCDEtDeltaPhiEB->SetLineColor(4);
  QCDEtDeltaPhiEB->SetLineWidth(2);
  QCDEtDeltaPhiEB->SetMarkerColor(4); 
  QCDEtDeltaPhiEB->SetMarkerStyle(22);
  QCDEtDeltaPhiEB->SetMarkerSize(0.875);
  //-----------------------------------------------------------------
  */
  mg3->Add(EleEtDeltaPhiEB);
  //mg3->Add(PhotonEtDeltaPhiEB);
  //mg3->Add(QCDEtDeltaPhiEB);
  mg3->Draw("AP");
  mg3->SetTitle("title;|#eta_{muons}|; Efficiency(N-1)");
  mg3->GetXaxis()->SetTitleOffset(1.7);
  mg3->GetYaxis()->SetTitleOffset(1.7);
  //mg3->GetXaxis()->SetRangeUser(0.0,800.0);
  //mg3->GetYaxis()->SetRangeUser(0.0,1.05);
  mg3->GetYaxis()->SetRangeUser(0.9,1.01);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  //TPaveText* tText2 = new TPaveText(0.40, 0.90, 0.90, 0.92, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  //Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t2 = tText2->AddText("DeltaPhi");
  //TText *t3 = tText2->AddText("NbOfPixelHits>0.0; CMSSW701 [1/fb]"); 
  tText2->Draw();
  //========================================================== 
  //TPaveText* tText2 = new TPaveText(0.20, 0.30, 0.60, 0.45, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  //Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t1 = tText2->AddText("Eff(N-1) = 0.9969 #pm 0.0002");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9684 #pm 0.0005");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8921 #pm 0.0057");

  //TText *t1 = tText2->AddText("Eff(N-1) = 0.9988 #pm 9.92665e-05");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9198 #pm 0.0008");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8316 #pm 0.0071");
  tText2->Draw();
  //========================================================== 
  TLegend *leg2 = new TLegend(0.60, 0.15, 0.80, 0.30);
  leg2->AddEntry(EleEtDeltaPhiEB,"Z' [M = 1 TeV/c^{2}]","p");
  //leg2->AddEntry(PhotonEtDeltaPhiEB,"#gamma + jets (EB)","p");
  //leg2->AddEntry(QCDEtDeltaPhiEB,"QCD (EB)","p");
  leg2->SetBorderSize(0.0);
  leg2->SetMargin(0.3);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(10);
  leg2->SetLineColor(0);
  //Float_t tsize2 = 0.03;
  leg2->SetTextSize(tsize2); 
  leg2->Draw();
  //======================================================================= 
  c3->Print("Eff_NbOfPixelHits_HEEP_ID_MC_CMSSW701_Eta.png","png");
  //c3->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.pdf","pdf");
  //c3->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.eps","eps");
  //=======================================================================
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  TCanvas *c4 = new TCanvas("c4","Et in EB",600,600);
  char textpro1111[100],textNDF1111[100],textRatio1111[100];           
  c4->cd(1);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  gPad->SetGridx();
  EleEtDeltaPhiEB     = new TGraphAsymmErrors;
  PhotonEtDeltaPhiEB  = new TGraphAsymmErrors;
  QCDEtDeltaPhiEB     = new TGraphAsymmErrors;
  TMultiGraph *mg4 = new TMultiGraph;
  //========================================================== 
  //                      DeltaPhi in EB                                                  
  //==========================================================
  TFile *f111 = new TFile("/lustre/cms/store/user/defilip/ZprimeAnalysis/80X/histos/histosZprimeMuMu/output_ZprimeToMuMu_M-5000_TuneCUETP8M1_13TeV-pythia8_tree.root","READ");
  //TH1F *EtaID; 
  f111->GetObject("EtaID",EtaID);
  TH1F *EtaEffptnumberOfMuonHits; f111->GetObject("EtaEffptnumberOfMuonHits",EtaEffptnumberOfMuonHits);
  //---------------------------------------------------
  cout<<"====================== NumberOfMuonHits =============================="<<endl;
  int nbEleTopEB4                             = EtaID->Integral();
  int nbEleBottomEB4                          = EtaEffptnumberOfMuonHits->Integral();
  float EffChargedHadEB4                      = (float)nbEleTopEB4/nbEleBottomEB4;
  cout<<"nbEleTopEB4                          = "<<nbEleTopEB4<<endl;
  cout<<"nbEleBottomEB4                       = "<<nbEleBottomEB4<<endl;
  cout<<"EffChargedHadEB4(electron)           = "<<EffChargedHadEB4<<endl;
  float StatErrorEB4                          = sqrt(nbEleBottomEB4 - nbEleTopEB4)/nbEleBottomEB4;
  cout<<"StatErrorEB4(electron)               = "<<StatErrorEB4<<endl;
  //---------------------------------------------------
  EleEtDeltaPhiEB->BayesDivide(EtaID,EtaEffptnumberOfMuonHits);
  EleEtDeltaPhiEB->SetLineStyle(0);
  EleEtDeltaPhiEB->SetLineColor(1);
  EleEtDeltaPhiEB->SetLineWidth(2);
  EleEtDeltaPhiEB->SetMarkerColor(1); 
  EleEtDeltaPhiEB->SetMarkerStyle(20);
  EleEtDeltaPhiEB->SetMarkerSize(0.875);
  /*
  //========================================================== 
  //                 DeltaPhi in EB                                                  
  //==========================================================
  TFile *f2 = new TFile("CMSSW701-Analyse_MC_GJet_pt_3000_V5_withPF.root","READ");
  PhotonEtDeltaPhiEB->BayesDivide(EtDeltaPhiEB2,EtDeltaPhiEB1);
  //---------------------------------------------------
  cout<<"====================== Gamma Jet EB =============================="<<endl;
  int nbEleTopEB2                             = EtDeltaPhiEB2->Integral();
  int nbEleBottomEB2                          = EtDeltaPhiEB1->Integral();
  float EffChargedHadEB3                      = (float)nbEleTopEB2/nbEleBottomEB2;
  cout<<"nbEleTopEB2                          = "<<nbEleTopEB2<<endl;
  cout<<"nbEleBottomEB2                       = "<<nbEleBottomEB2<<endl;
  cout<<"EffChargedHadEB3(photon)             = "<<EffChargedHadEB3<<endl;
  float StatErrorEB2                          = sqrt(nbEleBottomEB2 - nbEleTopEB2)/nbEleBottomEB2;
  cout<<"StatErrorEB2(photon)                 = "<<StatErrorEB2<<endl;
  //---------------------------------------------------
  PhotonEtDeltaPhiEB->SetLineStyle(0);
  PhotonEtDeltaPhiEB->SetLineColor(2);
  PhotonEtDeltaPhiEB->SetLineWidth(2);
  PhotonEtDeltaPhiEB->SetMarkerColor(2); 
  PhotonEtDeltaPhiEB->SetMarkerStyle(24);
  PhotonEtDeltaPhiEB->SetMarkerSize(0.875);
  //========================================================== 
  //                 QCD pfChargedHadron in EB                                                  
  //==========================================================
  TFile *f3 = new TFile("CMSSW701-Analyse_MC_QCD_15_pt_3000_V5_withPF.root","READ");
  QCDEtDeltaPhiEB->BayesDivide(EtDeltaPhiEB2,EtDeltaPhiEB1);
  //---------------------------------------------------
  cout<<"====================== QCD EB =============================="<<endl;
  int nbEleTopEB3                             = EtDeltaPhiEB2->Integral();
  int nbEleBottomEB3                          = EtDeltaPhiEB1->Integral();
  float EffChargedHadEB3                      = (float)nbEleTopEB3/nbEleBottomEB3;
  cout<<"nbEleTopEB3                          = "<<nbEleTopEB3<<endl;
  cout<<"nbEleBottomEB3                       = "<<nbEleBottomEB3<<endl;
  cout<<"EffChargedHadEB3(QCD)                = "<<EffChargedHadEB3<<endl;
  float StatErrorEB3                          = sqrt(nbEleBottomEB3 - nbEleTopEB3)/nbEleBottomEB3;
  cout<<"StatErrorEB3(QCD)                    = "<<StatErrorEB3<<endl;
  //---------------------------------------------------
  QCDEtDeltaPhiEB->SetLineStyle(0);
  QCDEtDeltaPhiEB->SetLineColor(4);
  QCDEtDeltaPhiEB->SetLineWidth(2);
  QCDEtDeltaPhiEB->SetMarkerColor(4); 
  QCDEtDeltaPhiEB->SetMarkerStyle(22);
  QCDEtDeltaPhiEB->SetMarkerSize(0.875);
  //-----------------------------------------------------------------
  */
  mg4->Add(EleEtDeltaPhiEB);
  //mg4->Add(PhotonEtDeltaPhiEB);
  //mg4->Add(QCDEtDeltaPhiEB);
  mg4->Draw("AP");
  mg4->SetTitle("title;|#eta_{muons}|; Efficiency(N-1)");
  mg4->GetXaxis()->SetTitleOffset(1.7);
  mg4->GetYaxis()->SetTitleOffset(1.7);
  //mg4->GetXaxis()->SetRangeUser(0.0,800.0);
  //mg4->GetYaxis()->SetRangeUser(0.0,1.05);
  mg4->GetYaxis()->SetRangeUser(0.9,1.01);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  //TPaveText* tText2 = new TPaveText(0.40, 0.90, 0.90, 0.92, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  //Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t2 = tText2->AddText("DeltaPhi");
  //  TText *t3 = tText2->AddText("NumberOfMuonHits>0.0; CMSSW701 [1/fb]"); 
  tText2->Draw();
  //========================================================== 
  //TPaveText* tText2 = new TPaveText(0.20, 0.30, 0.60, 0.45, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  //Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t1 = tText2->AddText("Eff(N-1) = 0.9969 #pm 0.0002");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9684 #pm 0.0005");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8921 #pm 0.0057");

  //TText *t1 = tText2->AddText("Eff(N-1) = 0.9982 #pm 0.00012");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9198 #pm 0.0008");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8316 #pm 0.0071");
  tText2->Draw();
  //========================================================== 
  TLegend *leg3 = new TLegend(0.60, 0.15, 0.80, 0.30);
  leg3->AddEntry(EleEtDeltaPhiEB,"Z' [M = 1 TeV/c^{2}]","p");
  //leg3->AddEntry(PhotonEtDeltaPhiEB,"#gamma + jets (EB)","p");
  //leg3->AddEntry(QCDEtDeltaPhiEB,"QCD (EB)","p");
  leg3->SetBorderSize(0.0);
  leg3->SetMargin(0.3);
  leg3->SetFillColor(0);
  leg3->SetFillStyle(10);
  leg3->SetLineColor(0);
  //Float_t tsize2 = 0.03;
  leg3->SetTextSize(tsize2); 
  leg3->Draw();
  //======================================================================= 
  c4->Print("Eff_NumberOfMuonHits_HEEP_ID_MC_CMSSW701_Eta.png","png");
  //c4->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.pdf","pdf");
  //c4->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.eps","eps");
  //=======================================================================
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  TCanvas *c5 = new TCanvas("c5","Et in EB",600,600);
  char textpro12[100],textNDF12[100],textRatio12[100];           
  c5->cd(1);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  gPad->SetGridx();
  EleEtDeltaPhiEB     = new TGraphAsymmErrors;
  PhotonEtDeltaPhiEB  = new TGraphAsymmErrors;
  QCDEtDeltaPhiEB     = new TGraphAsymmErrors;
  TMultiGraph *mg5 = new TMultiGraph;
  //========================================================== 
  //                      DeltaPhi in EB                                                  
  //==========================================================
  TFile *f1111 = new TFile("/lustre/cms/store/user/defilip/ZprimeAnalysis/80X/histos/histosZprimeMuMu/output_ZprimeToMuMu_M-5000_TuneCUETP8M1_13TeV-pythia8_tree.root","READ");
  //TH1F *EtaID; 
  f1111->GetObject("EtaID",EtaID);
  TH1F *EtaEffptnumberOfMatchedStations; f1111->GetObject("EtaEffptnumberOfMatchedStations",EtaEffptnumberOfMatchedStations);

  //---------------------------------------------------
  cout<<"====================== NumberOfMatchedStations =============================="<<endl;
  int nbEleTopEB5                             = EtaID->Integral();
  int nbEleBottomEB5                          = EtaEffptnumberOfMatchedStations->Integral();
  float EffChargedHadEB5                      = (float)nbEleTopEB5/nbEleBottomEB5;
  cout<<"nbEleTopEB5                          = "<<nbEleTopEB5<<endl;
  cout<<"nbEleBottomEB5                       = "<<nbEleBottomEB5<<endl;
  cout<<"EffChargedHadEB5(electron)           = "<<EffChargedHadEB5<<endl;
  float StatErrorEB5                          = sqrt(nbEleBottomEB5 - nbEleTopEB5)/nbEleBottomEB5;
  cout<<"StatErrorEB5(electron)               = "<<StatErrorEB5<<endl;
  //---------------------------------------------------
  EleEtDeltaPhiEB->BayesDivide(EtaID,EtaEffptnumberOfMatchedStations);
  EleEtDeltaPhiEB->SetLineStyle(0);
  EleEtDeltaPhiEB->SetLineColor(1);
  EleEtDeltaPhiEB->SetLineWidth(2);
  EleEtDeltaPhiEB->SetMarkerColor(1); 
  EleEtDeltaPhiEB->SetMarkerStyle(20);
  EleEtDeltaPhiEB->SetMarkerSize(0.875);
  /*
  //========================================================== 
  //                 DeltaPhi in EB                                                  
  //==========================================================
  TFile *f2 = new TFile("CMSSW701-Analyse_MC_GJet_pt_3000_V5_withPF.root","READ");
  PhotonEtDeltaPhiEB->BayesDivide(EtDeltaPhiEB2,EtDeltaPhiEB1);
  //---------------------------------------------------
  cout<<"====================== Gamma Jet EB =============================="<<endl;
  int nbEleTopEB2                             = EtDeltaPhiEB2->Integral();
  int nbEleBottomEB2                          = EtDeltaPhiEB1->Integral();
  float EffChargedHadEB3                      = (float)nbEleTopEB2/nbEleBottomEB2;
  cout<<"nbEleTopEB2                          = "<<nbEleTopEB2<<endl;
  cout<<"nbEleBottomEB2                       = "<<nbEleBottomEB2<<endl;
  cout<<"EffChargedHadEB3(photon)             = "<<EffChargedHadEB3<<endl;
  float StatErrorEB2                          = sqrt(nbEleBottomEB2 - nbEleTopEB2)/nbEleBottomEB2;
  cout<<"StatErrorEB2(photon)                 = "<<StatErrorEB2<<endl;
  //---------------------------------------------------
  PhotonEtDeltaPhiEB->SetLineStyle(0);
  PhotonEtDeltaPhiEB->SetLineColor(2);
  PhotonEtDeltaPhiEB->SetLineWidth(2);
  PhotonEtDeltaPhiEB->SetMarkerColor(2); 
  PhotonEtDeltaPhiEB->SetMarkerStyle(24);
  PhotonEtDeltaPhiEB->SetMarkerSize(0.875);
  //========================================================== 
  //                 QCD pfChargedHadron in EB                                                  
  //==========================================================
  TFile *f3 = new TFile("CMSSW701-Analyse_MC_QCD_15_pt_3000_V5_withPF.root","READ");
  QCDEtDeltaPhiEB->BayesDivide(EtDeltaPhiEB2,EtDeltaPhiEB1);
  //---------------------------------------------------
  cout<<"====================== QCD EB =============================="<<endl;
  int nbEleTopEB3                             = EtDeltaPhiEB2->Integral();
  int nbEleBottomEB3                          = EtDeltaPhiEB1->Integral();
  float EffChargedHadEB3                      = (float)nbEleTopEB3/nbEleBottomEB3;
  cout<<"nbEleTopEB3                          = "<<nbEleTopEB3<<endl;
  cout<<"nbEleBottomEB3                       = "<<nbEleBottomEB3<<endl;
  cout<<"EffChargedHadEB3(QCD)                = "<<EffChargedHadEB3<<endl;
  float StatErrorEB3                          = sqrt(nbEleBottomEB3 - nbEleTopEB3)/nbEleBottomEB3;
  cout<<"StatErrorEB3(QCD)                    = "<<StatErrorEB3<<endl;
  //---------------------------------------------------
  QCDEtDeltaPhiEB->SetLineStyle(0);
  QCDEtDeltaPhiEB->SetLineColor(4);
  QCDEtDeltaPhiEB->SetLineWidth(2);
  QCDEtDeltaPhiEB->SetMarkerColor(4); 
  QCDEtDeltaPhiEB->SetMarkerStyle(22);
  QCDEtDeltaPhiEB->SetMarkerSize(0.875);
  //-----------------------------------------------------------------
  */
  mg5->Add(EleEtDeltaPhiEB);
  //mg5->Add(PhotonEtDeltaPhiEB);
  //mg5->Add(QCDEtDeltaPhiEB);
  mg5->Draw("AP");
  mg5->SetTitle("title;|#eta_{muons}|; Efficiency(N-1)");
  mg5->GetXaxis()->SetTitleOffset(1.7);
  mg5->GetYaxis()->SetTitleOffset(1.7);
  //mg5->GetXaxis()->SetRangeUser(0.0,800.0);
  //mg5->GetYaxis()->SetRangeUser(0.0,1.05);
  mg5->GetYaxis()->SetRangeUser(0.9,1.01);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  //TPaveText* tText2 = new TPaveText(0.40, 0.90, 0.90, 0.92, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  //Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t2 = tText2->AddText("DeltaPhi");
  //TText *t3 = tText2->AddText("NumberOfMatchedStations>1.0; CMSSW701 [1/fb]"); 
  tText2->Draw();
  //========================================================== 
  //TPaveText* tText2 = new TPaveText(0.20, 0.30, 0.60, 0.45, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  //Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t1 = tText2->AddText("Eff(N-1) = 0.9969 #pm 0.0002");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9684 #pm 0.0005");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8921 #pm 0.0057");

  //TText *t1 = tText2->AddText("Eff(N-1) = 0.9937 #pm 0.0002");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9198 #pm 0.0008");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8316 #pm 0.0071");
  tText2->Draw();
  //========================================================== 
  TLegend *leg4 = new TLegend(0.60, 0.15, 0.80, 0.30);
  leg4->AddEntry(EleEtDeltaPhiEB,"Z' [M = 1 TeV/c^{2}]","p");
  //leg4->AddEntry(PhotonEtDeltaPhiEB,"#gamma + jets (EB)","p");
  //leg4->AddEntry(QCDEtDeltaPhiEB,"QCD (EB)","p");
  leg4->SetBorderSize(0.0);
  leg4->SetMargin(0.3);
  leg4->SetFillColor(0);
  leg4->SetFillStyle(10);
  leg4->SetLineColor(0);
  //Float_t tsize2 = 0.03;
  leg4->SetTextSize(tsize2); 
  leg4->Draw();
  //======================================================================= 
  c5->Print("Eff_NumberOfMatchedStations_HEEP_ID_MC_CMSSW701_Eta.png","png");
  //c5->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.pdf","pdf");
  //c5->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.eps","eps");
  //=======================================================================
  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  TCanvas *c6 = new TCanvas("c6","Et in EB",600,600);
  char textpro122[100],textNDF122[100],textRatio122[100];           
  c6->cd(1);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  gPad->SetGridx();
  EleEtDeltaPhiEB     = new TGraphAsymmErrors;
  PhotonEtDeltaPhiEB  = new TGraphAsymmErrors;
  QCDEtDeltaPhiEB     = new TGraphAsymmErrors;
  TMultiGraph *mg6 = new TMultiGraph;
  //========================================================== 
  //                      DeltaPhi in EB                                                  
  //==========================================================
  TFile *f11111 = new TFile("/lustre/cms/store/user/defilip/ZprimeAnalysis/80X/histos/histosZprimeMuMu/output_ZprimeToMuMu_M-5000_TuneCUETP8M1_13TeV-pythia8_tree.root","READ");
  //TH1F *EtaID; 
  f11111->GetObject("EtaID",EtaID);
  TH1F *EtaEffptTrackIso; f11111->GetObject("EtaEffptTrackIso",EtaEffptTrackIso);

  //---------------------------------------------------
  cout<<"====================== TrackIso =============================="<<endl;
  int nbEleTopEB6                             = EtaID->Integral();
  int nbEleBottomEB6                          = EtaEffptTrackIso->Integral();
  float EffChargedHadEB6                      = (float)nbEleTopEB6/nbEleBottomEB6;
  cout<<"nbEleTopEB6                          = "<<nbEleTopEB6<<endl;
  cout<<"nbEleBottomEB6                       = "<<nbEleBottomEB6<<endl;
  cout<<"EffChargedHadEB6(electron)           = "<<EffChargedHadEB6<<endl;
  float StatErrorEB6                          = sqrt(nbEleBottomEB6 - nbEleTopEB6)/nbEleBottomEB6;
  cout<<"StatErrorEB6(electron)               = "<<StatErrorEB6<<endl;
  //---------------------------------------------------
  EleEtDeltaPhiEB->BayesDivide(EtaID,EtaEffptTrackIso);
  EleEtDeltaPhiEB->SetLineStyle(0);
  EleEtDeltaPhiEB->SetLineColor(1);
  EleEtDeltaPhiEB->SetLineWidth(2);
  EleEtDeltaPhiEB->SetMarkerColor(1); 
  EleEtDeltaPhiEB->SetMarkerStyle(20);
  EleEtDeltaPhiEB->SetMarkerSize(0.875);
  /*
  //========================================================== 
  //                 DeltaPhi in EB                                                  
  //==========================================================
  TFile *f2 = new TFile("CMSSW701-Analyse_MC_GJet_pt_3000_V5_withPF.root","READ");
  PhotonEtDeltaPhiEB->BayesDivide(EtDeltaPhiEB2,EtDeltaPhiEB1);
  //---------------------------------------------------
  cout<<"====================== Gamma Jet EB =============================="<<endl;
  int nbEleTopEB2                             = EtDeltaPhiEB2->Integral();
  int nbEleBottomEB2                          = EtDeltaPhiEB1->Integral();
  float EffChargedHadEB3                      = (float)nbEleTopEB2/nbEleBottomEB2;
  cout<<"nbEleTopEB2                          = "<<nbEleTopEB2<<endl;
  cout<<"nbEleBottomEB2                       = "<<nbEleBottomEB2<<endl;
  cout<<"EffChargedHadEB3(photon)             = "<<EffChargedHadEB3<<endl;
  float StatErrorEB2                          = sqrt(nbEleBottomEB2 - nbEleTopEB2)/nbEleBottomEB2;
  cout<<"StatErrorEB2(photon)                 = "<<StatErrorEB2<<endl;
  //---------------------------------------------------
  PhotonEtDeltaPhiEB->SetLineStyle(0);
  PhotonEtDeltaPhiEB->SetLineColor(2);
  PhotonEtDeltaPhiEB->SetLineWidth(2);
  PhotonEtDeltaPhiEB->SetMarkerColor(2); 
  PhotonEtDeltaPhiEB->SetMarkerStyle(24);
  PhotonEtDeltaPhiEB->SetMarkerSize(0.875);
  //========================================================== 
  //                 QCD pfChargedHadron in EB                                                  
  //==========================================================
  TFile *f3 = new TFile("CMSSW701-Analyse_MC_QCD_15_pt_3000_V5_withPF.root","READ");
  QCDEtDeltaPhiEB->BayesDivide(EtDeltaPhiEB2,EtDeltaPhiEB1);
  //---------------------------------------------------
  cout<<"====================== QCD EB =============================="<<endl;
  int nbEleTopEB3                             = EtDeltaPhiEB2->Integral();
  int nbEleBottomEB3                          = EtDeltaPhiEB1->Integral();
  float EffChargedHadEB3                      = (float)nbEleTopEB3/nbEleBottomEB3;
  cout<<"nbEleTopEB3                          = "<<nbEleTopEB3<<endl;
  cout<<"nbEleBottomEB3                       = "<<nbEleBottomEB3<<endl;
  cout<<"EffChargedHadEB3(QCD)                = "<<EffChargedHadEB3<<endl;
  float StatErrorEB3                          = sqrt(nbEleBottomEB3 - nbEleTopEB3)/nbEleBottomEB3;
  cout<<"StatErrorEB3(QCD)                    = "<<StatErrorEB3<<endl;
  //---------------------------------------------------
  QCDEtDeltaPhiEB->SetLineStyle(0);
  QCDEtDeltaPhiEB->SetLineColor(4);
  QCDEtDeltaPhiEB->SetLineWidth(2);
  QCDEtDeltaPhiEB->SetMarkerColor(4); 
  QCDEtDeltaPhiEB->SetMarkerStyle(22);
  QCDEtDeltaPhiEB->SetMarkerSize(0.875);
  //-----------------------------------------------------------------
  */
  mg6->Add(EleEtDeltaPhiEB);
  //mg6->Add(PhotonEtDeltaPhiEB);
  //mg6->Add(QCDEtDeltaPhiEB);
  mg6->Draw("AP");
  mg6->SetTitle("title;|#eta_{muons}|; Efficiency(N-1)");
  mg6->GetXaxis()->SetTitleOffset(1.7);
  mg6->GetYaxis()->SetTitleOffset(1.7);
  //mg6->GetXaxis()->SetRangeUser(0.0,800.0);
  //mg6->GetYaxis()->SetRangeUser(0.0,1.05);
  mg6->GetYaxis()->SetRangeUser(0.9,1.01);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  //TPaveText* tText2 = new TPaveText(0.40, 0.90, 0.90, 0.92, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  //Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t2 = tText2->AddText("DeltaPhi");
  //TText *t3 = tText2->AddText("Track Iso>0.10; CMSSW701 [1/fb]"); 
  tText2->Draw();
  //========================================================== 
  //TPaveText* tText2 = new TPaveText(0.20, 0.30, 0.60, 0.45, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  //Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t1 = tText2->AddText("Eff(N-1) = 0.9969 #pm 0.0002");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9684 #pm 0.0005");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8921 #pm 0.0057");

  //TText *t1 = tText2->AddText("Eff(N-1) = 0.9952 #pm 0.0002");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9198 #pm 0.0008");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8316 #pm 0.0071");
  tText2->Draw();
  //========================================================== 
  TLegend *leg5 = new TLegend(0.60, 0.15, 0.80, 0.30);
  leg5->AddEntry(EleEtDeltaPhiEB,"Z' [M = 1 TeV/c^{2}]","p");
  //leg5->AddEntry(PhotonEtDeltaPhiEB,"#gamma + jets (EB)","p");
  //leg5->AddEntry(QCDEtDeltaPhiEB,"QCD (EB)","p");
  leg5->SetBorderSize(0.0);
  leg5->SetMargin(0.3);
  leg5->SetFillColor(0);
  leg5->SetFillStyle(10);
  leg5->SetLineColor(0);
  //Float_t tsize2 = 0.03;
  leg5->SetTextSize(tsize2); 
  leg5->Draw();
  //======================================================================= 
  c6->Print("Eff_Trackiso_HEEP_ID_MC_CMSSW701_Eta.png","png");
  //c6->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.pdf","pdf");
  //c6->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.eps","eps");
  //=======================================================================
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  TCanvas *c7 = new TCanvas("c7","Et in EB",600,600);
  char textpro123[100],textNDF123[100],textRatio123[100];           
  c7->cd(1);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  gPad->SetGridx();
  EleEtDeltaPhiEB     = new TGraphAsymmErrors;
  PhotonEtDeltaPhiEB  = new TGraphAsymmErrors;
  QCDEtDeltaPhiEB     = new TGraphAsymmErrors;
  mg6 = new TMultiGraph;
  //========================================================== 
  //                      DeltaPhi in EB                                                  
  //==========================================================
  TFile *f111111 = new TFile("/lustre/cms/store/user/defilip/ZprimeAnalysis/80X/histos/histosZprimeMuMu/output_ZprimeToMuMu_M-5000_TuneCUETP8M1_13TeV-pythia8_tree.root","READ");
  //TH1F *EtaID; 
  f111111->GetObject("EtaID",EtaID);
  TH1F *EtaEffptabsdsy; f111111->GetObject("EtaEffptabsdsy",EtaEffptabsdsy);

  //---------------------------------------------------
  cout<<"====================== dxy =============================="<<endl;
  int nbEleTopEB7                             = EtaID->Integral();
  int nbEleBottomEB7                          = EtaEffptabsdsy->Integral();
  float EffChargedHadEB7                      = (float)nbEleTopEB7/nbEleBottomEB7;
  cout<<"nbEleTopEB7                          = "<<nbEleTopEB7<<endl;
  cout<<"nbEleBottomEB7                       = "<<nbEleBottomEB7<<endl;
  cout<<"EffChargedHadEB7(electron)           = "<<EffChargedHadEB7<<endl;
  float StatErrorEB7                          = sqrt(nbEleBottomEB7 - nbEleTopEB7)/nbEleBottomEB7;
  cout<<"StatErrorEB7(electron)               = "<<StatErrorEB7<<endl;
  //---------------------------------------------------
  EleEtDeltaPhiEB->BayesDivide(EtaID,EtaEffptabsdsy);
  EleEtDeltaPhiEB->SetLineStyle(0);
  EleEtDeltaPhiEB->SetLineColor(1);
  EleEtDeltaPhiEB->SetLineWidth(2);
  EleEtDeltaPhiEB->SetMarkerColor(1); 
  EleEtDeltaPhiEB->SetMarkerStyle(20);
  EleEtDeltaPhiEB->SetMarkerSize(0.875);
  /*
  //========================================================== 
  //                 DeltaPhi in EB                                                  
  //==========================================================
  TFile *f2 = new TFile("CMSSW701-Analyse_MC_GJet_pt_3000_V5_withPF.root","READ");
  PhotonEtDeltaPhiEB->BayesDivide(EtDeltaPhiEB2,EtDeltaPhiEB1);
  //---------------------------------------------------
  cout<<"====================== Gamma Jet EB =============================="<<endl;
  int nbEleTopEB2                             = EtDeltaPhiEB2->Integral();
  int nbEleBottomEB2                          = EtDeltaPhiEB1->Integral();
  float EffChargedHadEB3                      = (float)nbEleTopEB2/nbEleBottomEB2;
  cout<<"nbEleTopEB2                          = "<<nbEleTopEB2<<endl;
  cout<<"nbEleBottomEB2                       = "<<nbEleBottomEB2<<endl;
  cout<<"EffChargedHadEB3(photon)             = "<<EffChargedHadEB3<<endl;
  float StatErrorEB2                          = sqrt(nbEleBottomEB2 - nbEleTopEB2)/nbEleBottomEB2;
  cout<<"StatErrorEB2(photon)                 = "<<StatErrorEB2<<endl;
  //---------------------------------------------------
  PhotonEtDeltaPhiEB->SetLineStyle(0);
  PhotonEtDeltaPhiEB->SetLineColor(2);
  PhotonEtDeltaPhiEB->SetLineWidth(2);
  PhotonEtDeltaPhiEB->SetMarkerColor(2); 
  PhotonEtDeltaPhiEB->SetMarkerStyle(24);
  PhotonEtDeltaPhiEB->SetMarkerSize(0.875);
  //========================================================== 
  //                 QCD pfChargedHadron in EB                                                  
  //==========================================================
  TFile *f3 = new TFile("CMSSW701-Analyse_MC_QCD_15_pt_3000_V5_withPF.root","READ");
  QCDEtDeltaPhiEB->BayesDivide(EtDeltaPhiEB2,EtDeltaPhiEB1);
  //---------------------------------------------------
  cout<<"====================== QCD EB =============================="<<endl;
  int nbEleTopEB3                             = EtDeltaPhiEB2->Integral();
  int nbEleBottomEB3                          = EtDeltaPhiEB1->Integral();
  float EffChargedHadEB3                      = (float)nbEleTopEB3/nbEleBottomEB3;
  cout<<"nbEleTopEB3                          = "<<nbEleTopEB3<<endl;
  cout<<"nbEleBottomEB3                       = "<<nbEleBottomEB3<<endl;
  cout<<"EffChargedHadEB3(QCD)                = "<<EffChargedHadEB3<<endl;
  float StatErrorEB3                          = sqrt(nbEleBottomEB3 - nbEleTopEB3)/nbEleBottomEB3;
  cout<<"StatErrorEB3(QCD)                    = "<<StatErrorEB3<<endl;
  //---------------------------------------------------
  QCDEtDeltaPhiEB->SetLineStyle(0);
  QCDEtDeltaPhiEB->SetLineColor(4);
  QCDEtDeltaPhiEB->SetLineWidth(2);
  QCDEtDeltaPhiEB->SetMarkerColor(4); 
  QCDEtDeltaPhiEB->SetMarkerStyle(22);
  QCDEtDeltaPhiEB->SetMarkerSize(0.875);
  //-----------------------------------------------------------------
  */
  mg6->Add(EleEtDeltaPhiEB);
  //mg6->Add(PhotonEtDeltaPhiEB);
  //mg6->Add(QCDEtDeltaPhiEB);
  mg6->Draw("AP");
  mg6->SetTitle("title;|#eta_{muons}|; Efficiency(N-1)");
  mg6->GetXaxis()->SetTitleOffset(1.7);
  mg6->GetYaxis()->SetTitleOffset(1.7);
  //mg6->GetXaxis()->SetRangeUser(0.0,800.0);
  //mg6->GetYaxis()->SetRangeUser(0.0,1.05);
  mg6->GetYaxis()->SetRangeUser(0.9,1.01);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  //TPaveText* tText2 = new TPaveText(0.40, 0.90, 0.90, 0.92, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  //Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t2 = tText2->AddText("DeltaPhi");
  //TText *t3 = tText2->AddText("|dxy|<0.2; CMSSW701 [1/fb]"); 
  tText2->Draw();
  //========================================================== 
  //TPaveText* tText2 = new TPaveText(0.20, 0.30, 0.60, 0.45, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  //Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t1 = tText2->AddText("Eff(N-1) = 0.9969 #pm 0.0002");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9684 #pm 0.0005");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8921 #pm 0.0057");

  //TText *t1 = tText2->AddText("Eff(N-1) = 1.0 #pm 0.0");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9198 #pm 0.0008");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8316 #pm 0.0071");
  tText2->Draw();
  //========================================================== 
  TLegend *leg6 = new TLegend(0.60, 0.15, 0.80, 0.30);
  leg6->AddEntry(EleEtDeltaPhiEB,"Z' [M = 1 TeV/c^{2}]","p");
  //leg6->AddEntry(PhotonEtDeltaPhiEB,"#gamma + jets (EB)","p");
  //leg6->AddEntry(QCDEtDeltaPhiEB,"QCD (EB)","p");
  leg6->SetBorderSize(0.0);
  leg6->SetMargin(0.3);
  leg6->SetFillColor(0);
  leg6->SetFillStyle(10);
  leg6->SetLineColor(0);
  //Float_t tsize2 = 0.03;
  leg6->SetTextSize(tsize2); 
  leg6->Draw();
  //======================================================================= 
  c7->Print("Eff_Dxy_HEEP_ID_MC_CMSSW701_Eta.png","png");
  //c7->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.pdf","pdf");
  //c7->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.eps","eps");
  //=======================================================================



  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  TCanvas *c99 = new TCanvas("c99","#delta pt/pt",700,700);
  char textpro13[100],textNDF13[100],textRatio13[100];           
  c99->cd(1);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  gPad->SetGridx();
  TGraphAsymmErrors *EtaEffForErrorOnPt     = new TGraphAsymmErrors;
  TGraphAsymmErrors *Effdxy                = new TGraphAsymmErrors;
  TGraphAsymmErrors *EtaEffForTrackIso      = new TGraphAsymmErrors;
  TGraphAsymmErrors *EtaEffForNbTrackLayers = new TGraphAsymmErrors;
  TGraphAsymmErrors *EtaEffForNbPixelHits   = new TGraphAsymmErrors;
  TGraphAsymmErrors *EtaEffForNbMuonHits    = new TGraphAsymmErrors;
  TGraphAsymmErrors *EtaEffForNbMatchedStations = new TGraphAsymmErrors;
  TMultiGraph *mg1000 = new TMultiGraph;
  //========================================================== 
  //                      DeltaPhi in EB                                                  
  //==========================================================
  TFile *f1111111 = new TFile("/lustre/cms/store/user/defilip/ZprimeAnalysis/80X/histos/histosZprimeMuMu/output_ZprimeToMuMu_M-5000_TuneCUETP8M1_13TeV-pythia8_tree.root","READ");
  //TH1F *EtaID; 
  f1111111->GetObject("EtaID",EtaID);
  //TH1F *EtaEffpterror; 
  f1111111->GetObject("EtaEffpterror",EtaEffpterror);
  //---------------------------------------------------
  cout<<"====================== dpt/pt =============================="<<endl;
  int nbEleTopEB10                             = EtaID->Integral();
  int nbEleBottomEB10                          = EtaEffpterror->Integral();
  float EffChargedHadEB10                      = (float)nbEleTopEB10/nbEleBottomEB10;
  cout<<"nbEleTopEB10                          = "<<nbEleTopEB10<<endl;
  cout<<"nbEleBottomEB10                       = "<<nbEleBottomEB10<<endl;
  cout<<"EffChargedHadEB10(electron)           = "<<EffChargedHadEB10<<endl;
  float StatErrorEB10                          = sqrt(nbEleBottomEB10 - nbEleTopEB10)/nbEleBottomEB10;
  cout<<"StatErrorEB10(electron)               = "<<StatErrorEB10<<endl;
  //---------------------------------------------------
  EtaEffForErrorOnPt->BayesDivide(EtaID,EtaEffpterror);
  EtaEffForErrorOnPt->SetLineStyle(0);
  EtaEffForErrorOnPt->SetLineColor(1);
  EtaEffForErrorOnPt->SetLineWidth(2);
  EtaEffForErrorOnPt->SetMarkerColor(1); 
  EtaEffForErrorOnPt->SetMarkerStyle(20);
  EtaEffForErrorOnPt->SetMarkerSize(0.875);
  //========================================================== 
  //                 Dxy                                                 
  //==========================================================
  cout<<"==================== |dxy| =============================="<<endl;
  int nbEleTopEB70                             = EtaID->Integral();
  int nbEleBottomEB70                          = EtaEffptabsdsy->Integral();
  float EffChargedHadEB70                      = (float)nbEleTopEB70/nbEleBottomEB70;
  cout<<"nbEleTopEB70                          = "<<nbEleTopEB70<<endl;
  cout<<"nbEleBottomEB70                       = "<<nbEleBottomEB70<<endl;
  cout<<"EffChargedHadEB70(electron)           = "<<EffChargedHadEB70<<endl;
  float StatErrorEB70                          = sqrt(nbEleBottomEB70 - nbEleTopEB70)/nbEleBottomEB70;
  cout<<"StatErrorEB70(electron)               = "<<StatErrorEB70<<endl;
  //---------------------------------------------------
  Effdxy->BayesDivide(EtaID,EtaEffptabsdsy);
  Effdxy->SetLineStyle(0);
  Effdxy->SetLineColor(2);
  Effdxy->SetLineWidth(2);
  Effdxy->SetMarkerColor(2); 
  Effdxy->SetMarkerStyle(20);
  Effdxy->SetMarkerSize(0.875);
  //========================================================== 
  //                 track isolation                                                 
  //==========================================================
  cout<<"====================== track iso =============================="<<endl;
  int nbEleTopEB20                             = EtaID->Integral();
  int nbEleBottomEB20                          = EtaEffptTrackIso->Integral();
  float EffChargedHadEB20                      = (float)nbEleTopEB20/nbEleBottomEB20;
  cout<<"nbEleTopEB20                          = "<<nbEleTopEB20<<endl;
  cout<<"nbEleBottomEB20                       = "<<nbEleBottomEB20<<endl;
  cout<<"EffChargedHadEB20(electron)           = "<<EffChargedHadEB20<<endl;
  float StatErrorEB20                          = sqrt(nbEleBottomEB20 - nbEleTopEB20)/nbEleBottomEB20;
  cout<<"StatErrorEB20(electron)               = "<<StatErrorEB20<<endl;
  //---------------------------------------------------
  EtaEffForTrackIso->BayesDivide(EtaID,EtaEffptTrackIso);
  EtaEffForTrackIso->SetLineStyle(0);
  EtaEffForTrackIso->SetLineWidth(2);
  EtaEffForTrackIso->SetLineColor(3);
  EtaEffForTrackIso->SetMarkerColor(3); 
  EtaEffForTrackIso->SetMarkerStyle(21);
  EtaEffForTrackIso->SetMarkerSize(0.875);
  //========================================================== 
  //                              NbTrackLayers                                 
  //==========================================================
  cout<<"====================== NbTrackLayers =============================="<<endl;
  //TH1F *EtaEffptnumberOftrackerLayers; f1->GetObject("EtaEffptnumberOftrackerLayers",EtaEffptnumberOftrackerLayers);
  int nbEleTopEB30                             = EtaID->Integral();
  int nbEleBottomEB30                          = EtaEffptnumberOftrackerLayers->Integral();
  float EffChargedHadEB30                      = (float)nbEleTopEB30/nbEleBottomEB30;
  cout<<"nbEleTopEB30                          = "<<nbEleTopEB30<<endl;
  cout<<"nbEleBottomEB30                       = "<<nbEleBottomEB30<<endl;
  cout<<"EffChargedHadEB30(electron)           = "<<EffChargedHadEB30<<endl;
  float StatErrorEB30                          = sqrt(nbEleBottomEB30 - nbEleTopEB30)/nbEleBottomEB30;
  cout<<"StatErrorEB30(electron)               = "<<StatErrorEB30<<endl;
  //---------------------------------------------------
  EtaEffForNbTrackLayers->BayesDivide(EtaID,EtaEffptnumberOftrackerLayers);
  EtaEffForNbTrackLayers->SetLineStyle(0);
  EtaEffForNbTrackLayers->SetLineColor(4);
  EtaEffForNbTrackLayers->SetLineWidth(2);
  EtaEffForNbTrackLayers->SetMarkerColor(4); 
  EtaEffForNbTrackLayers->SetMarkerStyle(22);
  EtaEffForNbTrackLayers->SetMarkerSize(0.875);
  //========================================================== 
  //                      NbPixelHits                                                  
  //==========================================================
  cout<<"==================== NbPixelHits =============================="<<endl;
  int nbEleTopEB40                             = EtaID->Integral();
  int nbEleBottomEB40                          = EtaEffptnumberOfPixelHits->Integral();
  float EffChargedHadEB40                      = (float)nbEleTopEB40/nbEleBottomEB40;
  cout<<"nbEleTopEB40                          = "<<nbEleTopEB40<<endl;
  cout<<"nbEleBottomEB40                       = "<<nbEleBottomEB40<<endl;
  cout<<"EffChargedHadEB40(electron)           = "<<EffChargedHadEB40<<endl;
  float StatErrorEB40                          = sqrt(nbEleBottomEB40 - nbEleTopEB40)/nbEleBottomEB40;
  cout<<"StatErrorEB40(electron)               = "<<StatErrorEB40<<endl;
  //---------------------------------------------------
  EtaEffForNbPixelHits->BayesDivide(EtaID,EtaEffptnumberOfPixelHits);
  EtaEffForNbPixelHits->SetLineStyle(0);
  EtaEffForNbPixelHits->SetLineColor(6);
  EtaEffForNbPixelHits->SetLineWidth(2);
  EtaEffForNbPixelHits->SetMarkerColor(6); 
  EtaEffForNbPixelHits->SetMarkerStyle(23);
  EtaEffForNbPixelHits->SetMarkerSize(0.875);
  //========================================================== 
  //                      NbMuonHits
  //==========================================================
  cout<<"==================== NbMuonHits =============================="<<endl;
  int nbEleTopEB50                             = EtaID->Integral();
  int nbEleBottomEB50                          = EtaEffptnumberOfMuonHits->Integral();
  float EffChargedHadEB50                      = (float)nbEleTopEB50/nbEleBottomEB50;
  cout<<"nbEleTopEB50                          = "<<nbEleTopEB50<<endl;
  cout<<"nbEleBottomEB50                       = "<<nbEleBottomEB50<<endl;
  cout<<"EffChargedHadEB50(electron)           = "<<EffChargedHadEB50<<endl;
  float StatErrorEB50                          = sqrt(nbEleBottomEB50 - nbEleTopEB50)/nbEleBottomEB50;
  cout<<"StatErrorEB50(electron)               = "<<StatErrorEB50<<endl;
  //---------------------------------------------------
  EtaEffForNbMuonHits->BayesDivide(EtaID,EtaEffptnumberOfMuonHits);
  EtaEffForNbMuonHits->SetLineStyle(0);
  EtaEffForNbMuonHits->SetLineColor(7);
  EtaEffForNbMuonHits->SetLineWidth(2);
  EtaEffForNbMuonHits->SetMarkerColor(7); 
  EtaEffForNbMuonHits->SetMarkerStyle(25);
  EtaEffForNbMuonHits->SetMarkerSize(0.875);
  //========================================================== 
  //                      NbMatchedStations
  //==========================================================
  cout<<"==================== NbMatchedStations =============================="<<endl;
  int nbEleTopEB60                             = EtaID->Integral();
  int nbEleBottomEB60                          = EtaEffptnumberOfMatchedStations->Integral();
  float EffChargedHadEB60                      = (float)nbEleTopEB60/nbEleBottomEB60;
  cout<<"nbEleTopEB60                          = "<<nbEleTopEB60<<endl;
  cout<<"nbEleBottomEB60                       = "<<nbEleBottomEB60<<endl;
  cout<<"EffChargedHadEB60(electron)           = "<<EffChargedHadEB60<<endl;
  float StatErrorEB60                          = sqrt(nbEleBottomEB60 - nbEleTopEB60)/nbEleBottomEB60;
  cout<<"StatErrorEB60(electron)               = "<<StatErrorEB60<<endl;
  //---------------------------------------------------
  EtaEffForNbMatchedStations->BayesDivide(EtaID,EtaEffptnumberOfMatchedStations);
  EtaEffForNbMatchedStations->SetLineStyle(0);
  EtaEffForNbMatchedStations->SetLineColor(9);
  EtaEffForNbMatchedStations->SetLineWidth(2);
  EtaEffForNbMatchedStations->SetMarkerColor(9); 
  EtaEffForNbMatchedStations->SetMarkerStyle(28);
  EtaEffForNbMatchedStations->SetMarkerSize(0.875);
  //----------------------------------------------------------------- 
  mg1000->Add(EtaEffForErrorOnPt);
  mg1000->Add(Effdxy);
  mg1000->Add(EtaEffForTrackIso);
  mg1000->Add(EtaEffForNbTrackLayers);
  mg1000->Add(EtaEffForNbPixelHits);
  mg1000->Add(EtaEffForNbMuonHits);
  mg1000->Add(EtaEffForNbMatchedStations);
  //mg1000->Add();
  //mg1000->Add();

  mg1000->Draw("AP");
  mg1000->SetTitle("title;pt^{BestMuonTrack} [GeV]; Efficiency(N-1)");
  mg1000->GetXaxis()->SetTitleOffset(1.7);
  mg1000->GetYaxis()->SetTitleOffset(1.7);
  mg1000->GetXaxis()->SetLabelSize(0.03);
  mg1000->GetYaxis()->SetLabelSize(0.03);
  //mg1000->GetXaxis()->SetRangeUser(0.0,800.0);
  //mg1000->GetYaxis()->SetRangeUser(0.0,1.05);
  mg1000->GetYaxis()->SetRangeUser(0.80,1.015);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  //TPaveText* tText2 = new TPaveText(0.40, 0.90, 0.90, 0.92, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  //Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t2 = tText2->AddText("DeltaPhi");
  //TText *t3 = tText2->AddText("Z'(M = 5 TeV); CMSSW720 [Phys14]"); 
  tText2->Draw();
  //========================================================== 
  //TPaveText* tText2 = new TPaveText(0.20, 0.20, 0.50, 0.60, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  //Float_t tsize = 0.025;
  tText2->SetTextSize(tsize); 
  //TText *t1 = tText2->AddText("Eff(|dxy|) =1.0 #pm 0.0");
  TText *t2 = tText2->AddText("Eff(NbTrackLayers) = 0.9997 #pm 0.000");
  //TText *t3 = tText2->AddText("Eff(NbMuonHits) = 0.9986 #pm 0.0000");
  TText *t4 = tText2->AddText("Eff(NbPixelHits) = 0.9980 #pm 0.0001");
  TText *t5 = tText2->AddText("Eff(track iso) = 0.9950 #pm 0.0002");
  TText *t6 = tText2->AddText("Eff(dpt/pt) = 0.9950 #pm 0.0002");
  TText *t7 = tText2->AddText("Eff(NbMatchedStations) = 0.9885 #pm 0.0003");
  tText2->Draw();
  //========================================================== 
  TLegend *leg7 = new TLegend(0.60, 0.20, 0.80, 0.50);
  leg7->AddEntry(Effdxy,"|dxy|","p");
  leg7->AddEntry(EtaEffForNbTrackLayers,"NbTrackLayers","p");  
  leg7->AddEntry(EtaEffForNbMuonHits,"NbMuonHits","p");
  leg7->AddEntry(EtaEffForNbPixelHits,"NbPixelHits","p");
  leg7->AddEntry(EtaEffForTrackIso,"track iso","p");
  leg7->AddEntry(EtaEffForErrorOnPt,"dpt/pt","p");
  leg7->AddEntry(EtaEffForNbMatchedStations,"NbMatchedStations","p");
  leg7->SetBorderSize(0.0);
  leg7->SetMargin(0.3);
  leg7->SetFillColor(0);
  leg7->SetFillStyle(10);
  leg7->SetLineColor(0);
  //Float_t tsize2 = 0.03;
  leg7->SetTextSize(tsize2); 
  leg7->Draw();
  //======================================================================= 
  c99->Print("Eff_All_Muon_ID_MC_CMSSW720_Et.png","png");
  //c99->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.pdf","pdf");
  //c99->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.eps","eps");
  //=======================================================================







  

}





