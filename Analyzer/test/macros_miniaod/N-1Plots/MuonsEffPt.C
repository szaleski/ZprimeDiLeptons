#include "TProfile.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TF2.h"
#include <string>
using std::string;
void MuonsEffPt(){
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
  TCanvas *c10 = new TCanvas("c10","#delta pt/pt",600,600);
  char textpro1[100],textNDF1[100],textRatio1[100];           
  c10->cd(1);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  gPad->SetGridx();
  TGraphAsymmErrors *EffErrorOnPtMuonBestTrack   = new TGraphAsymmErrors;
  TGraphAsymmErrors *EffErrorOnPtTunePMuonBestTrack   = new TGraphAsymmErrors;
  TGraphAsymmErrors *EffForErrorOnPt539   = new TGraphAsymmErrors;
  TMultiGraph *mg10 = new TMultiGraph;
  //========================================================== 
  //                      Dpt/pt cmssw720                                                  
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720_miniaod.root","READ");
  //---------------------------------------------------
  cout<<"========== dpt/pt MuonBestTrack ======================"<<endl;
  int nbEleTopEB101                             = PtID->Integral();
  int nbEleBottomEB101                          = PtEffpterror->Integral();
  float EffChargedHadEB101                      = (float)nbEleTopEB101/nbEleBottomEB101;
  cout<<"nbEleTopEB101                          = "<<nbEleTopEB101<<endl;
  cout<<"nbEleBottomEB101                       = "<<nbEleBottomEB101<<endl;
  cout<<"EffChargedHadEB101(electron)           = "<<EffChargedHadEB101<<endl;
  float StatErrorEB101                          = sqrt(nbEleBottomEB101 - nbEleTopEB101)/nbEleBottomEB101;
  cout<<"StatErrorEB101(electron)               = "<<StatErrorEB101<<endl;
  //---------------------------------------------------
  EffErrorOnPtMuonBestTrack->BayesDivide(PtID,PtEffpterror);
  EffErrorOnPtMuonBestTrack->SetLineStyle(0);
  EffErrorOnPtMuonBestTrack->SetLineColor(1);
  EffErrorOnPtMuonBestTrack->SetLineWidth(2);
  EffErrorOnPtMuonBestTrack->SetMarkerColor(1); 
  EffErrorOnPtMuonBestTrack->SetMarkerStyle(20);
  EffErrorOnPtMuonBestTrack->SetMarkerSize(0.875);
  //========================================================== 
  //                      Dpt/pt cmssw706                                                  
  //==========================================================
  TFile *f121 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  //---------------------------------------------------
  cout<<"========== dpt/pt tunePMuonBestTrack ======================"<<endl;
  int nbEleTopEB102                             = PtID->Integral();
  int nbEleBottomEB102                          = PtEffpterror->Integral();
  float EffChargedHadEB102                      = (float)nbEleTopEB102/nbEleBottomEB102;
  cout<<"nbEleTopEB102                          = "<<nbEleTopEB102<<endl;
  cout<<"nbEleBottomEB102                       = "<<nbEleBottomEB102<<endl;
  cout<<"EffChargedHadEB102(electron)           = "<<EffChargedHadEB102<<endl;
  float StatErrorEB102                          = sqrt(nbEleBottomEB102 - nbEleTopEB102)/nbEleBottomEB102;
  cout<<"StatErrorEB102(electron)               = "<<StatErrorEB102<<endl;
  //---------------------------------------------------
  EffErrorOnPtTunePMuonBestTrack->BayesDivide(PtID,PtEffpterror);
  EffErrorOnPtTunePMuonBestTrack->SetLineStyle(0);
  EffErrorOnPtTunePMuonBestTrack->SetLineColor(2);
  EffErrorOnPtTunePMuonBestTrack->SetLineWidth(2);
  EffErrorOnPtTunePMuonBestTrack->SetMarkerColor(2); 
  EffErrorOnPtTunePMuonBestTrack->SetMarkerStyle(22);
  EffErrorOnPtTunePMuonBestTrack->SetMarkerSize(0.875);

  //========================================================== 
  //                      Dpt/pt cmssw539                          
  //==========================================================
  TFile *f131 = new TFile("ZprimetoMuMu-MC-Mass3000-CMSSW539.root","READ");
  //---------------------------------------------------
  cout<<"====================== dpt/pt CMSSW539 =============================="<<endl;
  int nbEleTopEB103                             = PtID->Integral();
  int nbEleBottomEB103                          = PtEffpterror->Integral();
  float EffChargedHadEB103                      = (float)nbEleTopEB103/nbEleBottomEB103;
  cout<<"nbEleTopEB103                          = "<<nbEleTopEB103<<endl;
  cout<<"nbEleBottomEB103                       = "<<nbEleBottomEB103<<endl;
  cout<<"EffChargedHadEB103(electron)           = "<<EffChargedHadEB103<<endl;
  float StatErrorEB103                          = sqrt(nbEleBottomEB103 - nbEleTopEB103)/nbEleBottomEB103;
  cout<<"StatErrorEB103(electron)               = "<<StatErrorEB103<<endl;
  //---------------------------------------------------
  EffForErrorOnPt539->BayesDivide(PtID,PtEffpterror);
  EffForErrorOnPt539->SetLineStyle(0);
  EffForErrorOnPt539->SetLineColor(4);
  EffForErrorOnPt539->SetLineWidth(2);
  EffForErrorOnPt539->SetMarkerColor(4); 
  EffForErrorOnPt539->SetMarkerStyle(23);
  EffForErrorOnPt539->SetMarkerSize(0.875);
  //-----------------------------------------------------------------
  mg10->Add(EffErrorOnPtMuonBestTrack);
  mg10->Add(EffErrorOnPtTunePMuonBestTrack);
  mg10->Add(EffForErrorOnPt539);
  mg10->Draw("AP");
  mg10->SetTitle("title;pt [GeV]; Efficiency(N-1)");
  mg10->GetXaxis()->SetTitleOffset(1.7);
  mg10->GetYaxis()->SetTitleOffset(1.7);
  mg10->GetXaxis()->SetLabelSize(0.03);
  mg10->GetYaxis()->SetLabelSize(0.03);
  //mg10->GetXaxis()->SetRangeUser(0.0,800.0);
  //mg10->GetYaxis()->SetRangeUser(0.0,1.05);
  mg10->GetYaxis()->SetRangeUser(0.85,1.015);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText2 = new TPaveText(0.40, 0.90, 0.90, 0.92, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t2 = tText2->AddText("DeltaPhi");
  TText *t3 = tText2->AddText("#delta pt/pt<0.30"); 
  tText2->Draw();
  //========================================================== 
  TPaveText* tText1 = new TPaveText(0.40, 0.50, 0.50, 0.55, "brNDC");
  tText1->SetBorderSize(0);
  tText1->SetFillColor(0);
  tText1->SetFillStyle(0);
  Float_t tsize = 0.025;
  tText1->SetTextSize(tsize);
  tText1->SetTextColor(kBlack);
  TText *t1 = tText1->AddText("Eff(MuonBestTrack) = 0.9950 #pm 0.0002");
  tText1->Draw();
  //==========================================================
  TPaveText* tText2 = new TPaveText(0.40, 0.45, 0.50, 0.50, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(2);
  tText2->SetFillStyle(0);
  Float_t tsize = 0.025;
  tText2->SetTextSize(tsize); 
  tText2->SetTextColor(kRed);
  TText *t2 = tText2->AddText("Eff(tunePMuonBestTrack) = 0.9994 #pm 0.0001");
  tText2->Draw();
  //==========================================================
  TPaveText* tText3 = new TPaveText(0.40, 0.40, 0.50, 0.45, "brNDC");
  tText3->SetBorderSize(0);
  tText3->SetFillColor(4);
  tText3->SetFillStyle(0);
  Float_t tsize = 0.025;
  tText3->SetTextSize(tsize);
  tText3->SetTextColor(kBlue);
  TText *t3 = tText3->AddText("Eff(CMSSW539) = 0.9998 #pm 0.0000");
  tText3->Draw();
  //==========================================================
  TLegend *leg = new TLegend(0.60, 0.20, 0.80, 0.35);
  leg->AddEntry(EffErrorOnPtMuonBestTrack,"MuonBestTrack[Phys14]","p");
  leg->AddEntry(EffErrorOnPtTunePMuonBestTrack,"tunePMuonBestTrack[Phys14]","p");
  leg->AddEntry(EffForErrorOnPt539,"tuneP[CMSSW539]","p");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.03;
  leg->SetTextSize(tsize2); 
  leg->Draw();
  //======================================================================= 
  c10->Print("Eff_PtError_Muon_ID_MC_CMSSW720_Et.png","png");
  //c10->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.pdf","pdf");
  //c10->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.eps","eps");
  //=======================================================================
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  TCanvas *c20 = new TCanvas("c20","Track iso",600,600);
  char textpro1[100],textNDF1[100],textRatio1[100];           
  c20->cd(1);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  gPad->SetGridx();
  TGraphAsymmErrors *PtEffForTrackIso     = new TGraphAsymmErrors;
  TGraphAsymmErrors *PhotonEtDeltaPhiEB  = new TGraphAsymmErrors;
  TGraphAsymmErrors *QCDEtDeltaPhiEB     = new TGraphAsymmErrors;
  TMultiGraph *mg20 = new TMultiGraph;
  //========================================================== 
  //                      DeltaPhi in EB                                                  
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  //---------------------------------------------------
  cout<<"====================== track iso =============================="<<endl;
  int nbEleTopEB20                             = PtID->Integral();
  int nbEleBottomEB20                          = PtEffptTrackIso->Integral();
  float EffChargedHadEB20                      = (float)nbEleTopEB20/nbEleBottomEB20;
  cout<<"nbEleTopEB20                          = "<<nbEleTopEB20<<endl;
  cout<<"nbEleBottomEB20                       = "<<nbEleBottomEB20<<endl;
  cout<<"EffChargedHadEB20(electron)           = "<<EffChargedHadEB20<<endl;
  float StatErrorEB20                          = sqrt(nbEleBottomEB20 - nbEleTopEB20)/nbEleBottomEB20;
  cout<<"StatErrorEB20(electron)               = "<<StatErrorEB20<<endl;
  //---------------------------------------------------
  PtEffForTrackIso->BayesDivide(PtID,PtEffptTrackIso);
  PtEffForTrackIso->SetLineStyle(0);
  PtEffForTrackIso->SetLineColor(1);
  PtEffForTrackIso->SetLineWidth(2);
  PtEffForTrackIso->SetMarkerColor(1); 
  PtEffForTrackIso->SetMarkerStyle(20);
  PtEffForTrackIso->SetMarkerSize(0.875);
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
  mg20->Add(PtEffForTrackIso);
  //mg20->Add(PhotonEtDeltaPhiEB);
  //mg20->Add(QCDEtDeltaPhiEB);
  mg20->Draw("AP");
  mg20->SetTitle("title;pt [GeV]; Efficiency(N-1)");
  mg20->GetXaxis()->SetTitleOffset(1.7);
  mg20->GetYaxis()->SetTitleOffset(1.7);
  mg20->GetXaxis()->SetLabelSize(0.03);
  mg20->GetYaxis()->SetLabelSize(0.03);
  //mg20->GetXaxis()->SetRangeUser(0.0,800.0);
  //mg20->GetYaxis()->SetRangeUser(0.0,1.05);
  mg20->GetYaxis()->SetRangeUser(0.0,1.1);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText2 = new TPaveText(0.40, 0.90, 0.90, 0.92, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t2 = tText2->AddText("DeltaPhi");
  TText *t3 = tText2->AddText("TrackPtSum/pt<0.10; CMSSW706 [1/fb]"); 
  tText2->Draw();
  //========================================================== 
  TPaveText* tText2 = new TPaveText(0.20, 0.30, 0.60, 0.45, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t1 = tText2->AddText("Eff(N-1) = 0.9969 #pm 0.0002");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9684 #pm 0.0005");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8921 #pm 0.0057");

  TText *t1 = tText2->AddText("Eff(N-1) = 0.9949 #pm 0.0002");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9198 #pm 0.0008");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8316 #pm 0.0071");
  tText2->Draw();
  //========================================================== 
  TLegend *leg = new TLegend(0.60, 0.20, 0.80, 0.35);
  leg->AddEntry(PtEffForTrackIso,"Z' [M = 1 TeV/c^{2}]","p");
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
  c20->Print("Eff_TrackIso_Muon_ID_MC_CMSSW720_Et.png","png");
  //c20->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.pdf","pdf");
  //c20->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.eps","eps");
  //=======================================================================
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  TCanvas *c30 = new TCanvas("c30","NbTrackLayers",600,600);
  char textpro1[100],textNDF1[100],textRatio1[100];           
  c30->cd(1);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  gPad->SetGridx();
  TGraphAsymmErrors *PtEffForNbTrackLayers     = new TGraphAsymmErrors;
  TGraphAsymmErrors *PhotonEtDeltaPhiEB  = new TGraphAsymmErrors;
  TGraphAsymmErrors *QCDEtDeltaPhiEB     = new TGraphAsymmErrors;
  TMultiGraph *mg30 = new TMultiGraph;
  //========================================================== 
  //                      DeltaPhi in EB                                                  
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  //---------------------------------------------------
  cout<<"====================== NbTrackLayers =============================="<<endl;
  int nbEleTopEB30                             = PtID->Integral();
  int nbEleBottomEB30                          = PtEffptnumberOftrackerLayers->Integral();
  float EffChargedHadEB30                      = (float)nbEleTopEB30/nbEleBottomEB30;
  cout<<"nbEleTopEB30                          = "<<nbEleTopEB30<<endl;
  cout<<"nbEleBottomEB30                       = "<<nbEleBottomEB30<<endl;
  cout<<"EffChargedHadEB30(electron)           = "<<EffChargedHadEB30<<endl;
  float StatErrorEB30                          = sqrt(nbEleBottomEB30 - nbEleTopEB30)/nbEleBottomEB30;
  cout<<"StatErrorEB30(electron)               = "<<StatErrorEB30<<endl;
  //---------------------------------------------------
  PtEffForNbTrackLayers->BayesDivide(PtID,PtEffptnumberOftrackerLayers);
  PtEffForNbTrackLayers->SetLineStyle(0);
  PtEffForNbTrackLayers->SetLineColor(1);
  PtEffForNbTrackLayers->SetLineWidth(2);
  PtEffForNbTrackLayers->SetMarkerColor(1); 
  PtEffForNbTrackLayers->SetMarkerStyle(20);
  PtEffForNbTrackLayers->SetMarkerSize(0.875);
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
  mg30->Add(PtEffForNbTrackLayers);
  //mg30->Add(PhotonEtDeltaPhiEB);
  //mg30->Add(QCDEtDeltaPhiEB);
  mg30->Draw("AP");
  mg30->SetTitle("title;pt [GeV]; Efficiency(N-1)");
  mg30->GetXaxis()->SetTitleOffset(1.7);
  mg30->GetYaxis()->SetTitleOffset(1.7);
  mg30->GetXaxis()->SetLabelSize(0.03);
  mg30->GetYaxis()->SetLabelSize(0.03);
  //mg30->GetXaxis()->SetRangeUser(0.0,800.0);
  //mg30->GetYaxis()->SetRangeUser(0.0,1.05);
  mg30->GetYaxis()->SetRangeUser(0.0,1.1);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText2 = new TPaveText(0.40, 0.90, 0.90, 0.92, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t2 = tText2->AddText("DeltaPhi");
  TText *t3 = tText2->AddText("NbTrackerLayers>5; CMSSW706 [1/fb]"); 
  tText2->Draw();
  //========================================================== 
  TPaveText* tText2 = new TPaveText(0.20, 0.30, 0.60, 0.45, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t1 = tText2->AddText("Eff(N-1) = 0.9969 #pm 0.0002");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9684 #pm 0.0005");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8921 #pm 0.0057");

  TText *t1 = tText2->AddText("Eff(N-1) =0.9999 #pm 0.0000");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9198 #pm 0.0008");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8316 #pm 0.0071");
  tText2->Draw();
  //========================================================== 
  TLegend *leg = new TLegend(0.60, 0.20, 0.80, 0.35);
  leg->AddEntry(PtEffForNbTrackLayers,"Z' [M = 1 TeV/c^{2}]","p");
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
  c30->Print("Eff_NbTrackLayers_Muon_ID_MC_CMSSW720_Et.png","png");
  //c30->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.pdf","pdf");
  //c30->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.eps","eps");
  //=======================================================================
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  TCanvas *c40 = new TCanvas("c40","NbPixelHits",600,600);
  char textpro1[100],textNDF1[100],textRatio1[100];           
  c40->cd(1);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  gPad->SetGridx();
  TGraphAsymmErrors *PtEffForNbPixelHits     = new TGraphAsymmErrors;
  TGraphAsymmErrors *PhotonEtDeltaPhiEB  = new TGraphAsymmErrors;
  TGraphAsymmErrors *QCDEtDeltaPhiEB     = new TGraphAsymmErrors;
  TMultiGraph *mg40 = new TMultiGraph;
  //========================================================== 
  //                      DeltaPhi in EB                                                  
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  //---------------------------------------------------
  cout<<"==================== NbPixelHits =============================="<<endl;
  int nbEleTopEB40                             = PtID->Integral();
  int nbEleBottomEB40                          = PtEffptnumberOfPixelHits->Integral();
  float EffChargedHadEB40                      = (float)nbEleTopEB40/nbEleBottomEB40;
  cout<<"nbEleTopEB40                          = "<<nbEleTopEB40<<endl;
  cout<<"nbEleBottomEB40                       = "<<nbEleBottomEB40<<endl;
  cout<<"EffChargedHadEB40(electron)           = "<<EffChargedHadEB40<<endl;
  float StatErrorEB40                          = sqrt(nbEleBottomEB40 - nbEleTopEB40)/nbEleBottomEB40;
  cout<<"StatErrorEB40(electron)               = "<<StatErrorEB40<<endl;
  //---------------------------------------------------
  PtEffForNbPixelHits->BayesDivide(PtID,PtEffptnumberOfPixelHits);
  PtEffForNbPixelHits->SetLineStyle(0);
  PtEffForNbPixelHits->SetLineColor(1);
  PtEffForNbPixelHits->SetLineWidth(2);
  PtEffForNbPixelHits->SetMarkerColor(1); 
  PtEffForNbPixelHits->SetMarkerStyle(20);
  PtEffForNbPixelHits->SetMarkerSize(0.875);
  /*
  //========================================================== 
  //                 DeltaPhi in EB                                                  
  //==========================================================
  TFile *f2 = new TFile("CMSSW701-Analyse_MC_GJet_pt_4000_V5_withPF.root","READ");
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
  TFile *f3 = new TFile("CMSSW701-Analyse_MC_QCD_15_pt_4000_V5_withPF.root","READ");
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
  mg40->Add(PtEffForNbPixelHits);
  //mg40->Add(PhotonEtDeltaPhiEB);
  //mg40->Add(QCDEtDeltaPhiEB);
  mg40->Draw("AP");
  mg40->SetTitle("title;pt [GeV]; Efficiency(N-1)");
  mg40->GetXaxis()->SetTitleOffset(1.7);
  mg40->GetYaxis()->SetTitleOffset(1.7);
  mg40->GetXaxis()->SetLabelSize(0.03);
  mg40->GetYaxis()->SetLabelSize(0.03);
  //mg40->GetXaxis()->SetRangeUser(0.0,800.0);
  //mg40->GetYaxis()->SetRangeUser(0.0,1.05);
  mg40->GetYaxis()->SetRangeUser(0.0,1.1);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText2 = new TPaveText(0.40, 0.90, 0.90, 0.92, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t2 = tText2->AddText("DeltaPhi");
  TText *t3 = tText2->AddText("NbPixelHits>0; CMSSW706 [1/fb]"); 
  tText2->Draw();
  //========================================================== 
  TPaveText* tText2 = new TPaveText(0.20, 0.30, 0.60, 0.45, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t1 = tText2->AddText("Eff(N-1) = 0.9969 #pm 0.0002");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9684 #pm 0.0005");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8921 #pm 0.0057");

  TText *t1 = tText2->AddText("Eff(N-1) =0.9987 #pm 0.0000");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9198 #pm 0.0008");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8316 #pm 0.0071");
  tText2->Draw();
  //========================================================== 
  TLegend *leg = new TLegend(0.60, 0.20, 0.80, 0.35);
  leg->AddEntry(PtEffForNbPixelHits,"Z' [M = 1 TeV/c^{2}]","p");
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
  c40->Print("Eff_NbPixelHits_Muon_ID_MC_CMSSW720_Et.png","png");
  //c40->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.pdf","pdf");
  //c40->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.eps","eps");
  //=======================================================================
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  TCanvas *c50 = new TCanvas("c50","NbMuonHits",600,600);
  char textpro1[100],textNDF1[100],textRatio1[100];           
  c50->cd(1);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  gPad->SetGridx();
  TGraphAsymmErrors *PtEffForNbMuonHits     = new TGraphAsymmErrors;
  TGraphAsymmErrors *PhotonEtDeltaPhiEB  = new TGraphAsymmErrors;
  TGraphAsymmErrors *QCDEtDeltaPhiEB     = new TGraphAsymmErrors;
  TMultiGraph *mg50 = new TMultiGraph;
  //========================================================== 
  //                      DeltaPhi in EB                                                  
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  //---------------------------------------------------
  cout<<"==================== NbMuonHits =============================="<<endl;
  int nbEleTopEB50                             = PtID->Integral();
  int nbEleBottomEB50                          = PtEffptnumberOfMuonHits->Integral();
  float EffChargedHadEB50                      = (float)nbEleTopEB50/nbEleBottomEB50;
  cout<<"nbEleTopEB50                          = "<<nbEleTopEB50<<endl;
  cout<<"nbEleBottomEB50                       = "<<nbEleBottomEB50<<endl;
  cout<<"EffChargedHadEB50(electron)           = "<<EffChargedHadEB50<<endl;
  float StatErrorEB50                          = sqrt(nbEleBottomEB50 - nbEleTopEB50)/nbEleBottomEB50;
  cout<<"StatErrorEB50(electron)               = "<<StatErrorEB50<<endl;
  //---------------------------------------------------
  PtEffForNbMuonHits->BayesDivide(PtID,PtEffptnumberOfMuonHits);
  PtEffForNbMuonHits->SetLineStyle(0);
  PtEffForNbMuonHits->SetLineColor(1);
  PtEffForNbMuonHits->SetLineWidth(2);
  PtEffForNbMuonHits->SetMarkerColor(1); 
  PtEffForNbMuonHits->SetMarkerStyle(20);
  PtEffForNbMuonHits->SetMarkerSize(0.875);
  /*
  //========================================================== 
  //                 DeltaPhi in EB                                                  
  //==========================================================
  TFile *f2 = new TFile("CMSSW701-Analyse_MC_GJet_pt_5000_V5_withPF.root","READ");
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
  TFile *f3 = new TFile("CMSSW701-Analyse_MC_QCD_15_pt_5000_V5_withPF.root","READ");
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
  mg50->Add(PtEffForNbMuonHits);
  //mg50->Add(PhotonEtDeltaPhiEB);
  //mg50->Add(QCDEtDeltaPhiEB);
  mg50->Draw("AP");
  mg50->SetTitle("title;pt [GeV]; Efficiency(N-1)");
  mg50->GetXaxis()->SetTitleOffset(1.7);
  mg50->GetYaxis()->SetTitleOffset(1.7);
  mg50->GetXaxis()->SetLabelSize(0.03);
  mg50->GetYaxis()->SetLabelSize(0.03);
  //mg50->GetXaxis()->SetRangeUser(0.0,800.0);
  //mg50->GetYaxis()->SetRangeUser(0.0,1.05);
  mg50->GetYaxis()->SetRangeUser(0.0,1.1);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText2 = new TPaveText(0.40, 0.90, 0.90, 0.92, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t2 = tText2->AddText("DeltaPhi");
  TText *t3 = tText2->AddText("NbMuonHits>0; CMSSW706 [1/fb]"); 
  tText2->Draw();
  //========================================================== 
  TPaveText* tText2 = new TPaveText(0.20, 0.30, 0.60, 0.45, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t1 = tText2->AddText("Eff(N-1) = 0.9969 #pm 0.0002");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9684 #pm 0.0005");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8921 #pm 0.0057");

  TText *t1 = tText2->AddText("Eff(N-1) = 0.9983 #pm 0.0001");
  //TText *t2 = tText2->AddText("Eff(Photons) = 0.9198 #pm 0.0008");
  //TText *t3 = tText2->AddText("Eff(QCD) = 0.8316 #pm 0.0071");
  tText2->Draw();
  //========================================================== 
  TLegend *leg = new TLegend(0.60, 0.20, 0.80, 0.35);
  leg->AddEntry(PtEffForNbMuonHits,"Z' [M = 1 TeV/c^{2}]","p");
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
  c50->Print("Eff_NbMuonHits_Muon_ID_MC_CMSSW720_Et.png","png");
  //c50->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.pdf","pdf");
  //c50->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.eps","eps");
  //=======================================================================
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  TCanvas *c60 = new TCanvas("c60","NbMatchedStations",600,600);
  char textpro1[100],textNDF1[100],textRatio1[100];           
  c60->cd(1);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  gPad->SetGridx();
  TGraphAsymmErrors *PtEffForNbMatchedStations720  = new TGraphAsymmErrors;
  TGraphAsymmErrors *PtEffForNbMatchedStations706  = new TGraphAsymmErrors;
  TGraphAsymmErrors *PtEffForNbMatchedStations539  = new TGraphAsymmErrors;
  TMultiGraph *mg60 = new TMultiGraph;
  //========================================================== 
  //                      CMSSW720                                                  
  //==========================================================
  TFile *f221 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  //---------------------------------------------------
  cout<<"==================== NbMatchedStations =============================="<<endl;
  int nbEleTopEB601                             = PtID->Integral();
  int nbEleBottomEB601                          = PtEffptnumberOfMatchedStations->Integral();
  float EffChargedHadEB601                      = (float)nbEleTopEB601/nbEleBottomEB601;
  cout<<"nbEleTopEB601                          = "<<nbEleTopEB601<<endl;
  cout<<"nbEleBottomEB601                       = "<<nbEleBottomEB601<<endl;
  cout<<"EffChargedHadEB601(electron)           = "<<EffChargedHadEB601<<endl;
  float StatErrorEB601                          = sqrt(nbEleBottomEB601 - nbEleTopEB601)/nbEleBottomEB601;
  cout<<"StatErrorEB601(electron)               = "<<StatErrorEB601<<endl;
  //---------------------------------------------------
  PtEffForNbMatchedStations720->BayesDivide(PtID,PtEffptnumberOfMatchedStations);
  PtEffForNbMatchedStations720->SetLineStyle(0);
  PtEffForNbMatchedStations720->SetLineColor(1);
  PtEffForNbMatchedStations720->SetLineWidth(2);
  PtEffForNbMatchedStations720->SetMarkerColor(1); 
  PtEffForNbMatchedStations720->SetMarkerStyle(20);
  PtEffForNbMatchedStations720->SetMarkerSize(0.875);
  //========================================================== 
  //                      CMSSW706                                                  
  //==========================================================
  TFile *f222 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW706.root","READ");
  //---------------------------------------------------
  cout<<"==================== NbMatchedStations =============================="<<endl;
  int nbEleTopEB602                             = PtID->Integral();
  int nbEleBottomEB602                          = PtEffptnumberOfMatchedStations->Integral();
  float EffChargedHadEB602                      = (float)nbEleTopEB602/nbEleBottomEB602;
  cout<<"nbEleTopEB602                          = "<<nbEleTopEB602<<endl;
  cout<<"nbEleBottomEB602                       = "<<nbEleBottomEB602<<endl;
  cout<<"EffChargedHadEB602(electron)           = "<<EffChargedHadEB602<<endl;
  float StatErrorEB602                          = sqrt(nbEleBottomEB602 - nbEleTopEB602)/nbEleBottomEB602;
  cout<<"StatErrorEB602(electron)               = "<<StatErrorEB602<<endl;
  //---------------------------------------------------
  PtEffForNbMatchedStations706->BayesDivide(PtID,PtEffptnumberOfMatchedStations);
  PtEffForNbMatchedStations706->SetLineStyle(0);
  PtEffForNbMatchedStations706->SetLineColor(2);
  PtEffForNbMatchedStations706->SetLineWidth(2);
  PtEffForNbMatchedStations706->SetMarkerColor(2); 
  PtEffForNbMatchedStations706->SetMarkerStyle(22);
  PtEffForNbMatchedStations706->SetMarkerSize(0.875);
  //========================================================== 
  //                      CMSSW539                                                  
  //==========================================================
  TFile *f223 = new TFile("ZprimetoMuMu-MC-Mass3000-CMSSW539.root","READ");
  //---------------------------------------------------
  cout<<"==================== NbMatchedStations =============================="<<endl;
  int nbEleTopEB603                             = PtID->Integral();
  int nbEleBottomEB603                          = PtEffptnumberOfMatchedStations->Integral();
  float EffChargedHadEB603                      = (float)nbEleTopEB603/nbEleBottomEB603;
  cout<<"nbEleTopEB603                          = "<<nbEleTopEB603<<endl;
  cout<<"nbEleBottomEB603                       = "<<nbEleBottomEB603<<endl;
  cout<<"EffChargedHadEB603(electron)           = "<<EffChargedHadEB603<<endl;
  float StatErrorEB603                          = sqrt(nbEleBottomEB603 - nbEleTopEB603)/nbEleBottomEB603;
  cout<<"StatErrorEB603(electron)               = "<<StatErrorEB603<<endl;
  //---------------------------------------------------
  PtEffForNbMatchedStations539->BayesDivide(PtID,PtEffptnumberOfMatchedStations);
  PtEffForNbMatchedStations539->SetLineStyle(0);
  PtEffForNbMatchedStations539->SetLineColor(4);
  PtEffForNbMatchedStations539->SetLineWidth(2);
  PtEffForNbMatchedStations539->SetMarkerColor(4); 
  PtEffForNbMatchedStations539->SetMarkerStyle(23);
  PtEffForNbMatchedStations539->SetMarkerSize(0.875);
  //==========================================================================
  mg60->Add(PtEffForNbMatchedStations720);
  mg60->Add(PtEffForNbMatchedStations706);
  mg60->Add(PtEffForNbMatchedStations539);
  mg60->Draw("AP");
  mg60->SetTitle("title;pt [GeV]; Efficiency(N-1)");
  mg60->GetXaxis()->SetTitleOffset(1.7);
  mg60->GetYaxis()->SetTitleOffset(1.7);
  mg60->GetXaxis()->SetLabelSize(0.03);
  mg60->GetYaxis()->SetLabelSize(0.03);
  //mg60->GetXaxis()->SetRangeUser(0.0,800.0);
  //mg60->GetYaxis()->SetRangeUser(0.0,1.05);
  mg60->GetYaxis()->SetRangeUser(0.85,1.015);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText2 = new TPaveText(0.40, 0.90, 0.90, 0.92, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  TText *t3 = tText2->AddText("NbMatchedStations>1"); 
  tText2->Draw();
  //========================================================== 
  TPaveText* tText2 = new TPaveText(0.20, 0.30, 0.60, 0.45, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  TText *t1 = tText2->AddText("Eff(phys14) = 0.9885 #pm 0.0003");
  TText *t2 = tText2->AddText("Eff(CSA14) = 0.9883 #pm 0.0003");
  TText *t3 = tText2->AddText("Eff(CMSSW539) = 0.9849 #pm 0.0006");
  tText2->Draw();
  //========================================================== 
  TLegend *leg = new TLegend(0.60, 0.20, 0.80, 0.35);
  leg->AddEntry(PtEffForNbMatchedStations720,"CMSSW720[Phys14]","p");
  leg->AddEntry(PtEffForNbMatchedStations706,"CMSSW706[CSA14]","p");
  leg->AddEntry(PtEffForNbMatchedStations539,"CMSSW539","p");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.03;
  leg->SetTextSize(tsize2); 
  leg->Draw();
  //======================================================================= 
  c60->Print("Eff_NbMatchedStations_Muon_ID_MC_CMSSW720_Et.png","png");
  //c60->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.pdf","pdf");
  //c60->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.eps","eps");
  //=======================================================================
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  TCanvas *c70 = new TCanvas("c70","Dxy",600,600);
  char textpro1[100],textNDF1[100],textRatio1[100];           
  c70->cd(1);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  gPad->SetGridx();
  TGraphAsymmErrors *PtEffForDxy      = new TGraphAsymmErrors;
  TGraphAsymmErrors *PtEffForNewDxy   = new TGraphAsymmErrors;
  TGraphAsymmErrors *PtEffForTightDxy = new TGraphAsymmErrors;
  TMultiGraph *mg70 = new TMultiGraph;
  //========================================================== 
  //                          |dxy| (old)                                            
  //==========================================================
  TFile *f5555 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  //---------------------------------------------------
  cout<<"==================== |dxy| (old) =============================="<<endl;
  int nbEleTopEB70                             = PtID->Integral();
  int nbEleBottomEB70                          = PtEffptabsdsy->Integral();
  float EffChargedHadEB70                      = (float)nbEleTopEB70/nbEleBottomEB70;
  cout<<"nbEleTopEB70                          = "<<nbEleTopEB70<<endl;
  cout<<"nbEleBottomEB70                       = "<<nbEleBottomEB70<<endl;
  cout<<"EffChargedHadEB70(old)           = "<<EffChargedHadEB70<<endl;
  float StatErrorEB70                          = sqrt(nbEleBottomEB70 - nbEleTopEB70)/nbEleBottomEB70;
  cout<<"StatErrorEB70(old)               = "<<StatErrorEB70<<endl;
  //---------------------------------------------------
  PtEffForDxy->BayesDivide(PtID,PtEffptabsdsy);
  PtEffForDxy->SetLineStyle(0);
  PtEffForDxy->SetLineColor(1);
  PtEffForDxy->SetLineWidth(2);
  PtEffForDxy->SetMarkerColor(1); 
  PtEffForDxy->SetMarkerStyle(20);
  PtEffForDxy->SetMarkerSize(0.875);
  //========================================================== 
  //                          |dxy| (new)                                            
  //==========================================================
  TFile *f6666 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  //---------------------------------------------------
  cout<<"==================== |dxy| (new) =============================="<<endl;
  int nbEleTopEB77                             = PtNewID->Integral();
  int nbEleBottomEB77                          = PtEffptabsdsy->Integral();
  float EffChargedHadEB77                      = (float)nbEleTopEB77/nbEleBottomEB77;
  cout<<"nbEleTopEB77                          = "<<nbEleTopEB77<<endl;
  cout<<"nbEleBottomEB77                       = "<<nbEleBottomEB77<<endl;
  cout<<"EffChargedHadEB77(new)                = "<<EffChargedHadEB77<<endl;
  float StatErrorEB77                          = sqrt(nbEleBottomEB77 - nbEleTopEB77)/nbEleBottomEB77;
  cout<<"StatErrorEB77(new)                    = "<<StatErrorEB77<<endl;
  //---------------------------------------------------
  PtEffForNewDxy->BayesDivide(PtNewID,PtEffptabsdsy);
  PtEffForNewDxy->SetLineStyle(0);
  PtEffForNewDxy->SetLineWidth(2);
  PtEffForNewDxy->SetLineColor(2);
  PtEffForNewDxy->SetMarkerColor(2); 
  PtEffForNewDxy->SetMarkerStyle(21);
  PtEffForNewDxy->SetMarkerSize(0.875);
  //========================================================== 
  //                          |dxy| (tight)                                            
  //==========================================================
  TFile *f6666 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  //---------------------------------------------------
  cout<<"==================== |dxy| (tight) =============================="<<endl;
  int nbEleTopEB88                             = PtTightID->Integral();
  int nbEleBottomEB88                          = PtEffptabsdsy->Integral();
  float EffChargedHadEB88                      = (float)nbEleTopEB88/nbEleBottomEB88;
  cout<<"nbEleTopEB88                          = "<<nbEleTopEB88<<endl;
  cout<<"nbEleBottomEB88                       = "<<nbEleBottomEB88<<endl;
  cout<<"EffChargedHadEB88(new)                = "<<EffChargedHadEB88<<endl;
  float StatErrorEB88                          = sqrt(nbEleBottomEB88 - nbEleTopEB88)/nbEleBottomEB88;
  cout<<"StatErrorEB88(new)                    = "<<StatErrorEB88<<endl;
  //---------------------------------------------------
  PtEffForTightDxy->BayesDivide(PtTightID,PtEffptabsdsy);
  PtEffForTightDxy->SetLineStyle(0);
  PtEffForTightDxy->SetLineWidth(2);
  PtEffForTightDxy->SetLineColor(4);
  PtEffForTightDxy->SetMarkerColor(4); 
  PtEffForTightDxy->SetMarkerStyle(22);
  PtEffForTightDxy->SetMarkerSize(0.875);
  //==============================================================================                       
  mg70->Add(PtEffForDxy);
  mg70->Add(PtEffForNewDxy);
  mg70->Add(PtEffForTightDxy);
  mg70->Draw("AP");
  mg70->SetTitle("title;pt [GeV]; Efficiency(N-1)");
  mg70->GetXaxis()->SetTitleOffset(1.7);
  mg70->GetYaxis()->SetTitleOffset(1.7);
  mg70->GetXaxis()->SetLabelSize(0.03);
  mg70->GetYaxis()->SetLabelSize(0.03);
  mg70->GetYaxis()->SetRangeUser(0.98,1.002);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText2 = new TPaveText(0.40, 0.90, 0.90, 0.92, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  TText *t3 = tText2->AddText("Z'(M = 5 TeV); CMSSW720 [Phys14]"); 
  tText2->Draw();
  //========================================================== 
  TPaveText* tText2 = new TPaveText(0.20, 0.20, 0.50, 0.60, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  Float_t tsize = 0.030;
  tText2->SetTextSize(tsize); 
  TText *t1 = tText2->AddText("Eff(|dxy|<0.2) = 1.0 #pm 0.0");
  TText *t2 = tText2->AddText("Eff(|dxy|<0.02) = 0.99987 #pm 0.0");
  TText *t3 = tText2->AddText("Eff(|dxy|<0.01) = 0.99944 #pm 0.0");
  tText2->Draw();
  //========================================================== 
  TLegend *leg = new TLegend(0.60, 0.20, 0.80, 0.50);
  leg->AddEntry(PtEffForDxy,"|dxy|<0.2","p");
  leg->AddEntry(PtEffForNewDxy,"|dxy|<0.02","p");  
  leg->AddEntry(PtEffForTightDxy,"|dxy|<0.01","p");  
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.03;
  leg->SetTextSize(tsize2); 
  leg->Draw();
  //======================================================================= 
  c70->Print("Eff_Dxy_Muon_ID_MC_CMSSW720_Et.png","png");
  //c70->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.pdf","pdf");
  //c70->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.eps","eps");
  //=======================================================================
  //---------------------------------------------------------------------------
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  TCanvas *c99 = new TCanvas("c99","All vs Pt",600,600);
  char textpro1[100],textNDF1[100],textRatio1[100];           
  c99->cd(1);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  gPad->SetGridx();
  TGraphAsymmErrors *PtEffForErrorOnPt     = new TGraphAsymmErrors;
  TGraphAsymmErrors *Effdxy                = new TGraphAsymmErrors;
  TGraphAsymmErrors *PtEffForTrackIso      = new TGraphAsymmErrors;
  TGraphAsymmErrors *PtEffForNbTrackLayers = new TGraphAsymmErrors;
  TGraphAsymmErrors *PtEffForNbPixelHits   = new TGraphAsymmErrors;
  TGraphAsymmErrors *PtEffForNbMuonHits    = new TGraphAsymmErrors;
  TGraphAsymmErrors *PtEffForNbMatchedStations = new TGraphAsymmErrors;
  TMultiGraph *mg1000 = new TMultiGraph;
  //========================================================== 
  //                 track isolation                                                 
  //==========================================================
  TFile *f111 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720_miniaod.root","READ");
  cout<<"====================== track iso =============================="<<endl;
  int nbEleTopEB20                             = PtID->Integral();
  int nbEleBottomEB20                          = PtEffptTrackIso->Integral();
  float EffChargedHadEB20                      = (float)nbEleTopEB20/nbEleBottomEB20;
  cout<<"nbEleTopEB20                          = "<<nbEleTopEB20<<endl;
  cout<<"nbEleBottomEB20                       = "<<nbEleBottomEB20<<endl;
  cout<<"EffChargedHadEB20(electron)           = "<<EffChargedHadEB20<<endl;
  float StatErrorEB20                          = sqrt(nbEleBottomEB20 - nbEleTopEB20)/nbEleBottomEB20;
  cout<<"StatErrorEB20(electron)               = "<<StatErrorEB20<<endl;
  //---------------------------------------------------
  PtEffForTrackIso->BayesDivide(PtID,PtEffptTrackIso);
  PtEffForTrackIso->SetLineStyle(0);
  PtEffForTrackIso->SetLineWidth(2);
  PtEffForTrackIso->SetLineColor(3);
  PtEffForTrackIso->SetMarkerColor(3); 
  PtEffForTrackIso->SetMarkerStyle(21);
  PtEffForTrackIso->SetMarkerSize(0.875);
  //---------------------------------------------------
  //========================================================== 
  //                      DeltaPt/Pt                                                  
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  cout<<"====================== dpt/pt =============================="<<endl;
  int nbEleTopEB10                             = PtID->Integral();
  int nbEleBottomEB10                          = PtEffpterror->Integral();
  float EffChargedHadEB10                      = (float)nbEleTopEB10/nbEleBottomEB10;
  cout<<"nbEleTopEB10                          = "<<nbEleTopEB10<<endl;
  cout<<"nbEleBottomEB10                       = "<<nbEleBottomEB10<<endl;
  cout<<"EffChargedHadEB10(electron)           = "<<EffChargedHadEB10<<endl;
  float StatErrorEB10                          = sqrt(nbEleBottomEB10 - nbEleTopEB10)/nbEleBottomEB10;
  cout<<"StatErrorEB10(electron)               = "<<StatErrorEB10<<endl;
  //---------------------------------------------------
  PtEffForErrorOnPt->BayesDivide(PtID,PtEffpterror);
  PtEffForErrorOnPt->SetLineStyle(0);
  PtEffForErrorOnPt->SetLineColor(1);
  PtEffForErrorOnPt->SetLineWidth(2);
  PtEffForErrorOnPt->SetMarkerColor(1); 
  PtEffForErrorOnPt->SetMarkerStyle(20);
  PtEffForErrorOnPt->SetMarkerSize(0.875);
  //========================================================== 
  //                 Dxy                                                 
  //==========================================================
  cout<<"==================== |dxy| =============================="<<endl;
  int nbEleTopEB70                             = PtID->Integral();
  int nbEleBottomEB70                          = PtEffptabsdsy->Integral();
  float EffChargedHadEB70                      = (float)nbEleTopEB70/nbEleBottomEB70;
  cout<<"nbEleTopEB70                          = "<<nbEleTopEB70<<endl;
  cout<<"nbEleBottomEB70                       = "<<nbEleBottomEB70<<endl;
  cout<<"EffChargedHadEB70(electron)           = "<<EffChargedHadEB70<<endl;
  float StatErrorEB70                          = sqrt(nbEleBottomEB70 - nbEleTopEB70)/nbEleBottomEB70;
  cout<<"StatErrorEB70(electron)               = "<<StatErrorEB70<<endl;
  //---------------------------------------------------
  Effdxy->BayesDivide(PtID,PtEffptabsdsy);
  Effdxy->SetLineStyle(0);
  Effdxy->SetLineColor(2);
  Effdxy->SetLineWidth(2);
  Effdxy->SetMarkerColor(2); 
  Effdxy->SetMarkerStyle(20);
  Effdxy->SetMarkerSize(0.875);
  //========================================================== 
  //                              NbTrackLayers                                 
  //==========================================================
  cout<<"====================== NbTrackLayers =============================="<<endl;
  int nbEleTopEB30                             = PtID->Integral();
  int nbEleBottomEB30                          = PtEffptnumberOftrackerLayers->Integral();
  float EffChargedHadEB30                      = (float)nbEleTopEB30/nbEleBottomEB30;
  cout<<"nbEleTopEB30                          = "<<nbEleTopEB30<<endl;
  cout<<"nbEleBottomEB30                       = "<<nbEleBottomEB30<<endl;
  cout<<"EffChargedHadEB30(electron)           = "<<EffChargedHadEB30<<endl;
  float StatErrorEB30                          = sqrt(nbEleBottomEB30 - nbEleTopEB30)/nbEleBottomEB30;
  cout<<"StatErrorEB30(electron)               = "<<StatErrorEB30<<endl;
  //---------------------------------------------------
  PtEffForNbTrackLayers->BayesDivide(PtID,PtEffptnumberOftrackerLayers);
  PtEffForNbTrackLayers->SetLineStyle(0);
  PtEffForNbTrackLayers->SetLineColor(4);
  PtEffForNbTrackLayers->SetLineWidth(2);
  PtEffForNbTrackLayers->SetMarkerColor(4); 
  PtEffForNbTrackLayers->SetMarkerStyle(22);
  PtEffForNbTrackLayers->SetMarkerSize(0.875);
  //========================================================== 
  //                      NbPixelHits                                                  
  //==========================================================
  cout<<"==================== NbPixelHits =============================="<<endl;
  int nbEleTopEB40                             = PtID->Integral();
  int nbEleBottomEB40                          = PtEffptnumberOfPixelHits->Integral();
  float EffChargedHadEB40                      = (float)nbEleTopEB40/nbEleBottomEB40;
  cout<<"nbEleTopEB40                          = "<<nbEleTopEB40<<endl;
  cout<<"nbEleBottomEB40                       = "<<nbEleBottomEB40<<endl;
  cout<<"EffChargedHadEB40(electron)           = "<<EffChargedHadEB40<<endl;
  float StatErrorEB40                          = sqrt(nbEleBottomEB40 - nbEleTopEB40)/nbEleBottomEB40;
  cout<<"StatErrorEB40(electron)               = "<<StatErrorEB40<<endl;
  //---------------------------------------------------
  PtEffForNbPixelHits->BayesDivide(PtID,PtEffptnumberOfPixelHits);
  PtEffForNbPixelHits->SetLineStyle(0);
  PtEffForNbPixelHits->SetLineColor(6);
  PtEffForNbPixelHits->SetLineWidth(2);
  PtEffForNbPixelHits->SetMarkerColor(6); 
  PtEffForNbPixelHits->SetMarkerStyle(23);
  PtEffForNbPixelHits->SetMarkerSize(0.875);
  //========================================================== 
  //                      NbMuonHits
  //==========================================================
  cout<<"==================== NbMuonHits =============================="<<endl;
  int nbEleTopEB50                             = PtID->Integral();
  int nbEleBottomEB50                          = PtEffptnumberOfMuonHits->Integral();
  float EffChargedHadEB50                      = (float)nbEleTopEB50/nbEleBottomEB50;
  cout<<"nbEleTopEB50                          = "<<nbEleTopEB50<<endl;
  cout<<"nbEleBottomEB50                       = "<<nbEleBottomEB50<<endl;
  cout<<"EffChargedHadEB50(electron)           = "<<EffChargedHadEB50<<endl;
  float StatErrorEB50                          = sqrt(nbEleBottomEB50 - nbEleTopEB50)/nbEleBottomEB50;
  cout<<"StatErrorEB50(electron)               = "<<StatErrorEB50<<endl;
  //---------------------------------------------------
  PtEffForNbMuonHits->BayesDivide(PtID,PtEffptnumberOfMuonHits);
  PtEffForNbMuonHits->SetLineStyle(0);
  PtEffForNbMuonHits->SetLineColor(7);
  PtEffForNbMuonHits->SetLineWidth(2);
  PtEffForNbMuonHits->SetMarkerColor(7); 
  PtEffForNbMuonHits->SetMarkerStyle(25);
  PtEffForNbMuonHits->SetMarkerSize(0.875);
  //========================================================== 
  //                      NbMatchedStations
  //==========================================================
  cout<<"==================== NbMatchedStations =============================="<<endl;
  int nbEleTopEB60                             = PtID->Integral();
  int nbEleBottomEB60                          = PtEffptnumberOfMatchedStations->Integral();
  float EffChargedHadEB60                      = (float)nbEleTopEB60/nbEleBottomEB60;
  cout<<"nbEleTopEB60                          = "<<nbEleTopEB60<<endl;
  cout<<"nbEleBottomEB60                       = "<<nbEleBottomEB60<<endl;
  cout<<"EffChargedHadEB60(electron)           = "<<EffChargedHadEB60<<endl;
  float StatErrorEB60                          = sqrt(nbEleBottomEB60 - nbEleTopEB60)/nbEleBottomEB60;
  cout<<"StatErrorEB60(electron)               = "<<StatErrorEB60<<endl;
  //---------------------------------------------------
  PtEffForNbMatchedStations->BayesDivide(PtID,PtEffptnumberOfMatchedStations);
  PtEffForNbMatchedStations->SetLineStyle(0);
  PtEffForNbMatchedStations->SetLineColor(9);
  PtEffForNbMatchedStations->SetLineWidth(2);
  PtEffForNbMatchedStations->SetMarkerColor(9); 
  PtEffForNbMatchedStations->SetMarkerStyle(28);
  PtEffForNbMatchedStations->SetMarkerSize(0.875);
  //----------------------------------------------------------------- 
  mg1000->Add(PtEffForErrorOnPt);
  mg1000->Add(Effdxy);
  mg1000->Add(PtEffForTrackIso);
  mg1000->Add(PtEffForNbTrackLayers);
  mg1000->Add(PtEffForNbPixelHits);
  mg1000->Add(PtEffForNbMuonHits);
  mg1000->Add(PtEffForNbMatchedStations);
  //mg1000->Add();
  //mg1000->Add();

  mg1000->Draw("AP");
  mg1000->SetTitle("title;pt(tunePMuonBestTrack) [GeV]; Efficiency(N-1)");
  mg1000->GetXaxis()->SetTitleOffset(1.7);
  mg1000->GetYaxis()->SetTitleOffset(1.7);
  mg1000->GetXaxis()->SetLabelSize(0.03);
  mg1000->GetYaxis()->SetLabelSize(0.03);
  mg1000->GetXaxis()->SetRangeUser(0.0,3000.0);
  mg1000->GetYaxis()->SetRangeUser(0.85,1.015);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText2 = new TPaveText(0.40, 0.90, 0.90, 0.92, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t2 = tText2->AddText("DeltaPhi");
  TText *t3 = tText2->AddText("Z'(M = 5 TeV); CMSSW720 [Phys14]"); 
  tText2->Draw();
  //========================================================== 
  TPaveText* tText2 = new TPaveText(0.20, 0.20, 0.50, 0.60, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  Float_t tsize = 0.025;
  tText2->SetTextSize(tsize); 
  TText *t1 = tText2->AddText("Eff(|dxy|) =1.0 #pm 0.0");
  TText *t2 = tText2->AddText("Eff(NbTrackLayers) = 0.9997 #pm 0.000");
  TText *t3 = tText2->AddText("Eff(dpt/pt) = 0.9994 #pm 0.0001");
  TText *t4 = tText2->AddText("Eff(NbMuonHits) = 0.9986 #pm 0.0000");
  TText *t5 = tText2->AddText("Eff(NbPixelHits) = 0.9980 #pm 0.0001");
  TText *t6 = tText2->AddText("Eff(track iso) = 0.9950 #pm 0.0002");
  TText *t7 = tText2->AddText("Eff(NbMatchedStations) = 0.9885 #pm 0.0003");
  tText2->Draw();
  //========================================================== 
  TLegend *leg = new TLegend(0.60, 0.20, 0.80, 0.50);
  leg->AddEntry(Effdxy,"|dxy|","p");
  leg->AddEntry(PtEffForNbTrackLayers,"NbTrackLayers","p");  
  leg->AddEntry(PtEffForErrorOnPt,"dpt/pt","p");
  leg->AddEntry(PtEffForNbMuonHits,"NbMuonHits","p");
  leg->AddEntry(PtEffForNbPixelHits,"NbPixelHits","p");
  leg->AddEntry(PtEffForTrackIso,"track iso","p");
  leg->AddEntry(PtEffForNbMatchedStations,"NbMatchedStations","p");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.03;
  leg->SetTextSize(tsize2); 
  leg->Draw();
  //======================================================================= 
  c99->Print("Eff_All_Muon_ID_MC_CMSSW720_Et.png","png");
  //c99->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.pdf","pdf");
  //c99->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.eps","eps");
  //=======================================================================
  //==================================================================  
  //================================================================== 
  //==================================================================  
  //================================================================== 
  TCanvas *c333 = new TCanvas("c333","All vs Eta",600,600);
  char textpro1[100],textNDF1[100],textRatio1[100];           
  c333->cd(1);
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
  //                 track isolation                                                 
  //==========================================================
  TFile *f3123 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720_miniaod.root","READ");
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
  EtaEffForTrackIso->SetLineColor(3);
  EtaEffForTrackIso->SetLineWidth(2);
  EtaEffForTrackIso->SetMarkerColor(3); 
  EtaEffForTrackIso->SetMarkerStyle(21);
  EtaEffForTrackIso->SetMarkerSize(0.875);
  //========================================================== 
  //                      DeltaPt/Pt                                                  
  //==========================================================
  TFile *f11113 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
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
  //                              NbTrackLayers                                 
  //==========================================================
  cout<<"====================== NbTrackLayers =============================="<<endl;
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
  mg1000->SetTitle("title;|#eta|_{tunePMuonBestTrack}; Efficiency(N-1)");
  mg1000->GetXaxis()->SetTitleOffset(1.7);
  mg1000->GetYaxis()->SetTitleOffset(1.7);
  mg1000->GetXaxis()->SetLabelSize(0.03);
  mg1000->GetYaxis()->SetLabelSize(0.03);
  mg1000->GetXaxis()->SetRangeUser(0.0,2.3);
  mg1000->GetYaxis()->SetRangeUser(0.85,1.015);
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText2 = new TPaveText(0.40, 0.90, 0.90, 0.92, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  Float_t tsize = 0.03;
  tText2->SetTextSize(tsize); 
  //TText *t2 = tText2->AddText("DeltaPhi");
  TText *t3 = tText2->AddText("Z'(M = 5 TeV); CMSSW720 [Phys14]"); 
  tText2->Draw();
  //========================================================== 
  TPaveText* tText2 = new TPaveText(0.20, 0.20, 0.50, 0.60, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  Float_t tsize = 0.025;
  tText2->SetTextSize(tsize); 
  TText *t1 = tText2->AddText("Eff(|dxy|) =1.0 #pm 0.0");
  TText *t2 = tText2->AddText("Eff(NbTrackLayers) = 0.9997 #pm 0.000");
  TText *t3 = tText2->AddText("Eff(dpt/pt) = 0.9994 #pm 0.0001");
  TText *t4 = tText2->AddText("Eff(NbMuonHits) = 0.9986 #pm 0.0000");
  TText *t5 = tText2->AddText("Eff(NbPixelHits) = 0.9980 #pm 0.0001");
  TText *t6 = tText2->AddText("Eff(track iso) = 0.9950 #pm 0.0002");
  TText *t7 = tText2->AddText("Eff(NbMatchedStations) = 0.9885 #pm 0.0003");
  tText2->Draw();
  //========================================================== 
  TLegend *leg = new TLegend(0.60, 0.20, 0.80, 0.50);
  leg->AddEntry(Effdxy,"|dxy|","p");
  leg->AddEntry(EtaEffForNbTrackLayers,"NbTrackLayers","p");  
  leg->AddEntry(EtaEffForErrorOnPt,"dpt/pt","p");
  leg->AddEntry(EtaEffForNbMuonHits,"NbMuonHits","p");
  leg->AddEntry(EtaEffForNbPixelHits,"NbPixelHits","p");
  leg->AddEntry(EtaEffForTrackIso,"track iso","p");
  leg->AddEntry(EtaEffForNbMatchedStations,"NbMatchedStations","p");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.03;
  leg->SetTextSize(tsize2); 
  leg->Draw();
  //======================================================================= 
  c333->Print("Eff_All_Muon_ID_MC_CMSSW720_Eta.png","png");
  //c333->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.pdf","pdf");
  //c333->Print("PlotsDir4/Eff_DeltaPhi_HEEP_ID_MC_CMSSW701_Eta_EBEE.eps","eps");
  //=======================================================================






  





















}





