#include "TProfile.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TF2.h"
#include <string>
using std::string;
void Macro_MVA_MLP_Fit_EE(){
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
  float FitMin      = -0.07;
  float FitMax      = 0.07;
  //float FitMin      = -0.2;
  //float FitMax      = 0.2;
  float norm = 1.0;
  //==================================================================  
  //==================================================================  
  //================================================================== 
  //==================================================================
  TCanvas *c111 = new TCanvas("c111","MassResolution",600,600);
  c111->SetFillColor(0);
  TH1F *MassResolution1 = new TH1F("MassResolution1","",100.0,-0.5,0.5);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx();
  gPad->SetGridy();
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  MassResolution1->Add(MassRecoGenDif);
  MassResolution1->SetLineStyle(0);
  MassResolution1->SetLineColor(4);
  MassResolution1->SetLineWidth(2);
  MassResolution1->SetMarkerColor(4); 
  //MassResolution->SetMarkerStyle(20);
  MassResolution1->SetMarkerStyle(2);
  //MassResolution->SetMarkerSize(0.70);
  MassResolution1->SetTitle("");
  MassResolution1->GetXaxis()->SetTitle("(dil. mass - gen dil. mass)/(gen dil. mass)");
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
  TText *t3 = tText2->AddText("Phys14 (miniaod)");
  tText2->Draw(); 
  //==========================================================
  TLegend *leg = new TLegend(0.60, 0.60, 0.87, 0.80);
  leg->AddEntry(MassResolution1,"CMSSW720","f");
  //leg->AddEntry(MassResolution2,"CMSSW706","p");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.04;
  leg->SetTextSize(tsize2); 
  leg->Draw();  
  //======================================================================= 
  c111->Print("MC-Zprime5000-CMSSW720-MassResolution-fit.png","png");
  //c111->Print("PlotsDir/MC-DY400-MLP-Eta1-CMSSW532-500-E25-2000-EE.eps","eps");
  //=======================================================================
  //==================================================================  
  //==================================================================  
  //================================================================== 
  //==================================================================
  TCanvas *c1 = new TCanvas("c1","MassResolution",600,600);
  c1->SetFillColor(0);
  TH1F *MassResolution1 = new TH1F("MassResolution1","",100.0,-0.5,0.5);
  TH1F *MassResolution2 = new TH1F("MassResolution2","",100.0,-0.5,0.5);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx();
  gPad->SetGridy();
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  MassResolution1->Add(MassRecoGenDif);
  MassResolution1->SetLineStyle(0);
  MassResolution1->SetLineColor(4);
  MassResolution1->SetLineWidth(2);
  MassResolution1->SetMarkerColor(4); 
  //MassResolution->SetMarkerStyle(20);
  MassResolution1->SetMarkerStyle(2);
  //MassResolution->SetMarkerSize(0.70);
  MassResolution1->SetTitle("");
  MassResolution1->GetXaxis()->SetTitle("(dil. mass - gen dil. mass)/(gen dil. mass)");
  MassResolution1->GetYaxis()->SetTitle("A.U");
  MassResolution1->GetXaxis()->SetTitleOffset(1.3);
  MassResolution1->GetYaxis()->SetTitleOffset(2.0);
  MassResolution1->GetXaxis()->SetTitleSize(0.05);
  MassResolution1->DrawNormalized("",norm);
  //MassResolution->Sumw2();
  //MassResolution->GetXaxis()->SetRangeUser(0.95,1.05);
  //TF1* fn1 = new TF1("fn1","gaus",FitMin,FitMax);
  //fn1->SetLineColor(2);
  //MassResolution1->Fit("fn1","R");
  //MassResolution1->Draw("sames");
  //MassResolution->Draw("Dot,sames");
  //MassResolution->Sumw2();
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f10 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW706.root","READ");
  MassResolution2->Add(MassRecoGenDif);
  MassResolution2->SetLineStyle(0);
  MassResolution2->SetLineColor(2);
  MassResolution2->SetLineWidth(2);
  MassResolution2->SetMarkerColor(2); 
  MassResolution2->SetMarkerStyle(20);
  MassResolution2->SetMarkerSize(0.85);
  MassResolution2->DrawNormalized("sames",norm);
  //=================================================================
  /*
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
  */ 
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText2 = new TPaveText(0.20, 0.70, 0.50, 0.87, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("#sqrt{s} = 13 TeV, CMS, MC");
  TText *t2 = tText2->AddText("Z' #rightarrow #mu #mu [5 TeV]");
  TText *t3 = tText2->AddText("Phys14 (miniaod)");
  tText2->Draw(); 
  //==========================================================
  TLegend *leg = new TLegend(0.60, 0.60, 0.87, 0.80);
  leg->AddEntry(MassResolution1,"CMSSW720","f");
  leg->AddEntry(MassResolution2,"CMSSW706","p");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.04;
  leg->SetTextSize(tsize2); 
  leg->Draw();  
  //======================================================================= 
  c1->Print("MC-Zprime5000-CMSSW720-genMass-RecoMass-diff.png","png");
  //c1->Print("PlotsDir/MC-DY400-MLP-Eta1-CMSSW532-500-E25-2000-EE.eps","eps");
  //=======================================================================
  
  //==================================================================  
  //==================================================================  
  //================================================================== 
  //==================================================================
  TCanvas *c2 = new TCanvas("c2","RecoMass",600,600);
  c2->SetFillColor(0);
  TH1F *RecoMass1 = new TH1F("RecoMass1","",100.0,0.0,8000.0);
  TH1F *RecoMass2 = new TH1F("RecoMass2","",100.0,0.0,8000.0);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx();
  gPad->SetGridy();
  //gPad->SetLogy();
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  RecoMass1->Add(ZprimeRecomass);
  RecoMass1->SetLineStyle(0);
  RecoMass1->SetLineColor(1);
  RecoMass1->SetLineWidth(2);
  RecoMass1->SetMarkerColor(1); 
  RecoMass1->SetMarkerStyle(20);
  RecoMass1->SetMarkerSize(0.80);
  RecoMass1->SetTitle("");
  RecoMass1->GetXaxis()->SetTitle("M_{#mu#mu}");
  RecoMass1->GetYaxis()->SetTitle("A.U");
  RecoMass1->GetXaxis()->SetTitleOffset(1.3);
  RecoMass1->GetYaxis()->SetTitleOffset(2.0);
  RecoMass1->GetXaxis()->SetTitleSize(0.05);
  RecoMass1->GetXaxis()->SetLabelSize(0.03);
  RecoMass1->GetYaxis()->SetLabelSize(0.03);
  RecoMass1->DrawNormalized("",norm);
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f10 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW706.root","READ");
  RecoMass2->Add(ZprimeRecomass);
  RecoMass2->SetLineStyle(0);
  RecoMass2->SetLineColor(2);
  RecoMass2->SetLineWidth(2);
  RecoMass2->SetMarkerColor(2); 
  RecoMass2->SetMarkerStyle(20);
  RecoMass2->SetMarkerSize(0.85);
  RecoMass2->DrawNormalized("sames",norm);
  //=================================================================
  TPaveText* tText2 = new TPaveText(0.20, 0.70, 0.50, 0.87, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("#sqrt{s} = 13 TeV, CMS, MC");
  TText *t2 = tText2->AddText("Z' #rightarrow #mu #mu [5 TeV]");
  TText *t3 = tText2->AddText("Phys14 (miniaod)");
  tText2->Draw();
  //==========================================================
  TLegend *leg = new TLegend(0.60, 0.60, 0.87, 0.80);
  leg->AddEntry(RecoMass1,"CMSSW720","f");
  leg->AddEntry(RecoMass2,"CMSSW706","p");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.04;
  leg->SetTextSize(tsize2); 
  leg->Draw();  
  //======================================================================= 
  c2->Print("MC-Zprime5000-CMSSW720-RecoMuon-Mass.png","png");
  //c2->Print("PlotsDir/MC-DY400-MLP-Eta1-CMSSW532-500-E25-2000-EE.eps","eps");
  //=======================================================================

  //==================================================================  
  //==================================================================  
  //================================================================== 
  //==================================================================
  TCanvas *c3 = new TCanvas("c3","GenMass",600,600);
  c3->SetFillColor(0);
  TH1F *GenMass1 = new TH1F("GenMass1","",100.0,0.0,8000.0);
  TH1F *GenMass2 = new TH1F("GenMass2","",100.0,0.0,8000.0);
  TH1F *GenMass3 = new TH1F("GenMass3","",100.0,0.0,8000.0);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx();
  gPad->SetGridy();
  //gPad->SetLogy();
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  GenMass1->Add(ZprimeGenmass);
  GenMass1->SetLineStyle(0);
  GenMass1->SetLineColor(1);
  GenMass1->SetLineWidth(2);
  GenMass1->SetMarkerColor(1); 
  GenMass1->SetMarkerStyle(20);
  GenMass1->SetMarkerSize(0.80);
  GenMass1->SetTitle("");
  GenMass1->GetXaxis()->SetTitle("M_{#mu#mu}");
  GenMass1->GetYaxis()->SetTitle("A.U");
  GenMass1->GetXaxis()->SetTitleOffset(1.3);
  GenMass1->GetYaxis()->SetTitleOffset(2.0);
  GenMass1->GetXaxis()->SetTitleSize(0.05);
  GenMass1->GetXaxis()->SetLabelSize(0.03);
  GenMass1->GetYaxis()->SetLabelSize(0.03);
  GenMass1->DrawNormalized("",norm);
  //========================================================== 
  //                                                              
  //                   CMSSW539 M=3000                                    
  //==========================================================
  TFile *f11 = new TFile("ZprimetoMuMu-MC-Mass3000-CMSSW539.root","READ");
  GenMass3->Add(ZprimeGenmass);
  GenMass3->SetLineStyle(1);
  GenMass3->SetLineColor(4);
  GenMass3->SetLineWidth(2);
  GenMass3->SetMarkerColor(1);
  GenMass3->DrawNormalized("sames",norm);
  //========================================================== 
  //                                                              
  //                   CMSSW706 M=5000                                    
  //==========================================================
  TFile *f10 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW706.root","READ");
  GenMass2->Add(ZprimeGenmass);
  GenMass2->SetLineStyle(0);
  GenMass2->SetLineColor(2);
  GenMass2->SetLineWidth(2);
  GenMass2->SetMarkerColor(2); 
  GenMass2->SetMarkerStyle(20);
  GenMass2->SetMarkerSize(0.85);
  GenMass2->DrawNormalized("sames",norm);
  //=================================================================
  TPaveText* tText2 = new TPaveText(0.20, 0.70, 0.50, 0.87, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("#sqrt{s} = 13 TeV, CMS, MC");
  TText *t2 = tText2->AddText("Z' #rightarrow #mu #mu [5 TeV]");
  TText *t3 = tText2->AddText("Phys14 (miniaod)");
  tText2->Draw(); 
  //==========================================================
  TLegend *leg = new TLegend(0.60, 0.60, 0.87, 0.80);
  leg->AddEntry(GenMass1,"CMSSW720","f");
  leg->AddEntry(GenMass2,"CMSSW706","p");
  leg->AddEntry(GenMass3,"CMSSW539","p");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.04;
  leg->SetTextSize(tsize2); 
  leg->Draw();  
  //======================================================================= 
  c3->Print("MC-Zprime5000-CMSSW720-GenMuon-Mass.png","png");
  //c3->Print("PlotsDir/MC-DY400-MLP-Eta1-CMSSW532-500-E25-2000-EE.eps","eps");
  //=======================================================================



















  /*  
  //==================================================================  
  //==================================================================  
  //================================================================== 
  //==================================================================
  TCanvas *c3 = new TCanvas("c3","pt & En",600,600);
  c3->SetFillColor(0);
  TH1F *RecoPtMu1 = new TH1F("RecoPtMu1","",100,-8.0,8.0);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx();
  gPad->SetGridy();
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  RecoPtMu1->Add(LeadingRecoMuPt);
  RecoPtMu1->SetLineStyle(0);
  RecoPtMu1->SetLineColor(1);
  RecoPtMu1->SetLineWidth(2);
  RecoPtMu1->SetMarkerColor(1); 
  RecoPtMu1->SetMarkerStyle(20);
  RecoPtMu1->SetMarkerSize(0.80);
  RecoPtMu1->SetTitle("");
  RecoPtMu1->GetXaxis()->SetTitle("pt^{#mu1}");
  RecoPtMu1->GetYaxis()->SetTitle("Nb. events");
  RecoPtMu1->GetXaxis()->SetTitleOffset(1.3);
  RecoPtMu1->GetYaxis()->SetTitleOffset(2.0);
  RecoPtMu1->GetXaxis()->SetTitleSize(0.05);
  RecoPtMu1->GetXaxis()->SetLabelSize(0.03);
  RecoPtMu1->GetYaxis()->SetLabelSize(0.03);
  RecoPtMu1->Draw();
  //=================================================================
  gStyle->SetOptFit(11);
  //gStyle->SetOptStat(0);  
  TH1F *h=(TH1F*)gROOT->FindObject("RecoPtMu1");
  cout<<"histogramme "<<h<<endl;
  cout<<"Liste des fonctions histogramme "<<h->GetListOfFunctions()<<endl;
  cout<<"Objet histogramme "<<h->GetListOfFunctions()->FindObject("stats")<<endl;
  TPaveStats *st2 =(TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
  st2->SetX1NDC(0.60); //(0.16);
  st2->SetY1NDC(0.37); //(0.40);
  st2->SetX2NDC(0.88); //(0.60);
  st2->SetY2NDC(0.70); //(0.80);
  st2->SetLineColor(kBlack); st2->SetFillColor(0);
  st2->Draw(); 
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText2 = new TPaveText(0.40, 0.70, 0.70, 0.87, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("#sqrt{s} = 8 TeV, CMS, MC");
  TText *t2 = tText2->AddText("Z' #rightarrow #mu #mu [5 TeV]");
  TText *t3 = tText2->AddText("|#eta| < 2.40");
  tText2->Draw(); 
  //======================================================================= 
  c3->Print("MC-Zprime1000-CMSSW539-RecoMuon-pt1.png","png");
  //c2->Print("PlotsDir/MC-DY400-MLP-Eta1-CMSSW532-500-E25-2000-EE.eps","eps");
  //=======================================================================
  
  //==================================================================  
  //==================================================================  
  //================================================================== 
  //==================================================================
  TCanvas *c4 = new TCanvas("c4","pt & En",600,600);
  c4->SetFillColor(0);
  TH1F *RecoPtMu2 = new TH1F("RecoPtMu2","",100,-8.0,8.0);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx();
  gPad->SetGridy();
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  RecoPtMu2->Add(SecondLeadingMuPt);
  RecoPtMu2->SetLineStyle(0);
  RecoPtMu2->SetLineColor(1);
  RecoPtMu2->SetLineWidth(2);
  RecoPtMu2->SetMarkerColor(1); 
  RecoPtMu2->SetMarkerStyle(20);
  RecoPtMu2->SetMarkerSize(0.80);
  RecoPtMu2->SetTitle("");
  RecoPtMu2->GetXaxis()->SetTitle("pt^{#mu2}");
  RecoPtMu2->GetYaxis()->SetTitle("Nb. events");
  RecoPtMu2->GetXaxis()->SetTitleOffset(1.3);
  RecoPtMu2->GetYaxis()->SetTitleOffset(2.0);
  RecoPtMu2->GetXaxis()->SetTitleSize(0.05);
  RecoPtMu2->GetXaxis()->SetLabelSize(0.03);
  RecoPtMu2->GetYaxis()->SetLabelSize(0.03);
  RecoPtMu2->Draw();
  //=================================================================
  gStyle->SetOptFit(11);
  //gStyle->SetOptStat(0);  
  TH1F *h=(TH1F*)gROOT->FindObject("RecoPtMu2");
  cout<<"histogramme "<<h<<endl;
  cout<<"Liste des fonctions histogramme "<<h->GetListOfFunctions()<<endl;
  cout<<"Objet histogramme "<<h->GetListOfFunctions()->FindObject("stats")<<endl;
  TPaveStats *st2 =(TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
  st2->SetX1NDC(0.60); //(0.16);
  st2->SetY1NDC(0.37); //(0.40);
  st2->SetX2NDC(0.88); //(0.60);
  st2->SetY2NDC(0.70); //(0.80);
  st2->SetLineColor(kBlack); st2->SetFillColor(0);
  st2->Draw(); 
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText2 = new TPaveText(0.40, 0.70, 0.70, 0.87, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("#sqrt{s} = 8 TeV, CMS, MC");
  TText *t2 = tText2->AddText("Z' #rightarrow #mu #mu [5 TeV]");
  TText *t3 = tText2->AddText("|#eta| < 2.40");
  tText2->Draw(); 
  //======================================================================= 
  c4->Print("MC-Zprime1000-CMSSW539-RecoMuon-pt2.png","png");
  //c2->Print("PlotsDir/MC-DY400-MLP-Eta1-CMSSW532-500-E25-2000-EE.eps","eps");
  //=======================================================================
  
  //==================================================================  
  //==================================================================  
  //================================================================== 
  //==================================================================
  TCanvas *c5 = new TCanvas("c5","GenEta1",600,600);
  c5->SetFillColor(0);
  TH1F *GenEtaMu1 = new TH1F("GenEtaMu1","",100,-8.0,8.0);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx();
  gPad->SetGridy();
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  GenEtaMu1->Add(ZprimeGenEta1);
  GenEtaMu1->SetLineStyle(0);
  GenEtaMu1->SetLineColor(1);
  GenEtaMu1->SetLineWidth(2);
  GenEtaMu1->SetMarkerColor(1); 
  GenEtaMu1->SetMarkerStyle(20);
  GenEtaMu1->SetMarkerSize(0.80);
  GenEtaMu1->SetTitle("");
  GenEtaMu1->GetXaxis()->SetTitle("#eta^{#mu1}");
  GenEtaMu1->GetYaxis()->SetTitle("Nb. events");
  GenEtaMu1->GetXaxis()->SetTitleOffset(1.3);
  GenEtaMu1->GetYaxis()->SetTitleOffset(2.0);
  GenEtaMu1->GetXaxis()->SetTitleSize(0.05);
  GenEtaMu1->GetXaxis()->SetLabelSize(0.03);
  GenEtaMu1->GetYaxis()->SetLabelSize(0.03);
  GenEtaMu1->Draw();
  //=================================================================
  gStyle->SetOptFit(11);
  //gStyle->SetOptStat(0);  
  TH1F *h=(TH1F*)gROOT->FindObject("GenEtaMu1");
  cout<<"histogramme "<<h<<endl;
  cout<<"Liste des fonctions histogramme "<<h->GetListOfFunctions()<<endl;
  cout<<"Objet histogramme "<<h->GetListOfFunctions()->FindObject("stats")<<endl;
  TPaveStats *st2 =(TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
  st2->SetX1NDC(0.65); //(0.16);
  st2->SetY1NDC(0.37); //(0.40);
  st2->SetX2NDC(0.88); //(0.60);
  st2->SetY2NDC(0.70); //(0.80);
  st2->SetLineColor(kBlack); st2->SetFillColor(0);
  st2->Draw(); 
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText2 = new TPaveText(0.20, 0.70, 0.50, 0.87, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("#sqrt{s} = 8 TeV, CMS, MC");
  TText *t2 = tText2->AddText("Z' #rightarrow #mu #mu [5 TeV]");
  TText *t3 = tText2->AddText("Phys14 (miniaod)");
  //TText *t3 = tText2->AddText("|#eta| < 2.40");
  tText2->Draw(); 
  //======================================================================= 
  c5->Print("MC-Zprime5000-CMSSW720-GenMuon-eta1.png","png");
  //c5->Print("PlotsDir/MC-DY400-MLP-Eta1-CMSSW532-500-E25-2000-EE.eps","eps");
  //=======================================================================

  //==================================================================  
  //==================================================================  
  //================================================================== 
  //==================================================================
  TCanvas *c6 = new TCanvas("c6","GenEtaMu2",600,600);
  c6->SetFillColor(0);
  TH1F *GenEtaMu2 = new TH1F("GenEtaMu2","",100,-8.0,8.0);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx();
  gPad->SetGridy();
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  GenEtaMu2->Add(ZprimeGenEta2);
  GenEtaMu2->SetLineStyle(0);
  GenEtaMu2->SetLineColor(1);
  GenEtaMu2->SetLineWidth(2);
  GenEtaMu2->SetMarkerColor(1); 
  GenEtaMu2->SetMarkerStyle(20);
  GenEtaMu2->SetMarkerSize(0.80);
  GenEtaMu2->SetTitle("");
  GenEtaMu2->GetXaxis()->SetTitle("#eta^{#mu2}");
  GenEtaMu2->GetYaxis()->SetTitle("Nb. events");
  GenEtaMu2->GetXaxis()->SetTitleOffset(1.3);
  GenEtaMu2->GetYaxis()->SetTitleOffset(2.0);
  GenEtaMu2->GetXaxis()->SetTitleSize(0.05);
  GenEtaMu2->GetXaxis()->SetLabelSize(0.03);
  GenEtaMu2->GetYaxis()->SetLabelSize(0.03);
  GenEtaMu2->Draw();
  //=================================================================
  gStyle->SetOptFit(11);
  //gStyle->SetOptStat(0);  
  TH1F *h=(TH1F*)gROOT->FindObject("GenEtaMu2");
  cout<<"histogramme "<<h<<endl;
  cout<<"Liste des fonctions histogramme "<<h->GetListOfFunctions()<<endl;
  cout<<"Objet histogramme "<<h->GetListOfFunctions()->FindObject("stats")<<endl;
  TPaveStats *st2 =(TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
  st2->SetX1NDC(0.65); //(0.16);
  st2->SetY1NDC(0.37); //(0.40);
  st2->SetX2NDC(0.88); //(0.60);
  st2->SetY2NDC(0.70); //(0.80);
  st2->SetLineColor(kBlack); st2->SetFillColor(0);
  st2->Draw(); 
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText2 = new TPaveText(0.20, 0.70, 0.50, 0.87, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("#sqrt{s} = 8 TeV, CMS, MC");
  TText *t2 = tText2->AddText("Z' #rightarrow #mu #mu [5 TeV]");
  TText *t3 = tText2->AddText("Phys14 (miniaod)");
  //TText *t3 = tText2->AddText("|#eta| < 2.40");
  tText2->Draw(); 
  //======================================================================= 
  c6->Print("MC-Zprime5000-CMSSW720-GenMuon-eta2.png","png");
  //c6->Print("PlotsDir/MC-DY400-MLP-Eta1-CMSSW532-500-E25-2000-EE.eps","eps");
  //=======================================================================
  //==================================================================  
  //==================================================================  
  //================================================================== 
  //==================================================================
  
  TCanvas *c7 = new TCanvas("c7","pt & En",600,600);
  c7->SetFillColor(0);
  TH1F *RecoEnMu1 = new TH1F("RecoEnMu1","",100,0.0,2000.0);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx();
  gPad->SetGridy();
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  RecoEnMu1->Add(ZprimeRecoEn1);
  RecoEnMu1->SetLineStyle(0);
  RecoEnMu1->SetLineColor(1);
  RecoEnMu1->SetLineWidth(2);
  RecoEnMu1->SetMarkerColor(1); 
  RecoEnMu1->SetMarkerStyle(20);
  RecoEnMu1->SetMarkerSize(0.80);
  RecoEnMu1->SetTitle("");
  RecoEnMu1->GetXaxis()->SetTitle("En^{#mu1}");
  RecoEnMu1->GetYaxis()->SetTitle("Nb. events");
  RecoEnMu1->GetXaxis()->SetTitleOffset(1.3);
  RecoEnMu1->GetYaxis()->SetTitleOffset(2.0);
  RecoEnMu1->GetXaxis()->SetTitleSize(0.05);
  RecoEnMu1->GetXaxis()->SetLabelSize(0.03);
  RecoEnMu1->GetYaxis()->SetLabelSize(0.03);
  RecoEnMu1->Draw();
  //=================================================================
  gStyle->SetOptFit(11);
  //gStyle->SetOptStat(0);  
  TH1F *h=(TH1F*)gROOT->FindObject("RecoEnMu1");
  cout<<"histogramme "<<h<<endl;
  cout<<"Liste des fonctions histogramme "<<h->GetListOfFunctions()<<endl;
  cout<<"Objet histogramme "<<h->GetListOfFunctions()->FindObject("stats")<<endl;
  TPaveStats *st2 =(TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
  st2->SetX1NDC(0.60); //(0.16);
  st2->SetY1NDC(0.37); //(0.40);
  st2->SetX2NDC(0.88); //(0.60);
  st2->SetY2NDC(0.70); //(0.80);
  st2->SetLineColor(kBlack); st2->SetFillColor(0);
  st2->Draw(); 
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText2 = new TPaveText(0.50, 0.70, 0.80, 0.87, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("#sqrt{s} = 8 TeV, CMS, MC");
  TText *t2 = tText2->AddText("Z' #rightarrow #mu #mu [5 TeV]");
  //TText *t3 = tText2->AddText("|#eta| < 2.40");
  tText2->Draw(); 
  //======================================================================= 
  c7->Print("MC-Zprime1000-CMSSW539-RecoMuon-En1.png","png");
  //c7->Print("PlotsDir/MC-DY400-MLP-Eta1-CMSSW532-500-E25-2000-EE.eps","eps");
  //=======================================================================
  //==================================================================  
  //==================================================================  
  //================================================================== 
  //==================================================================
  TCanvas *c8 = new TCanvas("c8","pt & En",600,600);
  c8->SetFillColor(0);
  TH1F *RecoEnMu2 = new TH1F("RecoEnMu2","",100,0.0,2000.0);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx();
  gPad->SetGridy();
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  RecoEnMu2->Add(ZprimeRecoEn2);
  RecoEnMu2->SetLineStyle(0);
  RecoEnMu2->SetLineColor(1);
  RecoEnMu2->SetLineWidth(2);
  RecoEnMu2->SetMarkerColor(1); 
  RecoEnMu2->SetMarkerStyle(20);
  RecoEnMu2->SetMarkerSize(0.80);
  RecoEnMu2->SetTitle("");
  RecoEnMu2->GetXaxis()->SetTitle("En^{#mu2}");
  RecoEnMu2->GetYaxis()->SetTitle("Nb. events");
  RecoEnMu2->GetXaxis()->SetTitleOffset(1.3);
  RecoEnMu2->GetYaxis()->SetTitleOffset(2.0);
  RecoEnMu2->GetXaxis()->SetTitleSize(0.05);
  RecoEnMu2->GetXaxis()->SetLabelSize(0.03);
  RecoEnMu2->GetYaxis()->SetLabelSize(0.03);
  RecoEnMu2->Draw();
  //=================================================================
  gStyle->SetOptFit(11);
  //gStyle->SetOptStat(0);  
  TH1F *h=(TH1F*)gROOT->FindObject("RecoEnMu2");
  cout<<"histogramme "<<h<<endl;
  cout<<"Liste des fonctions histogramme "<<h->GetListOfFunctions()<<endl;
  cout<<"Objet histogramme "<<h->GetListOfFunctions()->FindObject("stats")<<endl;
  TPaveStats *st2 =(TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
  st2->SetX1NDC(0.60); //(0.16);
  st2->SetY1NDC(0.37); //(0.40);
  st2->SetX2NDC(0.88); //(0.60);
  st2->SetY2NDC(0.70); //(0.80);
  st2->SetLineColor(kBlack); st2->SetFillColor(0);
  st2->Draw(); 
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText2 = new TPaveText(0.50, 0.70, 0.80, 0.87, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("#sqrt{s} = 8 TeV, CMS, MC");
  TText *t2 = tText2->AddText("Z' #rightarrow #mu #mu [5 TeV]");
  //TText *t3 = tText2->AddText("|#eta| < 2.40");
  tText2->Draw(); 
  //======================================================================= 
  c8->Print("MC-Zprime1000-CMSSW539-RecoMuon-En2.png","png");
  //c8->Print("PlotsDir/MC-DY400-MLP-Eta1-CMSSW532-500-E25-2000-EE.eps","eps");
  //=======================================================================
  */



}
///////////////// Fit with 9 parameters function ////////////
Double_t FittingFunc5(Double_t *x,Double_t *par)
{
  double xx = x[0];
  Double_t function= par[0]*(par[1]+ par[2]*pow(cos(xx),2));
  
  //Double_t function= par[0]*(1+par[1]*x[0]+par[2]*pow(x[0],2)+par[3]*pow(x[0],3)+par[4]*pow(x[0],4))*(1+par[5]*x[1]+par[6]*pow(x[1],2)+par[7]*pow(x[1],3)+par[8]*pow(x[1],4));
  
  
  
  return function;
}


Double_t FittingFunc1(Double_t *x,Double_t *par)
{
  double xx = x[0];
  Double_t function= par[0]/pow((xx-par[1]),par[2]);
  return function;
}


Double_t FittingFunc2(Double_t *x,Double_t *par)
{
  double xx = x[0];
  Double_t function= par[0]/ pow((pow(xx,2) - par[1]*xx + par[2]),par[3]);
  return function;
}

Double_t FittingFunc4(Double_t *x,Double_t *par)
{
  double xx = x[0];
  Double_t function= exp( par[0] + par[1]*xx);
  
  //Double_t function= par[0]*(1+par[1]*x[0]+par[2]*pow(x[0],2)+par[3]*pow(x[0],3)+par[4]*pow(x[0],4))*(1+par[5]*x[1]+par[6]*pow(x[1],2)+par[7]*pow(x[1],3)+par[8]*pow(x[1],4));
  
  
  
  return function;
}





