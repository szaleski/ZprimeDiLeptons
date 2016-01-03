#include "TProfile.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TF2.h"
#include <string>
using std::string;
void Macro_Muon_ID_CMSSW539(){
  
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetStatFormat("5.3f");
  gStyle->SetFitFormat("5.3f");

  //gStyle->SetStatFormat("5.3f");
  //gStyle->SetFitFormat("5.3f");
  int BoxValue = 6611; //4680;  
  //int BoxValue = 11111111; //4680;  
  gStyle->SetOptFit(11);
  gStyle->SetOptDate(0);
  gStyle->SetOptTitle(0);
  //  gStyle->SetOptStat(BoxValue);
  gStyle->SetOptStat(kFALSE);
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
  float FitMin      = -0.075;
  float FitMax      = 0.08;
  float norm = 1.0;
  //==================================================================  
  //==================================================================  
  //================================================================== 
  //================================================================== 
  //================================================================== 
  //==================================================================
  TCanvas *c5 = new TCanvas("c5","[1]",600,600);
  c5->SetFillColor(0);
  TH1F *RecoEtaMu1 = new TH1F("RecoEtaMu1","",100,0.0,0.5);
  TH1F *RecoEtaMu2 = new TH1F("RecoEtaMu2","",100,0.0,0.5);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  RecoEtaMu1->Add(dPToverPT);
  RecoEtaMu1->SetLineStyle(0);
  RecoEtaMu1->SetLineColor(1);
  RecoEtaMu1->SetLineWidth(2);
  RecoEtaMu1->SetMarkerColor(1); 
  RecoEtaMu1->SetMarkerStyle(20);
  RecoEtaMu1->SetMarkerSize(0.80);
  RecoEtaMu1->SetTitle("");
  RecoEtaMu1->GetXaxis()->SetTitle("#delta pt^{BestTrack}/pt^{BestTrack}");
  RecoEtaMu1->GetYaxis()->SetTitle("A.U.");
  RecoEtaMu1->GetXaxis()->SetTitleOffset(1.3);
  RecoEtaMu1->GetYaxis()->SetTitleOffset(2.0);
  RecoEtaMu1->GetXaxis()->SetTitleSize(0.05);
  RecoEtaMu1->GetXaxis()->SetLabelSize(0.03);
  RecoEtaMu1->GetYaxis()->SetLabelSize(0.03);
  RecoEtaMu1->DrawNormalized("",norm);
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f10 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW706.root","READ");
  RecoEtaMu2->Add(dPToverPT);
  RecoEtaMu2->SetLineStyle(0);
  RecoEtaMu2->SetLineColor(2);
  RecoEtaMu2->SetLineWidth(2);
  RecoEtaMu2->SetMarkerColor(2); 
  RecoEtaMu2->SetMarkerStyle(20);
  RecoEtaMu2->SetMarkerSize(0.82);
  RecoEtaMu2->DrawNormalized("sames",norm);
  //=================================================================
  double LikelihoodCut    = 0.3 ;
  TLine *line = new TLine( LikelihoodCut , 0.0,  LikelihoodCut , 15500.0);
  line->SetLineColor(2);
  line->SetLineWidth(3);
  line->SetLineStyle(2);
  line->Draw();
  //=================================================================
  //=================================================================
  gStyle->SetOptFit(11);
  //gStyle->SetOptStat(0);  
  TH1F *h=(TH1F*)gROOT->FindObject("RecoEtaMu1");
  cout<<"histogramme "<<h<<endl;
  cout<<"Liste des fonctions histogramme "<<h->GetListOfFunctions()<<endl;
  cout<<"Objet histogramme "<<h->GetListOfFunctions()->FindObject("stats")<<endl;
  TPaveStats *st2 =(TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
  //st2->SetX1NDC(0.65); //(0.16);
  //st2->SetY1NDC(0.37); //(0.40);
  //st2->SetX2NDC(0.88); //(0.60);
  //st2->SetY2NDC(0.70); //(0.80);
  //st2->SetLineColor(kBlack); st2->SetFillColor(0);
  //st2->Draw(); 
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText2 = new TPaveText(0.25, 0.70, 0.55, 0.87, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("#sqrt{s} = 13 TeV, CMS, MC");
  TText *t2 = tText2->AddText("Z' #rightarrow #mu #mu [5 TeV]");
  //TText *t3 = tText2->AddText("|#eta| < 2.40");
  tText2->Draw(); 
  //==========================================================
  TLegend *leg = new TLegend(0.60, 0.60, 0.87, 0.80);
  leg->AddEntry(RecoEtaMu1,"CMSSW720","f");
  leg->AddEntry(RecoEtaMu2,"CMSSW706","p");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.04;
  leg->SetTextSize(tsize2); 
  leg->Draw();  
  //======================================================================= 
  c5->Print("MC-Zprime5000-CMSSW720-RecoMuon-errorPt.png","png");
  //c5->Print("PlotsDir/MC-DY400-MLP-Eta1-CMSSW532-500-E25-2000-EE.eps","eps");
  //=======================================================================

  //==================================================================  
  //==================================================================  
  //================================================================== 
  //==================================================================
  TCanvas *c6 = new TCanvas("c6","[2]",600,600);
  c6->SetFillColor(0);
  TH1F *RecoEtaMu22 = new TH1F("RecoEtaMu22","",60,0.0,60.0);
  TH1F *RecoEtaMu23 = new TH1F("RecoEtaMu23","",60,0.0,60.0);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f11 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  RecoEtaMu22->Add(numberOfValidMuonHits);
  RecoEtaMu22->SetLineStyle(0);
  RecoEtaMu22->SetLineColor(1);
  RecoEtaMu22->SetLineWidth(2);
  RecoEtaMu22->SetMarkerColor(1); 
  RecoEtaMu22->SetMarkerStyle(20);
  RecoEtaMu22->SetMarkerSize(0.80);
  RecoEtaMu22->SetTitle("");
  RecoEtaMu22->GetXaxis()->SetTitle("N(Mu hits)");
  RecoEtaMu22->GetYaxis()->SetTitle("A.U.");
  RecoEtaMu22->GetXaxis()->SetTitleOffset(1.3);
  RecoEtaMu22->GetYaxis()->SetTitleOffset(2.0);
  RecoEtaMu22->GetXaxis()->SetTitleSize(0.05);
  RecoEtaMu22->GetXaxis()->SetLabelSize(0.03);
  RecoEtaMu22->GetYaxis()->SetLabelSize(0.03);
  RecoEtaMu22->DrawNormalized("",norm);
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f112 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW706.root","READ");
  RecoEtaMu23->Add(numberOfValidMuonHits);
  RecoEtaMu23->SetLineStyle(0);
  RecoEtaMu23->SetLineColor(2);
  RecoEtaMu23->SetLineWidth(2);
  RecoEtaMu23->SetMarkerColor(2); 
  RecoEtaMu23->SetMarkerStyle(20);
  RecoEtaMu23->SetMarkerSize(0.82);
  RecoEtaMu23->DrawNormalized("sames",norm);
  //=================================================================
  double LikelihoodCut    = 1.0 ;
  TLine *line = new TLine( LikelihoodCut , 0.0,  LikelihoodCut , 3797.0);
  line->SetLineColor(2);
  line->SetLineWidth(3);
  line->SetLineStyle(2);
  line->Draw();
  //=================================================================
  TPaveText* tText2 = new TPaveText(0.30, 0.40, 0.60, 0.57, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("#sqrt{s} = 13 TeV, CMS, MC");
  TText *t2 = tText2->AddText("Z' #rightarrow #mu #mu [5 TeV]");
  //TText *t3 = tText2->AddText("|#eta| < 2.40");
  tText2->Draw(); 
  //==========================================================
  TLegend *leg = new TLegend(0.60, 0.60, 0.87, 0.80);
  leg->AddEntry(RecoEtaMu22,"CMSSW720","f");
  leg->AddEntry(RecoEtaMu23,"CMSSW706","p");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.04;
  leg->SetTextSize(tsize2); 
  leg->Draw();
  //======================================================================= 
  c6->Print("MC-Zprime5000-CMSSW720-RecoMuon-numberOfValidMuonHits.png","png");
  //c6->Print("PlotsDir/MC-DY400-MLP-Eta1-CMSSW532-500-E25-2000-EE.eps","eps");
  //=======================================================================
  //==================================================================  
  //==================================================================  
  //================================================================== 
  //==================================================================
  TCanvas *c7 = new TCanvas("c7","[3]",600,600);
  c7->SetFillColor(0);
  TH1F *RecoEnMu31 = new TH1F("RecoEnMu31","",10,0.0,10.0);
  TH1F *RecoEnMu32 = new TH1F("RecoEnMu32","",10,0.0,10.0);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  RecoEnMu31->Add(numberOfValidPixelHits);
  RecoEnMu31->SetLineStyle(0);
  RecoEnMu31->SetLineColor(1);
  RecoEnMu31->SetLineWidth(2);
  RecoEnMu31->SetMarkerColor(1); 
  RecoEnMu31->SetMarkerStyle(20);
  RecoEnMu31->SetMarkerSize(0.80);
  RecoEnMu31->SetTitle("");
  RecoEnMu31->GetXaxis()->SetTitle("N(Pixel hits)");
  RecoEnMu31->GetYaxis()->SetTitle("A.U.");
  RecoEnMu31->GetXaxis()->SetTitleOffset(1.3);
  RecoEnMu31->GetYaxis()->SetTitleOffset(2.0);
  RecoEnMu31->GetXaxis()->SetTitleSize(0.05);
  RecoEnMu31->GetXaxis()->SetLabelSize(0.03);
  RecoEnMu31->GetYaxis()->SetLabelSize(0.03);
  RecoEnMu31->DrawNormalized("",norm);
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f112 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW706.root","READ");
  RecoEnMu32->Add(numberOfValidPixelHits);
  RecoEnMu32->SetLineStyle(0);
  RecoEnMu32->SetLineColor(2);
  RecoEnMu32->SetLineWidth(2);
  RecoEnMu32->SetMarkerColor(2); 
  RecoEnMu32->SetMarkerStyle(20);
  RecoEnMu32->SetMarkerSize(0.82);
  RecoEnMu32->DrawNormalized("sames",norm);
  //=================================================================
  double LikelihoodCut    = 1.0 ;
  TLine *line = new TLine( LikelihoodCut , 0.0,  LikelihoodCut , 56736.0);
  line->SetLineColor(2);
  line->SetLineWidth(3);
  line->SetLineStyle(2);
  line->Draw();
  //=================================================================
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  TPaveText* tText2 = new TPaveText(0.50, 0.70, 0.80, 0.87, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("#sqrt{s} = 13 TeV, CMS, MC");
  TText *t2 = tText2->AddText("Z' #rightarrow #mu #mu [5 TeV]");
  //TText *t3 = tText2->AddText("|#eta| < 2.40");
  tText2->Draw(); 
  //==========================================================
  TLegend *leg = new TLegend(0.60, 0.60, 0.87, 0.80);
  leg->AddEntry(RecoEnMu31,"CMSSW720","f");
  leg->AddEntry(RecoEnMu32,"CMSSW706","p");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.04;
  leg->SetTextSize(tsize2); 
  leg->Draw();
  //======================================================================= 
  c7->Print("MC-Zprime5000-CMSSW720-RecoMuon-numberOfValidPixelHits.png","png");
  //c7->Print("PlotsDir/MC-DY400-MLP-Eta1-CMSSW532-500-E25-2000-EE.eps","eps");
  //=======================================================================
  //==================================================================  
  //==================================================================  
  //================================================================== 
  //==================================================================
  TCanvas *c8 = new TCanvas("c8","[4]",600,600);
  c8->SetFillColor(0);
  TH1F *RecoEnMu41 = new TH1F("RecoEnMu41","",10,0.0,10.0);
  TH1F *RecoEnMu42 = new TH1F("RecoEnMu42","",10,0.0,10.0);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  RecoEnMu41->Add(numberOfMatchedStations);
  RecoEnMu41->SetLineStyle(0);
  RecoEnMu41->SetLineColor(1);
  RecoEnMu41->SetLineWidth(2);
  RecoEnMu41->SetMarkerColor(1); 
  RecoEnMu41->SetMarkerStyle(20);
  RecoEnMu41->SetMarkerSize(0.80);
  RecoEnMu41->SetTitle("");
  RecoEnMu41->GetXaxis()->SetTitle("Number of Matches");
  RecoEnMu41->GetYaxis()->SetTitle("A.U.");
  RecoEnMu41->GetXaxis()->SetTitleOffset(1.3);
  RecoEnMu41->GetYaxis()->SetTitleOffset(2.0);
  RecoEnMu41->GetXaxis()->SetTitleSize(0.05);
  RecoEnMu41->GetXaxis()->SetLabelSize(0.03);
  RecoEnMu41->GetYaxis()->SetLabelSize(0.03);
  RecoEnMu41->DrawNormalized("",norm);
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f112 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW706.root","READ");
  RecoEnMu42->Add(numberOfMatchedStations);
  RecoEnMu42->SetLineStyle(0);
  RecoEnMu42->SetLineColor(2);
  RecoEnMu42->SetLineWidth(2);
  RecoEnMu42->SetMarkerColor(2); 
  RecoEnMu42->SetMarkerStyle(20);
  RecoEnMu42->SetMarkerSize(0.82);
  RecoEnMu42->DrawNormalized("sames",norm);
  //=================================================================
  double LikelihoodCut    = 2.0 ;
  TLine *line = new TLine( LikelihoodCut , 0.0,  LikelihoodCut , 39497.0);
  line->SetLineColor(2);
  line->SetLineWidth(3);
  line->SetLineStyle(2);
  line->Draw();
  //=================================================================
  //=================================================================
  TPaveText* tText2 = new TPaveText(0.50, 0.70, 0.80, 0.87, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("#sqrt{s} = 13 TeV, CMS, MC");
  TText *t2 = tText2->AddText("Z' #rightarrow #mu #mu [5 TeV]");
  //TText *t3 = tText2->AddText("|#eta| < 2.40");
  tText2->Draw(); 
  //==========================================================
  TLegend *leg = new TLegend(0.60, 0.60, 0.87, 0.80);
  leg->AddEntry(RecoEnMu41,"CMSSW720","f");
  leg->AddEntry(RecoEnMu42,"CMSSW706","p");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.04;
  leg->SetTextSize(tsize2); 
  leg->Draw();
  //======================================================================= 
  c8->Print("MC-Zprime5000-CMSSW720-RecoMuon-numberOfMatchedStations.png","png");
  //c8->Print("PlotsDir/MC-DY400-MLP-Eta1-CMSSW532-500-E25-2000-EE.eps","eps");
  //=======================================================================
  //==================================================================  
  //==================================================================  
  //================================================================== 
  //==================================================================
  TCanvas *c9 = new TCanvas("c9","[5]",600,600);
  c9->SetFillColor(0);
  TH1F *RecoPhiMu1 = new TH1F("RecoPhiMu1","",20,0.0,20.0);
  TH1F *RecoPhiMu2 = new TH1F("RecoPhiMu2","",20,0.0,20.0);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  RecoPhiMu1->Add(numberOftrackerLayersWithMeasurement);
  RecoPhiMu1->SetLineStyle(0);
  RecoPhiMu1->SetLineColor(1);
  RecoPhiMu1->SetLineWidth(2);
  RecoPhiMu1->SetMarkerColor(1); 
  RecoPhiMu1->SetMarkerStyle(20);
  RecoPhiMu1->SetMarkerSize(0.80);
  RecoPhiMu1->SetTitle("");
  RecoPhiMu1->GetXaxis()->SetTitle("N_{Layers}");
  RecoPhiMu1->GetYaxis()->SetTitle("A.U.");
  RecoPhiMu1->GetXaxis()->SetTitleOffset(1.3);
  RecoPhiMu1->GetYaxis()->SetTitleOffset(2.0);
  RecoPhiMu1->GetXaxis()->SetTitleSize(0.05);
  RecoPhiMu1->GetXaxis()->SetLabelSize(0.03);
  RecoPhiMu1->GetYaxis()->SetLabelSize(0.03);
  RecoPhiMu1->DrawNormalized("",norm);
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f112 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW706.root","READ");
  RecoPhiMu2->Add(numberOftrackerLayersWithMeasurement);
  RecoPhiMu2->SetLineStyle(0);
  RecoPhiMu2->SetLineColor(2);
  RecoPhiMu2->SetLineWidth(2);
  RecoPhiMu2->SetMarkerColor(2); 
  RecoPhiMu2->SetMarkerStyle(20);
  RecoPhiMu2->SetMarkerSize(0.82);
  RecoPhiMu2->DrawNormalized("sames",norm);
  //=================================================================
  double LikelihoodCut    = 6.0 ;
  TLine *line = new TLine( LikelihoodCut , 0.0,  LikelihoodCut , 39497.0);
  line->SetLineColor(2);
  line->SetLineWidth(3);
  line->SetLineStyle(2);
  line->Draw();
  //=================================================================
  //=================================================================
  TPaveText* tText2 = new TPaveText(0.20, 0.70, 0.50, 0.87, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("#sqrt{s} = 13 TeV, CMS, MC");
  TText *t2 = tText2->AddText("Z' #rightarrow #mu #mu [5 TeV]");
  //TText *t3 = tText2->AddText("|#eta| < 2.40");
  tText2->Draw(); 
  //==========================================================
  TLegend *leg = new TLegend(0.60, 0.60, 0.87, 0.80);
  leg->AddEntry(RecoPhiMu1,"CMSSW720","f");
  leg->AddEntry(RecoPhiMu2,"CMSSW706","p");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.04;
  leg->SetTextSize(tsize2); 
  leg->Draw();
  //======================================================================= 
  c9->Print("MC-Zprime5000-CMSSW720-RecoMuon-numberOftrackerLayers.png","png");
  //c9->Print("PlotsDir/MC-DY400-MLP-Eta1-CMSSW532-500-E25-2000-EE.eps","eps");
  //=======================================================================
  //==================================================================  
  //==================================================================  
  //================================================================== 
  //==================================================================
  TCanvas *c10 = new TCanvas("c10","[6]",600,600);
  c10->SetFillColor(0);
  TH1F *RecoPhiMu12 = new TH1F("RecoPhiMu12","",50,0.0,0.3);
  TH1F *RecoPhiMu13 = new TH1F("RecoPhiMu13","",50,0.0,0.3);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  RecoPhiMu12->Add(trackiso);
  RecoPhiMu12->SetLineStyle(0);
  RecoPhiMu12->SetLineColor(1);
  RecoPhiMu12->SetLineWidth(2);
  RecoPhiMu12->SetMarkerColor(1); 
  RecoPhiMu12->SetMarkerStyle(20);
  RecoPhiMu12->SetMarkerSize(0.80);
  RecoPhiMu12->SetTitle("");
  RecoPhiMu12->GetXaxis()->SetTitle("Tracker relative isolation");
  RecoPhiMu12->GetYaxis()->SetTitle("A.U.");
  RecoPhiMu12->GetXaxis()->SetTitleOffset(1.3);
  RecoPhiMu12->GetYaxis()->SetTitleOffset(2.0);
  RecoPhiMu12->GetXaxis()->SetTitleSize(0.05);
  RecoPhiMu12->GetXaxis()->SetLabelSize(0.03);
  RecoPhiMu12->GetYaxis()->SetLabelSize(0.03);
  RecoPhiMu12->DrawNormalized("",norm);
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f112 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW706.root","READ");
  RecoPhiMu13->Add(trackiso);
  RecoPhiMu13->SetLineStyle(0);
  RecoPhiMu13->SetLineColor(2);
  RecoPhiMu13->SetLineWidth(2);
  RecoPhiMu13->SetMarkerColor(2); 
  RecoPhiMu13->SetMarkerStyle(20);
  RecoPhiMu13->SetMarkerSize(0.82);
  RecoPhiMu13->DrawNormalized("sames",norm);
  //=================================================================
  double LikelihoodCut    = 0.1 ;
  TLine *line = new TLine( LikelihoodCut , 0.0,  LikelihoodCut , 77209.0);
  line->SetLineColor(2);
  line->SetLineWidth(3);
  line->SetLineStyle(2);
  line->Draw();
  //=================================================================
  //=================================================================
  TPaveText* tText2 = new TPaveText(0.20, 0.70, 0.50, 0.87, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("#sqrt{s} = 13 TeV, CMS, MC");
  TText *t2 = tText2->AddText("Z' #rightarrow #mu #mu [5 TeV]");
  //TText *t3 = tText2->AddText("|#eta| < 2.40");
  tText2->Draw();
  //==========================================================
  TLegend *leg = new TLegend(0.60, 0.60, 0.87, 0.80);
  leg->AddEntry(RecoPhiMu12,"CMSSW720","f");
  leg->AddEntry(RecoPhiMu13,"CMSSW706","p");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.04;
  leg->SetTextSize(tsize2); 
  leg->Draw();
  //======================================================================= 
  c10->Print("MC-Zprime5000-CMSSW720-RecoMuon-trackiso.png","png");
  //c10->Print("PlotsDir/MC-DY400-MLP-Eta1-CMSSW532-500-E25-2000-EE.eps","eps");
  //=======================================================================
  //==================================================================  
  //==================================================================  
  //================================================================== 
  //==================================================================
  TCanvas *c11 = new TCanvas("c11","[7]",600,600);
  c11->SetFillColor(0);
  TH1F *RecoPhiDif1 = new TH1F("RecoPhiDif1","",100,0.0,0.3);
  TH1F *RecoPhiDif2 = new TH1F("RecoPhiDif2","",100,0.0,0.3);
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f1 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW720.root","READ");
  RecoPhiDif1->Add(absdxy);
  RecoPhiDif1->SetLineStyle(0);
  RecoPhiDif1->SetLineColor(1);
  RecoPhiDif1->SetLineWidth(2);
  RecoPhiDif1->SetMarkerColor(1); 
  RecoPhiDif1->SetMarkerStyle(20);
  RecoPhiDif1->SetMarkerSize(0.80);
  RecoPhiDif1->SetTitle("");
  RecoPhiDif1->GetXaxis()->SetTitle("|dxy|");
  RecoPhiDif1->GetYaxis()->SetTitle("A.U.");
  RecoPhiDif1->GetXaxis()->SetTitleOffset(1.3);
  RecoPhiDif1->GetYaxis()->SetTitleOffset(2.0);
  RecoPhiDif1->GetXaxis()->SetTitleSize(0.05);
  RecoPhiDif1->GetXaxis()->SetLabelSize(0.03);
  RecoPhiDif1->GetYaxis()->SetLabelSize(0.03);
  RecoPhiDif1->DrawNormalized("",norm);
  //========================================================== 
  //                                                              
  //               MC Z -> ee                                    
  //==========================================================
  TFile *f112 = new TFile("ZprimetoMuMu-MC-Mass5000-CMSSW706.root","READ");
  RecoPhiDif2->Add(absdxy);
  RecoPhiDif2->SetLineStyle(0);
  RecoPhiDif2->SetLineColor(2);
  RecoPhiDif2->SetLineWidth(2);
  RecoPhiDif2->SetMarkerColor(2); 
  RecoPhiDif2->SetMarkerStyle(20);
  RecoPhiDif2->SetMarkerSize(0.82);
  RecoPhiDif2->DrawNormalized("sames",norm);
  //=================================================================
  double LikelihoodCut    = 0.2 ;
  TLine *line = new TLine( LikelihoodCut , 0.0,  LikelihoodCut , 1104.0);
  line->SetLineColor(2);
  line->SetLineWidth(3);
  line->SetLineStyle(2);
  line->Draw();
  //=================================================================
  //=================================================================
  TPaveText* tText2 = new TPaveText(0.30, 0.70, 0.60, 0.87, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("#sqrt{s} = 13 TeV, CMS, MC");
  TText *t2 = tText2->AddText("Z' #rightarrow #mu #mu [5 TeV]");
  //TText *t3 = tText2->AddText("|#eta| < 2.40");
  tText2->Draw(); 
  //==========================================================
  TLegend *leg = new TLegend(0.60, 0.60, 0.87, 0.80);
  leg->AddEntry(RecoPhiDif1,"CMSSW720","f");
  leg->AddEntry(RecoPhiDif2,"CMSSW706","p");
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.04;
  leg->SetTextSize(tsize2); 
  leg->Draw();
  //======================================================================= 
  c11->Print("MC-Zprime5000-CMSSW720-RecoMuon-dxy.png","png");
  //c11->Print("PlotsDir/MC-DY400-MLP-Eta1-CMSSW532-500-E25-2000-EE.eps","eps");
  //=======================================================================
  



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





