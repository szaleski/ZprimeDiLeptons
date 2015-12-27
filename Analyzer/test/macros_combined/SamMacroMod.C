#include <string>
#include <vector>
#include <sstream>
#include "TFile.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TPaveLabel.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TColor.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "Math/QuantFuncMathCore.h"
#include "TF1.h"
using std::string;

using std::string;
TH1* makeUnifiedPlot(const std::string& filename,bool isMuon,int histType,int regionCode);
//helper functions
TGraphAsymmErrors *makeDataGraph(TH1* dataHist,float normToBinWidth,bool xErrBar=true);
void normHistToBinWidth(TH1* hist,float binNormNr);
TH1* makeIntHist(const TH1* hist,bool intIsGreatThan=true);
//TGraphAsymmErrors *makeDataGraph  = new TGraphAsymmErrors;
//makeDataGraph(TH1* dataHist,float normToBinWidth,bool xErrBar=true);
//float getBkgErr(float mass,bool isMuon,int regionCode,TH1* zeeHistEBEB,TH1* zeeHistEBEE);
float getBkgErr(float mass);
void PaintOverflow(TH1 *h);
void SamMacroMod(){
  TCanvas *c1 = new TCanvas("c1", "c1",1000,700);
  TPad* spectrumPad=0;
  TPad* ratioPad=0;
  //c1->Divide(1,2);
  //gStyle->SetOptStat(111111);
  gStyle->SetOptFit(kFALSE);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetStatX(0.457589);
  gStyle->SetStatY(0.312937);
  gStyle->SetStatW(0.29241/2+0.0185);
  gStyle->SetStatH(0.169580+0.05);
  gStyle->SetStatFontSize(0.0402098);
  gStyle->SetStatFont(0.02);
  gStyle->SetFitFormat("5.2g");
  gStyle->SetStatBorderSize(0);
  gStyle->SetStatFontSize(0.040209);
  gStyle->SetStatFontSize(0.035209);
  c1->Range(1.592761,-5.173913,3.533814,6.006211);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetLogx(1);
  c1->SetLogy(1);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLeftMargin(0.13);
  c1->SetRightMargin(0.07);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);
  c1->SetTopMargin(0.085);
  c1->SetBottomMargin(0.11);
  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.002); // Upper and lower plot are joined
  pad1->SetTickx(1);
  pad1->SetTicky(1);
  pad1->SetLogx(1);
  pad1->SetLogy(1);
  pad1->SetGridx(0);
  pad1->SetGridy(0);
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  //========================================================== 
  //                                                              
  //               Get the histograms                                  
  //==========================================================
  TFile *file1 = new TFile("plots/h_ZprimeRecomass_DY_OF.root","READ");
  TH1* zeeHist = (TH1*) file1->Get("ZprimeRecomassBinWidthwtOverFlow");

  TFile *file2 = new TFile("plots/h_ZprimeRecomass_Tlike_OF.root","READ");
  TH1* ttbarHist = (TH1*) file2->Get("ZprimeRecomassBinWidthwtOverFlow");

  TFile *file3 = new TFile("plots/h_ZprimeRecomass_DiBoson_OF.root","READ");
  TH1* dibosonsBkgHist = (TH1*) file3->Get("ZprimeRecomassBinWidthwtOverFlow");

  TFile *file4 = new TFile("plots/h_ZprimeRecomass_data_OF.root","READ");
  TH1* dataHistTempbar = (TH1*) file4->Get("ZprimeRecomassBinWidthwtOverFlow");

  //TFile *file3 = new TFile("Wjets-25nsMC-OS-allbins-2673pb.root","READ");
  TFile *file5 = new TFile("plots/FR-Wjets-25nsMC-OS-BinWidth-2673pb_OF.root","READ");
  TH1* WjetsBkgHist = (TH1*) file5->Get("WjetsHistowtOverFlow");
  std::cout<<"nb.Wjets(MC)  = "<<WjetsBkgHist->Integral()<<endl;

  //TFile *file4 = new TFile("DiJets-Data-OS-2673pb-FR.root","READ");
  TFile *file6 = new TFile("plots/FR-DiJets-Data-OS-BinWidth-2673pb_OF.root","READ");
  TH1* jetBkgHist   = (TH1*) file6->Get("DataSubwtOverFlow");
  std::cout<<"nbQCD(dijets,Data)  = "<<jetBkgHist->Integral()<<endl;

  /*
  zeeHist->Rebin(20);
  ttbarHist->Rebin(20);
  dibosonsBkgHist->Rebin(20);
  dataHistTempbar->Rebin(20);
  WjetsBkgHist->Rebin(20);
  jetBkgHist->Rebin(20);
  */

  float binWidthNorm=-1;
  int zeeColour    =  TColor::GetColor("#99ccff");
  int jetBkgColour = TColor::GetColor("#ffff66"); 
  int ttbarColour  = TColor::GetColor("#ff6666");
  int bosonColour  = TColor::GetColorDark(3);
  int WjetsColour  = TColor::GetColorDark(5);
  int font = 42;
  //float xAxisMin = 60; //72
  //float xAxisMax = 1000.0;
  //float yAxisMin = 1e-4;
  //float yAxisMax = 1e3; 
  
  TGraphAsymmErrors* dataHist = makeDataGraph(dataHistTempbar,binWidthNorm,0);
  normHistToBinWidth(zeeHist,binWidthNorm);
  normHistToBinWidth(ttbarHist,binWidthNorm);
  normHistToBinWidth(dibosonsBkgHist,binWidthNorm);
  normHistToBinWidth(jetBkgHist,binWidthNorm);
  normHistToBinWidth(WjetsBkgHist,binWidthNorm);
  
  

  //gStyle->SetOptStat(111111);
  //PaintOverflow(zeeHist);
  //PaintOverflow(ttbarHist);
  //PaintOverflow(dibosonsBkgHist);
  //PaintOverflow(dataHistTempbar);
  //PaintOverflow(jetBkgHist);
  //PaintOverflow(WjetsBkgHist);
  

  THStack *axisHist = new THStack("axisHist","");

  zeeHist->SetFillColor(zeeColour);
  zeeHist->SetLineWidth(2);
  zeeHist->SetLineColor(1);
  zeeHist->SetTitle("");
  
  
  ttbarHist->SetFillColor(ttbarColour);
  ttbarHist->SetLineWidth(2); 
  ttbarHist->SetLineColor(1);
  
  dibosonsBkgHist->SetFillColor(bosonColour);
  dibosonsBkgHist->SetLineWidth(2); 
  dibosonsBkgHist->SetLineColor(1);
  
  WjetsBkgHist->SetFillColor(WjetsColour);
  //WjetsBkgHist->SetFillColor(jetBkgColour);
  WjetsBkgHist->SetLineWidth(2);
  WjetsBkgHist->SetLineColor(1);

  jetBkgHist->SetFillColor(jetBkgColour);
  jetBkgHist->SetLineWidth(2);
  jetBkgHist->SetLineColor(1);

  
  axisHist->Add(jetBkgHist,"histo");
  axisHist->Add(WjetsBkgHist,"histo");
  axisHist->Add(ttbarHist);
  axisHist->Add(dibosonsBkgHist);
  axisHist->Add(zeeHist);
  axisHist->Draw("histo");
  

  //axisHist->GetXaxis()->SetTitle("M(#mu^{+}#mu^{-}) [GeV]");
  //axisHist->GetYaxis()->SetTitle("Events / 20 GeV");
  axisHist->GetYaxis()->SetTitle("Events/GeV");
  //axisHist->GetXaxis()->SetTitleOffset(1.1);
  //axisHist->GetYaxis()->SetTitleOffset(1.1);
  axisHist->GetXaxis()->SetTitleSize(0.047);
  axisHist->GetYaxis()->SetTitleSize(0.047);
  //axisHist->GetXaxis()->SetLabelSize(0.040);
  axisHist->GetYaxis()->SetLabelSize(0.040);
  //axisHist->GetXaxis()->SetMoreLogLabels();
  //axisHist->GetXaxis()->SetNoExponent();
  axisHist->GetXaxis()->SetRangeUser(60.0,3000.0);
  axisHist->SetMinimum(0.01);
  axisHist->SetMaximum(20000.0);
  //axisHist->SetMaximum(2000.0);
  dataHist->SetMarkerStyle(20);
  dataHist->SetMarkerSize(0.9);
  dataHist->GetXaxis()->SetRange(5,83);
  dataHist->GetXaxis()->SetTitleSize(0.047);
  dataHist->GetXaxis()->SetTitleOffset(0.9);
  dataHist->GetYaxis()->SetTitleSize(0.047);
  dataHist->GetYaxis()->SetTitleOffset(1.2);
  dataHist->Draw("PZ");
  

  //==========================================================
  TLegend *leg = new TLegend(0.56741,0.48671,0.820536,0.83664,NULL,"brNDC"); //for lumi in plot
  leg->AddEntry(dataHist,"Data","PE");
  leg->AddEntry(zeeHist,"#gamma^{*}/Z#rightarrow#mu^{+}#mu^{-}","F");
  leg->AddEntry(ttbarHist,"t#bar{t}, Single Top","F");
  //leg->AddEntry(dibosonsBkgHist,"WW, WZ, ZZ","F");
  leg->AddEntry(dibosonsBkgHist,"Diboson, #tau#tau","F");
  leg->AddEntry(WjetsBkgHist,"W+jets (FR)","F");
  leg->AddEntry(jetBkgHist,"Di-Jets (data)","F");
  //leg->AddEntry(jetBkgHist,"Jets (FR)","F");
  leg->SetBorderSize(0);
  //leg->SetTextFont(62);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(19);
  leg->SetFillStyle(0);
  leg->SetTextFont(font);
  leg->SetTextSize(0.04);
  leg->Draw();
  //==========================================================
  TPaveText* tText1 = new TPaveText(0.75, 0.92, 0.87, 0.98, "brNDC");
  tText1->SetBorderSize(0);
  tText1->SetFillColor(0);
  tText1->SetFillStyle(0);
  TText *t1 = tText1->AddText("2.673 fb^{-1} (13 TeV)");
  tText1->SetTextSize(0.04);
  tText1->Draw(); 
  //==========================================================
  TPaveText* tText2 = new TPaveText(0.85, 0.86, 0.88, 0.87, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t2 = tText2->AddText("CMS");
  t2 = tText2->AddText("CMS");
  tText2->SetTextSize(0.04);
  tText2->Draw(); 
  //==========================================================
  TPaveText* tText3 = new TPaveText(0.80, 0.81, 0.85, 0.83, "brNDC");
  tText3->SetBorderSize(0);
  tText3->SetFillColor(0);
  tText3->SetFillStyle(0);
  TText *t3 = tText3->AddText("#it{Preliminary}");
  tText3->SetTextSize(0.04);
  tText3->Draw(); 

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
  c1->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.3);
  pad2->SetTickx(1);
  pad2->SetTicky(1);
  pad2->SetLogx(1);
  pad2->SetLogy(1);
  pad2->SetGridx(1);
  pad2->SetGridy(1);
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad
  //========================================================== 
  //                                                              
  //               Get the histograms                                  
  //==========================================================
  
  TFile *file1000 = new TFile("plots/h_ZprimeRecomass_Total_OF.root","READ");   
  TH1F* totalMC = (TH1*) file1000->Get("ZprimeRecomassBinWidthwtOverFlow");
  TH1F* totaldijet= (TH1*) file1000->Get("DataSubwtOverFlow");
  TH1F* totalwjet= (TH1*) file1000->Get("WjetsHistowtOverFlow");

  /*
  TH1F *total;
  total->Add(totalMC);
  total->Add(totaldijet);
  total->Add(totalwjet);
  */

  TFile *file2000 = new TFile("plots/h_ZprimeRecomass_data_OF.root","READ");
  TH1* dataHistTempbar  = (TH1*) file2000->Get("ZprimeRecomassBinWidthwtOverFlow");

  //dataHistTempbar->Divide(total);
  //dataHistTempbar->Draw();
 
  //TFile *file = new TFile("Data-DY-Dibosons-TTbarandTTbarLike-Wjets-Dijets-MC-OS-2673pb.root","READ");
  //TH1* dataHistTempbar = (TH1*) file1->Get("htotalHistoRatiowtOverFlow");
  
  float binWidthNorm=-1;
  int singletopColour    =  TColor::GetColor("#99ccff");
  int jetBkgColour = TColor::GetColor("#ffff66"); 
  int ttbarColour  = TColor::GetColor("#ff6666");
  int bosonColour  = TColor::GetColorDark(3);
  int wjetsColour  = TColor::GetColorDark(5);
  int DYtautauColour  = TColor::GetColorDark(4);
  int font = 42;
  //TGraphAsymmErrors* dataHist = makeDataGraph(dataHistTempbar,binWidthNorm,0);
  //normHistToBinWidth(singletopHist,binWidthNorm);
  dataHistTempbar->SetMarkerStyle(20);
  dataHistTempbar->SetMarkerColor(1);
  dataHistTempbar->SetMarkerSize(0.9);
  //dataHistTempbar->GetXaxis()->SetRange(5,83);

  dataHistTempbar->GetXaxis()->SetTitleOffset(1.1);
  dataHistTempbar->GetYaxis()->SetTitleOffset(0.4);
 
  dataHistTempbar->GetXaxis()->SetTitleSize(0.14);
  dataHistTempbar->GetYaxis()->SetTitleSize(0.12);
  
  dataHistTempbar->GetXaxis()->SetLabelSize(0.12);
  dataHistTempbar->GetYaxis()->SetLabelSize(0.08);  

  dataHistTempbar->GetXaxis()->SetTitle("M(#mu^{+}#mu^{-}) [GeV]");
  //dataHistTempbar->GetXaxis()->SetTitle("M(e#mu) [GeV]");
  dataHistTempbar->GetYaxis()->SetTitle("Data / MC");
  dataHistTempbar->GetXaxis()->SetMoreLogLabels();
  dataHistTempbar->GetXaxis()->SetNoExponent();
  dataHistTempbar->GetXaxis()->SetRangeUser(60.0,2500.0);
  dataHistTempbar->GetYaxis()->SetRangeUser(-1.0,3.0);
  dataHistTempbar->SetMarkerStyle(20);
  dataHistTempbar->SetMarkerSize(0.9);
  dataHistTempbar->Draw("");

  
  TF1* fn1 = new TF1("fn1","pol0",60,1300);
  fn1->SetLineColor(2);
  dataHistTempbar->Fit("fn1","R,smaes");
  //==========================================================
  TPaveText* tText10 = new TPaveText(0.70, 0.40, 0.88, 0.90, "brNDC");
  tText10->SetBorderSize(0);
  tText10->SetFillColor(0);
  tText10->SetFillStyle(0);
  TText *t1 = tText10->AddText("#chi^{2} / n.d.f = 53.91/54");
  TText *t2 = tText10->AddText("Prob = 0.478");
  TText *t3 = tText10->AddText("p0 = 1.009 #pm 0.005");
  tText10->SetTextSize(0.10);
  tText10->Draw(); 

  
  //=================================================================================== 
  c1->Print("Stack-DY-Spring15MCs-Data2015-mass-spectrum-MuMu-OS-2673pb.png","png");
  //c1->Print("Stack-DY-Spring15MCs-Data2015-mass-spectrum-MuMu-OS-2673pb.pdf","pdf");
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  TCanvas *c2 = new TCanvas("c2", "c2",1000,600);
  TPad* spectrumPad=0;
  TPad* ratioPad=0;
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0);
  gStyle->SetStatX(0.457589);
  gStyle->SetStatY(0.312937);
  gStyle->SetStatW(0.29241/2+0.0185);
  gStyle->SetStatH(0.169580+0.05);
  gStyle->SetStatFontSize(0.0402098);
  gStyle->SetStatFont(0.03);
  gStyle->SetFitFormat("5.2g");
  gStyle->SetStatBorderSize(0);
  gStyle->SetStatFontSize(0.040209);
  gStyle->SetStatFontSize(0.035209);
  c2->Range(1.592761,-5.173913,3.533814,6.006211);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetBorderSize(2);
  c2->SetLogx(1);
  c2->SetLogy(1);
  c2->SetTickx(1);
  c2->SetTicky(1);
  c2->SetLeftMargin(0.13);
  c2->SetRightMargin(0.07);
  c2->SetFrameBorderMode(0);
  c2->SetFrameBorderMode(0);
  c2->SetTopMargin(0.085);
  c2->SetBottomMargin(0.11);

  TFile *file1 = new TFile("plots/htotal_root_ZprimeRecomass_OF.root","READ");
  TH1* zeeHist = (TH1*) file1->Get("hfourlepbestmass_4l_afterSel_new_DYwtOverFlow");
 
  //TFile *file2 = new TFile("Data-Dibosons-TTbarandTTbarLike-MC-OS-2673pb.root","READ");
  TH1* ttbarHist = (TH1*) file1->Get("hfourlepbestmass_4l_afterSel_new_TlikewtOverFlow");
  TH1* dibosonsBkgHist = (TH1*) file1->Get("hfourlepbestmass_4l_afterSel_new_diBosonwtOverFlow");
  TH1* dataHistTempbar = (TH1*) file1->Get("htotaldatawtOverFlow");

  //TFile *file3 = new TFile("Wjets-25nsMC-OS-allbins-2673pb.root","READ");
  TH1* WjetsBkgHist = (TH1*) file1->Get("hfourlepbestmass_4l_afterSel_WJetsFromFakeRateFromMCwtOverFlow");
  std::cout<<"nb.Wjets(MC)  = "<<WjetsBkgHist->Integral()<<endl;

  //TFile *file4 = new TFile("DiJets-Data-OS-2673pb-FR.root","READ");
  TH1* jetBkgHist   = (TH1*) file1->Get("hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromDatawtOverFlow");
  std::cout<<"nbQCD(dijets,Data)  = "<<jetBkgHist->Integral()<<endl;

  /*
  zeeHist->Rebin(20);
  ttbarHist->Rebin(20);
  dibosonsBkgHist->Rebin(20);
  dataHistTempbar->Rebin(20);
  WjetsBkgHist->Rebin(20);
  jetBkgHist->Rebin(20);
  */

  float binWidthNorm2=-1;
  int zeeColour    =  TColor::GetColor("#99ccff");
  int jetBkgColour = TColor::GetColor("#ffff66"); 
  int ttbarColour  = TColor::GetColor("#ff6666");
  int bosonColour  = TColor::GetColorDark(3);
  int WjetsColour  = TColor::GetColorDark(5);
  int font = 42;
  float xAxisMin = 60; //72
  float xAxisMax = 3000.0;
  float yAxisMin = 1e-4;
  float yAxisMax = 1e4; 

  

  dataHistTempbar = makeIntHist(dataHistTempbar);
  zeeHist         = makeIntHist(zeeHist);
  ttbarHist       = makeIntHist(ttbarHist);
  dibosonsBkgHist = makeIntHist(dibosonsBkgHist);
  jetBkgHist      = makeIntHist(jetBkgHist);
  WjetsBkgHist    = makeIntHist(WjetsBkgHist);
  TGraphAsymmErrors* dataHist = makeDataGraph(dataHistTempbar,binWidthNorm2,0);
  THStack *axisHist2 = new THStack("axisHist2","");
  zeeHist->SetFillColor(zeeColour);
  zeeHist->SetLineWidth(2);
  zeeHist->SetLineColor(1);
  zeeHist->SetTitle("");
    
  ttbarHist->SetFillColor(ttbarColour);
  ttbarHist->SetLineWidth(2); 
  ttbarHist->SetLineColor(1);
    
  dibosonsBkgHist->SetFillColor(bosonColour);
  dibosonsBkgHist->SetLineWidth(2); 
  dibosonsBkgHist->SetLineColor(1);
   
  jetBkgHist->SetFillColor(jetBkgColour);
  jetBkgHist->SetLineWidth(2);
  jetBkgHist->SetLineColor(1);

  WjetsBkgHist->SetFillColor(WjetsColour);
  WjetsBkgHist->SetLineWidth(2);
  WjetsBkgHist->SetLineColor(1);
  
  
  axisHist2->Add(jetBkgHist,"histo");
  axisHist2->Add(WjetsBkgHist,"histo");
  axisHist2->Add(dibosonsBkgHist,"histo");
  axisHist2->Add(ttbarHist,"histo");
  axisHist2->Add(zeeHist,"histo");
  axisHist2->Draw();
  axisHist2->GetXaxis()->SetTitle("M(#mu^{+}#mu^{-}) [GeV]");
  axisHist2->GetYaxis()->SetTitle("Events #geq M(#mu^{+}#mu^{-}) [GeV]");
  axisHist2->GetXaxis()->SetTitleOffset(1.1);
  axisHist2->GetYaxis()->SetTitleOffset(1.1);
  axisHist2->GetXaxis()->SetTitleSize(0.047);
  axisHist2->GetYaxis()->SetTitleSize(0.047);
  axisHist2->GetXaxis()->SetLabelSize(0.040);
  axisHist2->GetYaxis()->SetLabelSize(0.040);
  axisHist2->GetXaxis()->SetMoreLogLabels();
  axisHist2->GetXaxis()->SetNoExponent();
  axisHist2->GetXaxis()->SetRangeUser(60.0,2500.0);
  axisHist2->SetMinimum(0.01);
  axisHist2->SetMaximum(60000.0);
  
  dataHist->SetMarkerStyle(20);
  dataHist->SetMarkerSize(0.9);
  dataHist->Draw("PZsames");
  //==========================================================
  TLegend *leg = new TLegend(0.56741,0.58671,0.820536,0.83664,NULL,"brNDC");
  leg->AddEntry(dataHist,"Data","PE");
  leg->AddEntry(zeeHist,"#gamma^{*}/Z#rightarrow#mu^{+}#mu^{-}","F");
  leg->AddEntry(ttbarHist,"t#bar{t}, Single Top","F");
  //leg->AddEntry(dibosonsBkgHist,"WW, WZ, ZZ","F");
  leg->AddEntry(dibosonsBkgHist,"Diboson, #tau#tau","F");
  leg->AddEntry(WjetsBkgHist,"W+jets (FR)","F");
  leg->AddEntry(jetBkgHist,"Di-Jets (data)","F");
  //leg->AddEntry(jetBkgHist,"Jets (FR)","F");
  leg->SetBorderSize(0);
  //leg->SetTextFont(62);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(19);
  leg->SetFillStyle(0);
  leg->SetTextFont(font);
  leg->SetTextSize(0.04);
  leg->Draw();
  //==========================================================
  TPaveText* tText1 = new TPaveText(0.75, 0.92, 0.87, 0.98, "brNDC");
  tText1->SetBorderSize(0);
  tText1->SetFillColor(0);
  tText1->SetFillStyle(0);
  TText *t1 = tText1->AddText("2.673 fb^{-1} (13 TeV)");
  tText1->SetTextSize(0.04);
  tText1->Draw(); 
  //==========================================================
  TPaveText* tText2 = new TPaveText(0.85, 0.86, 0.88, 0.87, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t1 = tText2->AddText("CMS");
  tText2->SetTextSize(0.04);
  tText2->Draw(); 
  //==========================================================
  TPaveText* tText3 = new TPaveText(0.80, 0.81, 0.85, 0.83, "brNDC");
  tText3->SetBorderSize(0);
  tText3->SetFillColor(0);
  tText3->SetFillStyle(0);
  TText *t1 = tText3->AddText("#it{Preliminary}");
  tText3->SetTextSize(0.04);
  tText3->Draw(); 
  //=================================================================================== 
  c2->Print("Stack-DY-Spring15MCs-Data2015-cumulative-spectrum-MuMu-OS-2673pb.png","png");
  //c2->Print("Stack-DY-Spring15MCs-Data2015-cumulative-spectrum-MuMu-OS-2673pb.pdf","pdf");

  

}


void fixStats()
{
  TH1* ratioHist = (TH1*) gROOT->FindObject("data");
  TPaveStats* stats = (TPaveStats*) ratioHist->FindObject("stats");
  if(stats)stats->SetTextSize(0.040209);
}

void normaliseToBinWidth(const std::string& inFilename,const std::string& outFilename,float binNormNr)
{
  TFile* inFile = new TFile(inFilename.c_str(),"READ");
  TFile* outFile = new TFile(outFilename.c_str(),"RECREATE");
  std::vector<TH1*> hists;
  hists.push_back((TH1*) inFile->Get("htotaldatawtOverFlow"));
  hists.push_back((TH1*) inFile->Get("hfourlepbestmass_4l_afterSel_new_DYwtOverFlow"));
  hists.push_back((TH1*) inFile->Get("hfourlepbestmass_4l_afterSel_new_diBosonwtOverFlow"));
  hists.push_back((TH1*) inFile->Get("hfourlepbestmass_4l_afterSel_new_TlikewtOverFlow"));
  hists.push_back((TH1*) inFile->Get("hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromDatawtOverFlow"));
  hists.push_back((TH1*) inFile->Get("hfourlepbestmass_4l_afterSel_WJetsFromFakeRateFromMCwtOverFlow"));

  const int NMBINS = hists[0]->GetNbinsX();
  const double MMIN = 0., MMAX = 3000.;
  double logMbins[NMBINS+1];

  for(size_t histNr=0;histNr<hists.size();histNr++){
    for(int binNr=1;binNr<=hists[histNr]->GetNbinsX();binNr++){
      logMbins[binNr] = exp(log(MMIN) + (log(MMAX)-log(MMIN))*binNR/NMBINS);
      cout << logMbins[binNr] << endl;
      binNormNr=logMbins[binNr]-logMbins[binNr-1];
      hists[histNr]->SetBinContent(binNr,hists[histNr]->GetBinContent(binNr)*binNormNr/hists[histNr]->GetBinWidth(binNr));
      hists[histNr]->SetBinError(binNr,hists[histNr]->GetBinError(binNr)*binNormNr/hists[histNr]->GetBinWidth(binNr));
    }
    hists[histNr]->SetDirectory(outFile);
  }
  outFile->Write();
 
  delete inFile;
  delete outFile;
}

void normHistToBinWidth(TH1* hist,float binNormNr)
{
  if(binNormNr>0){
    
    for(int binNr=1;binNr<=hist->GetNbinsX();binNr++){
      hist->SetBinContent(binNr,hist->GetBinContent(binNr)*binNormNr/hist->GetBinWidth(binNr));
      hist->SetBinError(binNr,hist->GetBinError(binNr)*binNormNr/hist->GetBinWidth(binNr));
    }
    
  }
  
}


TH1* makeIntHist(const TH1* hist,bool intIsGreatThan)
{
  TH1* cHist = (TH1*) hist->Clone("cHist");
  cHist->SetDirectory(0);
  cHist->SetName(hist->GetName());
  int maxBin = hist->GetNbinsX()+1;

  for(int binNr=0;binNr<=hist->GetNbinsX();binNr++){
    //if(hist->GetBinContent(binNr) == 0) continue;
    float nrEntries = intIsGreatThan ? hist->Integral(binNr,maxBin) : hist->Integral(0,binNr);
    cHist->SetBinContent(binNr,nrEntries);
  }
  return cHist;
    
}

void makeIntHists(const std::string& inFilename,const std::string& outFilename)
{
  TFile* inFile = new TFile(inFilename.c_str(),"READ");
  TFile* outFile = new TFile(outFilename.c_str(),"RECREATE");
  std::vector<TH1*> hists;
  hists.push_back((TH1*) inFile->Get("dataHist"));
  hists.push_back((TH1*) inFile->Get("zeeTTbarJetBkgHist"));
  hists.push_back((TH1*) inFile->Get("ttbarJetBkgHist"));
  hists.push_back((TH1*) inFile->Get("jetBkgHist"));
  for(size_t histNr=0;histNr<hists.size();histNr++){
    makeIntHist(hists[histNr])->SetDirectory(outFile);
  }
  outFile->Write();
 
  delete inFile;
  delete outFile;
}

TGraphAsymmErrors* makeDataGraph(TH1* dataHist,float normToBinWidth,bool xErrBars)
{
  std::vector<double> xPoint,yPoint,xErrLow,xErrHigh,yErrLow,yErrHigh;
  for(int binNr=1;binNr<=dataHist->GetNbinsX();binNr++){
    if(dataHist->GetBinContent(binNr)==0.0) continue;
    int nrData = dataHist->GetBinContent(binNr);
    float scale = 1;
    if(normToBinWidth>0) scale= normToBinWidth/dataHist->GetBinWidth(binNr);
    const double alpha = (1 - 0.6827)/2;
    const double beta  = (1 - 0.6827)/2;
    double dataLowBound=0;
    double dataHighBound=0;
    if(nrData!=0 || true){
      dataLowBound = nrData==0 ? 0 : 0.5*ROOT::Math::chisquared_quantile_c(1-alpha, 2*nrData);
      dataHighBound = 0.5*ROOT::Math::chisquared_quantile_c(beta, 2*(nrData+1));
    }
    double binCentre= dataHist->GetBinLowEdge(binNr)+0.5*dataHist->GetBinWidth(binNr);
    xPoint.push_back(binCentre);
    if(xErrBars){
      xErrLow.push_back(dataHist->GetBinWidth(binNr)*0.5);
      xErrHigh.push_back(dataHist->GetBinWidth(binNr)*0.5);
    }else{
      xErrHigh.push_back(0);
      xErrLow.push_back(0);
    }
    yPoint.push_back(nrData*scale);
    yErrLow.push_back((nrData-dataLowBound)*scale);
    yErrHigh.push_back((dataHighBound-nrData)*scale);
  }
  TGraphAsymmErrors* resultGraph = new TGraphAsymmErrors(xPoint.size(),&xPoint[0],&yPoint[0],&xErrLow[0],&xErrHigh[0],&yErrLow[0],&yErrHigh[0]);
  return resultGraph;
  
}


float getBkgErr(float mass)
{
  static bool  isInit=false;
  static TF1 ebErr2Func("ebErrFunc","pol4",0,3000);
  static TF1 eeErr2Func("eeErrFunc","pol4",0,3000);
  static TF1 pdfFunc("pdfFunc","pol3",0,3000);
  static TF1 ewkFunc("ewkFunc","pol2",0,3000);
  static TF1 muonErrFunc("muonEffFunc","pol2",0,3000);
  if(!isInit){
    ebErr2Func.SetParameters(19.2,6.79E-3,4.32E-5,9.81E-9,7.18E-12);
    eeErr2Func.SetParameters(21.2,6.79E-3,4.32E-5,9.81E-9,7.18E-12);
    pdfFunc.SetParameters(4.15,1.83E-3,2.68E-6);
    ewkFunc.SetParameters(-1,4.2E-3);
    muonErrFunc.SetParameters(1.045,1.043E-5,3.378E-8);
  }
  
  //if(isMuon){
  //float pdfErr = pdfFunc.Eval(mass);
  //    float ewkErr = ewkFunc.Eval(mass);
  return muonErrFunc.Eval(mass)-1;
  //}
}

void PaintOverflow(TH1 *h)
{
  // This function paint the histogram h with an extra bin for overflows
  
  char* name  = h->GetName();
  char* title = h->GetTitle();
  Int_t nx    = h->GetNbinsX()+1;
  Double_t x1 = h->GetBinLowEdge(1);
  Double_t bw = h->GetBinWidth(nx);
  Double_t x2 = h->GetBinLowEdge(nx)+bw;
  
  // Book a temporary histogram having ab extra bin for overflows
  TH1F *htmp = new TH1F(name, title, nx, x1, x2);
  
  // Fill the new hitogram including the extra bin for overflows
  for (Int_t i=1; i<=nx; i++) {
    htmp->Fill(htmp->GetBinCenter(i), h->GetBinContent(i));
  }
  
  // Fill the underflows
  htmp->Fill(x1-1, h->GetBinContent(0));
  
  // Restore the number of entries
  htmp->SetEntries(h->GetEntries());
  
  // Draw the temporary histogram
  htmp->Draw();
  /*
  TText *t = new TText(x2-bw/2,h->GetBinContent(nx),"Overflow");
  t->SetTextAngle(90);
  t->SetTextAlign(12);
  t->SetTextSize(0.03);;
  t->Draw();
  */
}




/*
  else{
  if(regionCode==0) return sqrt(ebErr2Func.Eval(mass))/100.;
  else if(regionCode==1) return sqrt(eeErr2Func.Eval(mass))/100.;
  else if(regionCode==-1){
  float bkgEBEB = zeeHistEBEB->GetBinContent(AnaFuncs::getBinNr(zeeHistEBEB,mass));
  float bkgEBEE = zeeHistEBEE->GetBinContent(AnaFuncs::getBinNr(zeeHistEBEE,mass));
  float bkgTot = bkgEBEB+bkgEBEE;
  float err = sqrt(ebErr2Func.Eval(mass))/100.*bkgEBEB + sqrt(eeErr2Func.Eval(mass))/100.*bkgEBEE;
  //      float err = sqrt(ebErr2Func.Eval(mass)/100./100.*bkgEBEB*bkgEBEB + eeErr2Func.Eval(mass)/100./100.*bkgEBEE*bkgEBEE);
  if(bkgTot!=0) err/=bkgTot;
  return err;
  }else return 0.;
  }
  
  }


/lustre/cms/store/user/defilip/ZprimeAnalysis/Data2015_ZprimeMuMu_13TeV_merged/SingleMuon_Run2015D-PromptReco-v4_MuonPhys_Nov13_25ns.root
/lustre/cms/store/user/defilip/ZprimeAnalysis/Data2015_ZprimeMuMu_13TeV_merged/SingleMuon_Run2015D-PromptReco-v3_MuonPhys_Nov13_25ns.root
[selgammal@cmssusy2 Data2015_ZprimeMuMu_13TeV_merged]$ 


*/
