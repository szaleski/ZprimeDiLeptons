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
TH1* makeUnifiedPlot(const std::string& filename,bool isMuon,int histType,int regionCode);
//helper functions
TGraphAsymmErrors *makeDataGraph(TH1* dataHist,float normToBinWidth,bool xErrBar=true);
void normHistToBinWidth(TH1* hist,float binNormNr);
TH1* makeIntHist(const TH1* hist,bool intIsGreatThan=true);
TGraphAsymmErrors *makeDataGraph  = new TGraphAsymmErrors;
makeDataGraph(TH1* dataHist,float normToBinWidth,bool xErrBar=true);
//float getBkgErr(float mass,bool isMuon,int regionCode,TH1* zeeHistEBEB,TH1* zeeHistEBEE);
float getBkgErr(float mass);
void PaintOverflow(TH1 *h);
//void PaintOverflow(THStack *h);
void SamMacroModBinWidth(){
  float value = 1.0;
  TCanvas *c1 = new TCanvas("c1", "c1",900,600);
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
  pad1->SetBottomMargin(0.008); // Upper and lower plot are joined
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
  TFile *file1 = new TFile("htotal_root_DYtoMuMu_OF.root","READ");
  TH1* zeeHist       = (TH1*) file1->Get("hMassDYAll7");
  TH1* hTTbarDiboson = (TH1*) file1->Get("ElectroWeakBG");
  TH1* jetBkgHist    = (TH1*) file1->Get("hDijetWjetsOF");
  std::cout<<"nb.AllJets(data)  = "<<jetBkgHist->Integral()<<endl;
  TH1* dataHistTempbar = (TH1*) file1->Get("ZprimeRecomassBinWidth");
  //zeeHist->Rebin(value);
  //hTTbarDiboson->Rebin(value);
  //dataHistTempbar->Rebin(value);
  //jetBkgHist->Rebin(value);
  float binWidthNorm=1;
  int zeeColour    =  TColor::GetColor("#99ccff");
  int jetBkgColour = TColor::GetColor("#ffff66"); 
  int hTTbarDibosonColor  = TColor::GetColor("#ff6666");
  int font = 42;
  //float xAxisMin = 60; //72
  //float xAxisMax = 1000.0;
  //float yAxisMin = 1e-4;
  //float yAxisMax = 1e3; 
  TGraphAsymmErrors* dataHist = makeDataGraph(dataHistTempbar,binWidthNorm,0);
  normHistToBinWidth(zeeHist,binWidthNorm);
  normHistToBinWidth(hTTbarDiboson,binWidthNorm);
  normHistToBinWidth(jetBkgHist,binWidthNorm);
  THStack *axisHist = new THStack("axisHist","");
  zeeHist->SetFillColor(zeeColour);
  zeeHist->SetLineWidth(2);
  zeeHist->SetLineColor(1);
  zeeHist->SetTitle("");
  
  hTTbarDiboson->SetFillColor(hTTbarDibosonColor);
  hTTbarDiboson->SetLineWidth(2); 
  hTTbarDiboson->SetLineColor(1);
  
  jetBkgHist->SetFillColor(jetBkgColour);
  jetBkgHist->SetLineWidth(2);
  jetBkgHist->SetLineColor(1);

  axisHist->Add(jetBkgHist,"histo");
  axisHist->Add(hTTbarDiboson,"histo");
  axisHist->Add(zeeHist,"histo");
  axisHist->Draw();

  //axisHist->GetYaxis()->SetTitle("Events / 20 GeV");
  axisHist->GetXaxis()->SetTitle("M(#mu^{+}#mu^{-}) [GeV]");
  axisHist->GetYaxis()->SetTitle("Events/GeV");
  //axisHist->GetXaxis()->SetTitleOffset(1.1);
  //axisHist->GetYaxis()->SetTitleOffset(1.1);
  axisHist->GetXaxis()->SetTitleSize(0.047);
  axisHist->GetYaxis()->SetTitleSize(0.047);
  axisHist->GetXaxis()->SetLabelSize(0.032);
  axisHist->GetYaxis()->SetLabelSize(0.032);
  axisHist->GetXaxis()->SetMoreLogLabels();
  axisHist->GetXaxis()->SetNoExponent();
  axisHist->GetXaxis()->SetRangeUser(70.0,3000.0);
  axisHist->SetMinimum(0.00001);
  axisHist->SetMaximum(80000.0);
  
  dataHist->SetMarkerStyle(20);
  dataHist->SetMarkerSize(0.9);
  dataHist->GetXaxis()->SetRange(5,83);
  dataHist->GetXaxis()->SetTitleSize(0.047);
  dataHist->GetXaxis()->SetTitleOffset(0.9);
  dataHist->GetYaxis()->SetTitleSize(0.047);
  dataHist->GetYaxis()->SetTitleOffset(1.2);
  dataHist->Draw("PZ");
  
  //==========================================================
  TLegend *leg = new TLegend(0.56741,0.62671,0.820536,0.87664,NULL,"brNDC"); //for lumi in plot
  leg->AddEntry(dataHist,"Data","PE");
  leg->AddEntry(zeeHist,"#gamma^{*}/Z#rightarrow#mu^{+}#mu^{-}","F");
  leg->AddEntry(hTTbarDiboson,"t#bar{t}, tW, WW, WZ, ZZ, #tau#tau","F");
  //leg->AddEntry(dibosonsBkgHist,"di-boson, #gamma^{*}/Z#rightarrow#tau^{+}#tau^{-}","F");
  //leg->AddEntry(WjetsBkgHist,"W+jets (FR)","F");
  leg->AddEntry(jetBkgHist,"Jets","F");
  leg->SetBorderSize(0);
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
  TText *t1 = tText1->AddText("2.8 fb^{-1} (13 TeV)");
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

  c1->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0.07);
  pad2->SetBottomMargin(0.35);
  //pad2->SetRightMargin(0.105);
  pad2->SetRightMargin(0.1015);
  pad2->SetTickx(1);
  pad2->SetTicky(1);
  pad2->SetLogx(1);
  pad2->SetLogy(0);
  pad2->SetGridx(1);
  pad2->SetGridy(1);
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad
  //========================================================== 
  //                                                              
  //               Get the histograms                                  
  //==========================================================
  TFile *file = new TFile("Data-DY-Dibosons-TTbarandTTbarLike-Wjets-Dijets-MC-OS-2800pb.root","READ");
  TH1* dataHistTempbar = (TH1*) file->Get("hDivideHisto");
  float binWidthNorm=1;
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
  dataHistTempbar->GetXaxis()->SetRangeUser(70.0,3000.0);
  dataHistTempbar->GetYaxis()->SetRangeUser(-1.0,3.0);
  //dataHistTempbar->GetYaxis()->SetRangeUser(0.0,3.0);
  dataHistTempbar->SetMarkerStyle(20);
  dataHistTempbar->SetMarkerSize(0.9);
  dataHistTempbar->Draw();
  TF1* fn1 = new TF1("fn1","pol0",70,1300);
  fn1->SetLineColor(2);
  dataHistTempbar->Fit("fn1","R,smaes");
  //==========================================================
  TPaveText* tText10 = new TPaveText(0.75, 0.40, 0.88, 0.90, "brNDC");
  tText10->SetBorderSize(0);
  tText10->SetFillColor(0);
  tText10->SetFillStyle(0);
  TText *t1 = tText10->AddText("#chi^{2} / n.d.f = 1.4");
  //TText *t2 = tText10->AddText("Prob = 0.478");
  TText *t3 = tText10->AddText("p0 = 0.990 #pm 0.005");
  tText10->SetTextSize(0.10);
  tText10->Draw(); 
  //=================================================================================== 
  c1->Print("Stack-DY-Spring15MCs-Data2015-mass-spectrum-MuMu-OS-2800pb_logx.C");
  c1->Print("Stack-DY-Spring15MCs-Data2015-mass-spectrum-MuMu-OS-2800pb_logx.png","png");
  //c1->Print("Stack-DY-Spring15MCs-Data2015-mass-spectrum-MuMu-OS-2673pb_logx.pdf","pdf");
  //=========================================================================





  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  TCanvas *c2 = new TCanvas("c2", "c2",900,600);
  TPad* spectrumPad=0;
  TPad* ratioPad=0;
  //c2->Divide(1,2);
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
  TFile *file1 = new TFile("htotal_root_DYtoMuMu_OF.root","READ");
  TH1* zeeHist       = (TH1*) file1->Get("hMassDYAll7");
  TH1* hTTbarDiboson = (TH1*) file1->Get("ElectroWeakBG");
  TH1* jetBkgHist    = (TH1*) file1->Get("hDijetWjetsOF");
  std::cout<<"nb.AllJets(data)  = "<<jetBkgHist->Integral()<<endl;
  TH1* dataHistTempbar = (TH1*) file1->Get("ZprimeRecomassBinWidth");

  //zeeHist->Rebin(value);
  //hTTbarDiboson->Rebin(value);
  //dataHistTempbar->Rebin(value);
  //jetBkgHist->Rebin(value);
  float binWidthNorm2=-1;
  int zeeColour    =  TColor::GetColor("#99ccff");
  int jetBkgColour = TColor::GetColor("#ffff66"); 
  int hTTbarDibosonColor  = TColor::GetColor("#ff6666");
  int font = 42;
  //float xAxisMin = 60; //72
  //float xAxisMax = 1000.0;
  //float yAxisMin = 1e-4;
  //float yAxisMax = 1e3; 
  dataHistTempbar = makeIntHist(dataHistTempbar);
  zeeHist         = makeIntHist(zeeHist);
  hTTbarDiboson   = makeIntHist(hTTbarDiboson);
  jetBkgHist      = makeIntHist(jetBkgHist);
  TGraphAsymmErrors* dataHist = makeDataGraph(dataHistTempbar,binWidthNorm2,0);
  THStack *axisHist2 = new THStack("axisHist2","");
  zeeHist->SetFillColor(zeeColour);
  zeeHist->SetLineWidth(2);
  zeeHist->SetLineColor(1);
  zeeHist->SetTitle("");
  
  hTTbarDiboson->SetFillColor(hTTbarDibosonColor);
  hTTbarDiboson->SetLineWidth(2); 
  hTTbarDiboson->SetLineColor(1);
  
  jetBkgHist->SetFillColor(jetBkgColour);
  jetBkgHist->SetLineWidth(2);
  jetBkgHist->SetLineColor(1);

  axisHist2->Add(jetBkgHist);
  axisHist2->Add(hTTbarDiboson);
  axisHist2->Add(zeeHist);
  axisHist2->Draw("histo");
  

  //axisHist2->Draw();
  axisHist2->GetXaxis()->SetTitle("M(#mu^{+}#mu^{-}) [GeV]");
  axisHist2->GetYaxis()->SetTitle("Events #geq M(#mu^{+}#mu^{-}) [GeV]");
 
  axisHist2->GetXaxis()->SetTitle("M(#mu^{+}#mu^{-}) [GeV]");
  axisHist2->GetYaxis()->SetTitle("Events/GeV");
  axisHist2->GetXaxis()->SetTitleSize(0.047);
  axisHist2->GetYaxis()->SetTitleSize(0.047);
  axisHist2->GetXaxis()->SetLabelSize(0.032);
  axisHist2->GetYaxis()->SetLabelSize(0.032);
  axisHist2->GetXaxis()->SetMoreLogLabels();
  axisHist2->GetXaxis()->SetNoExponent();
  axisHist2->GetXaxis()->SetRangeUser(70.0,3000.0);
  axisHist2->SetMinimum(0.0001);
  axisHist2->SetMaximum(60000.0);
  
  dataHist->SetMarkerStyle(20);
  dataHist->SetMarkerSize(0.9);
  dataHist->Draw("PZsames");
   //==========================================================
  TLegend *leg = new TLegend(0.56741,0.62671,0.820536,0.87664,NULL,"brNDC"); //for lumi in plot
  leg->AddEntry(dataHist,"Data","PE");
  leg->AddEntry(zeeHist,"#gamma^{*}/Z#rightarrow#mu^{+}#mu^{-}","F");
  leg->AddEntry(hTTbarDiboson,"t#bar{t}, tW, WW, WZ, ZZ, #tau#tau","F");
  //leg->AddEntry(dibosonsBkgHist,"di-boson, #gamma^{*}/Z#rightarrow#tau^{+}#tau^{-}","F");
  //leg->AddEntry(WjetsBkgHist,"W+jets (FR)","F");
  leg->AddEntry(jetBkgHist,"Jets","F");
  leg->SetBorderSize(0);
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
  TText *t1 = tText1->AddText("2.8 fb^{-1} (13 TeV)");
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
  //---------------------------------------------------------------------
  c2->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.35);
  //pad2->SetRightMargin(0.105);
  pad2->SetRightMargin(0.1015);
  pad2->SetTickx(1);
  pad2->SetTicky(1);
  pad2->SetLogx(1);
  pad2->SetLogy(0);
  pad2->SetGridx(1);
  pad2->SetGridy(1);
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad
  //========================================================== 
  //                                                              
  //               Get the histograms                                  
  //==========================================================
  TFile *file = new TFile("Data-MCs-cummulative-OS-2800pb.root","READ");
  TH1* dataHistTempbar10 = (TH1*) file->Get("hDivideHisto2");
  dataHistTempbar10->SetMarkerStyle(20);
  dataHistTempbar10->SetMarkerColor(1);
  dataHistTempbar10->SetMarkerSize(0.9);
  //dataHistTempbar10->GetXaxis()->SetRange(5,83);

  dataHistTempbar10->GetXaxis()->SetTitleOffset(1.1);
  dataHistTempbar10->GetYaxis()->SetTitleOffset(0.4);

 
  dataHistTempbar10->GetXaxis()->SetTitleSize(0.14);
  dataHistTempbar10->GetYaxis()->SetTitleSize(0.12);
  
  dataHistTempbar10->GetXaxis()->SetLabelSize(0.12);
  dataHistTempbar10->GetYaxis()->SetLabelSize(0.08);
  

  dataHistTempbar10->GetXaxis()->SetTitle("M(#mu^{+}#mu^{-}) [GeV]");
  //dataHistTempbar10->GetXaxis()->SetTitle("M(e#mu) [GeV]");
  dataHistTempbar10->GetYaxis()->SetTitle("Data / MC");
  dataHistTempbar10->GetXaxis()->SetMoreLogLabels();
  dataHistTempbar10->GetXaxis()->SetNoExponent();
  dataHistTempbar10->GetXaxis()->SetRangeUser(70.0,3000.0);
  dataHistTempbar10->GetYaxis()->SetRangeUser(-1.0,3.0);
  //dataHistTempbar10->GetYaxis()->SetRangeUser(0.0,3.0);
  dataHistTempbar10->SetMarkerStyle(20);
  dataHistTempbar10->SetMarkerSize(0.9);
  dataHistTempbar10->Draw();
  /*
  TF1* fn1 = new TF1("fn1","pol0",60,2300);
  fn1->SetLineColor(2);
  dataHistTempbar10->Fit("fn1","R,smaes");
  */
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  double LikelihoodCut    = 1.0 ;
  //TLine *line = new TLine( LikelihoodCut , 0.0,  LikelihoodCut , 2);
  TLine *line = new TLine( 70.0, LikelihoodCut ,3000,  LikelihoodCut);
  line->SetLineColor(2);
  line->SetLineWidth(3);
  line->SetLineStyle(2);
  line->Draw();
  //==========================================================
  TPaveText* tText10 = new TPaveText(0.70, 0.40, 0.88, 0.90, "brNDC");
  tText10->SetBorderSize(0);
  tText10->SetFillColor(0);
  tText10->SetFillStyle(0);
  TText *t1 = tText10->AddText("#chi^{2} / n.d.f = 53.91/54");
  TText *t2 = tText10->AddText("Prob = 0.478");
  TText *t3 = tText10->AddText("p0 = 1.009 #pm 0.005");
  tText10->SetTextSize(0.10);
  //tText10->Draw(); 
  //=================================================================================== 
  c2->Print("Stack-DY-Spring15MCs-Data2015-cumulative-spectrum-MuMu-OS-2800pb_logx.png","png");
  //c2->Print("Stack-DY-Spring15MCs-Data2015-cumulative-spectrum-MuMu-OS-2673pb_logx.pdf","pdf");
  

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
  hists.push_back((TH1*) inFile->Get("dataHist"));
  hists.push_back((TH1*) inFile->Get("zeeTTbarJetBkgHist"));
  hists.push_back((TH1*) inFile->Get("ttbarJetBkgHist"));
  hists.push_back((TH1*) inFile->Get("jetBkgHist"));
  
  for(size_t histNr=0;histNr<hists.size();histNr++){
    for(int binNr=1;binNr<=hists[histNr]->GetNbinsX();binNr++){
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
//void PaintOverflow(THStack *h)
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
  htmp->Draw("sames");
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
