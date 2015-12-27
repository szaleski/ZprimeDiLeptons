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

//#include "SHUtility/HistFuncs.hh"
//#include "Utility/AnaFuncs.hh"

// to load in .L zprimeUnifiedPlots.C+
//main function, muon option changes labels
//integral plot sets ranges and disables, enables normalising to a consistant bin width
//it assumes 4 histograms of format
//dataHist : data histogram, it will be make it with proper error bars
//zeeTTbarJetBkgHist : zee+ ttbar + jet backgrounds summed together
//ttbarJetBkgHist : ttbar + jet backgrounds summed together
//jetBkgHist : jet background
//histType = 0 (or really not 1,2,3) standard mass spectrum
//         = 1 integral mass spectrum
//         = 2 data-mc / mc plot
//         = 3 data-mc / mc plot (120-200 GeV)
//         = 4 linear mass spectrum 200 -2000 GeV
TH1* makeUnifiedPlot(const std::string& filename,bool isMuon,int histType,int regionCode);

//helper functions
void normHistToBinWidth(TH1* hist,float binNormNr);
TH1* makeIntHist(const TH1* hist,bool intIsGreatThan=true);
TGraphAsymmErrors* makeDataGraph(TH1* dataHist,float normToBinWidth,bool xErrBar=true);
float getBkgErr(float mass,bool isMuon,int regionCode,TH1* zeeHistEBEB,TH1* zeeHistEBEE);

TH1* makeUnifiedPlot(const std::string& filename,bool isMuon,int histType,int regionCode)
{
  bool isIntegral  = histType==1;
  bool isRatio = histType==2 || histType==3;
  bool isRatioControlRegion = histType==3;
  bool isLinear = histType==4;
  bool isRatioBelow = histType==5;

  int zeeColour=  TColor::GetColor("#99ccff");
  int jetBkgColour = TColor::GetColor("#ffff66"); 
  int ttbarColour = TColor::GetColor("#ff6666");
  int font = 42; //62
  float xAxisMin = 60; //72
  float xAxisMax = 2199.;
  float yAxisMin = 1e-4;
  float yAxisMax = 1e6;  

  if(isIntegral){
    yAxisMin = 1e-2;
    yAxisMax = 5e6;
    xAxisMin=50;
    xAxisMax=1999.;
  }
  if(isLinear){
    yAxisMin = 1e-3;
    yAxisMax = 5e3;
    xAxisMin=200;
    xAxisMax=1999.;
  }
 

  float binWidthNorm=1;
  if(isIntegral || isRatio || isLinear) binWidthNorm=-1;

  TCanvas *c1 = new TCanvas("c1", "c1",900,600);
  if(isRatioBelow) c1 =new TCanvas("c1","c1",900,900);
  TPad* spectrumPad=0;
  TPad* ratioPad=0;
  // TCanvas *c1 = new TCanvas("c1", "c1",5,24,900,600);
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0);
  gStyle->SetStatX(0.457589);
  gStyle->SetStatY(0.312937);
  gStyle->SetStatW(0.29241/2+0.0185);
  gStyle->SetStatH(0.169580+0.05);
  gStyle->SetStatFontSize(0.0402098);
  gStyle->SetStatFont(font);
  gStyle->SetFitFormat("5.2g");
  gStyle->SetStatBorderSize(0);
  gStyle->SetStatFontSize(0.040209);
  gStyle->SetStatFontSize(0.035209);
  c1->Range(1.592761,-5.173913,3.533814,6.006211);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetLogx(!isRatio && !isIntegral && !isLinear);
  c1->SetLogy(!isRatio);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLeftMargin(0.13);
  c1->SetRightMargin(0.07);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);

  c1->SetTopMargin(0.085);
  c1->SetBottomMargin(0.11);
  
  std::string regionStr;
  if(regionCode==0) regionStr="EBEB";
  else if(regionCode==1) regionStr="EBEE";
  else if(regionCode==2) regionStr="EEEE";

  TFile* file = new TFile(filename.c_str());
  TH1* dataHistTemp = (TH1*) file->Get(("dataHist"+regionStr).c_str());
  
  TH1* zeeHist = (TH1*) file->Get(("zeeTTbarJetBkgHist"+regionStr).c_str());
  if(zeeHist==0) zeeHist = (TH1*) file->Get(("zmumuTTbarJetBkgHist"+regionStr).c_str()); //not actually valid but done for completeness
  TH1* ttbarHist = (TH1*) file->Get(("ttbarJetBkgHist"+regionStr).c_str());
  TH1* jetBkgHist = (TH1*) file->Get(("jetBkgHist"+regionStr).c_str());

  
  TH1* zeeHistEBEB =  (TH1*) file->Get("zeeTTbarJetBkgHistEBEB");
  TH1* zeeHistEBEE =  (TH1*) file->Get("zeeTTbarJetBkgHistEBEE");
  

  //little hack to deal with muons not always having jet bkg
  //this should be harmless and its dy hist chosen so if its not it will be obiovus
  if(jetBkgHist==0) jetBkgHist = (TH1*) zeeHist->Clone("jetBkgHist");

 //  if(isMuon){
//     dataHistTemp = (TH1*) file->Get("data");
//     zeeHist = (TH1*) file->Get("zdy");
//     ttbarHist = (TH1*) file->Get("prompt");
//     jetBkgHist = (TH1*) ttbarHist->Clone("jetBkgHist");jetBkgHist->Reset();
//   }
    
  if(isIntegral){
    dataHistTemp = makeIntHist(dataHistTemp);
    zeeHist = makeIntHist(zeeHist);
    ttbarHist = makeIntHist(ttbarHist);
    jetBkgHist = makeIntHist(jetBkgHist);
  }

  TGraphAsymmErrors* dataHist = makeDataGraph(dataHistTemp,binWidthNorm,0);
  
  normHistToBinWidth(zeeHist,binWidthNorm);
  normHistToBinWidth(ttbarHist,binWidthNorm);
  normHistToBinWidth(jetBkgHist,binWidthNorm);

  if(isRatioBelow){
  
    c1->cd();
    spectrumPad = new TPad("spectrumPad", "newpad",0.01,0.33,0.99,0.99);
    spectrumPad->Draw(); 
    spectrumPad->cd();
    //    spectrumPad->SetTopMargin(0.2);
    // spectrumPad->SetBottomMargin(0.01);
    //spectrumPad->SetRightMargin(0.1);
  }


  TH1* axisHist = new TH1D("axisHist","",2100,0,2100);

  zeeHist->SetFillColor(zeeColour);
  zeeHist->SetLineWidth(2);
  zeeHist->SetLineColor(1);
  zeeHist->SetTitle("");
  if(isMuon) axisHist->GetXaxis()->SetTitle("m(#mu#mu) [GeV]");
  else axisHist->GetXaxis()->SetTitle("m(ee) [GeV]");
  axisHist->GetXaxis()->SetRange(5,84);
  axisHist->GetXaxis()->SetMoreLogLabels();
  axisHist->GetXaxis()->SetNoExponent();  
  axisHist->GetXaxis()->SetLabelSize(0.05); 
  axisHist->GetYaxis()->SetLabelSize(0.05);
  // axisHist->GetXaxis()->SetTitleSize(0.047);
  //axisHist->GetXaxis()->SetTitleOffset(1.1);  
  axisHist->GetXaxis()->SetTitleSize(0.05);
  axisHist->GetXaxis()->SetTitleOffset(1.0);
  if(!isIntegral) {
    if(!isLinear) axisHist->GetYaxis()->SetTitle(" Events / GeV");
    else {
      int binWidth = zeeHist->GetBinWidth(1);
      std::ostringstream yAxisTitle;
      yAxisTitle<<" Events / "<<binWidth<<" GeV";
      axisHist->GetYaxis()->SetTitle(yAxisTitle.str().c_str());
    }
  }else{
    if(isMuon) axisHist->GetYaxis()->SetTitle(" Events #geq m(#mu#mu)");
    else axisHist->GetYaxis()->SetTitle(" Events #geq m(ee)");
  }
  //axisHist->GetYaxis()->SetTitleSize(0.047);
  //axisHist->GetYaxis()->SetTitleOffset(1.2);
  axisHist->GetYaxis()->SetTitleSize(0.05);
  axisHist->GetYaxis()->SetTitleOffset(1.1);
  axisHist->GetYaxis()->SetRangeUser(yAxisMin,yAxisMax);
  axisHist->GetXaxis()->SetRangeUser(xAxisMin,xAxisMax);
  axisHist->GetYaxis()->SetTitleFont(font);
  axisHist->GetYaxis()->SetLabelFont(font);
  axisHist->GetXaxis()->SetTitleFont(font);
  axisHist->GetXaxis()->SetLabelFont(font);
  axisHist->Draw("HIST");
  
  zeeHist->Draw("SAME HIST");
  
  ttbarHist->SetFillColor(ttbarColour);
  ttbarHist->SetLineWidth(2); 
  ttbarHist->SetLineColor(1);
  ttbarHist->GetXaxis()->SetRange(5,83);
  ttbarHist->GetXaxis()->SetTitleSize(0.047);
  ttbarHist->GetXaxis()->SetTitleOffset(0.9);
  ttbarHist->GetYaxis()->SetTitleSize(0.047);
  ttbarHist->GetYaxis()->SetTitleOffset(1.2);
  ttbarHist->Draw(" SAME HIST");
 

  jetBkgHist->SetFillColor(jetBkgColour);
  jetBkgHist->SetLineWidth(2);
  jetBkgHist->SetLineColor(1);
  jetBkgHist->GetXaxis()->SetRange(5,83);
  jetBkgHist->GetXaxis()->SetTitleSize(0.047);
  jetBkgHist->GetXaxis()->SetTitleOffset(0.9);
  jetBkgHist->GetYaxis()->SetTitleSize(0.047);
  jetBkgHist->GetYaxis()->SetTitleOffset(1.2);
  jetBkgHist->Draw("SAME HIST");
  
 
  dataHist->SetMarkerStyle(20);
  dataHist->SetMarkerSize(1.1);
  dataHist->GetXaxis()->SetRange(5,83);
  dataHist->GetXaxis()->SetTitleSize(0.047);
  dataHist->GetXaxis()->SetTitleOffset(0.9);
  dataHist->GetYaxis()->SetTitleSize(0.047);
  dataHist->GetYaxis()->SetTitleOffset(1.2);
  
  dataHist->Draw("PZ");

 


  TLegend *leg = new TLegend(0.516741,0.578671,0.870536,0.835664,NULL,"brNDC"); //for lumi in plot
  //TLegend *leg = new TLegend(0.522321,0.620629,0.876116,0.877622,NULL,"brNDC"); //for lumi above plot
  leg->AddEntry(dataHist,"Data","PE");
  if(!isMuon) leg->AddEntry(zeeHist,"#gamma/Z#rightarrowe^{+}e^{-}","F");
  else leg->AddEntry(zeeHist,"#gamma/Z#rightarrow#mu^{+}#mu^{-}","F");
  leg->AddEntry(ttbarHist,"t#bar{t}, tW, WW, WZ, ZZ, #tau#tau","F");
  leg->AddEntry(jetBkgHist,"Jets (data)","F");
  //else leg->AddEntry(jetBkgHist,"jets","F");
  leg->SetBorderSize(0);
  //leg->SetTextFont(62);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(19);
  leg->SetFillStyle(0);
  leg->SetTextFont(font);
  leg->Draw();

  if(isRatioBelow){
    c1->cd();
    std::cout <<"ratio beloiw "<<std::endl;
    ratioPad = new TPad("ratioPad", "newpad",0.01,0.01,0.99,0.37);
    ratioPad->Draw();
    ratioPad->cd();
    ratioPad->SetTopMargin(0.05);
    ratioPad->SetBottomMargin(0.5);
    //    ratioPad->SetRightMargin(0.1);
    ratioPad->SetFillStyle(0);
    // TH1* ratioHist = (TH1*) dataHistTemp->Clone("ratioHist");
    //normHistToBinWidth(ratioHist,binWidthNorm);
    
    // ratioHist->Divide(zeeHist); 

    TGraphAsymmErrors* ratioHistErrs = new TGraphAsymmErrors(*dataHist);
    for(int i=0;i<ratioHistErrs->GetN();i++){
      float bkg = zeeHist->GetBinContent(i+1);
      ratioHistErrs->SetPoint(i,ratioHistErrs->GetX()[i],bkg!=0 ? ratioHistErrs->GetY()[i]/bkg-1 : 0.);
      ratioHistErrs->SetPointEYlow(i,bkg!=0 ? ratioHistErrs->GetEYlow()[i]/bkg : 0.);
      ratioHistErrs->SetPointEYhigh(i,bkg!=0 ? ratioHistErrs->GetEYhigh()[i]/bkg : 0.);   
      ratioHistErrs->SetPointEXlow(i,0);
      ratioHistErrs->SetPointEXhigh(i,0);
    }
    //for(int i=0;i<ratioHist->GetNbinsX()+2;i++) ratioHist->SetBinContent(i,ratioHist->GetBinContent(i)-1);

    //   TH1* ratioHist = HistFuncs::makeDataMinusBkgPlot(dataHistTemp,zeeHist,0);
    
    TH1* axisHistRatio = (TH1*) axisHist->Clone("axisHistRatio");
    float fontScale = spectrumPad->GetHNDC() / ratioPad->GetHNDC();
    axisHistRatio->GetYaxis()->SetTitleSize(0.05*fontScale); 
    axisHistRatio->GetYaxis()->SetLabelSize(0.05*fontScale);
    axisHistRatio->GetXaxis()->SetTitleSize(0.05*fontScale); 
    axisHistRatio->GetXaxis()->SetLabelSize(0.05*fontScale);
    axisHistRatio->GetYaxis()->SetRangeUser(-0.4,0.4);
    axisHistRatio->GetYaxis()->SetTitle("(data-bkg)/bkg");
    axisHistRatio->GetYaxis()->SetLabelOffset(0.015);
    axisHistRatio->GetYaxis()->SetTitleOffset(0.6);
    axisHistRatio->SetLineColor(0);
    axisHistRatio->SetMarkerColor(0);
    axisHistRatio->Draw();
    axisHist->GetXaxis()->SetTitle("");

    TH1* errBandHist  = (TH1*) axisHistRatio->Clone("errBandHist");
 
    for(int binNr=0;binNr<errBandHist->GetNbinsX()+2;binNr++){
     
      float mass = errBandHist->GetBinLowEdge(binNr)+errBandHist->GetBinWidth(binNr)/2; 
      float err = getBkgErr(mass,isMuon,regionCode,zeeHistEBEB,zeeHistEBEE);
      // float err = sqrt(ebErrFunc.Eval(mass))/100.;
      errBandHist->SetBinError(binNr,err);
    }
    errBandHist->Draw("SAME E2");
    errBandHist->SetFillColor(5);
    errBandHist->SetLineColor(0);
    errBandHist->SetMarkerColor(0);
   //  ratioHist->SetMarkerStyle(8);
//     ratioHist->SetMarkerSize(0.5);
//     ratioHist->SetLineColor(1);
//     ratioHist->SetLineWidth(1);
//     ratioHist->SetMarkerColor(1);
    ratioHistErrs->SetMarkerStyle(8);
    ratioHistErrs->SetMarkerSize(0.8);
    ratioHistErrs->SetLineColor(1);
    ratioHistErrs->SetLineWidth(1);
    ratioHistErrs->SetMarkerColor(1);
    ratioHistErrs->Draw("PZ");
    // ratioHist->Draw("SAME");
  
    TLine* zeroLine = new TLine(60,0,2099,0);
    zeroLine->SetLineStyle(3);
    zeroLine->SetLineColor(1);
    zeroLine->SetLineWidth(1);
    zeroLine->Draw();
    ratioPad->RedrawAxis();
    spectrumPad->cd();
    spectrumPad->RedrawAxis();
    spectrumPad->Draw();
  }


  if(isRatio){
    TH1* ratioHist;
    //= HistFuncs::makeDataMinusBkgPlot(dataHistTemp,zeeHist,10);
    ratioHist->SetTitle("");
    ratioHist->GetXaxis()->SetTitle(axisHist->GetXaxis()->GetTitle());
    ratioHist->GetYaxis()->SetTitle("(data-bkg)/bkg");
    ratioHist->GetYaxis()->SetTitleSize(0.05);
    ratioHist->GetYaxis()->SetTitleOffset(1.1); 
    float maxFit = ratioHist->GetBinLowEdge(ratioHist->GetNbinsX());
    float minFit=200;
    if(isRatioControlRegion){
      minFit=120;
      maxFit=200;  
      ratioHist->GetYaxis()->SetRangeUser(-0.2,0.2);
    }else ratioHist->GetYaxis()->SetRangeUser(-1,1);
    ratioHist->GetXaxis()->SetRangeUser(minFit,maxFit-0.001);
    ratioHist->GetYaxis()->SetTitleFont(font);
    ratioHist->GetYaxis()->SetLabelFont(font);
    ratioHist->GetXaxis()->SetTitleFont(font);
    ratioHist->GetXaxis()->SetLabelFont(font); 
    ratioHist->GetXaxis()->SetMoreLogLabels();
    ratioHist->GetXaxis()->SetNoExponent();   
    ratioHist->SetMarkerStyle(8);
    ratioHist->SetLineWidth(2);
    ratioHist->Fit("pol0","","",minFit,maxFit); 
    ratioHist->GetFunction("pol0")->SetLineColor(4);
    ratioHist->Draw();
    TPaveStats* stats = (TPaveStats*) ratioHist->FindObject("stats");
    std::cout <<"stats "<<stats<<std::endl;
    if(stats)stats->SetTextSize(0.040209);
   
    float lastBinContentRounded = static_cast<int>(ratioHist->GetBinContent(ratioHist->GetNbinsX())*100+0.5)/100.;
    float lastBinErrorRounded =  static_cast<int>(ratioHist->GetBinError(ratioHist->GetNbinsX())*100+0.5)/100.;
    std::ostringstream lastBinText;
    lastBinText<<"last bin ("<<maxFit<<" - 2000): "<<lastBinContentRounded<<" #pm "<<lastBinErrorRounded;
    TPaveLabel* lastBinLabel = new TPaveLabel(0.524554,0.15035,0.890625,0.237762,lastBinText.str().c_str(),"brNDC");
    lastBinLabel->SetFillColor(0);
    lastBinLabel->SetFillStyle(0);
    lastBinLabel->SetTextSize(0.48);
    lastBinLabel->SetTextFont(font);
    lastBinLabel->SetBorderSize(0);
    if(maxFit!=200) lastBinLabel->Draw();
  }
   
  //std::string lumiText("CMS Preliminary, 8 TeV, 19.7 fb^{-1}");
  //if(isMuon) lumiText ="CMS Preliminary, 8 TeV, 20.6 fb^{-1}";  
  std::string lumiText("19.7 fb^{-1} (8 TeV)");
  if(isMuon) lumiText ="20.6 fb^{-1} (8 TeV)";
 
  //TPaveLabel *pl = new TPaveLabel(0.534598,0.809441,0.877232,0.905594,lumiText.c_str(),"brNDC"); //for v6, in plot
  TPaveLabel *pl = new TPaveLabel(0.692,0.902,0.994,0.997,lumiText.c_str(),"brNDC"); //for paper, in plot
  
  pl->SetBorderSize(0);
  pl->SetFillColor(0);
  pl->SetFillStyle(0);
  pl->SetTextFont(font);
  pl->SetTextSize(0.44);
  pl->Draw();
  pl = new TPaveLabel(0.712,0.811,0.969,0.907,"CMS","brNDC"); //for paper, in plot

  if(isRatio) pl = new TPaveLabel(0.077,0.811,0.335,0.907,"CMS","brNDC"); //for paper, in plot
  pl->SetBorderSize(0);
  pl->SetFillColor(0);
  pl->SetFillStyle(0);
  pl->SetTextFont(61);
  pl->SetTextSize(0.44/0.75);
  pl->Draw();
  
  pl = new TPaveLabel(0.126,0.745,0.384,0.841,"unpublished","brNDC"); //for paper, in plot

  pl->SetBorderSize(0);
  pl->SetFillColor(0);
  pl->SetFillStyle(0);
  pl->SetTextFont(51);
  pl->SetTextSize(0.44/0.75);
  if(isRatio) pl->Draw();

  pl = new TPaveLabel(0.3490783,0.8007519,0.593318,0.8890977,"CMS Preliminary, #sqrt{s} = 8 TeV","brNDC");
  pl->SetBorderSize(0);
  pl->SetFillColor(0);
  pl->SetFillStyle(0);
  pl->SetTextFont(font);
  pl->SetTextSize(0.4680851);
  pl->SetTextSize(0.55);
  // pl->Draw();
 
  if(!isRatio) pl = new TPaveLabel(0.185268,0.812937,0.527902,0.907343,"barrel-barrel","brNDC"); //for v6, in plot
  else {
    pl = new TPaveLabel(0.170759,0.741259,0.512277,0.839161,"barrel-barrel","brNDC"); //for v6, in plot
    pl = new TPaveLabel(0.170759+0.56,0.741259,0.512277+0.56,0.839161,"barrel-barrel","brNDC"); //for v6, in plot
    pl->SetTextAlign(12);
  }  
  if(regionCode==1) pl->SetLabel("barrel-endcap");
  if(regionCode==2) pl->SetLabel("endcap-endcap");
  pl->SetBorderSize(0);
  pl->SetFillColor(0);
  pl->SetFillStyle(0);
  pl->SetTextFont(font);
  pl->SetTextSize(0.425926 ); 
  //pl->SetTextSize(0.45 );
  if(regionCode!=-1 && !isMuon) pl->Draw();
  
  if(isRatio){
    pl = new TPaveLabel(0.172991,0.802448,0.361607,0.898601,"signal region","brNDC");
    pl = new TPaveLabel(0.742,0.811,0.9308,0.905595,"signal region","brNDC");
    if(isRatioControlRegion) pl->SetLabel("control region");
    pl->SetBorderSize(0);
    pl->SetFillColor(0);
    pl->SetFillStyle(0);
    pl->SetTextFont(font);
    pl->SetTextAlign(12);
    pl->SetTextSize(0.425926 ); 
    pl->SetTextSize(0.45 );
    pl->Draw();
  }
  c1->RedrawAxis();
    
  return zeeHist;
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


float getBkgErr(float mass,bool isMuon,int regionCode,TH1* zeeHistEBEB,TH1* zeeHistEBEE)
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

  if(isMuon){
    //float pdfErr = pdfFunc.Eval(mass);
    //    float ewkErr = ewkFunc.Eval(mass);
    return muonErrFunc.Eval(mass)-1;
  }else{
    if(regionCode==0) return sqrt(ebErr2Func.Eval(mass))/100.;
    else if(regionCode==1) return sqrt(eeErr2Func.Eval(mass))/100.;
    else if(regionCode==-1){
      float bkgEBEB;
      //= zeeHistEBEB->GetBinContent(AnaFuncs::getBinNr(zeeHistEBEB,mass));
      float bkgEBEE;
      //= zeeHistEBEE->GetBinContent(AnaFuncs::getBinNr(zeeHistEBEE,mass));
      float bkgTot = bkgEBEB+bkgEBEE;
      float err = sqrt(ebErr2Func.Eval(mass))/100.*bkgEBEB + sqrt(eeErr2Func.Eval(mass))/100.*bkgEBEE;
			//      float err = sqrt(ebErr2Func.Eval(mass)/100./100.*bkgEBEB*bkgEBEB + eeErr2Func.Eval(mass)/100./100.*bkgEBEE*bkgEBEE);
      if(bkgTot!=0) err/=bkgTot;
      return err;
    }else return 0.;
  }
}

