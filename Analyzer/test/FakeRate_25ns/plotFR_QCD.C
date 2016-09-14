#include "TProfile.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TF2.h"
#include <string>
using std::string;

Double_t FittingFunc1(Double_t *x,Double_t *par)
{
  double y = x[0];
  Double_t function= par[0]+exp(par[1]+par[2]*y);
  //Double_t function= par[0]+exp(par[1]);
  return function;
}

Double_t FittingFunc2(Double_t *x,Double_t *par)
{
  double y = x[0];
  Double_t function= par[0]*(par[1]+exp(par[2]*y));

  //Double_t function= -log(par[0]*y+par[1]) + par[2];
  return function;
}

Double_t FittingFunc3(Double_t *x,Double_t *par)
{
  double y = x[0];
  Double_t function= par[0]+(par[1]/(y+par[2]));
  //Double_t function= par[0]+par[1]/pow(y+par[2],par[3]);
  return function;
}




void plotFR_QCD(){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetStatFormat("6.4f");
  gStyle->SetFitFormat("6.4f");
  int BoxValue = 11111111; //4680;  
  gStyle->SetOptFit(11);
  gStyle->SetOptDate(0);
  gStyle->SetOptTitle(0);
  //gStyle->SetOptStat(BoxValue);
  gStyle->SetOptStat(0);
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

  TPaveText* tText1 = new TPaveText(0.70, 0.90, 0.90, 0.95, "brNDC");
  tText1->SetBorderSize(0);
  tText1->SetFillColor(0);
  tText1->SetFillStyle(0);
  TText *t1 = tText1->AddText("(13 TeV)");
  tText1->SetTextSize(0.035);
  //tText1->Draw(); 
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TPaveText* tText2 = new TPaveText(0.2, 0.90, 0.4, 0.95, "brNDC");
  tText2->SetBorderSize(0);
  tText2->SetFillColor(0);
  tText2->SetFillStyle(0);
  TText *t2 = tText2->AddText("CMS Spring16 Simulation");
  tText2->SetTextSize(0.035);
  //tText2->Draw(); 


  TLegend *leg = new TLegend(0.17, 0.77, 0.4, 0.85);
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.03;
  leg->SetTextSize(tsize2); 

  
  TFile *_file0 = TFile::Open("ZprimetoMuMu-MC-CMSSW803_FR_QCD_25ns.root");
  

  TH1F *num; _file0->GetObject("h1_Num_Pt_w",num);
  TH1F *den; _file0->GetObject("h1_Den_Pt_w",den);

  TH1F *numMB; _file0->GetObject("h1_Num_Pt_Barrel_w",numMB);
  TH1F *denMB; _file0->GetObject("h1_Den_Pt_Barrel_w",denMB);

  TH1F *numME; _file0->GetObject("h1_Num_Pt_EndCap_w",numME);
  TH1F *denME; _file0->GetObject("h1_Den_Pt_EndCap_w",denME);

  
  double Nbins = num->GetNbinsX();
  cout << "Nbins=" << Nbins << endl;

  int nRebin=10;

  num->Rebin(nRebin);
  den->Rebin(nRebin);

  numMB->Rebin(nRebin);
  denMB->Rebin(nRebin);

  numME->Rebin(nRebin);
  denME->Rebin(nRebin);


  // 
  TH1F *fake=new TH1F("fake","fake",Nbins/nRebin,0.,2000.);
  TH1F *fakeMB=new TH1F("fakeMB","fakeMB",Nbins/nRebin,0.,2000.);
  TH1F *fakeME=new TH1F("fakeME","fakeME",Nbins/nRebin,0.,2000.);
   
  cout << "Bins=" << Nbins/nRebin << endl;

  int* arraysize = new int[1];
  arraysize[0]=Nbins/nRebin;
  Float_t xMB[arraysize[0]],yMB[arraysize[0]],exlMB[arraysize[0]],exhMB[arraysize[0]],eylMB[arraysize[0]],eyhMB[arraysize[0]];


  for (unsigned int i=1; i<=Nbins/nRebin;i++){
    // All Muon
    if (denMB->GetBinContent(i)>0.) {
      fakeMB->SetBinContent(i,double(numMB->GetBinContent(i)/denMB->GetBinContent(i)));

      xMB[i-1]=fakeMB->GetBinCenter(i);
      yMB[i-1]=fakeMB->GetBinContent(i);
      exlMB[i-1]=0.;
      exhMB[i-1]=0.;
      //eyhMB[i-1]=+0.5 + sqrt(fakeMB->GetBinContent(i)+0.25);
      //eylMB[i-1]=-0.5 + sqrt(fakeMB->GetBinContent(i)+0.25); 
      
      if (double(numMB->GetBinContent(i)/denMB->GetBinContent(i))<1.){
	fakeMB->SetBinError(i, sqrt( double(numMB->GetBinContent(i)/denMB->GetBinContent(i))*(1-double(numMB->GetBinContent(i)/denMB->GetBinContent(i)))/double(denMB->GetBinContent(i)) ) );}
      else fakeMB->SetBinError(i, sqrt( double(numMB->GetBinContent(i)/pow(denMB->GetBinContent(i),2)+1./pow(denMB->GetBinContent(i),3))));
      eyhMB[i-1]=fakeMB->GetBinError(i)/2.;
      eylMB[i-1]=fakeMB->GetBinError(i)/2.;      
    }
    else {
      fakeMB->SetBinContent(i,0.);
      fakeMB->SetBinError(i,0.);
      xMB[i-1] = 0.;
      eylMB[i-1] = 0.;
      eyhMB[i-1] = 0.;
      exlMB[i-1]=0.;
      exhMB[i-1]=0.;
    }

  }
 
  //TGraphAsymmErrors *grMB = new TGraphAsymmErrors(Nbins/nRebin,xMB,yMB,exlMB,exhMB,eylMB,eyhMB);
  TGraphAsymmErrors *grMB = new TGraphAsymmErrors();
  grMB->BayesDivide(numMB,denMB);

  //TGraphAsymmErrors *grMB = new TGraphAsymmErrors(numMB,denMB);
  grMB->SetMarkerColor(2);
  grMB->SetMarkerStyle(20);
  grMB->SetMarkerSize(0.95);

  //TF1* fn1 = new TF1("fn1","pol1",10,100);
  //fn1->SetLineColor(1);
  //grMB->Fit("fn1","R,smaes");

  //TF1* fn2 = new TF1("fn2",FittingFunc1,100,2000.0,3);
  //fn2->SetLineColor(1);
  //grMB->Fit("fn2","R");


  // ME
  Float_t xME[arraysize[0]],yME[arraysize[0]],exlME[arraysize[0]],exhME[arraysize[0]],eylME[arraysize[0]],eyhME[arraysize[0]];

  for (unsigned int i=1; i<=Nbins/nRebin;i++){
    // All Muon                                                                                                                                                                            
    if (denME->GetBinContent(i)>0.) {
      fakeME->SetBinContent(i,double(numME->GetBinContent(i)/denME->GetBinContent(i)));

      xME[i-1]=fakeME->GetBinCenter(i);
      yME[i-1]=fakeME->GetBinContent(i);
      exlME[i-1]=0.;
      exhME[i-1]=0.;
      //eyhME[i-1]=+0.5 + sqrt(fakeME->GetBinContent(i)+0.25);                                                                                                                               
      //eylME[i-1]=-0.5 + sqrt(fakeME->GetBinContent(i)+0.25);                                                                                

      if (double(numME->GetBinContent(i)/denME->GetBinContent(i))<1.){
        fakeME->SetBinError(i, sqrt( double(numME->GetBinContent(i)/denME->GetBinContent(i))*(1-double(numME->GetBinContent(i)/denME->GetBinContent(i)))/double(denME->GetBinContent(i)) ) );
      }
      else fakeME->SetBinError(i, sqrt( double(numME->GetBinContent(i)/pow(denME->GetBinContent(i),2)+1./pow(denME->GetBinContent(i),3))));
      eyhME[i-1]=fakeME->GetBinError(i)/2.;
      eylME[i-1]=fakeME->GetBinError(i)/2.;
    }
    else {
      fakeME->SetBinContent(i,0.);
      fakeME->SetBinError(i,0.);
      xME[i-1] = 0.;
      eylME[i-1] = 0.;
      eyhME[i-1] = 0.;
      exlME[i-1]=0.;
      exhME[i-1]=0.;
    }

  }


  //TGraphAsymmErrors *grME = new TGraphAsymmErrors(Nbins/nRebin,xME,yME,exlME,exhME,eylME,eyhME);  
  TGraphAsymmErrors *grME = new TGraphAsymmErrors();
  grME->BayesDivide(numME,denME);

  grME->SetMarkerColor(4);
  grME->SetMarkerStyle(20);
  grME->SetMarkerSize(0.95);
  

  
  TH2F *hframe=NULL,*hframe2=NULL; 
  hframe= new TH2F("hframe","hframe",500,50.,2000.,500,0.0001,2.0);
  hframe->GetYaxis()->SetTitle("Fake rate");
  hframe->GetXaxis()->SetTitle("p_{T} (GeV)");
  hframe->GetXaxis()->SetTitleOffset(1.7);
  hframe->GetYaxis()->SetTitleOffset(1.7);

  TCanvas *c1 = new TCanvas("c1","c1",800,600);

  c1->SetLogy(1);
  c1->SetLogx(1);
  //c1->SetGrid();
  
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  // gPad->SetGridy();
  //gPad->SetGridx();
    
  c1->SetLogy(1);
  c1->SetLogx(1);
  // c1->SetGrid();
  
  grMB->SetLineColor(1);
  grMB->GetXaxis()->SetTitleOffset(1.7);
  grMB->GetYaxis()->SetTitleOffset(1.7);
  grMB->GetXaxis()->SetLabelSize(0.03);
  grMB->GetYaxis()->SetLabelSize(0.03);
  grMB->GetXaxis()->SetTitle("p_{T} (GeV)");
  grMB->GetYaxis()->SetTitle("Fake rate");
  grMB->GetXaxis()->SetTitleOffset(1.5);
  grMB->GetYaxis()->SetTitleOffset(2.0);
  grMB->GetXaxis()->SetTitleSize(0.05);
  grMB->GetXaxis()->SetLabelSize(0.03);
  grMB->GetYaxis()->SetLabelSize(0.03);
  leg->AddEntry(grMB,"barrel","p");

  grME->SetLineColor(1);
  grME->GetXaxis()->SetTitleOffset(1.7);
  grME->GetYaxis()->SetTitleOffset(1.7);
  grME->GetXaxis()->SetLabelSize(0.03);
  grME->GetYaxis()->SetLabelSize(0.03);
  grME->GetXaxis()->SetTitle("p_{T} (GeV)");
  grME->GetYaxis()->SetTitle("Fake rate");
  grME->GetXaxis()->SetTitleOffset(1.5);
  grME->GetYaxis()->SetTitleOffset(2.0);
  grME->GetXaxis()->SetTitleSize(0.05);
  grME->GetXaxis()->SetLabelSize(0.03);
  grME->GetYaxis()->SetLabelSize(0.03);
  leg->AddEntry(grME,"endcap","p");

  hframe->Draw();
  grMB->Draw("Psame");
  grME->Draw("Psame");
  tText1->Draw("same");
  tText2->Draw("same");
  leg->Draw("same");
  
  c1->SaveAs("FR_QCD.png");
  c1->SaveAs("FR_QCD.pdf");
  c1->SaveAs("FR_QCD.eps");
  //c1->SaveAs("FR_QCD.root");

  TFile *fileX = new TFile("FR_QCD.root", "RECREATE");
  fileX->cd();
  numMB->Write();
  denMB->Write();
  numME->Write();
  denME->Write();
  //grMB->Write();
  //grME->Write();
  fileX->Write();
  fileX->Close();


}   


