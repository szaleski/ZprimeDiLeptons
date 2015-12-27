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




void plotFR_data(){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetStatFormat("6.4f");
  gStyle->SetFitFormat("6.4f");
  int BoxValue = 11111111; //4680;  
  gStyle->SetOptFit(1111);
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
  TText *t2 = tText2->AddText("CMS 2015C + 2015D Data");
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

  
  TFile *_file0 = TFile::Open("ZprimetoMuMu-MC-CMSSW745_FR_data.root");

  TH1F *num; _file0->GetObject("h1_Num_Pt_w",num);
  TH1F *den; _file0->GetObject("h1_Den_Pt_w",den);

  TH1F *numMB; _file0->GetObject("h1_Num_Pt_Barrel_w",numMB);
  TH1F *denMB; _file0->GetObject("h1_Den_Pt_Barrel_w",denMB);

  TH1F *numME; _file0->GetObject("h1_Num_Pt_EndCap_w",numME);
  TH1F *denME; _file0->GetObject("h1_Den_Pt_EndCap_w",denME);
  
  char htotal_noQCD_Num_Pt_w[300],htotal_noQCD_Den_Pt_w[300];
  sprintf(htotal_noQCD_Num_Pt_w,"htotal_noQCD_h1_Num_Pt_w.root");       
  sprintf(htotal_noQCD_Den_Pt_w,"htotal_noQCD_h1_Den_Pt_w.root");  

  TFile *file1=NULL;
  file1 = TFile::Open(htotal_noQCD_Num_Pt_w); 
  TH1F *num_noQCD = (TH1F*)file1->Get("fraction_noQCD");

  TFile *file2=NULL;
  file2 = TFile::Open(htotal_noQCD_Den_Pt_w); 
  TH1F *den_noQCD = (TH1F*)file2->Get("fraction_noQCD");
  
  double Nbins = num->GetNbinsX();
  cout << "Nbins=" << Nbins << endl;

  double Nbins_noQCD_Num = num_noQCD->GetNbinsX();
  cout << "Nbins_noQCD_Num=" << Nbins_noQCD_Num << endl;

  double Nbins_noQCD_Den = den_noQCD->GetNbinsX();
  cout << "Nbins_noQCD_Den=" << Nbins_noQCD_Den << endl;

  // Barrel
  char htotal_noQCD_Num_Pt_Barrel_w[300],htotal_noQCD_Den_Pt_Barrel_w[300];
  sprintf(htotal_noQCD_Num_Pt_Barrel_w,"htotal_noQCD_h1_Num_Pt_Barrel_w.root");       
  sprintf(htotal_noQCD_Den_Pt_Barrel_w,"htotal_noQCD_h1_Den_Pt_Barrel_w.root");  

  TFile *file3=NULL;
  file3 = TFile::Open(htotal_noQCD_Num_Pt_Barrel_w); 
  TH1F *num_noQCD_MB = (TH1F*)file3->Get("fraction_noQCD");

  TFile *file4=NULL;
  file4 = TFile::Open(htotal_noQCD_Den_Pt_Barrel_w); 
  TH1F *den_noQCD_MB = (TH1F*)file4->Get("fraction_noQCD");
  

  double Nbins_noQCD_Num_MB = num_noQCD_MB->GetNbinsX();
  cout << "Nbins_noQCD_Num_MB=" << Nbins_noQCD_Num_MB << endl;

  double Nbins_noQCD_Den_MB = den_noQCD_MB->GetNbinsX();
  cout << "Nbins_noQCD_Den_MB=" << Nbins_noQCD_Den_MB << endl;
  

  // Endcap
  char htotal_noQCD_Num_Pt_EndCap_w[300],htotal_noQCD_Den_Pt_EndCap_w[300];
  sprintf(htotal_noQCD_Num_Pt_EndCap_w,"htotal_noQCD_h1_Num_Pt_EndCap_w.root");       
  sprintf(htotal_noQCD_Den_Pt_EndCap_w,"htotal_noQCD_h1_Den_Pt_EndCap_w.root");  

  TFile *file5=NULL;
  file5 = TFile::Open(htotal_noQCD_Num_Pt_EndCap_w); 
  TH1F *num_noQCD_ME = (TH1F*)file5->Get("fraction_noQCD");

  TFile *file6=NULL;
  file6 = TFile::Open(htotal_noQCD_Den_Pt_EndCap_w); 
  TH1F *den_noQCD_ME = (TH1F*)file6->Get("fraction_noQCD");
  

  double Nbins_noQCD_Num_ME = num_noQCD_ME->GetNbinsX();
  cout << "Nbins_noQCD_Num_ME=" << Nbins_noQCD_Num_ME << endl;

  double Nbins_noQCD_Den_ME = den_noQCD_ME->GetNbinsX();
  cout << "Nbins_noQCD_Den_ME=" << Nbins_noQCD_Den_ME << endl;
  
  
  
  if (Nbins!=Nbins_noQCD_Num_MB || Nbins!=Nbins_noQCD_Den_MB){ 
    cout << "Different binning in data and noQCD histo...exiting" << endl;
    exit(1);
  }

  int nRebin=1;

  num->Rebin(nRebin);
  den->Rebin(nRebin);

  numMB->Rebin(nRebin);
  denMB->Rebin(nRebin);

  numME->Rebin(nRebin);
  denME->Rebin(nRebin);

  num_noQCD->Rebin(nRebin);
  den_noQCD->Rebin(nRebin);

  num_noQCD_MB->Rebin(nRebin);
  den_noQCD_MB->Rebin(nRebin);

  num_noQCD_ME->Rebin(nRebin);
  den_noQCD_ME->Rebin(nRebin);

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

  // Subtracking noQCD controbutions
  TH1F *denMB_noEWK=new TH1F("denMB_noEWK","denMB_noEWK",Nbins/nRebin,0.,2000.);
  TH1F *numMB_noEWK=new TH1F("numMB_noEWK","numMB_noEWK",Nbins/nRebin,0.,2000.);

  TH1F *tmp_denMB_noEWK=new TH1F("tmp_denMB_noEWK","tmp_denMB_noEWK",Nbins/nRebin,0.,2000.);
  TH1F *tmp_numMB_noEWK=new TH1F("tmp_numMB_noEWK","tmp_numMB_noEWK",Nbins/nRebin,0.,2000.);
  
  tmp_denMB_noEWK->Multiply(denMB,den_noQCD_MB,1.,1.);
  denMB_noEWK->Add(denMB,tmp_denMB_noEWK,1.,-1.);
  tmp_numMB_noEWK->Multiply(numMB,num_noQCD_MB,1.,1.);
  numMB_noEWK->Add(numMB,tmp_numMB_noEWK,1.,-1.);

    
  cout << denMB_noEWK->GetNbinsX() << " " << numMB_noEWK->GetNbinsX() << endl;
  cout << "Data MB= " << denMB->GetBinContent(1) << " " << numMB->GetBinContent(1) << endl;
  cout << "Histo noQCD= " << den_noQCD_MB->GetBinContent(1) << " " << num_noQCD_MB->GetBinContent(1) << endl;
  cout << "Rescaled= " << denMB_noEWK->GetBinContent(1) << " " << numMB_noEWK->GetBinContent(1) << endl;

  /*
  TCanvas *c0 = new TCanvas("c0","c0",800,600);
  c0->SetLogy(1);
  //denMB->SetMarkerStyle(29);
  //denMB->Draw("P");
  //den_noQCD_MB->SetLineColor(4);
  //den_noQCD_MB->Draw("same");
  //denMB_noEWK->SetFillStyle(1001);
  denMB_noEWK->SetMarkerStyle(29);
  denMB_noEWK->SetMarkerColor(4);
  //denMB_noEWK->Draw("");
  numMB_noEWK->SetMarkerStyle(21);
  numMB_noEWK->SetMarkerColor(2);
  numMB_noEWK->Draw("");
  */
  
  numMB_noEWK->Rebin(20);
  denMB_noEWK->Rebin(20);
  TGraphAsymmErrors *grMB_noEWK = new TGraphAsymmErrors();
  grMB_noEWK->BayesDivide(numMB_noEWK,denMB_noEWK);

  grMB_noEWK->SetMarkerColor(2);
  grMB_noEWK->SetMarkerStyle(20);
  grMB_noEWK->SetMarkerSize(0.95);

  //TF1* fn1 = new TF1("fn1","pol1",50.,900.);
  //fn1->SetLineColor(1);
  //grMB_noEWK->Fit("fn1","R,smaes");

  //TF1* fn2 = new TF1("fn2",FittingFunc1,0,200.0,3);
  //fn2->SetLineColor(1);
  //grMB_noEWK->Fit("fn2","R");

  //TF1* fn3 = new TF1("fn3",FittingFunc1,50.,900.,3);
  //fn3->SetLineColor(1);
  //grMB_noEWK->Fit("fn3","R");


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
  

  // Subtracking noQCD controbutions
  TH1F *denME_noEWK=new TH1F("denME_noEWK","denME_noEWK",Nbins/nRebin,0.,2000.);
  TH1F *numME_noEWK=new TH1F("numME_noEWK","numME_noEWK",Nbins/nRebin,0.,2000.);
  
  TH1F *tmp_denME_noEWK=new TH1F("tmp_denME_noEWK","tmp_denME_noEWK",Nbins/nRebin,0.,2000.);
  TH1F *tmp_numME_noEWK=new TH1F("tmp_numME_noEWK","tmp_numME_noEWK",Nbins/nRebin,0.,2000.);
  
  tmp_denME_noEWK->Multiply(denME,den_noQCD_ME,1.,1.);
  denME_noEWK->Add(denME,tmp_denME_noEWK,1.,-1.);
  tmp_numME_noEWK->Multiply(numME,num_noQCD_ME,1.,1.);
  numME_noEWK->Add(numME,tmp_numME_noEWK,1.,-1.);
    
  cout << denME_noEWK->GetNbinsX() << " " << numME_noEWK->GetNbinsX() << endl;
  cout << "Data ME= " << denME->GetBinContent(1) << " " << numME->GetBinContent(1) << endl;
  cout << "Histo noQCD= " << den_noQCD_ME->GetBinContent(1) << " " << num_noQCD_ME->GetBinContent(1) << endl;
  cout << "Rescaled= " << denME_noEWK->GetBinContent(1) << " " << numME_noEWK->GetBinContent(1) << endl;

  /*
  TCanvas *c0 = new TCanvas("c0","c0",800,600);
  c0->SetLogy(1);
  //denME->SetMarkerStyle(29);
  //denME->Draw("P");
  //den_noQCD_ME->SetLineColor(4);
  //den_noQCD_ME->Draw("same");
  //denME_noEWK->SetFillStyle(1001);
  denME_noEWK->SetMarkerStyle(29);
  denME_noEWK->SetMarkerColor(4);
  //denME_noEWK->Draw("");
  numME_noEWK->SetMarkerStyle(21);
  numME_noEWK->SetMarkerColor(2);
  numME_noEWK->Draw("");
  */

  numME_noEWK->Rebin(10);
  denME_noEWK->Rebin(10);
  TGraphAsymmErrors *grME_noEWK = new TGraphAsymmErrors();
  grME_noEWK->BayesDivide(numME_noEWK,denME_noEWK);

  grME_noEWK->SetMarkerColor(4);
  grME_noEWK->SetMarkerStyle(20);
  grME_noEWK->SetMarkerSize(0.95);

  
  //TF1* fn10 = new TF1("fn10","pol1",10.,100.);
  //fn10->SetLineColor(1);
  //grME_noEWK->Fit("fn10","R,smaes");

  //TF1* fn2 = new TF1("fn2",FittingFunc1,0,110.0,3);
  //fn2->SetLineColor(1);
  //grME_noEWK->Fit("fn2","R");

  //TF1* fn3 = new TF1("fn3",FittingFunc1,100.,900.,3);
  //fn3->SetLineColor(1);
  //grME_noEWK->Fit("fn3","R");

  
  TH2F *hframe=NULL,*hframe2=NULL; 
  hframe= new TH2F("hframe","hframe",500,50.,2000.,500,0.001,2.0);
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
  //grMB->Draw("Psame");
  //grME->Draw("Psame");
  grMB_noEWK->Draw("Psame");
  grME_noEWK->Draw("Psame");
  tText1->Draw("same");
  tText2->Draw("same");
  leg->Draw("same");
  
  c1->SaveAs("FR_data_25ns.png");
  c1->SaveAs("FR_data_25ns.pdf");
  c1->SaveAs("FR_data_25ns.eps");
  //c1->SaveAs("FR_data_25ns.root");

  TFile *fileX = new TFile("FR_data_25ns.root", "RECREATE");
  fileX->cd();
  numMB_noEWK->Write();
  denMB_noEWK->Write();
  numME_noEWK->Write();
  denME_noEWK->Write();
  //grMB_noEWK->Write();
  //grME_noEWK->Write();
  fileX->Write();
  fileX->Close();

}   


