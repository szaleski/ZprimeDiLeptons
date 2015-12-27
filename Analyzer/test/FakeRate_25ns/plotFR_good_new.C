{

  //#include "ZZStyle.C"

  //TStyle * style = getStyle("ZZ");
  //style->cd();
  //style->SetNdivisions(512, "X");
  //style->SetNdivisions(512, "Y");

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
  TText *t1 = tText2->AddText("CMS Spring15 Simulation");
  tText2->SetTextSize(0.035);
  //tText2->Draw(); 


  TLegend *leg = new TLegend(0.65, 0.20, 0.85, 0.35);
  leg->SetBorderSize(0.0);
  leg->SetMargin(0.3);
  leg->SetFillColor(0);
  leg->SetFillStyle(10);
  leg->SetLineColor(0);
  Float_t tsize2 = 0.035;
  leg->SetTextSize(tsize2); 

  
  TFile *_file0 = TFile::Open("ZprimetoMuMu-MC-CMSSW745_FR.root");

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

  TH1F *fake=new TH1F("fake","fake",Nbins/nRebin,0.,2000.);
  TH1F *fakeMB=new TH1F("fakeMB","fakeMB",Nbins/nRebin,0.,2000.);
  TH1F *fakeME=new TH1F("fakeME","fakeME",Nbins/nRebin,0.,2000.);
   
  cout << "Bins=" << Nbins/nRebin << endl;

  int* arraysize = new int[1];
  arraysize[0] =  Nbins/nRebin;
  Float_t x[arraysize[0]],y[arraysize[0]],exl[arraysize[0]],exh[arraysize[0]],eyl[arraysize[0]],eyh[arraysize[0]];
  delete [] arraysize;

  for (unsigned int i=1; i<=Nbins/nRebin;i++){
    // All Muon
    if (den->GetBinContent(i)>0.) {
      fake->SetBinContent(i,double(num->GetBinContent(i)/den->GetBinContent(i)));

      x[i-1]=fake->GetBinCenter(i);
      y[i-1]=fake->GetBinContent(i);
      eyh[i-1]=+0.5 + sqrt(fake->GetBinContent(i)+0.25);
      eyl[i-1]=-0.5 + sqrt(fake->GetBinContent(i)+0.25); 
      
      //if (double(num->GetBinContent(i)/den->GetBinContent(i))<1.){
      //fake->SetBinError(i, sqrt( double(num->GetBinContent(i)/den->GetBinContent(i))*(1-double(num->GetBinContent(i)/den->GetBinContent(i)))/double(den->GetBinContent(i)) ) );}
      //else fake->SetBinError(i, sqrt( double(num->GetBinContent(i)/pow(den->GetBinContent(i),2)+1./pow(den->GetBinContent(i),3))));
    }
    else {
      fake->SetBinContent(i,0.);
      fake->SetBinError(i,0.);
      x[i-1] = 0.;
      eyl[i-1] = 0.;
      eyh[i-1] = 0.;
    }

  }
  
 
  TGraphAsymmErrors *grMB = new TGraphAsymmErrors(numMB,denMB);
  grMB->SetMarkerColor(1);
  grMB->SetMarkerStyle(20);
  grMB->SetMarkerSize(0.95);
   
  for (unsigned int i=1; i<=Nbins/nRebin;i++){
    // Muon Barrel
    if (denMB->GetBinContent(i)>0.) { 
      fakeMB->SetBinContent(i,double(numMB->GetBinContent(i)/denMB->GetBinContent(i)));
      //if (double(numMB->GetBinContent(i)/denMB->GetBinContent(i))<1.){ 
      //fakeMB->SetBinError(i, sqrt( double(numMB->GetBinContent(i)/denMB->GetBinContent(i))*(1-double(numMB->GetBinContent(i)/denMB->GetBinContent(i)))/double(denMB->GetBinContent(i)) ) ;
      //}
      //else fakeMB->SetBinError(i, sqrt( double(numMB->GetBinContent(i)/pow(denMB->GetBinContent(i),2)+1./pow(denMB->GetBinContent(i),3))));
    }
    else {
      fakeMB->SetBinContent(i,0.);
      fakeMB->SetBinError(i,0.);
    }
  }
  
  
  TGraphAsymmErrors *grME = new TGraphAsymmErrors(numME,denME);
  grME->SetMarkerColor(1);
  grME->SetMarkerStyle(20);
  grME->SetMarkerSize(0.95);
  
  for (unsigned int i=1; i<=Nbins/nRebin;i++){
    // Muon EndCap
    if (denME->GetBinContent(i)>0.) {
      fakeME->SetBinContent(i,double(numME->GetBinContent(i)/denME->GetBinContent(i)));
      if (double(numME->GetBinContent(i)/denME->GetBinContent(i))<1.){ 
	fakeME->SetBinError(i, sqrt( double(numME->GetBinContent(i)/denME->GetBinContent(i))*(1-double(numME->GetBinContent(i)/denME->GetBinContent(i)))/double(denME->GetBinContent(i)) ) );
      }
      else fakeME->SetBinError(i, sqrt( double(numME->GetBinContent(i)/pow(denME->GetBinContent(i),2)+1./pow(denME->GetBinContent(i),3))));
    }		
    else {
      fakeME->SetBinContent(i,0.);
      fakeME->SetBinError(i,0.);
    }
    //cout << i << endl;
  } 
  

 TH2F *hframe=NULL,*hframe2=NULL; 
 hframe= new TH2F("hframe","hframe",500,0.,2000.,500,0.001,2.0);
  
  
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

 fake->GetXaxis()->SetTitleOffset(1.7);
 fake->GetYaxis()->SetTitleOffset(1.7);
 fake->GetXaxis()->SetLabelSize(0.03);
 fake->GetYaxis()->SetLabelSize(0.03);
 //fake->GetXaxis()->SetRangeUser(20.0,2000.0);
 fake->GetYaxis()->SetRangeUser(0.001,2.0);
 fake->SetMarkerColor(1);
 fake->SetMarkerStyle(24);
 fake->GetXaxis()->SetTitle("p_{T} (GeV)");
 fake->GetYaxis()->SetTitle("Fake rate (20 GeV bin)");
 fake->GetXaxis()->SetTitleOffset(1.5);
 fake->GetYaxis()->SetTitleOffset(2.0);
 fake->GetXaxis()->SetTitleSize(0.05);
 fake->GetXaxis()->SetLabelSize(0.03);
 fake->GetYaxis()->SetLabelSize(0.03);

 hframe->Draw();
 //gr->Draw("EPsame");
 //fake->Draw("EP1");
 //fake->Draw("same");
 
 tText1->Draw("same");
 tText2->Draw("same");

 c1->SaveAs("FR.png");
 c1->SaveAs("FR.pdf");
 c1->SaveAs("FR.eps");
 c1->SaveAs("FR.root");


 TCanvas *c2 = new TCanvas("c2","c2",800,600);

 c2->SetLogy(1);
 c2->SetLogx(1);
 // c2->SetGrid();

 gPad->SetTopMargin(0.12);
 gPad->SetLeftMargin(0.15);
 gPad->SetFillColor(0);
 gPad->SetTickx();
 gPad->SetTicky();
 // gPad->SetGridy();
 //gPad->SetGridx();


 fakeMB->GetXaxis()->SetTitleOffset(1.7);
 fakeMB->GetYaxis()->SetTitleOffset(1.7);
 fakeMB->GetXaxis()->SetLabelSize(0.03);
 fakeMB->GetYaxis()->SetLabelSize(0.03);
 //fake->GetXaxis()->SetRangeUser(20.0,2000.0);
 fakeMB->GetYaxis()->SetRangeUser(0.001,2.0);
 fakeMB->GetXaxis()->SetTitle("p_{T} (GeV)");
 fakeMB->GetYaxis()->SetTitle("Fake rate (10 GeV)");
 //fakeMB->GetXaxis()->SetTitleOffset(1.5);
 //fakeMB->GetYaxis()->SetTitleOffset(2.0);
 //fakeMB->GetXaxis()->SetTitleSize(0.05);
 //fakeMB->GetXaxis()->SetLabelSize(0.03);
 //fakeMB->GetYaxis()->SetLabelSize(0.03);

 //fakeMB->GetXaxis()->SetRangeUser(20.0,2000.0);
 fakeMB->GetYaxis()->SetRangeUser(0.001,2.0);
 fakeMB->SetMarkerColor(2);
 fakeMB->SetMarkerStyle(24);
 //fakeMB->Draw("EP");
 //fakeMB->Draw("same");
 grMB->Draw("EP");
 grMB->Draw("same");
 leg->AddEntry(grMB,"Muon barrel","p");


 
 //fakeME->GetXaxis()->SetRangeUser(20.0,2000.0);
 fakeME->GetYaxis()->SetRangeUser(0.001,2.0);
 fakeME->SetMarkerStyle(24);
 fakeME->SetMarkerColor(4);
 //fakeME->Draw("EP1");
 //fakeME->Draw("same");
 grME->Draw("EP");
 grME->Draw("same");
 leg->AddEntry(grME,"Muon endcap","p");

 tText1->Draw("same");
 tText2->Draw("same");
 leg->Draw("same");

 c2->SaveAs("plots_25ns/FR_MB_ME.png");
 c2->SaveAs("plots_25ns/FR_MB_ME.pdf");
 c2->SaveAs("plots_25ns/FR_MB_ME.eps");
 c2->SaveAs("plots_25ns/FR_MB_ME.root");
 
}

