{

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
  TText *t2 = tText2->AddText("CMS Spring15 Simulation");
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


  TFile *_file0 = TFile::Open("htotal_noQCD_h1_Den_Pt_w.root");
  TH1F *fraction_noQCD_Den; _file0->GetObject("fraction_noQCD",fraction_noQCD_Den);


  TFile *_file1 = TFile::Open("htotal_noQCD_h1_Num_Pt_w.root");
  TH1F *fraction_noQCD_Num; _file1->GetObject("fraction_noQCD",fraction_noQCD_Num);

  TH2F *hframe=NULL; 
  hframe= new TH2F("hframe","hframe",500,50.,1000.,500,0.,1.3);
  hframe->GetYaxis()->SetTitle("EWK fraction / 10 GeV");
  hframe->GetXaxis()->SetTitle("p_{T} (GeV)");
  hframe->GetXaxis()->SetTitleOffset(1.7);
  hframe->GetYaxis()->SetTitleOffset(1.7);

  TCanvas *c1 = new TCanvas("c1","c1",800,600);

  c1->SetLogy(0);
  c1->SetLogx(0);
  //c1->SetGrid();
  
  gPad->SetTopMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->SetFillColor(0);
  gPad->SetTickx();
  gPad->SetTicky();
  // gPad->SetGridy();
  //gPad->SetGridx();

 
  hframe->GetXaxis()->SetTitleOffset(1.7);
  hframe->GetYaxis()->SetTitleOffset(1.7);
  hframe->GetXaxis()->SetLabelSize(0.03);
  hframe->GetYaxis()->SetLabelSize(0.03);
  hframe->GetXaxis()->SetTitle("p_{T} (GeV)");
  hframe->GetYaxis()->SetTitle("EWK Fraction");
  hframe->GetXaxis()->SetTitleOffset(1.5);
  hframe->GetYaxis()->SetTitleOffset(2.0);
  hframe->GetXaxis()->SetTitleSize(0.05);
  hframe->GetXaxis()->SetLabelSize(0.03);
  hframe->GetYaxis()->SetLabelSize(0.03);

  fraction_noQCD_Den->SetMarkerColor(2);
  fraction_noQCD_Den->SetMarkerStyle(20);
  fraction_noQCD_Den->SetMarkerSize(0.95);
  leg->AddEntry(fraction_noQCD_Den,"denom.","p");
   
  fraction_noQCD_Num->SetMarkerColor(4);
  fraction_noQCD_Num->SetMarkerStyle(21);
  fraction_noQCD_Num->SetMarkerSize(0.95);
  leg->AddEntry(fraction_noQCD_Num,"num.","p");
  
  hframe->Draw();
  fraction_noQCD_Den->Draw("Psame");
  fraction_noQCD_Num->Draw("Psame");
  tText1->Draw("same");
  tText2->Draw("same");
  leg->Draw("same");

  c1->SaveAs("fraction_noQCD.png");
  c1->SaveAs("fraction_noQCD.pdf");
  c1->SaveAs("fraction_noQCD.eps");
  
}
