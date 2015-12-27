{
  
  //TFile *_file0 = TFile::Open("Results2673pb/Dibosons-TTbarandTTbarLike-MC-SS-2673pb.root");
  //TFile *_file0 = TFile::Open("Results2673pb/DY-MuMu-MC-OS-allbins-MC-2673pb.root");
  //TFile *_file0 = TFile::Open("Results2673pb/Data-Dibosons-TTbarandTTbarLike-MC-OS-2673pb.root");
  //TFile *_file0 = TFile::Open("Results2673pb/Wjets-25nsMC-OS-allbins-2673pb.root");
  //TFile *_file0 = TFile::Open("Results2673pb/DiJets-Data-OS-2673pb-FR.root");
  TFile *_file0 = TFile::Open("Results2673pb/SSplots/Data-SS-2673pb.root");
  //TH1F *Histo = (TH1F*)_file0->Get("DataSub");
  //TH1F *Histo = (TH1F*)_file0->Get("DiJetsbkgMC");  
  //TH1F *Histo = (TH1F*)_file0->Get("WjetsHisto");
  //TH1F *Histo = (TH1F*)_file0->Get("hMassDYAll6");
  //TH1F *Histo = (TH1F*)_file0->Get("TTbarAndTTbarlikeHisto");
  //TH1F *Histo = (TH1F*)_file0->Get("DiBosonBG");
  //TH1F *Histo = (TH1F*)_file0->Get("ZprimeRecomassSS"); 
  TH1F *Histo = (TH1F*)_file0->Get("DataHisto");
  //TH1F *Histo = (TH1F*)_file0->Get("AllBackgroundsFinal"); 
  //TH1F *Histo = (TH1F*)_file0->Get("hMassDY50GeV"); 
  //TH1F *Histo = (TH1F*)_file0->Get("WjetsHistofinal");    
  //TH1F *Histo = (TH1F*)_file0->Get("DiJetHisto");        

  TAxis *axis = Histo->GetXaxis();
  
  float xmin=120.,xmax=200.;
  //float xmin=200.,xmax=400.;
  //float xmin=400.,xmax=10000.;
  
  int bmin = axis->FindBin(xmin); 
  int bmax = axis->FindBin(xmax); 
  double integral = Histo->Integral(bmin,bmax);
  cout << "Integral before= " << integral << endl;
  integral-=Histo->GetBinContent(bmin)*(xmin-(axis->GetBinLowEdge(bmin)))/(axis->GetBinWidth(bmin));
  integral-=Histo->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-xmax)/(axis->GetBinWidth(bmax));
  cout << "Integral after= " << integral << endl;
  
  double entry=0;
  for (int i=bmin;i<bmax;i++){
    entry+=Histo->GetBinContent(i)*Histo->GetEntries()/Histo->Integral();
    // cout << "Entry= " << entry << endl;
  }
  double error=sqrt(entry)*Histo->Integral()/Histo->GetEntries();
  cout << "Entry final= "<< entry << " " << integral << " weighted statistical error= " << error << " weighted systematic error (20%)= " << integral*0.25<< endl;
}
