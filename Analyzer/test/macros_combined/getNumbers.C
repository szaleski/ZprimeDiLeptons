{

  //TFile *_file0 = TFile::Open("plots/h_ZprimeRecomass_data_OF.root");
  TFile *_file0 = TFile::Open("plots/h_ZprimeRecomass_data.root");
  //TFile *_file0 = TFile::Open("plots/h_ZprimeRecomass_Tlike.root");      
  //TFile *_file0 = TFile::Open("plots/h_ZprimeRecomass_DiBoson.root");  

 //TH1F *Histo = (TH1F*)_file0->Get("htotaldata");
 //TH1F *Histo = (TH1F*)_file0->Get("htotalHisto");
 
 TH1F *Histo = (TH1F*)_file0->Get("ZprimeRecomass");
 TAxis *axis = Histo->GetXaxis();

 //float xmin=120.,xmax=200.;
 //float xmin=200.,xmax=10000.;
 //float xmin=400.,xmax=10000.;
 float xmin=60.,xmax=120.;

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
 cout << "Entry final= "<< entry << " " << integral << " weighted statistical error= " << error << " weighted systematic error (20%)= " << integral*0.2<< endl;
}
