{

  TFile *fDiJets = TFile::Open("plots/FR-estimate-DiJets-Data-OS-1GeVBin-2673.root");

  TH1F *hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData= new TH1F("hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData", "Mass of four leptons after fullselecti
on", 2500, 0.,2500. );
  hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData = (TH1F*)fDiJets->Get("DataSub");
  TH1F *hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData_new = new TH1F("hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData_new", "Mass of di-jets by FR from da
ta", 10000, 0.,10000.);
      
  int mbins=1;
  for (int nbins=1;nbins<=hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData->GetNbinsX(); nbins++){
    if (hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData->GetBinCenter(nbins)>0.5 && hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData->GetBinCenter(nbins)<999
	9.5){
      hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData_new->SetBinContent(mbins,double(hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData->GetBinContent(nbins))
										 );
      mbins++;
    }
  }

}
