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
#include "TH1.h"

using namespace std;

class getBinning{
  
public: 
  getBinning();
  
};

getBinning::getBinning(){
  // const int NMBINS = 100;
  // const double MMIN = 60., MMAX = 3000.;
  // double logMbins[NMBINS+1];
  // for (int ibin = 0; ibin <= NMBINS; ibin++) {
  //  logMbins[ibin] = exp(log(MMIN) + (log(MMAX)-log(MMIN))*ibin/NMBINS);
  //  cout << logMbins[ibin] << endl;
  // }
  
  TFile* inFile = new TFile("plots/htotal_root_ZprimeRecomass_OF.root","READ");
  TFile* outFile = new TFile("plots/htotal_root_ZprimeRecomass_binned.root","RECREATE");
  
  std::vector<TH1*> hists;
  
  hists.push_back((TH1*) inFile->Get("htotaldatawtOverFlow"));
  hists.push_back((TH1*) inFile->Get("hfourlepbestmass_4l_afterSel_new_DYwtOverFlow"));
  hists.push_back((TH1*) inFile->Get("hfourlepbestmass_4l_afterSel_new_diBosonwtOverFlow"));
  hists.push_back((TH1*) inFile->Get("hfourlepbestmass_4l_afterSel_new_TlikewtOverFlow"));
  hists.push_back((TH1*) inFile->Get("hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromDatawtOverFlow"));
  hists.push_back((TH1*) inFile->Get("hfourlepbestmass_4l_afterSel_WJetsFromFakeRateFromMCwtOverFlow"));
  
  const int NMBINS = hists[0]->GetNbinsX();
  const double MMIN = 60., MMAX = 3001.;
  double logMbins[NMBINS+1];
  float binNormNr=0.;
  
  for(size_t histNr=0;histNr<hists.size();histNr++){    
    for(int binNr=1;binNr<=hists[histNr]->GetNbinsX();binNr++){
      logMbins[binNr] = exp(log(MMIN) + (log(MMAX)-log(MMIN))*binNr/NMBINS);
      cout << logMbins[binNr] << endl;
      if (binNr>1) binNormNr=logMbins[binNr]-logMbins[binNr-1];
      else binNormNr=logMbins[binNr];
      hists[histNr]->SetBinContent(binNr,hists[histNr]->GetBinContent(binNr)*binNormNr/hists[histNr]->GetBinWidth(binNr));
      hists[histNr]->SetBinError(binNr,hists[histNr]->GetBinError(binNr)*binNormNr/hists[histNr]->GetBinWidth(binNr));
    }
    hists[histNr]->SetDirectory(outFile);    
  }

  outFile->Write();
  
  delete inFile;
  delete outFile;
  
  
}
