#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TROOT.h>
#include <string>
#include <iostream>
#include <TSystem.h>
#include <TH2.h>
#include "TChain.h"
#include <stdlib.h>
#include <TDCacheFile.h>
#include "ZprimeMuMu_FR.h"

#endif

using namespace std;

int main(int argc, char ** argv){

  cout << "This is " << argv[1] << endl;
  string sampletype=argv[1];

  ifstream fdata;
  int nlines;

  string dataconf="";
  string mcconf="";
  
  if (sampletype.find("sig") < 10) {
    fdata.open(argv[2]);
    nlines = atoi(argv[3]);
    mcconf="Fall15";
  }
  else if (sampletype.find("bkg") < 10) {
    fdata.open(argv[4]);
    nlines = atoi(argv[5]);
    mcconf="Fall15";
  }
  else if (sampletype.find("data") < 10) {
    fdata.open(argv[6]);
    nlines = atoi(argv[7]);
    dataconf="2015";
  }
  
  string samples[nlines];
  float ninput[nlines];
  float nhlt[nlines];
  float nskim[nlines];
  float xsection[nlines];

  for(int i=0;i<nlines;i++){
    fdata >> samples[i] >> ninput[i] >> nhlt[i] >> nskim[i] >> xsection[i];
    cout << "Sample=" << samples[i] << " Ninput=" << ninput[i] << " NHLT=" << nhlt[i] << " NSkim=" << nskim[i] << " Xsection(pb)=" << xsection[i] << endl;
  }

  //
  float lumifb=0.;

  if (mcconf.find("Fall15")<5) lumifb=2.906; // 2015B+2015C+2015D

  string site=argv[8];
  //string site="Bari";
  cout << "Site is " << site.c_str() << " MC conf.= " << mcconf.c_str() << " data conf.= " << dataconf.c_str() << endl;

  // Run on data
 
  for(int i=0;i<nlines;i++){
    
    string name;
    if (mcconf.find("Fall15")<10) name= "CMSSW763_Data2015_ZprimeMuMu_13TeV_"+samples[i]+".root";
    if (dataconf.find("2015")<10) name=samples[i]+".root";
    
    TString dirInput;
    if (site.find("CERN")<5){
      dirInput="/castor/cern.ch/user/n/ndefilip/Paper/MCFall11";    // to run at CERN
    }
    else if (site.find("DESY")<5){
      dirInput="/nfs/dust/test/cmsdas/school16/ZllExercise/"; //to run at DESY
    }
    else if (site.find("FNAL")<5 && mcconf.find("Spring15_combined")<5){
      dirInput="/eos/uscms/store/user/cmsdas/2016/LONG_EXERCISES/ZprimeDiLeptons/Spring15_25ns_merged";
    }
    else if (site.find("FNAL")<5 && dataconf.find("2015")<5){
      dirInput="/eos/uscms/store/user/cmsdas/2016/LONG_EXERCISES/ZprimeDiLeptons/Data2015_ZprimeMuMu_13TeV_merged";
    }
    else if (mcconf.find("Fall15")<5){
      dirInput="/lustre/cms/store/user/defilip/ZprimeAnalysis/Fall15_merged";
    }
    else if (dataconf.find("2015")<5){
      dirInput="/lustre/cms/store/user/defilip/ZprimeAnalysis/Data2015rereco_ZprimeMuMu_13TeV_merged";
    }
    
    TString File=name;
    
    Char_t namechar[300];
    sprintf(namechar,"%s/%s",dirInput.Data(),File.Data());
    
    float weight= lumifb*(xsection[i]*1000.*nskim[i]/ninput[i])/nskim[i];
    cout << "weight is " << weight << endl;
    
    TFile *file3;
    TTree *tree3;
    
    file3 = TFile::Open(namechar);
    cout << "Read file with name: " << namechar << endl;
    tree3 = (TTree*)file3->Get("tree");
    
    cout << "Read file with name: " << namechar << " " << tree3->GetEntries() << endl;
    ZprimeMuMu_FR b(namechar,tree3,weight,dataconf,mcconf);
    b.Loop();
    
    delete tree3;
    file3 -> Close();
    
  }

  return 0; 

}

