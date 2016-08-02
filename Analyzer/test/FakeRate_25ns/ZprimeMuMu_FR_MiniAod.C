#define ZprimeMuMu_FR_MiniAod_cxx
#include "ZprimeMuMu_FR_MiniAod.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TFile.h"
#include <math.h>
#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include <time.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
    
using namespace std;

const double PI = 3.141592654; 


void ZprimeMuMu_FR_MiniAod::Loop()
{

  TString inputfile=name; 
  cout << "Name of the input file is= " << inputfile.Data() << endl;

  if (fChain == 0) return;

  if( DATA_type=="2016") weight=1;

  time_t start,end;
  double dif;
  time (&start);

  TFile *file = new TFile("CMSSW803-Analyse_ZprimeToMuMu_13TeV.root", "RECREATE");
  file->cd();

  TH1F *histo_EVENT              = new TH1F("histo_EVENT", "histo_EVENT",600,0,600.);
  
  TH1D* h1_Den_Nloose            = new TH1D("h1_Den_Nloose","h1_Den_Nloose",5,-0.5,4.5);
  TH1D* h1_Num_Ngood             = new TH1D("h1_Num_Ngood","h1_Num_Ngood",5,-0.5,4.5);
  TH1F* h1_Den_Pt                = new TH1F("h1_Den_Pt","h1_Den_Pt",1000,0.,2000.);
  TH1F* h1_Num_Pt                = new TH1F("h1_Num_Pt","h1_Num_Pt",1000,0.,2000.);
  TH1F* h1_Den_Pt_w              = new TH1F("h1_Den_Pt_w","h1_Den_Pt_w",1000,0.,2000.);
  TH1F* h1_Num_Pt_w              = new TH1F("h1_Num_Pt_w","h1_Num_Pt_w",1000,0.,2000.);

  TH1F* h1_Den_Pt_Barrel         = new TH1F("h1_Den_Pt_Barrel","h1_Den_Pt_Barrel",1000,0.,2000.);
  TH1F* h1_Num_Pt_Barrel         = new TH1F("h1_Num_Pt_Barrel","h1_Num_Pt_Barrel",1000,0.,2000.);
  TH1F* h1_Den_Pt_EndCap         = new TH1F("h1_Den_Pt_EndCap","h1_Den_Pt_EndCap",1000,0.,2000.);
  TH1F* h1_Num_Pt_EndCap         = new TH1F("h1_Num_Pt_EndCap","h1_Num_Pt_EndCap",1000,0.,2000.);
  
  TH1F* h1_Den_Pt_Barrel_w         = new TH1F("h1_Den_Pt_Barrel_w","h1_Den_Pt_Barrel_w",1000,0.,2000.);
  TH1F* h1_Num_Pt_Barrel_w         = new TH1F("h1_Num_Pt_Barrel_w","h1_Num_Pt_Barrel_w",1000,0.,2000.);
  TH1F* h1_Den_Pt_EndCap_w         = new TH1F("h1_Den_Pt_EndCap_w","h1_Den_Pt_EndCap_w",1000,0.,2000.);
  TH1F* h1_Num_Pt_EndCap_w         = new TH1F("h1_Num_Pt_EndCap_w","h1_Num_Pt_EndCap_w",1000,0.,2000.);
  

  TH1D* h1_NJets                 = new TH1D("h1_NJets","h1_NJets",5,-0.5,4.5);
  TH1F* h1_mW_T                  = new TH1F("h1_mW_T","h1_mW_T",1000,0.,1000.);
  TH1F* h1_MET                   = new TH1F("h1_MET","h1_MET",1000,0.,1000.);
  TH1F* h1_METSign               = new TH1F("h1_METSign","h1_METSign",1000,0.,10.);
  TH1F* h1_METSign_LH            = new TH1F("h1_METSign_LH","h1_METSign_LH",1000,0.,1000.);
  
  TH2F* h2_Den_Pt_Eta            = new TH2F("h2_Den_Pt_Eta","h2_Den_Pt_Eta",500,0.,500.,25,0.,2.5);
  TH2F* h2_Num_Pt_Eta            = new TH2F("h2_Num_Pt_Eta","h2_Num_Pt_Eta",500,0.,500.,25,0.,2.5);    

  int event=0;
  deltaRcut=0.2;
  RecoHLTMatchingDeltaRcut = 0.20;

  double muon_mass = 0.1056583;
  float gMassInv=-999,gMaxMassInv=-999.;
  int igfirstmu=-10,igsecondmu=-10;
  float MassInv=-999,MaxMassInv=-999.;
  int ifirstmu=-10,isecondmu=-10;
  double Zmass = 91.188;
  float mW_T = -999;
  // Book txt file for candidate events
  Char_t txtOUT[500];
  //sprintf(txtOUT,"%s_txt.txt",datasetName.Data());
  sprintf(txtOUT,"CMSSW803-Analyse_ZprimeToMuMu_13TeV_cand.txt");
  cout << "Opening a txt file with candidate events " << txtOUT << endl;
  //ofstream output_txt;
  output_txt.open(txtOUT);
  output_txt << "CANDIDATES Events:" << endl;

  Char_t outform[20000];
  sprintf (outform,"run: lumi: event: pTmu1: etamu1:");
  output_txt  << outform << endl;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;   
    //if( inputfile.Contains("amcatnlo")) {      
    //  cout << "Reweighting sample of amcatnlo with weight= " << MC_weighting->at(0)<< endl;
    //  weight=weight*MC_weighting->at(0);
    //}
    //if (!(event_runNo==254227 && event_evtNo==123724167 && event_lumi==124)) continue;
    cout<<"Analyzing entry: "<<jentry<<" Run: "<<event_runNo<<" Event: "<<event_evtNo<<" LumiSection: "<<event_lumi<<endl;
    event++;
    histo_EVENT->Fill(event);
    bool firstGenMu  = SelectFirstGenMu(genET1,genPhi1,genEta1,genEn1,genID1,genStat1,flag1);
    bool secondGenMu = SelectSecondGenMu(flag1,genET1,genET2,genPhi2,genEta2,genEn2,genID2,genStat2);

    bool firstMuFinal  = SelectFirstMuon(PtRecTunePMuBestTrack1,EnRecMu1,EtaRecMu1,PhiRecMu1,ChargeRecMu1,flagmu1,
					 pxRecMu1,pyRecMu1,pzRecMu1,pRecMu1,dxyRecMu1,PtRecTunePMu1,
					 PtRecMuBestTrack1);     
    bool secondMuFinal = SelectSecondMuon(ChargeRecMu1,flagmu1,PtRecTunePMuBestTrack1,EtaRecMu1,PhiRecMu1,
					  PtRecTunePMuBestTrack2,EnRecMu2,
					  EtaRecMu2,PhiRecMu2,ChargeRecMu2,pxRecMu2,pyRecMu2,pzRecMu2,pRecMu2,dxyRecMu2,
					  PtRecTunePMu2,PtRecMuBestTrack2);
    //if (!(firstMuFinal && secondMuFinal)) continue; // DY
    //if (!(firstMuFinal)) continue; // W+jets
    //if (secondMuFinal) continue;

    // Select the gen muons with the highest invariant mass
    gMassInv=999,gMaxMassInv=999.;
    igfirstmu=-10,igsecondmu=-10;

    //if (iGen->size()<2) continue;
    for(unsigned i=0; i<iGen->size(); i++){
      if( fabs(idGen->at(i)) != 13 ) continue;
      //if( fabs(etaGen->at(i))>2.4 ) continue;
      for(unsigned j=0; j<iGen->size(); j++){
	if( fabs(idGen->at(j)) != 13 ) continue;
	//if( fabs(etaGen->at(j))>2.4 ) continue;
	if (chargeGen->at(i)*chargeGen->at(j)>0.) continue;
	gMassInv=-999.;
	gMassInv = Mass(ptGen->at(i),etaGen->at(i),phiGen->at(i),EnergyGen->at(i),
			ptGen->at(j),etaGen->at(j),phiGen->at(j),EnergyGen->at(j));
	gMassInv=fabs(gMassInv-Zmass);
	if (gMassInv<gMaxMassInv && (i==flag1 || j==flag1)) {
	  igfirstmu=i;
	  igsecondmu=j;
	  gMaxMassInv=gMassInv;
	  cout << "Gen leptons for Max Inv Mass /pT/eta)= " 
	       << ptGen->at(i) << " " << etaGen->at(i) << " " 
	       << ptGen->at(j) << " " << etaGen->at(j) << " "
	       << " with Inv mass -mZ= " << gMassInv << endl;
	}
      }
    }
    cout << "Generated MaxMassInv - mZ= " << gMaxMassInv << endl;

    //    if (!(igfirstmu>=0 && igsecondmu>=0)) continue;    
    //if (!(GenToRecoMatchMu(ptGen->at(igfirstmu),etaGen->at(igfirstmu), phiGen->at(igfirstmu))>=0 && 
    //	  GenToRecoMatchMu(ptGen->at(igsecondmu),etaGen->at(igsecondmu), phiGen->at(igsecondmu))>=0)) continue; // DY 


    // Select the reco muons with the highest invariant mass
    MassInv=-999,MaxMassInv=-999.;
    ifirstmu=-10,isecondmu=-10;

    for(unsigned i=0; i<Mu_nbMuon->size(); i++){
      //if( fabs(Mu_ptTunePMuonBestTrack->at(i))>2.4 ) continue;
      if( Mu_isGlobalMuon->at(i) == 1 &&
	  //Mu_ptTunePMuonBestTrack->at(i) > 45.0 &&
	  Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
	  (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
	  Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
	  Mu_numberOfValidPixelHits->at(i) > 0 &&
	  Mu_numberOfValidMuonHits->at(i) > 0 &&
	  Mu_numberOfMatchedStations->at(i) > 1 &&
	  Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ) { 
	
	for(unsigned j=0; j!=i && j<Mu_nbMuon->size(); j++){
	  //if( fabs(Mu_ptTunePMuonBestTrack->at(j))>2.4 ) continue;
	  if( Mu_isGlobalMuon->at(j) == 1 &&
	      //Mu_ptTunePMuonBestTrack->at(j) > 45.0 &&
	      Mu_absdxyTunePMuonBestTrack->at(j) < 0.2 &&
	      (Mu_trackiso->at(j)/Mu_ptInnerTrack->at(j)) < 0.10  &&
	      Mu_numberOftrackerLayersWithMeasurement->at(j) > 5 && 
	      Mu_numberOfValidPixelHits->at(j) > 0 &&
	      Mu_numberOfValidMuonHits->at(j) > 0 &&
	      Mu_numberOfMatchedStations->at(j) > 1 &&
	      Mu_dPToverPTTunePMuonBestTrack->at(j) < 0.3 ) { 
	    cout << "Charge= " << Mu_chargeTunePMuonBestTrack->at(i)<< " " << Mu_chargeTunePMuonBestTrack->at(j) << endl;
	    if (Mu_chargeTunePMuonBestTrack->at(i)*Mu_chargeTunePMuonBestTrack->at(j)>0.) continue;
	    float enMu1=sqrt((Mu_pTunePMuonBestTrack->at(i)*Mu_pTunePMuonBestTrack->at(i))-muon_mass*muon_mass);
	    float enMu2=sqrt((Mu_pTunePMuonBestTrack->at(j)*Mu_pTunePMuonBestTrack->at(j))-muon_mass*muon_mass);
	    MassInv = Mass(Mu_ptTunePMuonBestTrack->at(i),Mu_etaTunePMuonBestTrack->at(i),Mu_phiTunePMuonBestTrack->at(i),enMu1,
			   Mu_ptTunePMuonBestTrack->at(j),Mu_etaTunePMuonBestTrack->at(j),Mu_phiTunePMuonBestTrack->at(j),enMu2);
	    if (MassInv>MaxMassInv) {
	      ifirstmu=i;
	      isecondmu=j;
	      MaxMassInv=MassInv;
	    }
	  }
	}
      }
    }
    cout << "MaxMassInv= " << MaxMassInv << endl;
    //if (MaxMassInv<0) {
    //  gMaxMassInv
    //  continue;
    //}

    if (firstMuFinal && secondMuFinal){
      MassInv = Mass(PtRecTunePMuBestTrack1,EtaRecMu1,PhiRecMu1,EnRecMu1,
		     PtRecTunePMuBestTrack2,EtaRecMu2,PhiRecMu2,EnRecMu2);
      //if (MassInv>0.) continue; 
      //if (MassInv>=60. || MassInv<=120.) continue; // Z constraint
    }
      //if (gMassInv<60. || gMassInv>120.) continue; // Z constraint
    //cout << "MassInv= "<< MassInv << endl;

    bool secondMuSSFinal = SelectSecondSSMuon(ChargeRecMu1,flagmu1,PtRecTunePMuBestTrack1,PtRecTunePMuBestTrack3,EnRecMu3,
					  EtaRecMu3,PhiRecMu3,ChargeRecMu3,pxRecMu3,pyRecMu3,pzRecMu3,pRecMu3,dxyRecMu3,
					  PtRecTunePMu3,PtRecMuBestTrack3);

    if (firstMuFinal && secondMuSSFinal){
      MassInv = Mass(PtRecTunePMuBestTrack1,EtaRecMu1,PhiRecMu1,EnRecMu1,
		     PtRecTunePMuBestTrack3,EtaRecMu3,PhiRecMu3,EnRecMu3);
      //if (MassInv>=0.) continue; 
      //if (MassInv>=12) continue; 
      //if (MassInv>=60. || MassInv<=120.) continue; // Z constraint
   }
 

    //bool GenRecoMatch1 = GenRecoMatchMu1(EtaRecMu1,PhiRecMu1,mPtGen1,mPhiGen1,mEtaGen1,mEnGen1,mGenFlag1);
    //bool GenRecoMatch2 = GenRecoMatchMu2(mGenFlag1,EtaRecMu2,PhiRecMu2,mPtGen2,mPhiGen2,mEtaGen2,mEnGen2);
  
    //if (!(firstGenMu)) continue;
    //if (!(GenToRecoMatchMu(genET1,genEta1, genPhi1)>=0 && GenToRecoMatchMu(genET2,genEta2, genPhi2)>=0)) continue; // DY 
    //if (!(GenToRecoMatchMu(genET1,genEta1, genPhi1)>=0)) continue; // W+jets
   
    
    //cout << "GenToreco1 " << GenToRecoMatchMu(genET1,genEta1, genPhi1) 
    //	 << "  GenToReco2= "<< GenToRecoMatchMu(genET2,genEta2, genPhi2)<< endl;

    //cout << " Run=  " <<event_runNo
    //	 << " Lumi= " << event_lumi
    //	 << " Event= " << event_evtNo << endl;
      
    //if (!firstMuFinal || !secondMuFinal) continue;
    //if (!firstMuFinal) continue;
    //cout << "First selected muon= " << PtRecTunePMuBestTrack1 << " Second selected muon= " << PtRecTunePMuBestTrack2 << endl;
    //if (!(firstGenMu && secondGenMu && GenRecoMatch1 && GenRecoMatch2)) continue;
    
    //int N_loose = 0;    
    //int N_good = 0;
    
    int iLl[8]= { -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1};
    int iL[8]= { -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1};
    

    // MET cut to reduce the W+jets
    //if (PFMet_pt >20.) continue;   
    //if (METSign>10.) continue;
    
    h1_MET->Fill(PFMet_pt_cor,weight);
    if (PFMet_sumEt_cor>0.) h1_METSign->Fill(PFMet_pt_cor/PFMet_sumEt_cor,weight);
    h1_METSign_LH->Fill(METSign,weight);
    
    // Number of jets
    // int N_jets=0;
    // for (int j=0;j<jetpt->size();j++){
    //   if (jetpt->at(j)>30. && fabs(jeteta->at(j)<2.4)) N_jets++;
    // }
    // h1_NJets->Fill(N_jets,weight);
    //if (N_jets<2) continue;
    
    for(int i=0;i<Mu_nbMuon->size();i++){

      int N_loose=  0;    
      int N_good = 0;
  

      //      if (i>2) continue;
      cout << "Muon with pT= "<< Mu_ptTunePMuonBestTrack->at(i) << " and eta=" <<  fabs(Mu_etaTunePMuonBestTrack->at(i)) << endl;

      //if(i==GenToRecoMatchMu(genET1,genEta1, genPhi1) || 
      //	 i==GenToRecoMatchMu(genET2,genEta2, genPhi2)) continue; // DY
      //if(i==GenToRecoMatchMu(genET1,genEta1, genPhi1)) continue; // W+jets
      //if (MaxMassInv>=0){
      //	if (i==ifirstmu || i==isecondmu) {
      //  cout << "Discarding muon with pT=" << Mu_ptTunePMuonBestTrack->at(i) << endl;
      //  continue;	  
      //}
      //}
      //else {
      //	cout << "Looking for muons from generated Invmass" << endl;
      //	if (gMaxMassInv>=0){
      //  if(i==GenToRecoMatchMu(ptGen->at(igfirstmu),etaGen->at(igfirstmu), phiGen->at(igfirstmu)) || 
      //     i==GenToRecoMatchMu(ptGen->at(igsecondmu),etaGen->at(igsecondmu), phiGen->at(igsecondmu))) continue; // 
      //}
      //else break;
      //}
      
      //if (Mu_ptTunePMuonBestTrack->at(i)==PtRecTunePMuBestTrack1 || Mu_ptTunePMuonBestTrack->at(i)==PtRecTunePMuBestTrack2) continue;

      if( fabs(Mu_etaTunePMuonBestTrack->at(i)) >= 2.4 ) continue;

      // NO HLT martching
      // bool RecoMuon1MatchingWithHLT = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(i),Mu_phiTunePMuonBestTrack->at(i));
      // if(RecoMuon1MatchingWithHLT==false) continue;


      //if(i==GenToRecoMatchMu(ptGen->at(igfirstmu),etaGen->at(igfirstmu), phiGen->at(igfirstmu)) || 
      // i==GenToRecoMatchMu(ptGen->at(igsecondmu),etaGen->at(igsecondmu), phiGen->at(igsecondmu))) continue; // 

      //if(Mu_ptTunePMuonBestTrack->at(i) == PtRecTunePMuBestTrack1) continue;
      //if(Mu_ptTunePMuonBestTrack->at(i) == PtRecTunePMuBestTrack2) continue;  
      //if(Mu_ptTunePMuonBestTrack->at(i) == PtRecTunePMuBestTrack3) continue;

      //if (Mu_ptTunePMuonBestTrack->at(i)<100) continue;
      //if (Mu_ptcocktail->at(i)<20) continue;

      cout << "Muon with pT= "<< Mu_ptTunePMuonBestTrack->at(i) << endl;
      

      cout << "Pfmet= " << PFMet_pt_cor << endl;

      if(Mu_ptTunePMuonBestTrack->at(i) == PtRecTunePMuBestTrack1) {
	float dp=std::abs(PFMet_phi_cor-Mu_phiTunePMuonBestTrack->at(i));
	if (dp>PI) dp-=float(2*PI);
	float mW_T=sqrt(2*Mu_ptTunePMuonBestTrack->at(i)*PFMet_pt_cor*(1-cos(dp)));
	h1_mW_T->Fill(mW_T,weight);
	if (mW_T<35 )  continue; // against W+jets
      }
      
      if(
	 Mu_isGlobalMuon->at(i) == 1  && Mu_isTrackerMuon->at(i) == 1 && 
	 // fabs(Mu_etaCocktail->at(i))<2.5 &&
	 //Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
	 //Mu_numberOfValidPixelHits->at(i) > 0 
	 //Mu_numberOfValidMuonHits->at(i) > 0
	 // Mu_numberOfMatchedStations->at(i) > 1 &&
	 // Mu_dPToverPTcocktail->at(i) < 0.3 
	 //){
	 // Mu_isGlobalMuon->at(i) == 1 &&	   
	 Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
	 Mu_absdzTunePMuonBestTrack->at(i) < 1. && 
	 Mu_numberOftrackerLayersWithMeasurement->at(i) > 5.  &&
	 Mu_numberOfValidPixelHits->at(i) > 0. 
	 
	 )
	{	
	  iLl[ N_loose ] = i ;
	  ++N_loose ;
	  
	  //if (fabs(Mu_etaTunePMuonBestTrack->at(i))<1.2){
	  //if (fabs(Mu_etaCocktail->at(i))<1.2){
	  h1_Den_Nloose->Fill(N_loose);
	  cout << "Filling Denominator with pT= " << Mu_ptTunePMuonBestTrack->at(i) << endl;
	  h1_Den_Pt->Fill(Mu_ptTunePMuonBestTrack->at(i));
	  h1_Den_Pt_w->Fill(Mu_ptTunePMuonBestTrack->at(i),weight);
	  //h1_Den_Pt->Fill(Mu_ptcocktail->at(i),1.);
	  h2_Den_Pt_Eta->Fill(Mu_ptTunePMuonBestTrack->at(i),fabs(Mu_etaTunePMuonBestTrack->at(i)));
	  if (fabs(Mu_etaTunePMuonBestTrack->at(i))<=1.2) {
	    h1_Den_Pt_Barrel->Fill(Mu_ptTunePMuonBestTrack->at(i));
	    h1_Den_Pt_Barrel_w->Fill(Mu_ptTunePMuonBestTrack->at(i),weight);
	  }
	  if (fabs(Mu_etaTunePMuonBestTrack->at(i))>1.2) {
	    h1_Den_Pt_EndCap->Fill(Mu_ptTunePMuonBestTrack->at(i));
	    h1_Den_Pt_EndCap_w->Fill(Mu_ptTunePMuonBestTrack->at(i),weight);
	  }
	  
	  
	  // Passing ID and Isolation	
	  if(
	     Mu_isGlobalMuon->at(i) == 1 &&
	     Mu_isTrackerMuon->at(i) == 1 &&
	     Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
	     Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
	     Mu_numberOfValidPixelHits->at(i) > 0 &&
	     Mu_numberOfValidMuonHits->at(i) > 0 &&
	     Mu_numberOfMatchedStations->at(i) > 1 &&
	     Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 &&
	     Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i) < 0.10){
	    
	    // Passing isolation
	    //if((Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10){
	    //if((Mu_isPF->at(i)>0 && Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.6){ // 4L definition
	    iL[ N_good ] = i ;
	    ++N_good ;
	    //if (fabs(Mu_etaTunePMuonBestTrack->at(i))<1.2){
	    //if (fabs(Mu_etaCocktail->at(i))<1.2){
	    h1_Num_Ngood->Fill(N_good);
	    cout << "Filling Numerator with pT= " << Mu_ptTunePMuonBestTrack->at(i) << endl;
	    h1_Num_Pt->Fill(Mu_ptTunePMuonBestTrack->at(i));
	    h1_Num_Pt_w->Fill(Mu_ptTunePMuonBestTrack->at(i),weight);
	    //h1_Num_Pt->Fill(Mu_ptcocktail->at(i),1.);
	    h2_Num_Pt_Eta->Fill(Mu_ptTunePMuonBestTrack->at(i),fabs(Mu_etaTunePMuonBestTrack->at(i)));	  
	    if (fabs(Mu_etaTunePMuonBestTrack->at(i))<=1.2) 
	      {
		h1_Num_Pt_Barrel->Fill(Mu_ptTunePMuonBestTrack->at(i));
		h1_Num_Pt_Barrel_w->Fill(Mu_ptTunePMuonBestTrack->at(i),weight);
	      }
	    if (fabs(Mu_etaTunePMuonBestTrack->at(i))>1.2) {
	      h1_Num_Pt_EndCap->Fill(Mu_ptTunePMuonBestTrack->at(i));
	      h1_Num_Pt_EndCap_w->Fill(Mu_ptTunePMuonBestTrack->at(i),weight); 
	    }
	    
	    Char_t outformat[20000];
	    sprintf (outformat,"%d %d %d %4.2f %4.2f",
		     event_runNo,event_lumi,event_evtNo,Mu_ptTunePMuonBestTrack->at(i),Mu_etaTunePMuonBestTrack->at(i));
	    
	    output_txt  << outformat << endl;
	    
	    
	  }
	  else
	    cout << "Not passing the selection " << endl;    
	  //if (fabs(Mu_etaTunePMuonBestTrack->at(i))<2.1) 
	  //h1_Num_Pt->Fill(Mu_ptTunePMuonBestTrack->at(i));
	  //if (fabs(Mu_etaCocktail->at(i))<2.1) h1_Num_Pt->Fill(Mu_ptcocktail->at(i),0.);
	}
    }
  }
  
  file->Write();
  file->Close();
  output_txt.close();
  //========================================================================
  time (&end);
  dif = difftime (end,start);
  printf ("It took you %.2lf minutes to run your program.\n", (dif/60.0) );
}

//================================================================================== 
//                                                                                 =
//                                    Start methods                                =             
//                                                                                 =
//==================================================================================
//===================== Methode to calculate the mass ========================
float ZprimeMuMu_FR_MiniAod::Mass(float Pt1,float Eta1,float Phi1,float En1,
		       float Pt2,float Eta2,float Phi2,float En2){
  float MuMuMass = 0.0;
  TLorentzVector Mu1;
  TLorentzVector Mu2;
  Mu1.SetPtEtaPhiE(Pt1,Eta1,Phi1,En1);
  Mu2.SetPtEtaPhiE(Pt2,Eta2,Phi2,En2);
  MuMuMass = (Mu1 + Mu2).M();
  return MuMuMass;
}
void ZprimeMuMu_FR_MiniAod::PrintEventInformation(int runNumber, int lumiNumber, int eventNumber,
				       float vtxChi2, float vtxMass, float CosmicRejection)
{
  if(event_runNo == runNumber && event_lumi == lumiNumber && event_evtNo == eventNumber)
    {
      output_txt << event_runNo
		 << "        " << event_lumi
		 << "        " << event_evtNo
		 << "        " << vtxChi2
		 << "        " << vtxMass << endl;
      for(unsigned i=0; i<Mu_nbMuon->size(); i++){
	//if( fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.5 ){cout<<"[1] eta="<<Mu_etaTunePMuonBestTrack->at(i)<<endl;}
	cout<<"[0] phi="<<Mu_phiTunePMuonBestTrack->at(i)<<endl;
	cout<<"[1] eta="<<Mu_etaTunePMuonBestTrack->at(i)<<endl;
	if(Mu_isGlobalMuon->at(i) == 1) {cout<<"[2] isGlobal="<<Mu_isGlobalMuon->at(i)<<endl;}
	if(Mu_ptTunePMuonBestTrack->at(i) > 45.0) {cout<<"[3] ptcocktail="<<Mu_ptTunePMuonBestTrack->at(i)<<endl;}
	if(Mu_absdxyTunePMuonBestTrack->at(i) < 0.2) {cout<<"[4] absdxy="<<Mu_absdxyTunePMuonBestTrack->at(i)<<endl;}
	if(Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i) < 0.10) {cout<<"[5] trackiso="<<Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)<<endl;}
	if(Mu_numberOftrackerLayersWithMeasurement->at(i) > 5) {cout<<"[6] nbTrackerLayer="<<Mu_numberOftrackerLayersWithMeasurement->at(i)<<endl;}
	if(Mu_numberOfValidPixelHits->at(i) > 0) {cout<<"[7] nbPixelHits="<<Mu_numberOfValidPixelHits->at(i)<<endl;}
	if(Mu_numberOfValidMuonHits->at(i) > 0) {cout<<"[8] nbMuonHits="<<Mu_numberOfValidMuonHits->at(i)<<endl;}
	if(Mu_numberOfMatchedStations->at(i) > 1) {cout<<"[9] nbStation="<<Mu_numberOfMatchedStations->at(i)<<endl;}
	if(Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3) {cout<<"[10] DeltaPterror="<<Mu_dPToverPTTunePMuonBestTrack->at(i)<<endl;}
	cout<<"[11] Charge="<<Mu_chargeTunePMuonBestTrack->at(i)<<endl;
      }
      cout<<"[000] vtxMassMu="<<vtxMass<<endl;
      cout<<"[000] vtxChi2Mu="<<vtxChi2<<endl;
      cout<<"[000] CosAngle="<<CosmicRejection<<endl;
    }
}

void ZprimeMuMu_FR_MiniAod::PickThehighestMass(float &vtxHighestMass,float &vtxHighestChi2,int EvtNb)
{

  float Massinv  = -10.0;
  unsigned iflag = -10;
  vtxHighestMass = -10.0;
  vtxHighestChi2 = 1000.0;
  int NbMu = 0;
  int Nb = 0;
  for(unsigned i=0; i<Mu_vtxMass->size(); i++)
    {
      Nb++;
      if(Mu_vtxNormChi2->at(i)> 10) continue;
      if(Mu_vtxMass->at(i)>Massinv){
	Massinv = Mu_vtxMass->at(i);
	iflag  = i;
	NbMu++;
      }
    }
  if( NbMu > 0 ){
    vtxHighestMass = Mu_vtxMass->at(iflag);
    vtxHighestChi2 = Mu_vtxNormChi2->at(iflag);
  }
}
double ZprimeMuMu_FR_MiniAod::ThreeDangle(float pxMu1,float pyMu1,float pzMu1,float pMu1,
			       float pxMu2,float pyMu2,float pzMu2,float pMu2)
{
  TVector3 Mu1(pxMu1,pyMu1,pzMu1);
  TVector3 Mu2(pxMu2,pyMu2,pzMu2);
  double cos_angle = Mu1.Dot(Mu2) / pMu1 / pMu2;
  return cos_angle;
}
/*
double ZprimeMuMu_FR_MiniAod::ThreeDangle(float pxMu1,float pyMu1,float pzMu1,float pMu1,
			       float pxMu2,float pyMu2,float pzMu2,float pMu2)
{
  double angle = acos((pxMu1*pxMu2 + pyMu1*pyMu2 + pzMu1*pzMu2)/(pMu1*pMu2));
  double alpha = PI - angle;
  return alpha;
}
*/
bool ZprimeMuMu_FR_MiniAod::SelectFirstGenMu(float &ETMu1,float &PhiSCMu1,
				  float &EtaSCMu1,float &EnMu1,
				  int &IDele1,int &Statele1,
				  unsigned &GenFlag1){
  int NbHEEPele = 0;
  int iflag = -10;
  ETMu1 = 0.0;
  for(unsigned i=0; i<iGen->size(); i++){
    //printf ("BS[%d] ID = %d status = %d \n",i,idGen->at(i),statusGen->at(i));
    if( fabs(idGen->at(i)) != 13 ) continue;
    if( statusGen->at(i) != 1 )  continue;
    if( ptGen->at(i) > ETMu1) {
      ETMu1 = ptGen->at(i);
      iflag  = i;
      NbHEEPele ++;
    }
    else continue;
  }
  if(NbHEEPele>0) {
    GenFlag1       = iflag;
    ETMu1          = ptGen->at(iflag);
    PhiSCMu1       = phiGen->at(iflag);
    EtaSCMu1       = etaGen->at(iflag);
    EnMu1          = EnergyGen->at(iflag);
    IDele1         = idGen->at(iflag);
    Statele1       = statusGen->at(iflag);
    cout << "First generated muon= " << ETMu1 << " " << PhiSCMu1 << " " << EtaSCMu1 << endl;
    return true;
  }         
  else return false;
}
//============================ Method to select second Gen Mu ========================
bool ZprimeMuMu_FR_MiniAod::SelectSecondGenMu(unsigned GenFlag1,float ETMu1,float &ETMu2,float &PhiSCMu2,
				   float &EtaSCMu2,float &EnMu2,int &IDele2,int &Statele2){
  int NbHEEPele = 0;
  int iflag = -10;
  ETMu2 = 0.0;
  for(unsigned i=0; i<iGen->size(); i++){
    if( fabs(idGen->at(i)) != 13 ) continue;
    if( statusGen->at(i) != 1 )  continue;
    if(i == GenFlag1) continue;
    if( fabs(ptGen->at(i) - ETMu1) <0.00001 ) continue;
    if( ptGen->at(i) > ETMu2) {
      ETMu2 = ptGen->at(i);
      iflag  = i;
      NbHEEPele ++;
    }
    else continue;
  }
  if(NbHEEPele>0) {
    ETMu2      = ptGen->at(iflag);
    PhiSCMu2   = phiGen->at(iflag);
    EtaSCMu2   = etaGen->at(iflag);
    EnMu2      = EnergyGen->at(iflag);
    IDele2     = idGen->at(iflag);
    Statele2   = statusGen->at(iflag); 
    cout << "Second generated muon= " << ETMu2 << " " << PhiSCMu2 << " " << EtaSCMu2 << endl;
    return true;
  }
  else return false;
}

//----------------------------------------------------
//                                                   -
//       Part for Gen & Reco Matching                -
//                                                   -  
//----------------------------------------------------
int ZprimeMuMu_FR_MiniAod::GenToRecoMatchMu(float pTgen, float etagen, float phigen){

  int NbrecomuMatched1=0;
  int recomuflag1=-10;

  for(int i=0;i<Mu_nbMuon->size();i++){
    float deltaEta = etagen - Mu_etaTunePMuonBestTrack->at(i);
    float deltaPhi = phigen - Mu_phiTunePMuonBestTrack->at(i);
    float deltaR   = sqrt(pow(deltaEta,2)+pow(deltaPhi,2));
    float deltaPt_over_Pt  = (pTgen-Mu_ptTunePMuonBestTrack->at(i))/pTgen;
    cout << "Dr=" << deltaR << " " << deltaPt_over_Pt << " " << Mu_ptTunePMuonBestTrack->at(i) << " " << pTgen << endl;
    if(!(fabs(deltaR)<=deltaRcut && fabs(deltaPt_over_Pt)<=0.3)) continue;
    recomuflag1=i;
    NbrecomuMatched1++;
  }
  if (recomuflag1>-10) cout << "Gen Muon matched to reco muon= "<< recomuflag1 << " " << Mu_ptTunePMuonBestTrack->at(recomuflag1) << endl;

  return recomuflag1;
}

bool ZprimeMuMu_FR_MiniAod::GenRecoMatchMu1(float RecoEta1,float RecoPhi1,
				 float &ETMu1,float &PhiSCMu1,
				 float &EtaSCMu1,float &EnMu1,
				 unsigned &GenFlag1){
  int NbHEEPele = 0;
  unsigned iflag = -10;
  for(unsigned i=0; i<iGen->size(); i++){
    //printf ("BS[%d] genEta = %f  genPhi = %f recoEta = %f recoPhi = %f deltaR = %f\n",
    //    i,etaGen->at(i),phiGen->at(i),RecoEta1,RecoPhi1,deltaR);
    if( fabs(idGen->at(i)) != 13 ) continue;
    if( statusGen->at(i) != 1 )  continue;
    float deltaEta = RecoEta1 - etaGen->at(i);
    float deltaPhi = RecoPhi1 - phiGen->at(i);
    float deltaR   = sqrt(pow(deltaEta,2)+pow(deltaPhi,2));
    if(fabs(deltaR)>deltaRcut) continue;
    iflag  = i;
    NbHEEPele ++;
  }
  if(NbHEEPele > 0) {
    GenFlag1       = iflag;
    ETMu1          = ptGen->at(iflag);
    PhiSCMu1       = phiGen->at(iflag);
    EtaSCMu1       = etaGen->at(iflag);
    EnMu1          = EnergyGen->at(iflag);
    //    printf ("AS[%d] genEta = %f  genPhi = %f recoEta = %f recoPhi = %f\n",iflag,etaGen->at(iflag),phiGen->at(iflag),RecoEta1,RecoPhi1);
    return true;         
  }
  else return false;
}

bool ZprimeMuMu_FR_MiniAod::GenRecoMatchMu2(unsigned GenFlag1,float RecoEta2,float RecoPhi2,
				 float &ETMu2,float &PhiSCMu2,
				 float &EtaSCMu2,float &EnMu2){
  int NbHEEPele = 0;
  unsigned iflag = -10;
  for(unsigned i=0; i<iGen->size(); i++){
    //printf ("BS[%d] genEta = %f  genPhi = %f recoEta = %f recoPhi = %f deltaR = %f\n",
    //    i,etaGen->at(i),phiGen->at(i),RecoEta1,RecoPhi1,deltaR);
    if( fabs(idGen->at(i)) != 13 ) continue;
    if( statusGen->at(i) != 1 )  continue;
    if( i ==  GenFlag1)  continue;
    float deltaEta2 = RecoEta2 - etaGen->at(i);
    float deltaPhi2 = RecoPhi2 - phiGen->at(i);
    float deltaR2   = sqrt(pow(deltaEta2,2)+pow(deltaPhi2,2));
    if(fabs(deltaR2)>deltaRcut) continue;
    iflag  = i;
    NbHEEPele ++;
  }
  if(NbHEEPele > 0) {
    ETMu2          = ptGen->at(iflag);
    PhiSCMu2       = phiGen->at(iflag);
    EtaSCMu2       = etaGen->at(iflag);
    EnMu2          = EnergyGen->at(iflag);
    //    printf ("AS[%d] genEta = %f  genPhi = %f recoEta = %f recoPhi = %f\n",iflag,etaGen->at(iflag),phiGen->at(iflag),RecoEta1,RecoPhi1);
    return true;         
  }
  else return false;
}


bool ZprimeMuMu_FR_MiniAod::SelectFirstMuon(float &pTmuon1,float &Enmuon1,float &Etamuon1,
				    float &Phimuon1,int &ChargeMu1,unsigned &FlagMu1,
				    float &pxmuon1,float &pymuon1,float &pzmuon1,
				    float &pmuon1,float &dxymuon1,float &pTmuon1tuneP,
				    float &pTmuonBestTrack1){
  int NbHEEPele = 0;
  unsigned iflag = -10;
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if( Mu_isGlobalMuon->at(i) == 1 &&
	Mu_ptTunePMuonBestTrack->at(i) > 45.0 &&
	Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
	(Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
	Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
	Mu_numberOfValidPixelHits->at(i) > 0 &&
	Mu_numberOfValidMuonHits->at(i) > 0 &&
	Mu_numberOfMatchedStations->at(i) > 1 &&
	Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ) { 
      iflag  = i;
      NbHEEPele ++;
    }
    else continue;
  }
  if( NbHEEPele > 0 ){
    FlagMu1             = iflag;
    pTmuon1             = Mu_ptTunePMuonBestTrack->at(iflag);
    Enmuon1             = Mu_en->at(iflag);
    Etamuon1            = Mu_etaTunePMuonBestTrack->at(iflag);
    Phimuon1            = Mu_phiTunePMuonBestTrack->at(iflag);
    ChargeMu1           = Mu_chargeTunePMuonBestTrack->at(iflag);
    pxmuon1             = Mu_pxTunePMuonBestTrack->at(iflag);
    pymuon1             = Mu_pyTunePMuonBestTrack->at(iflag);
    pzmuon1             = Mu_pzTunePMuonBestTrack->at(iflag);
    pmuon1              = Mu_pTunePMuonBestTrack->at(iflag);
    dxymuon1            = Mu_absdxyTunePMuonBestTrack->at(iflag);
    pTmuon1tuneP        = Mu_ptTunePMuonBestTrack->at(iflag);
    pTmuonBestTrack1    = Mu_ptTunePMuonBestTrack->at(iflag);
    return true;  
  }
  else return false;
}
//============================ Method to select second high pt muon ========================
bool ZprimeMuMu_FR_MiniAod::SelectSecondMuon(int ChargeMu1,unsigned FlagMu1,float pTmuon1,float Etamuon1,float Phimuon1,
				     float &pTmuon2,float &Enmuon2,
				     float &Etamuon2,float &Phimuon2,int &ChargeMu2,float &pxmuon2,
				     float &pymuon2,float &pzmuon2,float &pmuon2,float &dxymuon2,
				     float &pTmuon2tuneP,float &pTmuonBestTrack2)
{
  int NbHEEPele = 0;
  unsigned iflag = -10;
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(i == FlagMu1) continue;
    if(Mu_ptTunePMuonBestTrack->at(i) == pTmuon1) continue;
    if(ChargeMu1*Mu_chargeTunePMuonBestTrack->at(i)>0) continue;
    if( Mu_isGlobalMuon->at(i) == 1 &&
	Mu_ptTunePMuonBestTrack->at(i) > 45.0 &&
	Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
	(Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
	Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 && 
	Mu_numberOfValidPixelHits->at(i) > 0 &&
	Mu_numberOfValidMuonHits->at(i) > 0 &&
	Mu_numberOfMatchedStations->at(i) > 1 &&
	Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ) { 
      iflag  = i;
      NbHEEPele ++;
    }
    else continue;
  }
  if( NbHEEPele > 0 ){
    pTmuon2          = Mu_ptTunePMuonBestTrack->at(iflag);
    Enmuon2          = Mu_en->at(iflag);
    Etamuon2         = Mu_etaTunePMuonBestTrack->at(iflag);
    Phimuon2         = Mu_phiTunePMuonBestTrack->at(iflag);
    ChargeMu2        = Mu_chargeTunePMuonBestTrack->at(iflag);
    pxmuon2          = Mu_pxTunePMuonBestTrack->at(iflag);
    pymuon2          = Mu_pyTunePMuonBestTrack->at(iflag);
    pzmuon2          = Mu_pzTunePMuonBestTrack->at(iflag);
    pmuon2           = Mu_pTunePMuonBestTrack->at(iflag);
    dxymuon2         = Mu_absdxyTunePMuonBestTrack->at(iflag);
    pTmuon2tuneP     = Mu_ptTunePMuonBestTrack->at(iflag);
    pTmuonBestTrack2 = Mu_ptTunePMuonBestTrack->at(iflag);
    return true;  
  }
  else return false;
}
bool ZprimeMuMu_FR_MiniAod::SelectSecondSSMuon(int ChargeMu1,unsigned FlagMu1,float pTmuon1,float &pTmuon2,float &Enmuon2,
                                     float &Etamuon2,float &Phimuon2,int &ChargeMu2,float &pxmuon2,
                                     float &pymuon2,float &pzmuon2,float &pmuon2,float &dxymuon2,
                                     float &pTmuon2tuneP,float &pTmuonBestTrack2)
{
  int NbHEEPele = 0;
  unsigned iflag = -10;
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(i == FlagMu1) continue;
    if(Mu_ptTunePMuonBestTrack->at(i) == pTmuon1) continue;
    if(ChargeMu1*Mu_chargeTunePMuonBestTrack->at(i)<0) continue;
    if( Mu_isGlobalMuon->at(i) == 1 &&
	Mu_ptTunePMuonBestTrack->at(i) > 45.0 &&
	Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
	(Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
	Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
	Mu_numberOfValidPixelHits->at(i) > 0 &&
	Mu_numberOfValidMuonHits->at(i) > 0 &&
	Mu_numberOfMatchedStations->at(i) > 1 &&
	Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ) {
      iflag  = i;
      NbHEEPele ++;
    }
    else continue;
  }
  if( NbHEEPele > 0 ){
    pTmuon2          = Mu_ptTunePMuonBestTrack->at(iflag);
    Enmuon2          = Mu_en->at(iflag);
    Etamuon2         = Mu_etaTunePMuonBestTrack->at(iflag);
    Phimuon2         = Mu_phiTunePMuonBestTrack->at(iflag);
    ChargeMu2        = Mu_chargeTunePMuonBestTrack->at(iflag);
    pxmuon2          = Mu_pxTunePMuonBestTrack->at(iflag);
    pymuon2          = Mu_pyTunePMuonBestTrack->at(iflag);
    pzmuon2          = Mu_pzTunePMuonBestTrack->at(iflag);
    pmuon2           = Mu_pTunePMuonBestTrack->at(iflag);
    dxymuon2         = Mu_absdxyTunePMuonBestTrack->at(iflag);
    pTmuon2tuneP     = Mu_ptTunePMuonBestTrack->at(iflag);
    pTmuonBestTrack2 = Mu_ptTunePMuonBestTrack->at(iflag);
    return true;
  }
  else return false;
}

bool ZprimeMuMu_FR_MiniAod::RecoHLTMuonMatching(float RecoEta,float RecoPhi){
  int nbMatch = 0;
  float deltaR   = -10000.0;
  //for(unsigned i=0; i<MuHLTMatch_nbMuonMatchHLT->size(); i++){
  // for(unsigned i=0; i<MuHLTMatch_nbHLTMuonMatchReco->size(); i++){
  for(unsigned i=0; i<HLTObj_nbObj->size(); i++){
    if( HLTObj_collection->at(i) == "HLT_Mu50_v1" || 
	HLTObj_collection->at(i) == "HLT_Mu50_v2" ||
	HLTObj_collection->at(i) == "HLT_Mu50_v3" || 
        HLTObj_collection->at(i) == "HLT_Mu50_v4" ||
	HLTObj_collection->at(i) == "HLT_Mu50_v5" || 
        HLTObj_collection->at(i) == "HLT_Mu50_v6" ||
        HLTObj_collection->at(i) == "HLT_Mu50_v7" ||
        HLTObj_collection->at(i) == "HLT_Mu50_v8" ||
        HLTObj_collection->at(i) == "HLT_Mu50_v9" ||
        HLTObj_collection->at(i) == "HLT_Mu50_v10" ){
	deltaR   = delR(HLTObj_eta->at(i),HLTObj_phi->at(i),RecoEta,RecoPhi);
    //    deltaR   = delR(MuHLTMatch_Trigger_eta->at(i),MuHLTMatch_Trigger_phi->at(i),RecoEta,RecoPhi);
    //printf ("HLT_Eta = %f  HLT_Phi = %f recoEta = %f recoPhi = %f DelR_trigger = %f\n",MuHLTMatch_Trigger_eta->at(i),MuHLTMatch_Trigger_phi->at(i),RecoEta,RecoPhi,deltaR);
    if(fabs(deltaR)>RecoHLTMatchingDeltaRcut) continue;
    nbMatch++;
  }
  }
  if(nbMatch>0) return true;
  else return false;
}


float ZprimeMuMu_FR_MiniAod::delR(float eta1,float phi1,float eta2,float phi2){
  float mpi=3.14;
  float dp=std::abs(phi1-phi2);
  if (dp>mpi) dp-=float(2*mpi);
  return sqrt((eta1-eta2)*(eta1-eta2) + dp*dp);
}
