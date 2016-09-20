//=========================================================================
//          Analysis code for Z' boson to Mu Mu analysis                  =  
//                  Tag and Proble code for Muon ID                       = 
//                        To run over MC                                  =
//                        Counting method                                 =  
//                    Author: Sherif Elgammal                             =
//                         1/1/2015                                       =
//=========================================================================
#define TagProbeMuon_cxx
#include "TagProbeMuon.h"
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
using namespace std;
void TagProbeMuon::Loop()
{
  time_t start,end;
  double dif;
  time (&start);
  FILE * pFile;
  pFile = fopen ("myfile.txt","w");
  binEta = 50;
  EtaMin = -3.0;
  EtaMax = 3.0;
  binMass = 60.0;
  MassMin = 60.0;
  MassMax = 120.0;
  MassCutMin = 60.0;
  MassCutMax = 120.0;
  ptMin = 0.0; 
  ptMax = 200.0; 
  RecoGenMatchingDeltaRcut = 0.15;
  RecoHLTMatchingDeltaRcut = 0.20;
  etaCut = 2.4;
  TagPtCut = 20.0; 
  TagProbeEtCut = 20.0;
  nbTP = 0;
  nbTT = 0;
  nbTF = 0;
  //nbPP = 0;
  Eff  = 0.0;
  nbTP1  = 0;
  nbTP2  = 0;
  nbTP3  = 0;
  nbTP4  = 0;
  nbTP5  = 0;
  nbTP6  = 0;
  nbTP7  = 0;
  nbTP8  = 0;
  nbTP9  = 0;
  nbTP10 = 0;
  nbTP11 = 0;
  nbTP12 = 0;
  nbTP13 = 0;
  nbTP14 = 0;
  nbTP15 = 0;
  nbTP16 = 0;
  nbTP17 = 0;
  
  nbTT1  = 0;
  nbTT2  = 0;
  nbTT3  = 0;
  nbTT4  = 0;
  nbTT5  = 0;
  nbTT6  = 0;
  nbTT7  = 0;
  nbTT8  = 0;
  nbTT9  = 0;
  nbTT10 = 0;
  nbTT11 = 0;
  nbTT12 = 0;
  nbTT13 = 0;
  nbTT14 = 0;
  nbTT15 = 0;
  nbTT16 = 0;
  nbTT17 = 0;
  //================================================================================== 
  //                                                                                 =
  //            Start the histograms for muons from Z                             =
  //                                                                                 =
  //==================================================================================
  h1_invMassTP_  = new TH1F("invMassTP","",binMass,MassCutMin,MassCutMax);
  h1_invMassTF_  = new TH1F("invMassTF","",binMass,MassCutMin,MassCutMax);
  h1_invMassTT_  = new TH1F("invMassTT","",binMass,MassCutMin,MassCutMax);
  h1_PtProbeTT_  = new TH1F("PtProbeTT","",10,0.0,100.0);
  h1_PtProbeTP_  = new TH1F("PtProbeTP","",10,0.0,100.0);
  h1_PtProbeTF_  = new TH1F("PtProbeTF","",10,0.0,100.0);
  h1_EtaProbeTT_ = new TH1F("EtaProbeTT","",binEta,EtaMin,EtaMax);
  h1_EtaProbeTP_ = new TH1F("EtaProbeTP","",binEta,EtaMin,EtaMax);
  h1_EtaProbeTF_ = new TH1F("EtaProbeTF","",binEta,EtaMin,EtaMax);
  h1_ZpeakMassTP_  = new TH1F("ZpeakMassTP","",12,60.0,120.0);
  h1_ZpeakMassTF_  = new TH1F("ZpeakMassTF","",12,60.0,120.0);
  h1_ZpeakMassTT_  = new TH1F("ZpeakMassTT","",12,60.0,120.0);
  h1_PhiProbeTT_ = new TH1F("PhiProbeTT","",20,-4.0,4.0);
  h1_PhiProbeTP_ = new TH1F("PhiProbeTP","",20,-4.0,4.0);
  h1_PhiProbeTF_ = new TH1F("PhiProbeTF","",20,-4.0,4.0);
  h1_PVProbeTT_  = new TH1F("PVProbeTT","",15,0.0,50.0);
  h1_PVProbeTP_  = new TH1F("PVProbeTP","",15,0.0,50.0);
  h1_PVProbeTF_  = new TH1F("PVProbeTF","",15,0.0,50.0);
  h1_NBofTP_     = new TH1F("NBofTP","",100,0.0,10.0);
  h1_NBofTF_     = new TH1F("NBofTF","",100,0.0,10.0);
  h1_NBofTT_     = new TH1F("NBofTT","",100,0.0,10.0);

  h1_invMassPP_    = new TH1F("invMassPP","",binMass,MassCutMin,MassCutMax);
  h1_ZpeakMassPP_  = new TH1F("ZpeakMassPP","",12,60.0,120.0);
  h1_PtProbePP_    = new TH1F("PtProbePP","",100,0.0,500.0);
  h1_EtaProbePP_   = new TH1F("EtaProbePP","",binEta,EtaMin,EtaMax);
  h1_PhiProbePP_   = new TH1F("PhiProbePP","",20,-4.0,4.0);
  h1_PVProbePP_    = new TH1F("PVProbePP","",15,0.0,50.0);
  h1_NBofPP_       = new TH1F("NBofPP","",100,0.0,10.0);


  //==================================================================================
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    /*
      cout<<"=======> jentry = "<<jentry<< 
      "=======> Evt = "<<event_evtNo<< 
      "=======> Run = "<<event_runNo<< 
      "=======> Lumi = "<<event_lumi<< 
      "=======> bunch = "<<event_bunch<<endl;
    */
    //=====================================================
    //                                                    =
    //  Calling methods to get the generated information  =
    //                                                    =
    //=====================================================
    //bool firstGenMu  = SelectFirstGenMu(genET1,genPhi1,genEta1,genEn1,genID1,genStat1,flag1);
    //bool secondGenMu = SelectSecondGenMu(flag1,genET1,genET2,genPhi2,genEta2,genEn2,genID2,genStat2);
    //if(firstGenMu == 0 || secondGenMu == 0) continue;
    //MassGen = Mass(genET1,genEta1,genPhi1,genEn1,genET2,genEta2,genPhi2,genEn2);
    //acceptance conditions
    //if(fabs(genEta1)>etaCut || fabs(genEta2)>etaCut) continue;
    //if(MassGen<60.0) continue;
    TagandProbeResults();
  }
  /*
  cout<<"nbTP = "<< nbTP << endl;
  cout<<"nbTT = "<< nbTT << endl;
  cout<<"nbTF = "<< nbTF << endl;
  cout<<"nbPP = "<< nbPP << endl;
  //Eff = (float)(nbTT)/(nbTP);
  //Eff = (float)(nbTP)/(nbTP+nbTF);
  */
  cout<<"float Eff1    = "<< (float) nbTT1/nbTP1 << ";" << endl;
  cout<<"float Eff2    = "<< (float) nbTT2/nbTP2 << ";" <<endl;
  cout<<"float Eff3    = "<< (float) nbTT3/nbTP3 << ";" << endl;
  cout<<"float Eff4    = "<< (float) nbTT4/nbTP4 << ";" << endl;
  cout<<"float Eff5    = "<< (float) nbTT5/nbTP5 << ";" << endl;
  cout<<"float Eff6    = "<< (float) nbTT6/nbTP6 << ";" << endl;
  cout<<"float Eff7    = "<< (float) nbTT7/nbTP7 << ";" << endl;
  cout<<"float Eff8    = "<< (float) nbTT8/nbTP8 << ";" << endl;
  cout<<"float Eff9    = "<< (float) nbTT9/nbTP9 << ";" << endl;
  cout<<"float Eff10   = "<< (float) nbTT10/nbTP10 << ";" << endl;
  cout<<"float Eff11   = "<< (float) nbTT11/nbTP11 << ";" << endl;
  cout<<"float Eff12   = "<< (float) nbTT12/nbTP12 << ";" << endl;
  cout<<"float Eff13   = "<< (float) nbTT13/nbTP13 << ";" << endl;
  cout<<"[0.0,0.9]   = "<< (float) nbTT14/nbTP14 << endl;
  cout<<"[0.9,1.2]   = "<< (float) nbTT15/nbTP15 << endl;
  cout<<"[1.2,2.1]   = "<< (float) nbTT16/nbTP16 << endl;
  cout<<"[2.1,2.4]   = "<< (float) nbTT17/nbTP17 << endl;



  cout<<"float error1     = "<< sqrt(nbTP1-nbTT1)/nbTP1 << ";" << endl;
  cout<<"float error2     = "<< sqrt(nbTP2-nbTT2)/nbTP2 << ";" << endl;
  cout<<"float error3     = "<< sqrt(nbTP3-nbTT3)/nbTP3 << ";" << endl;
  cout<<"float error4     = "<< sqrt(nbTP4-nbTT4)/nbTP4 << ";" << endl;
  cout<<"float error5     = "<< sqrt(nbTP5-nbTT5)/nbTP5 << ";" << endl;
  cout<<"float error6     = "<< sqrt(nbTP6-nbTT6)/nbTP6 << ";" << endl;
  cout<<"float error7     = "<< sqrt(nbTP7-nbTT7)/nbTP7 << ";" << endl;
  cout<<"float error8     = "<< sqrt(nbTP8-nbTT8)/nbTP8 << ";" << endl;
  cout<<"float error9     = "<< sqrt(nbTP9-nbTT9)/nbTP9 << ";" << endl;
  cout<<"float error10    = "<< sqrt(nbTP10-nbTT10)/nbTP10 << ";" << endl;
  cout<<"float error11    = "<< sqrt(nbTP11-nbTT11)/nbTP11 << ";" << endl;
  cout<<"float error12    = "<< sqrt(nbTP12-nbTT12)/nbTP12 << ";" << endl;
  cout<<"float error13    = "<< sqrt(nbTP13-nbTT13)/nbTP13 << ";" << endl;
  cout<<"error[0.0,0.9]    = "<< sqrt(nbTP14-nbTT14)/nbTP14 << ";" << endl;
  cout<<"error[0.9,1.2]    = "<< sqrt(nbTP15-nbTT15)/nbTP15 << ";" << endl;
  cout<<"error[1.2,2.1]    = "<< sqrt(nbTP16-nbTT16)/nbTP16 << ";" << endl;
  cout<<"error[2.1,2.4]    = "<< sqrt(nbTP17-nbTT17)/nbTP17 << ";" << endl;

  cout<<"**************************************"<<endl;
  printf ("* Tag & Probe Eff at Z pole is %.3lf *\n",Eff);
  cout<<"**************************************"<<endl;
  cout<<"************** The T&P exercise is ****************"<<endl;
  cout<<"Try to write a Macro to plot Eff vs eta of the probe muon!"<<endl;
  cout<<"***************************************************"<<endl;
  fclose (pFile);
  //================================================================================== 
  //                                                                                 =
  //                         writing Histograms to a file                            =             
  //                                                                                 =
  //==================================================================================
  //TFile *output = new TFile("Analyse_TagAndProbe_CMSSW_8_0_13_ZPRIMEMuMu_13TeV-DataG-V1-JSON-tree.root","recreate");
  TFile *output = new TFile("Analyse_TagAndProbe_CMSSW8012_MC_reHLT_DYtoMuMu120to200_13TeV_tree.root","recreate");

  output->cd();
  /** Write the histograms for PF Jets in general */
  h1_invMassTP_->Write();
  h1_invMassTF_->Write();
  h1_invMassTT_->Write();
  h1_PtProbeTT_->Write();
  h1_PtProbeTP_->Write();
  h1_PtProbeTF_->Write();
  h1_EtaProbeTT_->Write();
  h1_EtaProbeTP_->Write();
  h1_EtaProbeTF_->Write();
  h1_ZpeakMassTP_->Write();
  h1_ZpeakMassTF_->Write();
  h1_ZpeakMassTT_->Write();
  h1_PhiProbeTT_->Write();
  h1_PhiProbeTP_->Write();
  h1_PhiProbeTF_->Write();
  h1_PVProbeTT_->Write();
  h1_PVProbeTP_->Write();
  h1_PVProbeTF_->Write();
  h1_NBofTP_->Write();
  h1_NBofTF_->Write();
  h1_NBofTT_->Write();

  h1_invMassPP_->Write();
  h1_ZpeakMassPP_->Write();
  h1_PtProbePP_->Write();
  h1_EtaProbePP_->Write();
  h1_PhiProbePP_->Write();
  h1_PVProbePP_->Write(); 
  h1_NBofPP_->Write();

  output->Close();
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
//========================== Method to select Tag & Probe passing Probe events ======================= 
void TagProbeMuon::SelectTagandProbePassingProbe(){
  //find Tag first
  int nbTnP = 0;
  for(unsigned mu=0; mu<Mu_nbMuon->size(); mu++){
    if(	Mu_isTrackerMuon->at(mu)==1 &&
	Mu_isGlobalMuon->at(mu) == 1 &&
	Mu_ptTunePMuonBestTrack->at(mu) > TagPtCut &&
	Mu_absdxyTunePMuonBestTrack->at(mu) < 0.2 &&
	//(Mu_trackiso->at(mu)/Mu_ptInnerTrack->at(mu)) < 0.1  //&&
	Mu_numberOftrackerLayersWithMeasurement->at(mu) > 5 && 
	Mu_numberOfValidPixelHits->at(mu) > 0 &&
	Mu_numberOfValidMuonHits->at(mu) > 0 &&
	Mu_numberOfMatchedStations->at(mu) > 1 &&
	Mu_dPToverPTTunePMuonBestTrack->at(mu) < 0.3
	){
      //bool GenRecoMatch = GenRecoMatchMuon(Mu_etaTunePMuonBestTrack->at(mu),Mu_phiTunePMuonBestTrack->at(mu));
      //if(GenRecoMatch == 0) continue;
      //bool RecoMuonMatchingWithHLT = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(mu),Mu_phiCocktail->at(mu));
      //if(RecoMuonMatchingWithHLT == 0) continue;
      //find Probe passing Probe
      for(unsigned probe=0; probe<Mu_nbMuon->size(); probe++){
	if( probe==mu ) continue;
	if(Mu_chargeTunePMuonBestTrack->at(mu)*Mu_chargeTunePMuonBestTrack->at(probe)>0) continue;
	if (Mu_isTrackerMuon->at(probe)==1 &&
	    Mu_ptTunePMuonBestTrack->at(probe) > TagProbeEtCut ) {
	  nbTnP++;
	  if(nbTnP>1) continue;
	  float MassTP = Mass(Mu_ptTunePMuonBestTrack->at(mu),Mu_etaTunePMuonBestTrack->at(mu),Mu_phiTunePMuonBestTrack->at(mu),Mu_en->at(mu),
			      Mu_ptTunePMuonBestTrack->at(probe),Mu_etaTunePMuonBestTrack->at(probe),Mu_phiTunePMuonBestTrack->at(probe),Mu_en->at(probe));
	  if( MassTP >MassCutMin && MassTP <MassCutMax )
	    {
	      h1_invMassTP_->Fill(MassTP);
	      h1_ZpeakMassTP_->Fill(MassTP);
	      //Start bins in eta of the probe
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > -2.4 && Mu_etaTunePMuonBestTrack->at(probe) < -2.1 )
		{
		  nbTP1++;
		  h1_PtProbeTP_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTP_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTP_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTP_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > -2.1 && Mu_etaTunePMuonBestTrack->at(probe) < -1.6 )
		{
		  nbTP2++;
		  h1_PtProbeTP_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTP_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTP_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTP_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > -1.6 && Mu_etaTunePMuonBestTrack->at(probe) < -1.2 )
		{
		  nbTP3++;
		  h1_PtProbeTP_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTP_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTP_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTP_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > -1.2 && Mu_etaTunePMuonBestTrack->at(probe) < -0.9 ) 
		{
		  nbTP4++;
		  h1_PtProbeTP_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTP_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTP_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTP_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > -0.9 && Mu_etaTunePMuonBestTrack->at(probe) < -0.3 ) 
		{
		  nbTP5++;
		  h1_PtProbeTP_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTP_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTP_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTP_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > -0.3 && Mu_etaTunePMuonBestTrack->at(probe) < -0.2 ) 
		{
		  nbTP6++;
		  h1_PtProbeTP_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTP_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTP_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTP_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > -0.2 && Mu_etaTunePMuonBestTrack->at(probe) < 0.2 ) 
		{
		  nbTP7++;
		  h1_PtProbeTP_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTP_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTP_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTP_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > 0.2 && Mu_etaTunePMuonBestTrack->at(probe) < 0.3 ) 
		{
		  nbTP8++;
		  h1_PtProbeTP_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTP_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTP_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTP_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > 0.3 && Mu_etaTunePMuonBestTrack->at(probe) < 0.9 ) 
		{
		  nbTP9++;
		  h1_PtProbeTP_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTP_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTP_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTP_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > 0.9 && Mu_etaTunePMuonBestTrack->at(probe) < 1.2 )
		{
		  nbTP10++;
		  h1_PtProbeTP_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTP_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTP_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTP_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > 1.2 && Mu_etaTunePMuonBestTrack->at(probe) < 1.6 )
		{
		  nbTP11++;
		  h1_PtProbeTP_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTP_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTP_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTP_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > 1.6 && Mu_etaTunePMuonBestTrack->at(probe) < 2.1 )
		{
		  nbTP12++;
		  h1_PtProbeTP_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTP_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTP_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTP_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > 2.1 && Mu_etaTunePMuonBestTrack->at(probe) < 2.4 )
		{
		  nbTP13++;
		  h1_PtProbeTP_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTP_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTP_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTP_->Fill(Mu_nbofpv->at(probe));
		}
	      //new 3 eta bins
	      if ( fabs(Mu_etaTunePMuonBestTrack->at(probe)) > 0.0 && fabs(Mu_etaTunePMuonBestTrack->at(probe)) < 0.9 )
		{
		  nbTP14++;
		}
	      if ( fabs(Mu_etaTunePMuonBestTrack->at(probe)) > 0.9 && fabs(Mu_etaTunePMuonBestTrack->at(probe)) < 1.2 )
		{
		  nbTP15++;
		}
	      if ( fabs(Mu_etaTunePMuonBestTrack->at(probe)) > 1.2 && fabs(Mu_etaTunePMuonBestTrack->at(probe)) < 2.1 )
		{
		  nbTP16++;
		}
	      if ( fabs(Mu_etaTunePMuonBestTrack->at(probe)) > 2.1 && fabs(Mu_etaTunePMuonBestTrack->at(probe)) < 2.4 )
		{
		  nbTP17++;
		}
	      nbTP++;
	      h1_NBofTP_->Fill(nbTP);
	    }
	}
      }
    }
  }
}
//========================== Method to select Tag & Probe passing Tag events =======================
void TagProbeMuon::SelectTagandProbePassingTag(){
  //find Tag first
  int nbTnT = 0;
  for(unsigned mu=0; mu<Mu_nbMuon->size(); mu++){
    if( Mu_isTrackerMuon->at(mu)==1 &&
	Mu_isGlobalMuon->at(mu) == 1 &&
	Mu_ptTunePMuonBestTrack->at(mu) > TagPtCut &&
	Mu_absdxyTunePMuonBestTrack->at(mu) < 0.2 &&
	//(Mu_trackiso->at(mu)/Mu_ptInnerTrack->at(mu)) < 0.1  //&&
	Mu_numberOftrackerLayersWithMeasurement->at(mu) > 5 && 
	Mu_numberOfValidPixelHits->at(mu) > 0 &&
	Mu_numberOfValidMuonHits->at(mu) > 0 &&
	Mu_numberOfMatchedStations->at(mu) > 1 &&
	Mu_dPToverPTTunePMuonBestTrack->at(mu) < 0.3 
	) {
      //bool GenRecoMatch1 = GenRecoMatchMuon(Mu_etaTunePMuonBestTrack->at(mu),Mu_phiCocktail->at(mu));
      //if( GenRecoMatch1 == 0 ) continue; 
      //find Proble passing Tag
      for(unsigned probe=0; probe<Mu_nbMuon->size(); probe++){
	if( probe==mu ) continue;
	if(Mu_chargeTunePMuonBestTrack->at(mu)*Mu_chargeTunePMuonBestTrack->at(probe)>0) continue;
	if( Mu_isTrackerMuon->at(probe)==1 &&
	    Mu_isGlobalMuon->at(probe) == 1 &&
	    Mu_ptTunePMuonBestTrack->at(probe) > TagPtCut &&
	    Mu_absdxyTunePMuonBestTrack->at(probe) < 0.2 &&
	    //(Mu_trackiso->at(probe)/Mu_ptInnerTrack->at(probe)) < 0.1  //&&
	    Mu_numberOftrackerLayersWithMeasurement->at(probe) > 5 && 
	    Mu_numberOfValidPixelHits->at(probe) > 0 &&
	    Mu_numberOfValidMuonHits->at(probe) > 0 &&
	    Mu_numberOfMatchedStations->at(probe) > 1 &&
	    Mu_dPToverPTTunePMuonBestTrack->at(probe) < 0.3 
	    ) {
	  //bool GenRecoMatch2 = GenRecoMatchMuon(Mu_etaTunePMuonBestTrack->at(probe),Mu_phiCocktail->at(probe));
	  //if( GenRecoMatch2 == 0 ) continue; 
	  nbTnT++;
	  if(nbTnT>1) continue;
	  float MassTT = Mass(Mu_ptTunePMuonBestTrack->at(mu),Mu_etaTunePMuonBestTrack->at(mu),Mu_phiTunePMuonBestTrack->at(mu),Mu_en->at(mu),
			      Mu_ptTunePMuonBestTrack->at(probe),Mu_etaTunePMuonBestTrack->at(probe),Mu_phiTunePMuonBestTrack->at(probe),Mu_en->at(probe));
	  if( MassTT >MassCutMin && MassTT <MassCutMax )
	    {
	      h1_invMassTT_ ->Fill(MassTT);
	      h1_ZpeakMassTT_->Fill(MassTT);
	      
	      //Start bins in eta of the probe
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > -2.4 && Mu_etaTunePMuonBestTrack->at(probe) < -2.1 )
		{
		  nbTT1++;
		  h1_PtProbeTT_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTT_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTT_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTT_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > -2.1 && Mu_etaTunePMuonBestTrack->at(probe) < -1.6 )
		{
		  nbTT2++;
		  h1_PtProbeTT_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTT_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTT_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTT_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > -1.6 && Mu_etaTunePMuonBestTrack->at(probe) < -1.2 )
		{
		  nbTT3++;
		  h1_PtProbeTT_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTT_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTT_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTT_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > -1.2 && Mu_etaTunePMuonBestTrack->at(probe) < -0.9 ) 
		{
		  nbTT4++;
		  h1_PtProbeTT_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTT_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTT_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTT_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > -0.9 && Mu_etaTunePMuonBestTrack->at(probe) < -0.3 ) 
		{
		  nbTT5++;
		  h1_PtProbeTT_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTT_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTT_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTT_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > -0.3 && Mu_etaTunePMuonBestTrack->at(probe) < -0.2 ) 
		{
		  nbTT6++;
		  h1_PtProbeTT_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTT_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTT_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTT_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > -0.2 && Mu_etaTunePMuonBestTrack->at(probe) < 0.2 ) 
		{
		  nbTT7++;
		  h1_PtProbeTT_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTT_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTT_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTT_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > 0.2 && Mu_etaTunePMuonBestTrack->at(probe) < 0.3 ) 
		{
		  nbTT8++;
		  h1_PtProbeTT_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTT_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTT_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTT_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > 0.3 && Mu_etaTunePMuonBestTrack->at(probe) < 0.9 ) 
		{
		  nbTT9++;
		  h1_PtProbeTT_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTT_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTT_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTT_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > 0.9 && Mu_etaTunePMuonBestTrack->at(probe) < 1.2 )
		{
		  nbTT10++;
		  h1_PtProbeTT_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTT_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTT_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTT_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > 1.2 && Mu_etaTunePMuonBestTrack->at(probe) < 1.6 )
		{
		  nbTT11++;
		  h1_PtProbeTT_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTT_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTT_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTT_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > 1.6 && Mu_etaTunePMuonBestTrack->at(probe) < 2.1 )
		{
		  nbTT12++;
		  h1_PtProbeTT_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTT_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTT_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTT_->Fill(Mu_nbofpv->at(probe));
		}
	      if ( Mu_etaTunePMuonBestTrack->at(probe) > 2.1 && Mu_etaTunePMuonBestTrack->at(probe) < 2.4 )
		{
		  nbTT13++;
		  h1_PtProbeTT_->Fill(Mu_ptTunePMuonBestTrack->at(probe));
		  h1_EtaProbeTT_->Fill(Mu_etaTunePMuonBestTrack->at(probe));
		  h1_PhiProbeTT_->Fill(Mu_phiTunePMuonBestTrack->at(probe));
		  h1_PVProbeTT_->Fill(Mu_nbofpv->at(probe));
		}
	      //new 3 eta bins
	      if ( fabs(Mu_etaTunePMuonBestTrack->at(probe)) > 0.0 && fabs(Mu_etaTunePMuonBestTrack->at(probe)) < 0.9 )
		{
		  nbTT14++;
		}
	      if ( fabs(Mu_etaTunePMuonBestTrack->at(probe)) > 0.9 && fabs(Mu_etaTunePMuonBestTrack->at(probe)) < 1.2 )
		{
		  nbTT15++;
		}
	      if ( fabs(Mu_etaTunePMuonBestTrack->at(probe)) > 1.2 && fabs(Mu_etaTunePMuonBestTrack->at(probe)) < 2.1 )
		{
		  nbTT16++;
		}
	      if ( fabs(Mu_etaTunePMuonBestTrack->at(probe)) > 2.1 && fabs(Mu_etaTunePMuonBestTrack->at(probe)) < 2.4 )
		{
		  nbTT17++;
		}
	      nbTT++;
	      h1_NBofTT_->Fill(nbTT);
	    }
	}
      }
    }
  }
}

void TagProbeMuon::TagandProbeResults(){
  SelectTagandProbePassingProbe();
  SelectTagandProbePassingTag();
}
//----------------------------------------------------
//                                                   -
//       Part for Gen & Reco Matching                -
//                                                   -  
//----------------------------------------------------
//========================== Method to select firt Gen Mu =======================
bool TagProbeMuon::SelectFirstGenMu(float &ETMu1,float &PhiSCMu1,
				    float &EtaSCMu1,float &EnMu1,
				    int &IDMu1,int &StatMu1,
				    unsigned &GenFlag1){
  int NbMuons = 0;
  int iflag = -10;
  ETMu1 = 0.0;
  for(unsigned i=0; i<iGen->size(); i++){
    if( fabs(idGen->at(i)) != 13 ) continue;
    if( statusGen->at(i) != 1 )  continue;
    if( ptGen->at(i) > ETMu1) {
      ETMu1 = ptGen->at(i);
      iflag  = i;
      NbMuons ++;
    }
    else continue;
  }
  if(NbMuons>0) {
    GenFlag1       = iflag;
    ETMu1          = ptGen->at(iflag);
    PhiSCMu1       = phiGen->at(iflag);
    EtaSCMu1       = etaGen->at(iflag);
    EnMu1          = EnergyGen->at(iflag);
    IDMu1         = idGen->at(iflag);
    StatMu1       = statusGen->at(iflag);
    return true;
  }         
  else return false;
}
//============================ Method to select second Gen Mu ========================
bool TagProbeMuon::SelectSecondGenMu(unsigned GenFlag1,float ETMu1,float &ETMu2,float &PhiSCMu2,
				     float &EtaSCMu2,float &EnMu2,int &IDMu2,int &StatMu2){
  int NbMuons = 0;
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
      NbMuons ++;
    }
    else continue;
  }
  if(NbMuons>0) {
    ETMu2      = ptGen->at(iflag);
    PhiSCMu2   = phiGen->at(iflag);
    EtaSCMu2   = etaGen->at(iflag);
    EnMu2      = EnergyGen->at(iflag);
    IDMu2     = idGen->at(iflag);
    StatMu2   = statusGen->at(iflag); 
    return true;
  }
  else return false;
}

bool TagProbeMuon::GenRecoMatchMuon(float RecoEta1,float RecoPhi1){
  int NbMuons = 0;
  unsigned iflag = -10;
  for(unsigned i=0; i<iGen->size(); i++){
    if( fabs(idGen->at(i)) != 13 ) continue;
    if( statusGen->at(i) != 1 )  continue;
    float deltaEta = RecoEta1 - etaGen->at(i);
    float deltaPhi = RecoPhi1 - phiGen->at(i);
    float deltaR   = sqrt(pow(deltaEta,2)+pow(deltaPhi,2));
    if(fabs(deltaR)>RecoGenMatchingDeltaRcut) continue;
    iflag  = i;
    NbMuons ++;
  }
  if(NbMuons > 0) {
    return true;         
  }
  else return false;
}
//===================== Methode to calculate the mass ========================
float TagProbeMuon::Mass(float Et1,float Eta1,float Phi1,float En1,
			 float Et2,float Eta2,float Phi2,float En2){
  float MassMuMu = 0.0;
  TLorentzVector Mu1;
  TLorentzVector Mu2;
  Mu1.SetPtEtaPhiE(Et1,Eta1,Phi1,En1);
  Mu2.SetPtEtaPhiE(Et2,Eta2,Phi2,En2);
  MassMuMu = (Mu1 + Mu2).M();
  return MassMuMu;
}

float TagProbeMuon::delR(float eta1,float phi1,float eta2,float phi2){
   float mpi=3.14;
   float dp=std::abs(phi1-phi2);
   if (dp>mpi) dp-=float(2*mpi);
   return sqrt((eta1-eta2)*(eta1-eta2) + dp*dp);
}

//----------------------------------------------------
//                                                   -
//       Part for HLT & Reco Matching                -
//                                                   -  
//----------------------------------------------------
bool TagProbeMuon::RecoHLTMuonMatching(float RecoEta,float RecoPhi){
  int nbMatch = 0;
  float deltaR   = -10000.0;
  for(unsigned i=0; i<HLTObj_nbObj->size(); i++){
    //cout<<"[before]triggerName"<<HLTObj_collection->at(i) <<endl;
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
      //cout<<"[after]triggerName"<<HLTObj_collection->at(i) <<endl;
	deltaR   = delR(HLTObj_eta->at(i),HLTObj_phi->at(i),RecoEta,RecoPhi);
        //printf ("HLT_Eta = %f  HLT_Phi = %f recoEta = %f recoPhi = %f DelR_trigger = %f\n",HLTObj_eta->at(i),HLTObj_phi->at(i),RecoEta,RecoPhi,deltaR);
      if(fabs(deltaR)>RecoHLTMatchingDeltaRcut) continue;
      nbMatch++;
    }
  }
  if(nbMatch>0) return true;
  else return false;
}



