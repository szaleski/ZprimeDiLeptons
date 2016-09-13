//==============================================================
//          Analysis code for Z' boson to Mu Mu analysis       =  
//          In this code we select the high pt di-muons events = 
//                  To run over MINIAOD MC                     = 
//                  Author:  Sherif Elgammal                   = 
//                                                             = 
//                       17/5/2016                             = 
//==============================================================
#ifndef ZprimeMuMu_FR_MiniAod_h
#define ZprimeMuMu_FR_MiniAod_h
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
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <fstream>
#include <sstream>
#include "TVector3.h"
// Header file for the classes stored in the TTree if any.
#include <vector>
using namespace std;

// Fixed size dimensions of array or collections stored in the TTree if any.

class ZprimeMuMu_FR_MiniAod {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           event_runNo;
   Int_t           event_evtNo;
   Int_t           event_lumi;
   Int_t           event_bunch;
   vector<int>     *iGen;
   vector<int>     *idGen;
   vector<int>     *statusGen;
   vector<float>   *ptGen;
   vector<float>   *etaGen;
   vector<float>   *phiGen;
   vector<int>     *chargeGen;
   vector<float>   *EnergyGen;
   vector<float>   *pxGen;
   vector<float>   *pyGen;
   vector<float>   *pzGen;
   vector<int>     *Mu_nbMuon;
   vector<bool>    *Mu_isGlobalMuon;
   vector<bool>    *Mu_isPF;
   vector<bool>    *Mu_isTrackerMuon;
   vector<float>   *Mu_et;
   vector<float>   *Mu_en;
   vector<float>   *Mu_pt;
   vector<float>   *Mu_eta;
   vector<float>   *Mu_phi;
   vector<float>   *Mu_charge;
   vector<float>   *Mu_ptTunePMuonBestTrack;
   vector<float>   *Mu_pxTunePMuonBestTrack;
   vector<float>   *Mu_pyTunePMuonBestTrack;
   vector<float>   *Mu_pzTunePMuonBestTrack;
   vector<float>   *Mu_pTunePMuonBestTrack;
   vector<float>   *Mu_etaTunePMuonBestTrack;
   vector<float>   *Mu_phiTunePMuonBestTrack;
   vector<float>   *Mu_thetaTunePMuonBestTrack;
   vector<float>   *Mu_chargeTunePMuonBestTrack;
   vector<float>   *Mu_dPToverPTTunePMuonBestTrack;
   vector<float>   *Mu_absdxyTunePMuonBestTrack;
   vector<float>   *Mu_absdzTunePMuonBestTrack;
   vector<float>   *Mu_ptInnerTrack;
   vector<float>   *Mu_normalizedChi2;
   vector<float>   *Mu_vtxMass;
   vector<float>   *Mu_vtxNormChi2;
   vector<float>   *Mu_vtxMassLept;
   vector<int>     *Mu_numberOfMatchedStations;
   vector<int>     *Mu_numberOfValidPixelHits;
   vector<int>     *Mu_numberOfValidMuonHits;
   vector<int>     *Mu_numberOftrackerLayersWithMeasurement;
   vector<float>   *Mu_emIso;
   vector<float>   *Mu_hadIso;
   vector<float>   *Mu_trackiso;
   vector<float>   *Mu_pfSumChargedHadronPt;
   vector<float>   *Mu_pfSumNeutralHadronEt;
   vector<float>   *Mu_PFSumPhotonEt;
   vector<float>   *Mu_pfSumPUPt;
   vector<int>     *Mu_nbofpv;
   vector<int>     *HLT_nb;
   vector<string>  *HLT_name;
   vector<bool>    *HLT_isaccept;
   vector<int>     *HLTObj_nbObj;
   vector<float>   *HLTObj_pt;
   vector<float>   *HLTObj_eta;
   vector<float>   *HLTObj_phi;
   vector<string>  *HLTObj_collection;
   Double_t        GenMet_pt;
   Double_t        PFMet_pt_cor;
   Double_t        PFMet_phi_cor;
   Double_t        PFMet_en_cor;
   Double_t        PFMet_px_cor;
   Double_t        PFMet_py_cor;
   Double_t        PFMet_pz_cor;
   Double_t        PFMet_sumEt_cor;
   Double_t        CaloMet_pt;
   Double_t        CaloMet_phi;
   Double_t        CaloMet_sumEt;
   Double_t        METSign;
   Double_t        PFMet_shiftedPt_JetEnUp;
   Double_t        PFMet_shiftedPt_JetEnDown;

   // List of branches
   TBranch        *b_event_runNo;   //!
   TBranch        *b_event_evtNo;   //!
   TBranch        *b_event_lumi;   //!
   TBranch        *b_event_bunch;   //!
   TBranch        *b_iGen;   //!
   TBranch        *b_idGen;   //!
   TBranch        *b_statusGen;   //!
   TBranch        *b_ptGen;   //!
   TBranch        *b_etaGen;   //!
   TBranch        *b_phiGen;   //!
   TBranch        *b_chargeGen;   //!
   TBranch        *b_EnergyGen;   //!
   TBranch        *b_pxGen;   //!
   TBranch        *b_pyGen;   //!
   TBranch        *b_pzGen;   //!
   TBranch        *b_Mu_nbMuon;   //!
   TBranch        *b_Mu_isGlobalMuon;   //!
   TBranch        *b_Mu_isPF;   //!
   TBranch        *b_Mu_isTrackerMuon;   //!
   TBranch        *b_Mu_et;   //!
   TBranch        *b_Mu_en;   //!
   TBranch        *b_Mu_pt;   //!
   TBranch        *b_Mu_eta;   //!
   TBranch        *b_Mu_phi;   //!
   TBranch        *b_Mu_charge;   //!
   TBranch        *b_Mu_ptTunePMuonBestTrack;   //!
   TBranch        *b_Mu_pxTunePMuonBestTrack;   //!
   TBranch        *b_Mu_pyTunePMuonBestTrack;   //!
   TBranch        *b_Mu_pzTunePMuonBestTrack;   //!
   TBranch        *b_Mu_pTunePMuonBestTrack;   //!
   TBranch        *b_Mu_etaTunePMuonBestTrack;   //!
   TBranch        *b_Mu_phiTunePMuonBestTrack;   //!
   TBranch        *b_Mu_thetaTunePMuonBestTrack;   //!
   TBranch        *b_Mu_chargeTunePMuonBestTrack;   //!
   TBranch        *b_Mu_dPToverPTTunePMuonBestTrack;   //!
   TBranch        *b_Mu_absdxyTunePMuonBestTrack;   //!
   TBranch        *b_Mu_absdzTunePMuonBestTrack;   //!
   TBranch        *b_Mu_ptInnerTrack;   //!
   TBranch        *b_Mu_normalizedChi2;   //!
   TBranch        *b_Mu_vtxMass;   //!
   TBranch        *b_Mu_vtxNormChi2;   //!
   TBranch        *b_Mu_vtxMassLept;   //!
   TBranch        *b_Mu_numberOfMatchedStations;   //!
   TBranch        *b_Mu_numberOfValidPixelHits;   //!
   TBranch        *b_Mu_numberOfValidMuonHits;   //!
   TBranch        *b_Mu_numberOftrackerLayersWithMeasurement;   //!
   TBranch        *b_Mu_emIso;   //!
   TBranch        *b_Mu_hadIso;   //!
   TBranch        *b_Mu_trackiso;   //!
   TBranch        *b_Mu_pfSumChargedHadronPt;   //!
   TBranch        *b_Mu_pfSumNeutralHadronEt;   //!
   TBranch        *b_Mu_PFSumPhotonEt;   //!
   TBranch        *b_Mu_pfSumPUPt;   //!
   TBranch        *b_Mu_nbofpv;   //!
   TBranch        *b_HLT_nb;   //!
   TBranch        *b_HLT_name;   //!
   TBranch        *b_HLT_isaccept;   //!
   TBranch        *b_HLTObj_nbObj;   //!
   TBranch        *b_HLTObj_pt;   //!
   TBranch        *b_HLTObj_eta;   //!
   TBranch        *b_HLTObj_phi;   //!
   TBranch        *b_HLTObj_collection;   //!
   TBranch        *b_GenMet_pt;   //!
   TBranch        *b_PFMet_pt_cor;   //!
   TBranch        *b_PFMet_phi_cor;   //!
   TBranch        *b_PFMet_en_cor;   //!
   TBranch        *b_PFMet_px_cor;   //!
   TBranch        *b_PFMet_py_cor;   //!
   TBranch        *b_PFMet_pz_cor;   //!
   TBranch        *b_PFMet_sumEt_cor;   //!
   TBranch        *b_CaloMet_pt;   //!
   TBranch        *b_CaloMet_phi;   //!
   TBranch        *b_CaloMet_sumEt;   //!
   TBranch        *b_METSign;   //!
   TBranch        *b_PFMet_shiftedPt_JetEnUp;   //!
   TBranch        *b_PFMet_shiftedPt_JetEnDown;   //!
   
   ZprimeMuMu_FR_MiniAod(Char_t namechar_[300],TTree *tree=0,Double_t weight_=1.,string DATA_type_="DATA",string MC_type_="MC");
			virtual ~ZprimeMuMu_FR_MiniAod();
			Double_t weight;
			Char_t name[300];
			string DATA_type,MC_type;
			
			virtual Int_t    Cut(Long64_t entry);
			virtual Int_t    GetEntry(Long64_t entry);
			virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void PrintEventInformation(int runNumber, int lumiNumber, int eventNumber,
			      float vtxChi2, float vtxMass, float CosmicRejection);
   bool SelectFirstMuon(float &pTmuon1,float &Enmuon1,float &Etamuon1,
			float &Phimuon1,int &ChargeMu1,unsigned &FlagMu1,
			float &pxmuon1,float &pymuon1,float &pzmuon1,
			float &pmuon1,float &dxymuon1,float &pTmuon1tuneP,
			float &pTmuonBestTrack1);
   bool SelectSecondMuon(int ChargeMu1,unsigned FlagMu1,float pTmuon1,float Etamuon1,float Phimuon1,
			 float &pTmuon2,float &Enmuon2,
					      float &Etamuon2,float &Phimuon2,int &ChargeMu2,float &pxmuon2,
			 float &pymuon2,float &pzmuon2,float &pmuon2,float &dxymuon2,
			 float &pTmuon2tuneP,float &pTmuonBestTrack2);
   bool SelectSecondSSMuon(int ChargeMu1,unsigned FlagMu1,float pTmuon1,float &pTmuon3,float &Enmuon3,
                           float &Etamuon3,float &Phimuon3,int &ChargeMu3,float &pxmuon3,
                           float &pymuon3,float &pzmuon3,float &pmuon3,float &dxymuon3,
                           float &pTmuon3tuneP,float &pTmuonBestTrack3);

   float Mass(float Pt1,float Eta1,float Phi1,float En1,
	      float Pt2,float Eta2,float Phi2,float En2);
   void PlotRecoInfo(float CosmicMuonRejec,float vertexMassMu,float MassGenerated,
		     float PtTunePMuBestTrack,float PtTunePMu,float PtMuBestTrack,
		     float PtGenerated,float etaMu1,
		     float PtTunePMuBestTrack2,float PtTunePMu2,float PtMuBestTrack2,
		     float PtGenerated2,float etaMu2);
   void PickThehighestMass(float &vtxHighestMass,float &vtxHighestChi2,int EvtNb);
   double ThreeDangle(float pxMu1,float pyMu1,float pzMu1,float pMu1,
                      float pxMu2,float pyMu2,float pzMu2,float pMu2);
   
   bool SelectFirstGenMu(float &ETMu1,float &PhiSCMu1,
                         float &EtaSCMu1,float &EnMu1,
			 int &IDele1,int &Statele1,
			 unsigned &GenFlag1);
   
   bool SelectSecondGenMu(unsigned GenFlag1,float ETMu1,float &ETMu2,float &PhiSCMu2,
			  float &EtaSCMu2,float &EnMu2,int &IDele2,int &Statele2);
   /*
     bool GenRecoMatchMu1(float RecoEta1,float RecoPhi1,
   float &ETMu1,float &PhiSCMu1,
   float &EtaSCMu1,float &EnMu1,
   unsigned &GenFlag1);
   */
   //   bool GenRecoMatchMu(float RecoEta1,float RecoPhi1);
			/*
			  bool GenRecoMatchMu2(unsigned GenFlag1,float RecoEta2,float RecoPhi2,
			  float &ETMu2,float &PhiSCMu2,
			  float &EtaSCMu2,float &EnMu2);
			*/
bool GenRecoMatchMu1(float RecoEta1,float RecoPhi1,
                        float &ETMu1,float &PhiSCMu1,
                        float &EtaSCMu1,float &EnMu1,
                        unsigned &GenFlag1);
   
   bool GenRecoMatchMu2(unsigned GenFlag1,float RecoEta2,float RecoPhi2,
                        float &ETMu2,float &PhiSCMu2,
                        float &EtaSCMu2,float &EnMu2);
   
   int GenToRecoMatchMu(float pTgen, float etagen, float phigen);

   void PlotGenInfo(float ZprimeGenMass,float EtaGenMu1,float EtaGenMu2,float PtGenMu1,
		    float PtGenMu2,float EnGenMu1,float EnGenMu2);
   
   bool RecoHLTMuonMatching(float RecoEta,float RecoPhi);
   //bool RecoHLTMuonMatching(float hltEta,float hltPhi,float RecoEta,float RecoPhi);
   bool isPassHLT();
			
   void MuonPassingID();
   void PlotPterror();
   void PlotNbTrackLayers();
   void PlotNBValidPixelHits();
   void PlotNbValidMuonHits();
   void PlotNbMatchedStations();
   void PlotTrackiso();
   void PlotAbsDxy();
   void PlotPtTuneP();
   void plotAllHighPtMuonsID();
   void MuonPassingNewID();
   void MuonPassingTightID();
   void CosThetaCollinSoper(float Et1,float Eta1,float Phi1,float En1,
			    float Et2,float Eta2,float Phi2,float En2,
			    float ChargeEle1,float RecoMass);
   float delR(float eta1,float phi1,float eta2,float phi2);
   //================================================================================
   float HLT_pt,HLT_eta,HLT_phi;
   float ptEffCut;
   float PtDYTRecMu1,PtDYTRecMu2,PtRecTunePMu1,PtRecTunePMu2,PtRecTunePMu3,
    PtRecMuBestTrack1,PtRecMuBestTrack2,PtRecMuBestTrack3;
   float RecoHLTMatchingDeltaRcut,deltaRcut,minMassCut,maxMassCut;
   float vtxChi2Mu,vtxMassMu;
   float mPtGen1,mPhiGen1,mEtaGen1,mEnGen1;
   unsigned mGenFlag1;
   float mPtGen2,mPhiGen2,mEtaGen2,mEnGen2;
   int ChargeRecMu1,ChargeRecMu2,ChargeRecMu3;
   unsigned flagmu1;
   unsigned flag1;
   float PtRecTunePMuBestTrack1,EnRecMu1,EtaRecMu1,PhiRecMu1;
   float PtRecTunePMuBestTrack2,EnRecMu2,EtaRecMu2,PhiRecMu2;
   float PtRecTunePMuBestTrack3,EnRecMu3,EtaRecMu3,PhiRecMu3;
   float pxRecMu1,pyRecMu1,pzRecMu1,pRecMu1,dxyRecMu1;
   float pxRecMu2,pyRecMu2,pzRecMu2,pRecMu2,dxyRecMu2;
   float pxRecMu3,pyRecMu3,pzRecMu3,pRecMu3,dxyRecMu3;
   float genET1,genPhi1,genEta1,genEn1;
   int genID1,genStat1;
   float genET2,genPhi2,genEta2,genEn2;
   int genID2,genStat2;
   float MassGen,RecoMass;
   int NbGen,NbReco;
   int nbTP,nbTT,nbTF;
   float TagProbeEtaCut;
   float Eff;
   float MassCutMin,MassCutMax;
   float MassResolution;
   float EtaCut;
   TH1F* h1_ZprimeRecomassBeforeTrigger_;
   TH1F* h1_ZprimeRecomass_;
   TH1F* h1_ZprimeRecomasslogscale_;                                                              
   TH1F* h1_ZprimeRecomass60to120_;
   TH1F* h1_ZprimeRecomassAbove400GeV_;
   TH1F* h1_ZprimeRecomassAbove1000GeV_;
   TH1F* h1_3Dangle_;
   TH1F* h1_DxyDiff_;
   TH1F* h1_ZprimeGenmass_;
   TH1F* h1_ZprimeGenEta1_;
   TH1F* h1_ZprimeGenEta2_;
   TH1F* h1_ZprimeGenPt1_;
   TH1F* h1_ZprimeGenPt2_;
   TH1F* h1_ZprimeGenEn1_;
   TH1F* h1_ZprimeGenEn2_;
   TH1F* h1_MassRecoGenDif_;
   TH1F* h1_dPToverPT_;
   TH1F* h1_normalizedChi2_;
   TH1F* h1_numberOftrackerLayersWithMeasurement_;
   TH1F* h1_numberOfValidPixelHits_;
   TH1F* h1_numberOfValidMuonHits_;
   TH1F* h1_numberOfMatchedStations_;
   TH1F* h1_trackiso_;
   TH1F* h1_absdxy_;
   TH1F* h1_PtEffpterror_;
   TH1F* h1_PtEffptnumberOftrackerLayers_;
   TH1F* h1_PtEffptnumberOfPixelHits_;
   TH1F* h1_PtEffptnumberOfMuonHits_;
   TH1F* h1_PtEffptnumberOfMatchedStations_;
   TH1F* h1_PtEffptTrackIso_;
   TH1F* h1_PtEffptabsdsy_;
   TH1F* h1_PtEffpfSumChargedHadron_;
   TH1F* h1_PtEffpfSumNeutralHadron_;
   TH1F* h1_PtEffpfPhotonIso_;
   TH1F* h1_FracpfSumChargedHadron_;
   TH1F* h1_FracpfSumNeutralHadron_;
   TH1F* h1_FracpfPhotonIso_;
   TH1F* h1_EtaEffpterror_;
   TH1F* h1_EtaEffptnumberOftrackerLayers_;
   TH1F* h1_EtaEffptnumberOfPixelHits_;
   TH1F* h1_EtaEffptnumberOfMuonHits_;
   TH1F* h1_EtaEffptnumberOfMatchedStations_;
   TH1F* h1_EtaEffptTrackIso_;
   TH1F* h1_EtaEffptabsdsy_;
   TH1F* h1_nbPVID_;
   TH1F* h1_PtID_;
   TH1F* h1_EtaID_;
   TH1F* h1_nbPVNewID_;
   TH1F* h1_PtNewID_;
   TH1F* h1_EtaNewID_;
   TH1F* h1_nbPVTightID_;
   TH1F* h1_PtTightID_;
   TH1F* h1_EtaTightID_;
   ofstream output_txt; 
   TH1F* h1_PtResolutionMBT_;
   TH1F* h1_PtResolutionTunePMBT_;
   TH1F* h1_PtResolutiontuneP_;
   TH1F* h1_MassResultionEBEB1_;
   TH1F* h1_MassResultionEBEB2_;
   TH1F* h1_MassResultionEBEB3_;
   TH1F* h1_MassResultionEBEB4_;
   TH1F* h1_MassResultionEBEB5_;
   TH1F* h1_MassResultionEBEB6_;
   TH1F* h1_MassResultionEBEB7_;
   TH1F* h1_MassResultionEBEB8_;
   TH1F* h1_MassResultionEBEB9_;
   TH1F* h1_MassResultionEBEB10_; 
   TH1F* h1_MassGenInAccep_;
   TH1F* h1_MassRecoInAccep_;
   TH1F* h1_CosAngleCollinSoperCorrect60Mass120_;
   TH1F* h1_CosAngleCollinSoperCorrect120Mass300_;
   TH1F* h1_CosAngleCollinSoperCorrect300Mass700_;
   TH1F* h1_CosAngleCollinSoperCorrect700Mass3000_;
   TH1F* h1_CosAngleCollinSoperCorrect4900Mass5100_;
   TH1F* h1_absCosAngleCollinSoperCorrect4500Mass5500_;
   TH1F* h1_ptHistoBefor_;
   TH1F* h1_ptHistoPassingVtxChi2Mu_;
   TH1F* h1_ptHistoPassingCosmicRejec_;
   TH1F* h1_ptHistoPassingHLT_;
   TH1F* h1_etaHistoBefor_;
   TH1F* h1_etaHistoPassingVtxChi2Mu_;
   TH1F* h1_etaHistoPassingCosmicRejec_;
   TH1F* h1_etaHistoPassingHLT_;
   TH1F* h1_3DangleHisto1_;
   TH1F* h1_3DangleHisto2_;
   TH1F* h1_Fail3DangleHistoMass_;
   TH1F* h1_Fail3DangleHistoPhi_;
   TH1F* h_ptHistoFRDum_;
   TH1F* h_ptHistoFRNum_;
   TH1F* h1_ZprimeRecomass50_;
   TH1F* h1_ZprimeRecomassBB_;
   TH1F* h1_ZprimeRecomassEE_;
   TH1F* h1_ZprimeRecomassBE_;
   TH1F* h1_ZprimeRecomass20_;      
   TH2F* h2_ZprimeRecomassNewbin_;
   TH1F *h1_ZprimeRecomassBinWidth_;
   int NbFireHLT;
   TH1F *h_NbFireHLT;
   TH1F* h1_ptAfterTrigger_;
   TH1F* h1_ptBeforeTrigger_;
   TH1F* h1_PtTuneP_;
};

#endif
#ifdef ZprimeMuMu_FR_MiniAod_cxx
ZprimeMuMu_FR_MiniAod::ZprimeMuMu_FR_MiniAod(Char_t namechar_[300], TTree *tree,Double_t weight_,string DATA_type_,string MC_type_) : fChain(0) 
{
 sprintf(name,"%s",namechar_);
  cout << "Name is= " << name << endl;
  weight = weight_;
  DATA_type = DATA_type_;
  MC_type = MC_type_;
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/lustre/cms/store/user/selgammal/ZprimeMuMu/MCs/CMSSW803_MC_TTbar_13TeV_pattuple.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/lustre/cms/store/user/selgammal/ZprimeMuMu/MCs/CMSSW803_MC_TTbar_13TeV_pattuple.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

ZprimeMuMu_FR_MiniAod::~ZprimeMuMu_FR_MiniAod()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ZprimeMuMu_FR_MiniAod::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ZprimeMuMu_FR_MiniAod::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ZprimeMuMu_FR_MiniAod::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   iGen = 0;
   idGen = 0;
   statusGen = 0;
   ptGen = 0;
   etaGen = 0;
   phiGen = 0;
   chargeGen = 0;
   EnergyGen = 0;
   pxGen = 0;
   pyGen = 0;
   pzGen = 0;
   Mu_nbMuon = 0;
   Mu_isGlobalMuon = 0;
   Mu_isPF = 0;
   Mu_isTrackerMuon = 0;
   Mu_et = 0;
   Mu_en = 0;
   Mu_pt = 0;
   Mu_eta = 0;
   Mu_phi = 0;
   Mu_charge = 0;
   Mu_ptTunePMuonBestTrack = 0;
   Mu_pxTunePMuonBestTrack = 0;
   Mu_pyTunePMuonBestTrack = 0;
   Mu_pzTunePMuonBestTrack = 0;
   Mu_pTunePMuonBestTrack = 0;
   Mu_etaTunePMuonBestTrack = 0;
   Mu_phiTunePMuonBestTrack = 0;
   Mu_thetaTunePMuonBestTrack = 0;
   Mu_chargeTunePMuonBestTrack = 0;
   Mu_dPToverPTTunePMuonBestTrack = 0;
   Mu_absdxyTunePMuonBestTrack = 0;
   Mu_absdzTunePMuonBestTrack = 0;
   Mu_ptInnerTrack = 0;
   Mu_normalizedChi2 = 0;
   Mu_vtxMass = 0;
   Mu_vtxNormChi2 = 0;
   Mu_vtxMassLept = 0;
   Mu_numberOfMatchedStations = 0;
   Mu_numberOfValidPixelHits = 0;
   Mu_numberOfValidMuonHits = 0;
   Mu_numberOftrackerLayersWithMeasurement = 0;
   Mu_emIso = 0;
   Mu_hadIso = 0;
   Mu_trackiso = 0;
   Mu_pfSumChargedHadronPt = 0;
   Mu_pfSumNeutralHadronEt = 0;
   Mu_PFSumPhotonEt = 0;
   Mu_pfSumPUPt = 0;
   Mu_nbofpv = 0;
   HLT_nb = 0;
   HLT_name = 0;
   HLT_isaccept = 0;
   HLTObj_nbObj = 0;
   HLTObj_pt = 0;
   HLTObj_eta = 0;
   HLTObj_phi = 0;
   HLTObj_collection = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event_runNo", &event_runNo, &b_event_runNo);
   fChain->SetBranchAddress("event_evtNo", &event_evtNo, &b_event_evtNo);
   fChain->SetBranchAddress("event_lumi", &event_lumi, &b_event_lumi);
   fChain->SetBranchAddress("event_bunch", &event_bunch, &b_event_bunch);
   fChain->SetBranchAddress("iGen", &iGen, &b_iGen);
   fChain->SetBranchAddress("idGen", &idGen, &b_idGen);
   fChain->SetBranchAddress("statusGen", &statusGen, &b_statusGen);
   fChain->SetBranchAddress("ptGen", &ptGen, &b_ptGen);
   fChain->SetBranchAddress("etaGen", &etaGen, &b_etaGen);
   fChain->SetBranchAddress("phiGen", &phiGen, &b_phiGen);
   fChain->SetBranchAddress("chargeGen", &chargeGen, &b_chargeGen);
   fChain->SetBranchAddress("EnergyGen", &EnergyGen, &b_EnergyGen);
   fChain->SetBranchAddress("pxGen", &pxGen, &b_pxGen);
   fChain->SetBranchAddress("pyGen", &pyGen, &b_pyGen);
   fChain->SetBranchAddress("pzGen", &pzGen, &b_pzGen);
   fChain->SetBranchAddress("Mu_nbMuon", &Mu_nbMuon, &b_Mu_nbMuon);
   fChain->SetBranchAddress("Mu_isGlobalMuon", &Mu_isGlobalMuon, &b_Mu_isGlobalMuon);
   fChain->SetBranchAddress("Mu_isPF", &Mu_isPF, &b_Mu_isPF);
   fChain->SetBranchAddress("Mu_isTrackerMuon", &Mu_isTrackerMuon, &b_Mu_isTrackerMuon);
   fChain->SetBranchAddress("Mu_et", &Mu_et, &b_Mu_et);
   fChain->SetBranchAddress("Mu_en", &Mu_en, &b_Mu_en);
   fChain->SetBranchAddress("Mu_pt", &Mu_pt, &b_Mu_pt);
   fChain->SetBranchAddress("Mu_eta", &Mu_eta, &b_Mu_eta);
   fChain->SetBranchAddress("Mu_phi", &Mu_phi, &b_Mu_phi);
   fChain->SetBranchAddress("Mu_charge", &Mu_charge, &b_Mu_charge);
   fChain->SetBranchAddress("Mu_ptTunePMuonBestTrack", &Mu_ptTunePMuonBestTrack, &b_Mu_ptTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_pxTunePMuonBestTrack", &Mu_pxTunePMuonBestTrack, &b_Mu_pxTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_pyTunePMuonBestTrack", &Mu_pyTunePMuonBestTrack, &b_Mu_pyTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_pzTunePMuonBestTrack", &Mu_pzTunePMuonBestTrack, &b_Mu_pzTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_pTunePMuonBestTrack", &Mu_pTunePMuonBestTrack, &b_Mu_pTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_etaTunePMuonBestTrack", &Mu_etaTunePMuonBestTrack, &b_Mu_etaTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_phiTunePMuonBestTrack", &Mu_phiTunePMuonBestTrack, &b_Mu_phiTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_thetaTunePMuonBestTrack", &Mu_thetaTunePMuonBestTrack, &b_Mu_thetaTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_chargeTunePMuonBestTrack", &Mu_chargeTunePMuonBestTrack, &b_Mu_chargeTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_dPToverPTTunePMuonBestTrack", &Mu_dPToverPTTunePMuonBestTrack, &b_Mu_dPToverPTTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_absdxyTunePMuonBestTrack", &Mu_absdxyTunePMuonBestTrack, &b_Mu_absdxyTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_absdzTunePMuonBestTrack", &Mu_absdzTunePMuonBestTrack, &b_Mu_absdzTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_ptInnerTrack", &Mu_ptInnerTrack, &b_Mu_ptInnerTrack);
   fChain->SetBranchAddress("Mu_normalizedChi2", &Mu_normalizedChi2, &b_Mu_normalizedChi2);
   fChain->SetBranchAddress("Mu_vtxMass", &Mu_vtxMass, &b_Mu_vtxMass);
   fChain->SetBranchAddress("Mu_vtxNormChi2", &Mu_vtxNormChi2, &b_Mu_vtxNormChi2);
   fChain->SetBranchAddress("Mu_vtxMassLept", &Mu_vtxMassLept, &b_Mu_vtxMassLept);
   fChain->SetBranchAddress("Mu_numberOfMatchedStations", &Mu_numberOfMatchedStations, &b_Mu_numberOfMatchedStations);
   fChain->SetBranchAddress("Mu_numberOfValidPixelHits", &Mu_numberOfValidPixelHits, &b_Mu_numberOfValidPixelHits);
   fChain->SetBranchAddress("Mu_numberOfValidMuonHits", &Mu_numberOfValidMuonHits, &b_Mu_numberOfValidMuonHits);
   fChain->SetBranchAddress("Mu_numberOftrackerLayersWithMeasurement", &Mu_numberOftrackerLayersWithMeasurement, &b_Mu_numberOftrackerLayersWithMeasurement);
   fChain->SetBranchAddress("Mu_emIso", &Mu_emIso, &b_Mu_emIso);
   fChain->SetBranchAddress("Mu_hadIso", &Mu_hadIso, &b_Mu_hadIso);
   fChain->SetBranchAddress("Mu_trackiso", &Mu_trackiso, &b_Mu_trackiso);
   fChain->SetBranchAddress("Mu_pfSumChargedHadronPt", &Mu_pfSumChargedHadronPt, &b_Mu_pfSumChargedHadronPt);
   fChain->SetBranchAddress("Mu_pfSumNeutralHadronEt", &Mu_pfSumNeutralHadronEt, &b_Mu_pfSumNeutralHadronEt);
   fChain->SetBranchAddress("Mu_PFSumPhotonEt", &Mu_PFSumPhotonEt, &b_Mu_PFSumPhotonEt);
   fChain->SetBranchAddress("Mu_pfSumPUPt", &Mu_pfSumPUPt, &b_Mu_pfSumPUPt);
   fChain->SetBranchAddress("Mu_nbofpv", &Mu_nbofpv, &b_Mu_nbofpv);
   fChain->SetBranchAddress("HLT_nb", &HLT_nb, &b_HLT_nb);
   fChain->SetBranchAddress("HLT_name", &HLT_name, &b_HLT_name);
   fChain->SetBranchAddress("HLT_isaccept", &HLT_isaccept, &b_HLT_isaccept);
   fChain->SetBranchAddress("HLTObj_nbObj", &HLTObj_nbObj, &b_HLTObj_nbObj);
   fChain->SetBranchAddress("HLTObj_pt", &HLTObj_pt, &b_HLTObj_pt);
   fChain->SetBranchAddress("HLTObj_eta", &HLTObj_eta, &b_HLTObj_eta);
   fChain->SetBranchAddress("HLTObj_phi", &HLTObj_phi, &b_HLTObj_phi);
   fChain->SetBranchAddress("HLTObj_collection", &HLTObj_collection, &b_HLTObj_collection);
   fChain->SetBranchAddress("GenMet_pt", &GenMet_pt, &b_GenMet_pt);
   fChain->SetBranchAddress("PFMet_pt_cor", &PFMet_pt_cor, &b_PFMet_pt_cor);
   fChain->SetBranchAddress("PFMet_phi_cor", &PFMet_phi_cor, &b_PFMet_phi_cor);
   fChain->SetBranchAddress("PFMet_en_cor", &PFMet_en_cor, &b_PFMet_en_cor);
   fChain->SetBranchAddress("PFMet_px_cor", &PFMet_px_cor, &b_PFMet_px_cor);
   fChain->SetBranchAddress("PFMet_py_cor", &PFMet_py_cor, &b_PFMet_py_cor);
   fChain->SetBranchAddress("PFMet_pz_cor", &PFMet_pz_cor, &b_PFMet_pz_cor);
   fChain->SetBranchAddress("PFMet_sumEt_cor", &PFMet_sumEt_cor, &b_PFMet_sumEt_cor);
   fChain->SetBranchAddress("CaloMet_pt", &CaloMet_pt, &b_CaloMet_pt);
   fChain->SetBranchAddress("CaloMet_phi", &CaloMet_phi, &b_CaloMet_phi);
   fChain->SetBranchAddress("CaloMet_sumEt", &CaloMet_sumEt, &b_CaloMet_sumEt);
   fChain->SetBranchAddress("METSign", &METSign, &b_METSign);
   fChain->SetBranchAddress("PFMet_shiftedPt_JetEnUp", &PFMet_shiftedPt_JetEnUp, &b_PFMet_shiftedPt_JetEnUp);
   fChain->SetBranchAddress("PFMet_shiftedPt_JetEnDown", &PFMet_shiftedPt_JetEnDown, &b_PFMet_shiftedPt_JetEnDown);
   Notify();
}

Bool_t ZprimeMuMu_FR_MiniAod::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ZprimeMuMu_FR_MiniAod::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ZprimeMuMu_FR_MiniAod::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ZprimeMuMu_FR_MiniAod_cxx
