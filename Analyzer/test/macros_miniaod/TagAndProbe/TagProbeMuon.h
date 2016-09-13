//=========================================================================
//          Analysis code for Z' boson to Mu Mu analysis                  =  
//                  Tag and Proble code for Muon ID                       = 
//                        To run over MC                                  =
//                        Counting method                                 =  
//                    Author: Sherif Elgammal                             =
//                         1/1/2015                                       =
//=========================================================================
#ifndef TagProbeEB_h
#define TagProbeEB_h
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
// Fixed size dimensions of array or collections stored in the TTree if any.

class TagProbeMuon {
public :

      TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // save a tree
   TTree          *tree_; // pointer to tree before cuts
   TTree          *tree_twomuons_;   // pointer to tree after 2mu cut
   TTree          *tree_vtxmass_;   // pointer to tree after vtxmass cut
   TTree          *tree_hlt_;   // pointer to tree after hlt cut
   TTree          *tree_vtxchi2mu_;  // pointer to tree after vtxchi2mu cut

   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           event_runNo;
   Int_t           event_evtNo;
   Int_t           event_lumi;
   Int_t           event_bunch;
   vector<int>     *iGenJet;
   vector<int>     *idGenJet;
   vector<int>     *statusGenJet;
   vector<int>     *chargeGenJet;
   vector<float>   *ptGenJet;
   vector<float>   *etaGenJet;
   vector<float>   *phiGenJet;
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
   vector<bool>    *Mu_isTightMuon;
   vector<bool>    *Mu_isLooseMuon;
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
   vector<float>   *Mu_pxInnerTrack;
   vector<float>   *Mu_pyInnerTrack;
   vector<float>   *Mu_pzInnerTrack;
   vector<float>   *Mu_pInnerTrack;
   vector<float>   *Mu_etaInnerTrack;
   vector<float>   *Mu_phiInnerTrack;
   vector<float>   *Mu_thetaInnerTrack;
   vector<float>   *Mu_chargeInnerTrack;
   vector<float>   *Mu_dPToverPTInnerTrack;
   vector<float>   *Mu_absdxyInnerTrack;
   vector<float>   *Mu_absdzInnerTrack;
   vector<float>   *Mu_normalizedChi2;
   vector<float>   *Mu_absdxy;
   vector<float>   *Mu_absdz;
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
   vector<float>   *Mu_patDeltaBeta;
   vector<int>     *HLT_nb;
   vector<string>  *HLT_name;
   vector<bool>    *HLT_isaccept;
   vector<int>     *HLTObj_nbObj;
   vector<float>   *HLTObj_pt;
   vector<float>   *HLTObj_eta;
   vector<float>   *HLTObj_phi;
   vector<string>  *HLTObj_collection;
   Double_t        GenMet_pt;
   Double_t        PFMet_et_cor;
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
   Double_t        PFMet_shiftedPt_JetEnUp;
   Double_t        PFMet_shiftedPt_JetEnDown;
   vector<int>     *jet_nb;
   vector<float>   *jet_charge;
   vector<float>   *jet_et;
   vector<float>   *jet_pt;
   vector<float>   *jet_eta;
   vector<float>   *jet_phi;
   vector<float>   *jet_en;
   vector<float>   *jet_theta;
   vector<float>   *jet_beta;
   vector<float>   *jet_pileup_mva_disc;
   vector<int>     *Nb_bDiscriminators;
   vector<float>   *jet_btag_pt;
   vector<float>   *jet_btag_eta;
   vector<float>   *jet_btag_phi;
   vector<int>     *jet_btag_flavor;
   vector<float>   *jet_btag_pfCSVv2IVF_discriminator;
   vector<int>     *Nb_taus;
   vector<float>   *Tau_pt;
   vector<float>   *Tau_eta;
   vector<float>   *Tau_phi;
   vector<int>     *Tau_id;
   vector<float>   *Tau_LooseCombinedIsolationDeltaBetaCorr3Hits;
   Int_t           pfphoton_size;
   vector<float>   *pfphoton_pt;
   vector<float>   *pfphoton_eta;
   vector<float>   *pfphoton_phi;
   vector<float>   *pfphoton_theta;
   Int_t           num_PU_vertices;
   Int_t           PU_BunchCrossing;
   Int_t           num_PU_gen_vertices;
   Float_t         Rho;
   vector<float>   *MC_weighting;

   // List of branches
   TBranch        *b_event_runNo;   //!
   TBranch        *b_event_evtNo;   //!
   TBranch        *b_event_lumi;   //!
   TBranch        *b_event_bunch;   //!
   TBranch        *b_iGenJet;   //!
   TBranch        *b_idGenJet;   //!
   TBranch        *b_statusGenJet;   //!
   TBranch        *b_chargeGenJet;   //!
   TBranch        *b_ptGenJet;   //!
   TBranch        *b_etaGenJet;   //!
   TBranch        *b_phiGenJet;   //!
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
   TBranch        *b_Mu_isTightMuon;   //!
   TBranch        *b_Mu_isLooseMuon;   //!
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
   TBranch        *b_Mu_pxInnerTrack;   //!
   TBranch        *b_Mu_pyInnerTrack;   //!
   TBranch        *b_Mu_pzInnerTrack;   //!
   TBranch        *b_Mu_pInnerTrack;   //!
   TBranch        *b_Mu_etaInnerTrack;   //!
   TBranch        *b_Mu_phiInnerTrack;   //!
   TBranch        *b_Mu_thetaInnerTrack;   //!
   TBranch        *b_Mu_chargeInnerTrack;   //!
   TBranch        *b_Mu_dPToverPTInnerTrack;   //!
   TBranch        *b_Mu_absdxyInnerTrack;   //!
   TBranch        *b_Mu_absdzInnerTrack;   //!
   TBranch        *b_Mu_normalizedChi2;   //!
   TBranch        *b_Mu_absdxy;   //!
   TBranch        *b_Mu_absdz;   //!
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
   TBranch        *b_Mu_patDeltaBeta;   //!
   TBranch        *b_HLT_nb;   //!
   TBranch        *b_HLT_name;   //!
   TBranch        *b_HLT_isaccept;   //!
   TBranch        *b_HLTObj_nbObj;   //!
   TBranch        *b_HLTObj_pt;   //!
   TBranch        *b_HLTObj_eta;   //!
   TBranch        *b_HLTObj_phi;   //!
   TBranch        *b_HLTObj_collection;   //!
   TBranch        *b_GenMet_pt;   //!
   TBranch        *b_PFMet_et_cor;   //!
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
   TBranch        *b_PFMet_shiftedPt_JetEnUp;   //!
   TBranch        *b_PFMet_shiftedPt_JetEnDown;   //!
   TBranch        *b_jet_nb;   //!
   TBranch        *b_jet_charge;   //!
   TBranch        *b_jet_et;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_en;   //!
   TBranch        *b_jet_theta;   //!
   TBranch        *b_jet_beta;   //!
   TBranch        *b_jet_pileup_mva_disc;   //!
   TBranch        *b_Nb_bDiscriminators;   //!
   TBranch        *b_jet_btag_pt;   //!
   TBranch        *b_jet_btag_eta;   //!
   TBranch        *b_jet_btag_phi;   //!
   TBranch        *b_jet_btag_flavor;   //!
   TBranch        *b_jet_btag_pfCSVv2IVF_discriminator;   //!
   TBranch        *b_Nb_taus;   //!
   TBranch        *b_Tau_pt;   //!
   TBranch        *b_Tau_eta;   //!
   TBranch        *b_Tau_phi;   //!
   TBranch        *b_Tau_id;   //!
   TBranch        *b_Tau_LooseCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_pfphoton_size;   //!
   TBranch        *b_pfphoton_pt;   //!
   TBranch        *b_pfphoton_eta;   //!
   TBranch        *b_pfphoton_phi;   //!
   TBranch        *b_pfphoton_theta;   //!
   TBranch        *b_num_PU_vertices;   //!
   TBranch        *b_PU_BunchCrossing;   //!
   TBranch        *b_num_PU_gen_vertices;   //!
   TBranch        *b_Rho;   //!
   TBranch        *b_MC_weighting;   //!


   TagProbeMuon(TTree *tree=0);
   virtual ~TagProbeMuon();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   
   //================================================================================== 
   //                                                                                 =
   //                                    Start methods                                =             
   //                                                                                 =
   //==================================================================================
   //========================== Method to select first HEEP Ele =======================
   float Mass(float Et1,float Eta1,float Phi1,float En1,
	      float Et2,float Eta2,float Phi2,float En2);
   void SelectTagandProbePassingProbe();
   void SelectTagandProbePassingTag();
   void TagandProbeResults();
   bool SelectFirstGenMu(float &ETMu1,float &PhiSCMu1,
			 float &EtaSCMu1,float &EnMu1,
			 int &IDMu1,int &StatMu1,
			 unsigned &GenFlag1);
   bool SelectSecondGenMu(unsigned GenFlag1,float ETMu1,float &ETMu2,float &PhiSCMu2,
			  float &EtaSCMu2,float &EnMu2,int &IDMu2,int &StatMu2);
   bool GenRecoMatchMuon(float RecoEta1,float RecoPhi1);
   bool RecoHLTMuonMatching(float RecoEta,float RecoPhi);
   float delR(float eta1,float phi1,float eta2,float phi2);
   //============== start Histos ==============
   int nbTP1,
     nbTP2,
     nbTP3,
     nbTP4,
     nbTP5,
     nbTP6,
     nbTP7,
     nbTP8,
     nbTP9,  
     nbTP10,
     nbTP11,
     nbTP12,
     nbTP13,
     nbTP14,
     nbTP15,
     nbTP16,
     nbTP17;

   int nbTT1,
     nbTT2,
     nbTT3,
     nbTT4,
     nbTT5,
     nbTT6,
     nbTT7,
     nbTT8,
     nbTT9,  
     nbTT10,
     nbTT11,
     nbTT12,
     nbTT13,
     nbTT14,
     nbTT15,
     nbTT16,
     nbTT17;

   int NewnbTP,NewnbTF;
   float ptMin,ptMax;
   float binMass,MassCutMin,MassCutMax,TagProbeEtCut,RecoGenMatchingDeltaRcut;
   float ET1               ;
   float PhiSC1            ;
   float EtaSC1            ;
   int EleNb1              ;
   //float En1               ;
   float distEle1          ;
   float dcotEle1          ;
   float radiusEle1        ;
   float TrackIsoEle1      ;
   float deltaEtaEle1      ;
   float deltaPhiEle1      ;
   float hoeEle1           ;
   float ecalIsoEle1       ;
   float hcalIsoEle1       ;
   float sigmaIetaIetaEle1 ;
   float ChargeEle1        ;
   int eleClassEle1        ;
   float hcalDepth2IsoEle1 ;
   float thetaEleSC1       ;
   float etaTrakEle1       ;
   float phiTrakEle1       ;
   float thetaTrackEle1    ;
   int eelectronFromConv1  ;
   float NewHoeBcEle1      ;
   float NewHcal1BcIsoEle1 ;
   float NewEcalPlusHcal1BcIsoEle1;
   float eLxyFromConvEle1;
   float eProbFromConvEle1;
   int eNHitsBeforeVtxMaxFromConvEle1;
   //================================
   float ET2               ;
   float PhiSC2            ;
   float EtaSC2            ;
   int EleNb2              ;
   //float En2               ;
   float distEle2          ;
   float dcotEle2          ;
   float radiusEle2        ;
   float TrackIsoEle2      ;
   float deltaEtaEle2      ;
   float deltaPhiEle2      ;
   float hoeEle2           ;
   float ecalIsoEle2       ;
   float hcalIsoEle2       ;
   float sigmaIetaIetaEle2 ;
   float ChargeEle2        ;
   int eleClassEle2        ;
   float hcalDepth2IsoEle2 ;
   float thetaEleSC2       ;
   float etaTrakEle2       ;
   float phiTrakEle2       ;
   float thetaTrackEle2    ;
   int eelectronFromConv2  ;
   float NewHoeBcEle2      ;
   float NewHcal1BcIsoEle2 ;
   float NewEcalPlusHcal1BcIsoEle2;
   float eLxyMinFromConvEle2;
   //================================
   float HeepMassPass          ;
   float HeepMassFail          ;
   float ptjet1            ;
   float etajet1           ;
   float phijet1           ;
   float enjet1            ;
   float ptjet2            ;
   float etajet2           ;
   float phijet2           ;
   float enjet2            ;
   int nbofPV              ;
   int nbdof               ;
   float rhoPV,xPV,yPV,zPV ;
   float etTag1,etaTag1,phiTag1,enTag1;
   float etProbe1,etaProbe1,phiProbe1,enProbe1,deltaEtaProbe1,deltaPhiProbe1,hoeProbe1;
   float etTag2,etaTag2,phiTag2,enTag2;
   float etProbe2,etaProbe2,phiProbe2,enProbe2,deltaEtaProbe2,deltaPhiProbe2,hoeProbe2;
   float etTagfail1,etaTagfail1,phiTagfail1,enTagfail1;
   float etEle3,phiTrEle3,etaTrEle3,enEle3;
   float etEle4,phiTrEle4,etaTrEle4,enEle4;
   float etEle5,phiTrEle5,etaTrEle5,enEle5;
   float MassMin,MassMax;
   int nbTP,nbTT,nbTF,nbPP;
   float Eff,error;
   float genET1,genPhi1,genEta1,genEn1;
   int genID1,genStat1;
   float genET2,genPhi2,genEta2,genEn2;
   int genID2,genStat2;
   unsigned flagmu1;
   unsigned flag1;
   float MassGen;
   float etaCut,TagPtCut;
   float RecoHLTMatchingDeltaRcut;
   int binEta;
   float EtaMin;
   float EtaMax;
   //-----------------------------------------------
   
   /** Write the histograms for electrons in general */
   TH1F* h1_EtTagTag_;
   TH1F* h1_EtTagProble_;
   TH1F* h1_EtTagProbleFail_;

   TH1F* h1_EtTag1_;   
   TH1F* h1_EtTag2_;
   TH1F* h1_MassPass_;
   TH1F* h1_MassFail_;

   TH1F* h1_HeepMassTP_;
   TH1F* h1_HeepMassTT_;
   TH1F* h1_HeepMassTF_;
   
   TH1F* h1_NbHEEPele_;
   TH1F* h1_NbHEEPele2_;
   
   TH1F* h1_NbeleFail_;

   TH1F* h1_EtProbeFail_;
   TH1F* h1_EtaProbeFail_;
   
   
   TH1F* h1_ETEleTag1_;
   TH1F* h1_ETEleTag2_;

   TH1F* h1_deltaEtaFail_;
   TH1F* h1_deltaPhiFail_;
   TH1F* h1_hoeFail_;


   TH1F* h1_EtProbePass_;
   TH1F* h1_EtaProbePass_;
   TH1F* h1_deltaEtaProbePass_;
   TH1F* h1_deltaPhiProbePass_;
   TH1F* h1_hoeProbePass_;


   TH1F* h1_invMassTP_;
   TH1F* h1_invMassTF_;
   TH1F* h1_invMassTT_;

   TH1F* h1_PtProbeTT_;
   TH1F* h1_PtProbeTP_;
   TH1F* h1_PtProbeTF_;
   TH1F* h1_EtaProbeTT_;
   TH1F* h1_EtaProbeTP_;
   TH1F* h1_EtaProbeTF_;



   ///fake rate hiso
   TH1F* h1_EtFakeRateEB1_;
   TH1F* h1_EtaFakeRateEB1_;
   TH1F* h1_NbPvFakeRateEB1_;
   TH1F* h1_HoeFakeRateEB1_;
   TH1F* h1_SigmaFakeRateEB1_;
   TH1F* h1_NbhitsFakeRateEB1_;
   TH1F* h1_EtFakeRateEE1_;
   TH1F* h1_EtaFakeRateEE1_;
   TH1F* h1_NbPvFakeRateEE1_;
   TH1F* h1_HoeFakeRateEE1_;
   TH1F* h1_SigmaFakeRateEE1_;
   TH1F* h1_NbhitsFakeRateEE1_;
   TH1F* h1_HeepEtEB_;
   TH1F* h1_ZpeakMassTP_;
   TH1F* h1_ZpeakMassTF_;
   TH1F* h1_ZpeakMassTT_;
   TH1F* h1_PhiProbeTT_;
   TH1F* h1_PhiProbeTP_;
   TH1F* h1_PhiProbeTF_;
   TH1F* h1_PVProbeTT_;
   TH1F* h1_PVProbeTP_;
   TH1F* h1_PVProbeTF_;
   TH1F* h1_NBofTP_;
   TH1F* h1_NBofTF_;
   TH1F* h1_NBofTT_;

   TH1F* h1_invMassPP_;
   TH1F* h1_ZpeakMassPP_;
   TH1F* h1_PtProbePP_;
   TH1F* h1_EtaProbePP_;
   TH1F* h1_PhiProbePP_;
   TH1F* h1_PVProbePP_;
   TH1F* h1_NBofPP_;

};

#endif

#ifdef TagProbeMuon_cxx
TagProbeMuon::TagProbeMuon(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/lustre/cms/store/user/defilip/ZprimeAnalysis/Spring16_ZprimeMuMu_13TeV_merged_HLT/CMSSW8012_MC_reHLT_DYtoMuMu120to200_13TeV_tree.root");
     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/lustre/cms/store/user/defilip/ZprimeAnalysis/Data2016_ZprimeMuMu_13TeV_merged_HLT/CMSSW_8_0_13_ZPRIMEMuMu_13TeV-DataG-V1-JSON-tree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("CMSSW8012_MC_reHLT_DYtoMuMu120to200_13TeV_tree.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

TagProbeMuon::~TagProbeMuon()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TagProbeMuon::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TagProbeMuon::LoadTree(Long64_t entry)
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

void TagProbeMuon::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   iGenJet = 0;
   idGenJet = 0;
   statusGenJet = 0;
   chargeGenJet = 0;
   ptGenJet = 0;
   etaGenJet = 0;
   phiGenJet = 0;
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
   Mu_isTightMuon = 0;
   Mu_isLooseMuon = 0;
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
   Mu_pxInnerTrack = 0;
   Mu_pyInnerTrack = 0;
   Mu_pzInnerTrack = 0;
   Mu_pInnerTrack = 0;
   Mu_etaInnerTrack = 0;
   Mu_phiInnerTrack = 0;
   Mu_thetaInnerTrack = 0;
   Mu_chargeInnerTrack = 0;
   Mu_dPToverPTInnerTrack = 0;
   Mu_absdxyInnerTrack = 0;
   Mu_absdzInnerTrack = 0;
   Mu_normalizedChi2 = 0;
   Mu_absdxy = 0;
   Mu_absdz = 0;
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
   Mu_patDeltaBeta = 0;
   HLT_nb = 0;
   HLT_name = 0;
   HLT_isaccept = 0;
   HLTObj_nbObj = 0;
   HLTObj_pt = 0;
   HLTObj_eta = 0;
   HLTObj_phi = 0;
   HLTObj_collection = 0;
   jet_nb = 0;
   jet_charge = 0;
   jet_et = 0;
   jet_pt = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_en = 0;
   jet_theta = 0;
   jet_beta = 0;
   jet_pileup_mva_disc = 0;
   Nb_bDiscriminators = 0;
   jet_btag_pt = 0;
   jet_btag_eta = 0;
   jet_btag_phi = 0;
   jet_btag_flavor = 0;
   jet_btag_pfCSVv2IVF_discriminator = 0;
   Nb_taus = 0;
   Tau_pt = 0;
   Tau_eta = 0;
   Tau_phi = 0;
   Tau_id = 0;
   Tau_LooseCombinedIsolationDeltaBetaCorr3Hits = 0;
   pfphoton_pt = 0;
   pfphoton_eta = 0;
   pfphoton_phi = 0;
   pfphoton_theta = 0;
   MC_weighting = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event_runNo", &event_runNo, &b_event_runNo);
   fChain->SetBranchAddress("event_evtNo", &event_evtNo, &b_event_evtNo);
   fChain->SetBranchAddress("event_lumi", &event_lumi, &b_event_lumi);
   fChain->SetBranchAddress("event_bunch", &event_bunch, &b_event_bunch);
   fChain->SetBranchAddress("iGenJet", &iGenJet, &b_iGenJet);
   fChain->SetBranchAddress("idGenJet", &idGenJet, &b_idGenJet);
   fChain->SetBranchAddress("statusGenJet", &statusGenJet, &b_statusGenJet);
   fChain->SetBranchAddress("chargeGenJet", &chargeGenJet, &b_chargeGenJet);
   fChain->SetBranchAddress("ptGenJet", &ptGenJet, &b_ptGenJet);
   fChain->SetBranchAddress("etaGenJet", &etaGenJet, &b_etaGenJet);
   fChain->SetBranchAddress("phiGenJet", &phiGenJet, &b_phiGenJet);
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
   fChain->SetBranchAddress("Mu_isTightMuon", &Mu_isTightMuon, &b_Mu_isTightMuon);
   fChain->SetBranchAddress("Mu_isLooseMuon", &Mu_isLooseMuon, &b_Mu_isLooseMuon);
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
   fChain->SetBranchAddress("Mu_pxInnerTrack", &Mu_pxInnerTrack, &b_Mu_pxInnerTrack);
   fChain->SetBranchAddress("Mu_pyInnerTrack", &Mu_pyInnerTrack, &b_Mu_pyInnerTrack);
   fChain->SetBranchAddress("Mu_pzInnerTrack", &Mu_pzInnerTrack, &b_Mu_pzInnerTrack);
   fChain->SetBranchAddress("Mu_pInnerTrack", &Mu_pInnerTrack, &b_Mu_pInnerTrack);
   fChain->SetBranchAddress("Mu_etaInnerTrack", &Mu_etaInnerTrack, &b_Mu_etaInnerTrack);
   fChain->SetBranchAddress("Mu_phiInnerTrack", &Mu_phiInnerTrack, &b_Mu_phiInnerTrack);
   fChain->SetBranchAddress("Mu_thetaInnerTrack", &Mu_thetaInnerTrack, &b_Mu_thetaInnerTrack);
   fChain->SetBranchAddress("Mu_chargeInnerTrack", &Mu_chargeInnerTrack, &b_Mu_chargeInnerTrack);
   fChain->SetBranchAddress("Mu_dPToverPTInnerTrack", &Mu_dPToverPTInnerTrack, &b_Mu_dPToverPTInnerTrack);
   fChain->SetBranchAddress("Mu_absdxyInnerTrack", &Mu_absdxyInnerTrack, &b_Mu_absdxyInnerTrack);
   fChain->SetBranchAddress("Mu_absdzInnerTrack", &Mu_absdzInnerTrack, &b_Mu_absdzInnerTrack);
   fChain->SetBranchAddress("Mu_normalizedChi2", &Mu_normalizedChi2, &b_Mu_normalizedChi2);
   fChain->SetBranchAddress("Mu_absdxy", &Mu_absdxy, &b_Mu_absdxy);
   fChain->SetBranchAddress("Mu_absdz", &Mu_absdz, &b_Mu_absdz);
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
   fChain->SetBranchAddress("Mu_patDeltaBeta", &Mu_patDeltaBeta, &b_Mu_patDeltaBeta);
   fChain->SetBranchAddress("HLT_nb", &HLT_nb, &b_HLT_nb);
   fChain->SetBranchAddress("HLT_name", &HLT_name, &b_HLT_name);
   fChain->SetBranchAddress("HLT_isaccept", &HLT_isaccept, &b_HLT_isaccept);
   fChain->SetBranchAddress("HLTObj_nbObj", &HLTObj_nbObj, &b_HLTObj_nbObj);
   fChain->SetBranchAddress("HLTObj_pt", &HLTObj_pt, &b_HLTObj_pt);
   fChain->SetBranchAddress("HLTObj_eta", &HLTObj_eta, &b_HLTObj_eta);
   fChain->SetBranchAddress("HLTObj_phi", &HLTObj_phi, &b_HLTObj_phi);
   fChain->SetBranchAddress("HLTObj_collection", &HLTObj_collection, &b_HLTObj_collection);
   fChain->SetBranchAddress("GenMet_pt", &GenMet_pt, &b_GenMet_pt);
   fChain->SetBranchAddress("PFMet_et_cor", &PFMet_et_cor, &b_PFMet_et_cor);
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
   fChain->SetBranchAddress("PFMet_shiftedPt_JetEnUp", &PFMet_shiftedPt_JetEnUp, &b_PFMet_shiftedPt_JetEnUp);
   fChain->SetBranchAddress("PFMet_shiftedPt_JetEnDown", &PFMet_shiftedPt_JetEnDown, &b_PFMet_shiftedPt_JetEnDown);
   fChain->SetBranchAddress("jet_nb", &jet_nb, &b_jet_nb);
   fChain->SetBranchAddress("jet_charge", &jet_charge, &b_jet_charge);
   fChain->SetBranchAddress("jet_et", &jet_et, &b_jet_et);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_en", &jet_en, &b_jet_en);
   fChain->SetBranchAddress("jet_theta", &jet_theta, &b_jet_theta);
   fChain->SetBranchAddress("jet_beta", &jet_beta, &b_jet_beta);
   fChain->SetBranchAddress("jet_pileup_mva_disc", &jet_pileup_mva_disc, &b_jet_pileup_mva_disc);
   fChain->SetBranchAddress("Nb_bDiscriminators", &Nb_bDiscriminators, &b_Nb_bDiscriminators);
   fChain->SetBranchAddress("jet_btag_pt", &jet_btag_pt, &b_jet_btag_pt);
   fChain->SetBranchAddress("jet_btag_eta", &jet_btag_eta, &b_jet_btag_eta);
   fChain->SetBranchAddress("jet_btag_phi", &jet_btag_phi, &b_jet_btag_phi);
   fChain->SetBranchAddress("jet_btag_flavor", &jet_btag_flavor, &b_jet_btag_flavor);
   fChain->SetBranchAddress("jet_btag_pfCSVv2IVF_discriminator", &jet_btag_pfCSVv2IVF_discriminator, &b_jet_btag_pfCSVv2IVF_discriminator);
   fChain->SetBranchAddress("Nb_taus", &Nb_taus, &b_Nb_taus);
   fChain->SetBranchAddress("Tau_pt", &Tau_pt, &b_Tau_pt);
   fChain->SetBranchAddress("Tau_eta", &Tau_eta, &b_Tau_eta);
   fChain->SetBranchAddress("Tau_phi", &Tau_phi, &b_Tau_phi);
   fChain->SetBranchAddress("Tau_id", &Tau_id, &b_Tau_id);
   fChain->SetBranchAddress("Tau_LooseCombinedIsolationDeltaBetaCorr3Hits", &Tau_LooseCombinedIsolationDeltaBetaCorr3Hits, &b_Tau_LooseCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("pfphoton_size", &pfphoton_size, &b_pfphoton_size);
   fChain->SetBranchAddress("pfphoton_pt", &pfphoton_pt, &b_pfphoton_pt);
   fChain->SetBranchAddress("pfphoton_eta", &pfphoton_eta, &b_pfphoton_eta);
   fChain->SetBranchAddress("pfphoton_phi", &pfphoton_phi, &b_pfphoton_phi);
   fChain->SetBranchAddress("pfphoton_theta", &pfphoton_theta, &b_pfphoton_theta);
   fChain->SetBranchAddress("num_PU_vertices", &num_PU_vertices, &b_num_PU_vertices);
   fChain->SetBranchAddress("PU_BunchCrossing", &PU_BunchCrossing, &b_PU_BunchCrossing);
   fChain->SetBranchAddress("num_PU_gen_vertices", &num_PU_gen_vertices, &b_num_PU_gen_vertices);
   fChain->SetBranchAddress("Rho", &Rho, &b_Rho);
   fChain->SetBranchAddress("MC_weighting", &MC_weighting, &b_MC_weighting);
   Notify();
}

Bool_t TagProbeMuon::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TagProbeMuon::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TagProbeMuon::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TagProbeMuon_cxx
