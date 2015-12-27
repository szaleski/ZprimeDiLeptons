//=========================================================================  
//      Make the root tree for Z' boson to Mu Mu analysis                 =  
//                                                                        =
//                         CMSSW_7_4_5                                    =
//                                                                        =
//       Written by Sherif Elgammal & Nicola De Filippis                  =
//                                                                        =
//                         25/12/2015                                     =
//=========================================================================
#ifndef Analyzer_MaketreeMuons_h
#define Analyzer_MaketreeMuons_h
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include <FWCore/Framework/interface/ESHandle.h>
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
//====================== Get PixelMatchGsfElectronCollection =========
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
//====================== Get PixelMatchGsfElectronCollection =========
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"
#include "CLHEP/Geometry/Transform3D.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
//#include "DQMOffline/JetMET/interface/CaloMETAnalyzer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/METReco/interface/SpecificCaloMETData.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/CorrMETData.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
//==================== start a part for Jet ================================
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
//==================== end a part for Jet ================================
//==================== start a part for photons ==========================
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
//#include "tmp/CaloGeometryBuilder/interface/CaloGeometryBuilder.h"
#include <DataFormats/MuonReco/interface/Muon.h>
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Common/interface/RefCore.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/TypeID.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
 
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "CommonTools/TriggerUtils/interface/GenericTriggerEventFlag.h"
 
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h" 
#include "DataFormats/MuonReco/interface/MuonEnergy.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

// Trigger
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerFilter.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "CommonTools/ParticleFlow/interface/PFMETAlgo.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
//====================== end a part for photons ==========================
#include <vector>
#include <set>
#include <stdio.h>
#include "TProfile.h"
#include "TProfile2D.h"
//#include <iostream.h>
//#include <ostream.h>
//#include <fstream.h>
//#include <vector.h>
#include "TFile.h"
#include <math.h>
#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "Math/LorentzVector.h"
#include "TTree.h"
#include "TMath.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
class TFile;

//
// class declaration
//


class MaketreeMuons : public edm::EDAnalyzer {
   public:
      explicit MaketreeMuons( const edm::ParameterSet& );  //constructor
      ~MaketreeMuons();  //distructor
      virtual void analyze( const edm::Event&, const edm::EventSetup& );
      virtual void beginJob();
      virtual void endJob();
      bool PrimaryVertex(const reco::VertexCollection &vertices);
      void IntialValues();
      void MuonTree(const edm::Event& evt,const edm::EventSetup& es,
		    const reco::MuonCollection* muons,const reco::Vertex &vertex,
		    const reco::VertexCollection &verticess);
      void GenParticleTree(const edm::Event& evt);
      void PrimaryVertexTree(const reco::VertexCollection &vertices);
      void JetsTree(const reco::PFJetCollection &pfjets);
      bool IsMuMatchedToHLTMu ( const reco::Muon &mu, std::vector<reco::Particle> HLTMu , std::vector<std::string> HLTMuNames, double DR, double DPtRel);
      void TriggerMatchingTree(const edm::Event& iEvent,const reco::MuonCollection* muons,const reco::Vertex &vertex);
      float delR(float eta1,float phi1,float eta2,float phi2);
      bool GenParticleRecoMatching(const edm::Event& evt,float eta2,float phi2);
      void ComputeMuonMassVtx(const edm::Event& evt,const edm::EventSetup& es,
			      const reco::MuonCollection* muons,const reco::Vertex &vertex);
      void pfMETtree(const edm::Event& evt);
      void fillPU(const edm::Event& iEvent);
      void fillRho(const edm::Event& evt);
      void EventsReWeighting(const edm::Event& evt);
 private:

      TTree* mytree;
      //====================================================
      //
      //    Create variable for Nb of event,run,lumi
      //
      //====================================================
      TH1F* h1_DaughterID_;

      GlobalPoint eventVertexPosition_;
      int Run;
      int Event;
      int lumi;
      int bunch;
      int event_runNo;
      int event_evtNo;
      int event_lumi;
      int event_bunch;
      //==================================================
      //
      //           Create vectors for Muons variables
      //
      //==================================================
      int TotalNbEvents;
      int NbEventsPass;
      int nbTk;
      std::vector<int> NbEventsPassTrigger; 
      std::vector<int> NbEventsPassTriggerandPVcond;
      int NbMuons;
      std::vector<int> Mu_nbMuon;
      std::vector<float> EtaMuon;
      std::vector<bool> Mu_isGlobalMuon;
      std::vector<bool> Mu_isPF;
      std::vector<bool> Mu_isGoodMuon;
      std::vector<bool> Mu_isTrackerMuon;      
      std::vector<float> Mu_etaCocktail;
      std::vector<float> Mu_ptCocktail;
      std::vector<float> Mu_et;
      std::vector<float> Mu_phiCocktail;	
      std::vector<float> Mu_en;
      std::vector<float> Mu_thetaCocktail;
      std::vector<float> Mu_normalizedChi2;
      std::vector<float> Mu_absdxy;
      std::vector<float> Mu_absdz;
      std::vector<float> Mu_trackiso;
      std::vector<int> Mu_numberOfMatchedStations;
      std::vector<int> Mu_numberOfValidPixelHits;
      std::vector<int> Mu_numberOftrackerLayersWithMeasurement; 
      std::vector<int> Mu_innerTK_numberOfValidPixelHits;
      std::vector<int> Mu_innerTK_numberOfValidMuonHits;
      std::vector<float> Mu_pfSumChargedHadronPt;
      std::vector<float> Mu_pfSumNeutralHadronEt;
      std::vector<float> Mu_PFSumPhotonEt;
      std::vector<float> Mu_pfSumPUPt;
      std::vector<float> Mu_calEnergy;
      std::vector<int> Mu_numberOfValidMuonHits;
      std::vector<int> Mu_chargeCocktail;
      std::vector<float> Mu_emIso;
      std::vector<float> Mu_hadIso;	
      std::vector<float> Mu_VTXnormalizedChi2;
      std::vector<float> Mu_vtxMass;
      std::vector<float> Mu_vtxMassLept;
      std::vector<float> Mu_vtxTotalChi2;
      std::vector<float> Mu_vtxNormChi2;
      std::vector<float> Mu_ptcocktail;
      std::vector<float> Mu_dPToverPTcocktail;
      std::vector<int> Mu_nbofpv;
      std::vector<float> Mu_pxCocktail;
      std::vector<float> Mu_pyCocktail;
      std::vector<float> Mu_pzCocktail;
      std::vector<float> Mu_pCocktail;
      std::vector<float> Mu_ptInnerTrack;
      std::vector<float> Mu_ptTunePMuonBestTrack;
      std::vector<float> Mu_dPToverPTTunePMuonBestTrack;
      std::vector<float> Mu_pxTunePMuonBestTrack;
      std::vector<float> Mu_pyTunePMuonBestTrack;
      std::vector<float> Mu_pzTunePMuonBestTrack;
      std::vector<float> Mu_pTunePMuonBestTrack; 
      std::vector<float> Mu_etaTunePMuonBestTrack;
      std::vector<float> Mu_phiTunePMuonBestTrack;
      std::vector<float> Mu_thetaTunePMuonBestTrack;
      std::vector<float> Mu_chargeTunePMuonBestTrack;
      std::vector<float> Mu_absdxyTunePMuonBestTrack;
      std::vector<float> Mu_ptDYTTrack;
      std::vector<float> Mu_pxDYTTrack;
      std::vector<float> Mu_pyDYTTrack;
      std::vector<float> Mu_pzDYTTrack;
      std::vector<float> Mu_pDYTTrack; 
      std::vector<float> Mu_etaDYTTrack;
      std::vector<float> Mu_phiDYTTrack;
      std::vector<float> Mu_thetaDYTTrack;
      std::vector<float> Mu_chargeDYTTrack;    
      std::vector<float> Mu_absdxyDYTTrack;
      std::vector<float> Mu_dPToverPTDYTTrack;
      std::vector<float> Mu_ptPickyTrack;
      std::vector<float> Mu_pxPickyTrack;
      std::vector<float> Mu_pyPickyTrack;
      std::vector<float> Mu_pzPickyTrack;
      std::vector<float> Mu_pPickyTrack; 
      std::vector<float> Mu_etaPickyTrack;
      std::vector<float> Mu_phiPickyTrack;
      std::vector<float> Mu_thetaPickyTrack;
      std::vector<float> Mu_chargePickyTrack;
      std::vector<float> Mu_absdxyPickyTrack;
      std::vector<float> Mu_dPToverPTPickyTrack;
      std::vector<float> Mu_ptMuonBestTrack;
      std::vector<float> Mu_dPToverPTMuonBestTrack;
      std::vector<float> Mu_pxMuonBestTrack;
      std::vector<float> Mu_pyMuonBestTrack;
      std::vector<float> Mu_pzMuonBestTrack;
      std::vector<float> Mu_pMuonBestTrack; 
      std::vector<float> Mu_etaMuonBestTrack;
      std::vector<float> Mu_phiMuonBestTrack;
      std::vector<float> Mu_thetaMuonBestTrack;
      std::vector<float> Mu_chargeMuonBestTrack;
      std::vector<float> Mu_absdxyMuonBestTrack;
      //===================================================  
      //
      //    Create vectors for gen particles variables
      //
      //===================================================
      int value2_;
      std::vector<int> iGen;
      std::vector<int> idGen;
      std::vector<int> statusGen; 
      std::vector<float> ptGen;
      std::vector<float> etaGen;
      std::vector<float> phiGen;
      std::vector<float> massGen;
      std::vector<int> chargeGen;
      std::vector<float> EnergyGen;
      std::vector<float> pxGen;
      std::vector<float> pyGen;
      std::vector<float> pzGen;
      std::vector<float> vxGen;
      std::vector<float> vyGen;
      std::vector<float> vzGen;
      std::vector<int> NbOfDaughters;
      std::vector<float> McZmass;
      std::vector<float> McZpt;
      std::vector<float> McZpx;
      std::vector<float> McZpy;
      std::vector<float> McZpz;
      std::vector<float> McZen;
      //====================================================
      //
      //   Create vectors for Pimary Vertice variables
      //
      //====================================================
      int value3_;
      std::vector<int> nbPv;
      std::vector<int> Nbdof;
      std::vector<float> PositionX;
      std::vector<float> PositionY;
      std::vector<float> PositionZ;
      std::vector<float> PositionRho;
      //====================================================
      //
      //   Create vectors for Jets variables
      //
      //====================================================
      int jetnumber;
      std::vector<float> jetcharge;
      std::vector<float> jetet;
      std::vector<float> jetpt;
      std::vector<float> jeteta;
      std::vector<float> jetphi;

      //=============================================================
      //
      //           Create Branchs for Muons match HLT variables
      //
      //=============================================================
      std::vector<int> MuHLTMatch_nbMuonMatchHLT;
      std::vector<float> MuHLTMatch_pt;
      std::vector<float> MuHLTMatch_eta;
      std::vector<float> MuHLTMatch_phi;
      std::vector<int> MuHLTMatch_nbHLTMuonMatchReco;
      std::vector<float> MuHLTMatch_Trigger_pt;
      std::vector<float> MuHLTMatch_Trigger_eta;
      std::vector<float> MuHLTMatch_Trigger_phi;
      //=============================================================
      //                   
      //           Create Branchs for PF MET
      //
      //=============================================================
      float PFMet_pt;
      float PFMet_eta;
      float PFMet_phi;
      float PFMet_en;
      float PFMet_px;
      float PFMet_py;
      float PFMet_pz;
      float PFMet_sumEt; 
      double METSign;
      //=============================================================
      //                   
      //           Create Branchs for PileUp tree
      //
      //=============================================================
      int num_PU_vertices;
      int PU_BunchCrossing;
      //=============================================================
      //                   
      //           Create Branch for Rho
      //
      //=============================================================
      std::vector<float> Rho;
      //=============================================================
      //                   
      //           Create Branch for events reweighting
      //
      //=============================================================
      std::vector<float> MC_weighting;
      //=============================================================
      double maxAbsZ_;
      double maxd0_;
      int minNdof_;
      int NbGoodPv_;
      std::string outputFile_; // output file
      edm::InputTag vertexSrc;
      edm::InputTag pfMuons_;
      edm::InputTag pfMuonToken_;
      edm::InputTag srcSelectedMuons_;  
      edm::InputTag srcPFCandidates_; 
      edm::InputTag theRecoLabel_;
      edm::InputTag PatElectronsrc_;
      edm::InputTag genParticlesColl_;
      edm::InputTag tokentracks_;
      edm::InputTag globalMuons_;
      edm::InputTag token_globalMuons;
      edm::InputTag thePFMETCollectionToken_;
      edm::InputTag theMETSignificance_;
      edm::InputTag thejetsTag_;
      edm::InputTag PileupSrc_;
      edm::InputTag rhoIsoInputTag_;
      edm::InputTag genEventInfo_;
      /// HLT TriggerResults EDProduct
      edm::InputTag inputTag_;
      /// HLT trigger names
      edm::TriggerNames triggerNames_;
      
      /// number of HLT trigger paths requested in configuration
      unsigned int n_;
      bool firstevent_;
      
      // PG and FRC 06-07-11 try to reduce printout!
      bool debug;
      
      /// list of required HLT triggers by HLT name
      std::vector<std::string > HLTPathsByName_;
      /// list of required HLT triggers by HLT index
      std::vector<unsigned int> HLTPathsByIndex_;
      edm::InputTag triggerEvent;
      std::vector<std::string> triggerFilter;
      std::vector<std::string> triggerFilter_asym;
      //===============================
      // root file to store histograms
      TFile*  rootFile_;
      // min and max of energy histograms
      
      
      

};


#endif



