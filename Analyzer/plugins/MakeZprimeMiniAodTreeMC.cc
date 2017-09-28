//======================================================
// Make the tree for Z' boson to Di-lepton analysis    =
// electron & muon common tree                         =
// Author: Sherif Elgammal                             =
// CMSSW: 8_0_25                                       =
// Data Format: MINIAOD                                =
// Date: 04/06/2017                                    =
//======================================================
//these header files give us easy to use shortcuts for which cut
//corresponds to which which cutnr
//this is fixed for a given ID (and can be different for each ID)
//hence its hard coded
//also these headerfiles are intentionally completely standalone 
//so you can easily include them in your analysis if you find them
//useful
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "HEEP/VID/interface/CutNrs.h"
#include "HEEP/VID/interface/VIDCutCodes.h"
//new by Sherif
#include "TTree.h"
#include "TFile.h"
using namespace std;
using namespace edm;
using namespace reco;
using namespace pat;
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
// Ecal includes
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"    
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/Common/interface/RefCore.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "RecoTracker/DeDx/interface/DeDxTools.h"
// Jets
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "TH1.h"
// New by Sam
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/ValueMap.h"
//#include "SHarper/HEEPAnalyzer/interface/HEEPCutCodes.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//a struct to count the number of electrons passing/failing
//see https://twiki.cern.ch/twiki/bin/view/CMSPublic/FWMultithreadedFrameworkStreamModuleInterface
//if you dont understand why I'm doing this
namespace{
struct NrPassFail {
NrPassFail():nrPass(0),nrFail(0){}
  mutable std::atomic<int> nrPass;
  mutable std::atomic<int> nrFail;
};
}
//
// class declaration
//
class MakeZprimeMiniAodTreeMC : public edm::EDAnalyzer {
public:
  explicit MakeZprimeMiniAodTreeMC(const edm::ParameterSet&);
  ~MakeZprimeMiniAodTreeMC();
  static std::unique_ptr<NrPassFail> initializeGlobalCache(const edm::ParameterSet&) {
    return std::make_unique<NrPassFail>();
  }
  static void globalEndJob(const NrPassFail* nrPassFail) {
    std::cout <<"nr eles pass "<<nrPassFail->nrPass<<" / "<<nrPassFail->nrPass+nrPassFail->nrFail<<std::endl;
  }
  void PatElectronTree(const edm::Event& iEvent,const edm::EventSetup& es);  
  void TriggerMatchingTree(const edm::Event& iEvent,const edm::EventSetup& es);
  bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);
  void accessGenInfo(const edm::Event& iEvent,const edm::EventSetup& es);
  void EBecalRecHitsTree(const edm::Event& iEvent,const edm::EventSetup& es);
  void EEecalRecHitsTree(const edm::Event& iEvent,const edm::EventSetup& es);
  void PatMuonTree(const edm::Event& evt,const edm::EventSetup& es);
  void SuperClusterTree(const edm::Event& iEvent,const edm::EventSetup& es);
  void ComputeMuonMassVtx(const edm::Event& evt,const edm::EventSetup& es);
  void ComputeMuonMassVtx30GeV(const edm::Event& evt,const edm::EventSetup& es);
  void PrimaryVertexTree(const edm::Event& iEvent,const edm::EventSetup& es);
  bool PrimaryVertex(const reco::VertexCollection &vertices);
  void fillMET(const edm::Event& iEvent);
  void GenJetTree(const edm::Event& iEvent);
  void JetsTree(const edm::Event& iEvent,const edm::EventSetup& es);
  void EventsReWeighting(const edm::Event& evt);
  void ParticleFlowPhotonTree(const edm::Event& iEvent,const edm::EventSetup& es);
  void fillPU(const edm::Event& iEvent);
  void fillRho(const edm::Event& evt);
  void BtaggingTree(const edm::Event& iEvent);
  void TauTree(const edm::Event& iEvent);
  
private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void beginJob() override;
  virtual void endJob() override;
  TTree* mytree;  
  edm::EDGetTokenT<edm::View<pat::Electron> > eleToken_;
  // ----------member data ---------------------------   //
  edm::EDGetTokenT<double> rhoToken;
  edm::EDGetTokenT<GenEventInfoProduct> genInfoProductToken;
  edm::EDGetTokenT<reco::GenJetCollection> EDMGenJetsToken_; 
  edm::EDGetTokenT<reco::SuperClusterCollection> scProducer_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<pat::TauCollection> tauToken_;
  edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
  edm::EDGetTokenT<pat::JetCollection> jetToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
  edm::EDGetTokenT<EcalRecHitCollection> tok_EB_;
  edm::EDGetTokenT<EcalRecHitCollection> tok_EE_;
  edm::EDGetTokenT<EcalRecHitCollection> ecalRechitEBToken_;
  edm::EDGetTokenT<EcalRecHitCollection> ecalRechitEEToken_;
  edm::EDGetTokenT<double> rhoIsoInputTag_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > PileupSrc_;
  const std::vector<std::string> bDiscriminators_;
  int BosonID_;
  int ParticleID1_;
  int ParticleID2_;
  int ParticleID3_;
  int ParticleStatus_;
  double maxAbsZ_;
  double maxd0_;
  int minNdof_;
  int NbGoodPv_;
  //=============================================================
  TFile*  rootFile_;
  std::string outputFile_; // output file
  //=============================================================
  //
  //           Create Branchs for Nb of event,run,lumi
  //
  //=============================================================
  unsigned int Run;
  unsigned int Event;
  unsigned int lumi;
  unsigned int bunch;
  //==================================================
  //
  //      Create vectors for Electrons variables
  //
  //==================================================
  std::vector<int> Ele_rawId;
  std::vector<int> Ele_nbElectrons;
  std::vector<bool> Ele_isEcalDrivenSeed;
  std::vector<bool> Ele_isPassConversionVeto;
  std::vector<int> Ele_charge;
  std::vector<int> Ele_nbOfMissingHits;
  std::vector<int> Ele_nbVtx;
  std::vector<float> Ele_Et;
  std::vector<float> Ele_EtFromCaloEn;
  std::vector<float> Ele_pt;
  std::vector<float> Ele_thetaSC;
  std::vector<float> Ele_etaSC;
  std::vector<float> Ele_phiSC;
  std::vector<float> Ele_energySC;
  std::vector<float> Ele_preshowerEnergySC;
  std::vector<float> Ele_thetaTrack;
  std::vector<float> Ele_etaTrack;
  std::vector<float> Ele_phiTrack;
  std::vector<float> Ele_hadronicOverEm;
  std::vector<float> Ele_deltaEtaInSeedCluster;
  std::vector<float> Ele_deltaPhiInSeedCluster;      
  std::vector<float> Ele_deltaEtaInSC;
  std::vector<float> Ele_deltaPhiInSC;
  std::vector<float> Ele_sigmaIetaIeta;
  std::vector<float> Ele_e2x5Max;
  std::vector<float> Ele_e1x5;
  std::vector<float> Ele_frac51;
  std::vector<float> Ele_frac15;
  std::vector<float> Ele_e5x5;
  std::vector<float> Ele_e3x3;
  std::vector<float> Ele_e2x5MaxOver5x5;
  std::vector<float> Ele_e1x5Over5x5;
  std::vector<float> Ele_sigmaIetaIetaFull5x5;
  std::vector<float> Ele_e2x5MaxFull5x5;
  std::vector<float> Ele_e1x5Full5x5;
  std::vector<float> Ele_e5x5Full5x5;
  std::vector<float> Ele_e2x5MaxOver5x5Full5x5;
  std::vector<float> Ele_e1x5Over5x5Full5x5;
  std::vector<float> Ele_e2x5Right;
  std::vector<float> Ele_e2x5Left;
  std::vector<float> Ele_e2x5Top;
  std::vector<float> Ele_e2x5Bottom;
  std::vector<float> Ele_eMax;
  std::vector<float> Ele_eRight;
  std::vector<float> Ele_eLeft;
  std::vector<float> Ele_eTop;
  std::vector<float> Ele_eBottom;
  std::vector<float> Ele_dxy;
  std::vector<float> Ele_dz;
  std::vector<float> Ele_rhoIso;
  std::vector<float> Ele_fbrem;
  std::vector<float> Ele_EoverP;
  std::vector<float> Ele_Xposition;
  std::vector<float> Ele_Yposition;
  //------------- detector isolation -------------------------
  std::vector<float> Ele_EcalPlusHcald1iso;
  std::vector<float> Ele_dr03TkSumPt;
  std::vector<float> Ele_dr03TkSumPt_corrected;
  std::vector<float> Ele_dr03EcalRecHitSumEt;
  std::vector<float> Ele_dr03HcalDepth1TowerSumEt;
  std::vector<float> Ele_dr03HcalDepth1TowerSumEtBc;
  std::vector<float> Ele_hcalDepth1OverEcal;
  std::vector<float> Ele_hcalDepth2OverEcal;
  std::vector<float> Ele_dr03HcalDepth2TowerSumEt;
  std::vector<float> Ele_hcalDepth2TowerSumEtNoVeto;
  std::vector<float> Ele_hcalDepth1TowerSumEtNoVeto;                       
  //------------- PF isolation from pat::ele -------------------------
  std::vector<float> Ele_pfSumPhotonEt;
  std::vector<float> Ele_pfSumChargedHadronPt;
  std::vector<float> Ele_pfSumNeutralHadronEt;
  std::vector<float> Ele_pfSumPUPt;
  std::vector<float> Ele_pfDeltaBeta;
  std::vector<float> Ele_x;
  std::vector<float> Ele_y;
  std::vector<float> Ele_z;
  std::vector<float> Ele_zTrackPositionAtVtx;
  std::vector<int> Ele_ieta;
  std::vector<float> Ele_phiWidth;
  std::vector<float> Ele_etaWidth;
  std::vector<float> Ele_nrSatCrys;
  std::vector<bool>  Ele_isPassHeepID;
  //=============================================================
  //
  //           Create Branchs for Muons match HLT variables
  //
  //=============================================================
  std::vector<int>    HLT_nb;
  std::vector<string> HLT_name;
  std::vector<bool>   HLT_isaccept;
  std::vector<int>    HLTObj_nbObj;
  std::vector<float>  HLTObj_pt;
  std::vector<float>  HLTObj_eta;
  std::vector<float>  HLTObj_phi;
  std::vector<string> HLTObj_collection;
   //=============================================================
  //
  //           Create Branchs for Pimary Vertice variables
  //
  //=============================================================
  std::vector<int> nbPv;
  std::vector<int> Nbdof;
  std::vector<float> PositionRho;
  std::vector<float> PositionX;
  std::vector<float> PositionY;
  std::vector<float> PositionZ;
  //====================================================
  //
  // Create vectors for Jets variables
  //
  //====================================================
  std::vector<int>   jet_nb;
  std::vector<float> jet_charge;
  std::vector<float> jet_et;
  std::vector<float> jet_pt;
  std::vector<float> jet_eta;
  std::vector<float> jet_phi;
  std::vector<float> jet_en;
  std::vector<float> jet_theta;
  std::vector<float> jet_beta;
  std::vector<float> jet_pileup_mva_disc;
  //====================================================
  //
  // Create vectors for BTagging variables
  //
  //====================================================
  std::vector<int>   Nb_bDiscriminators;
  std::vector<int>   jet_btag_flavor;
  std::vector<float> jet_btag_pfCSVv2IVF_discriminator;
  std::vector<float> jet_btag_pt;
  std::vector<float> jet_btag_eta;
  std::vector<float> jet_btag_phi;
  //====================================================
  //
  // Create vectors for Taus variables
  //
  //====================================================
  std::vector<int>   Nb_taus;
  std::vector<float> Tau_pt;
  std::vector<float> Tau_eta;
  std::vector<float> Tau_phi;
  std::vector<int>   Tau_id;
  std::vector<float> Tau_LooseCombinedIsolationDeltaBetaCorr3Hits;
  //=============================================================
  //                   
  //           Create Branchs for PileUp tree
  //
  //=============================================================
  int num_PU_vertices;
  int PU_BunchCrossing;
  int num_PU_gen_vertices;
  //=============================================================
  //                   
  //           Create Branch for Rho
  //
  //=============================================================
  float Rho;
  //=============================================================
  //                   
  //           Create Branch for events reweighting
  //
  //=============================================================
  std::vector<float> MC_weighting;
  //===================================================
  //
  //    Create vectors for gen particles variables
  //
  //===================================================
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
  //==================================================
  //
  //           Create vectors for Muons variables
  //
  //==================================================
  std::vector<int> Mu_nbMuon;
  std::vector<bool> Mu_passOldMatchedStationsCut;
  std::vector<bool> Mu_passNewMatchedStationsCut;
  std::vector<bool> Mu_isTightMuon;
  std::vector<bool> Mu_isLooseMuon;
  std::vector<bool> Mu_isGlobalMuon;
  std::vector<bool> Mu_isHighPtMuon;
  //std::vector<bool> Mu_isMuonsCleaned;
  std::vector<bool> Mu_isPF;
  std::vector<bool> Mu_isTrackerMuon;
  std::vector<float> Mu_en;
  std::vector<float> Mu_pt;
  std::vector<float> Mu_eta;
  std::vector<float> Mu_phi;
  std::vector<float> Mu_et;
  std::vector<float> Mu_charge;      
  std::vector<float> Mu_normalizedChi2;
  std::vector<float> Mu_trackiso;
  std::vector<int> Mu_numberOfMatchedStations;
  std::vector<int> Mu_numberOfValidPixelHits;
  std::vector<int> Mu_numberOftrackerLayersWithMeasurement;
  std::vector<float> Mu_pfSumChargedHadronPt;
  std::vector<float> Mu_pfSumNeutralHadronEt;
  std::vector<float> Mu_PFSumPhotonEt;
  std::vector<float> Mu_pfSumPUPt;
  std::vector<float> Mu_calEnergy;
  std::vector<int> Mu_numberOfValidMuonHits;
  std::vector<float> Mu_emIso;
  std::vector<float> Mu_hadIso;
  std::vector<float> Mu_VTXnormalizedChi2;
  std::vector<float> Mu_vtxTotalChi2;
  std::vector<float> Mu_vtxMass;
  std::vector<float> Mu_vtxMassLept;
  std::vector<float> Mu_vtxNormChi2;
  std::vector<float> Mu_vtxMass30GeV;
  std::vector<float> Mu_vtxNormChi30GeV;
  std::vector<float> Mu_vtxMass30GeVLept;
  std::vector<int> Mu_nbofpv;
  std::vector<float> Mu_pxTunePMuonBestTrack;
  std::vector<float> Mu_pyTunePMuonBestTrack;
  std::vector<float> Mu_pzTunePMuonBestTrack;
  std::vector<float> Mu_pTunePMuonBestTrack;
  std::vector<float> Mu_etaTunePMuonBestTrack;
  std::vector<float> Mu_ptTunePMuonBestTrack;
  std::vector<float> Mu_phiTunePMuonBestTrack;
  std::vector<float> Mu_thetaTunePMuonBestTrack;
  std::vector<float> Mu_chargeTunePMuonBestTrack;
  std::vector<float> Mu_dPToverPTTunePMuonBestTrack;
  std::vector<float> Mu_absdxyTunePMuonBestTrack;
  std::vector<float> Mu_absdzTunePMuonBestTrack;
  std::vector<float> Mu_ptInnerTrack;
  std::vector<float> Mu_pxInnerTrack;
  std::vector<float> Mu_pyInnerTrack;
  std::vector<float> Mu_pzInnerTrack;
  std::vector<float> Mu_pInnerTrack;
  std::vector<float> Mu_etaInnerTrack;
  std::vector<float> Mu_phiInnerTrack;
  std::vector<float> Mu_thetaInnerTrack;
  std::vector<float> Mu_chargeInnerTrack;
  std::vector<float> Mu_dPToverPTInnerTrack;
  std::vector<float> Mu_absdxyInnerTrack;
  std::vector<float> Mu_absdzInnerTrack;
  std::vector<float> Mu_absdxy;
  std::vector<float> Mu_absdz;
  std::vector<float> Mu_patDeltaBeta;
  std::vector<int> Mu_stationMask;
  std::vector<int> Mu_numberOfMatchedRPCLayers;
  //=============================================================
  //                   
  //           Create Branches for PF MET
  //
  //=============================================================
  double GenMet_pt;
  //The default type1 corrected MET
  double PFMet_et_cor;
  double PFMet_pt_cor;
  double PFMet_phi_cor;
  double PFMet_en_cor;
  double PFMet_px_cor;
  double PFMet_py_cor;
  double PFMet_pz_cor;
  double PFMet_sumEt_cor;
  //The Raw PF Met (un-corrected MET)
  double PFMet_pt_uncor;
  double PFMet_phi_uncor;
  double PFMet_sumEt_uncor;
  //The raw calo ETmiss  
  double CaloMet_pt;
  double CaloMet_phi;
  double CaloMet_sumEt;
  //double METSign;
  double PFMet_shiftedPt_JetEnUp;
  double PFMet_shiftedPt_JetEnDown;
  //Type I MET Uncertainties
  double pfMETshiftedPtJetEnUp;
  double pfMETshiftedPtJetEnDn;
  double pfMETshiftedPtEleEnUp;
  double pfMETshiftedPtEleEnDn;
  double pfMETshiftedPtMuEnUp; 
  double pfMETshiftedPtMuEnDn; 
  double pfMETshiftedPtJetResUp;
  double pfMETshiftedPtJetResDn;
  double pfMETshiftedPtUnclEnUp;
  double pfMETshiftedPtUnclEnDn;
  double pfMETshiftedPtPhoEnUp;
  double pfMETshiftedPtPhoEnDn;
  //double pfMETshiftedPtJetEnUpSmear;
  //double pfMETshiftedPtJetEnDnSmear;
  //double pfMETUncertaintySize;      
  //double pfMETFullUncertaintySize;
  //===================================================
  //
  //    Create vectors for gen Jets variables
  //
  //===================================================
  std::vector<int> iGenJet;
  std::vector<int> idGenJet;
  std::vector<int> statusGenJet;
  std::vector<int> chargeGenJet;
  std::vector<float> ptGenJet;
  std::vector<float> etaGenJet;
  std::vector<float> phiGenJet;
  //===================================================
  //
  //    Create vectors for photons variables
  //
  //===================================================
  int pfphoton_size;
  std::vector<float> pfphoton_pt;
  std::vector<float> pfphoton_eta;
  std::vector<float> pfphoton_phi;
  std::vector<float> pfphoton_theta;
  //std::vector<float> pfphoton_PFchHad;
  //std::vector<float> pfphoton_PFneuHad;
  //std::vector<float> pfphoton_PFphoton;
  //std::vector<float> pfphoton_PFPUchAllPart;
  //std::vector<float> pfphoton_PFX_rho;
};

MakeZprimeMiniAodTreeMC::MakeZprimeMiniAodTreeMC(const edm::ParameterSet& iConfig):
  eleToken_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("eles"))),
  rhoToken(consumes <double> (edm::InputTag(std::string("fixedGridRhoFastjetAll")))),
  genInfoProductToken(consumes <GenEventInfoProduct> (edm::InputTag(std::string("generator")))),
  EDMGenJetsToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("JetSource"))),
  scProducer_(consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("scProducer"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
  packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
  tok_EB_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBrecHitCollectionLabel"))),
  tok_EE_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EErecHitCollectionLabel"))),
  ecalRechitEBToken_(consumes<EcalRecHitCollection>(iConfig.getParameter< edm::InputTag >("ecalRechitEB"))),
  ecalRechitEEToken_(consumes<EcalRecHitCollection>(iConfig.getParameter< edm::InputTag >("ecalRechitEE"))),
  rhoIsoInputTag_(consumes<double>(iConfig.getParameter< edm::InputTag >("rhoIsoInputTag"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter< edm::InputTag >("pfCands"))),
  PileupSrc_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("PileupSrc"))),
  bDiscriminators_(iConfig.getParameter<std::vector<std::string> >("bDiscriminators")),
  outputFile_(iConfig.getParameter<std::string>("outputFile"))
{
  rootFile_   = TFile::Open(outputFile_.c_str(),"RECREATE"); // open output file to store histograms  
  BosonID_            = iConfig.getParameter<int>("GenBosonID");
  ParticleID1_        = iConfig.getParameter<int>("ParticleID1");
  ParticleID2_        = iConfig.getParameter<int>("ParticleID2");
  ParticleID3_        = iConfig.getParameter<int>("ParticleID3");
  ParticleStatus_     = iConfig.getParameter<int>("ParticleStatus");
  maxAbsZ_            = iConfig.getParameter<double>("maxAbsZ");
  maxd0_              = iConfig.getParameter<double>("maxd0");
  minNdof_            = iConfig.getParameter<int>("minndof");
  NbGoodPv_           = iConfig.getParameter<int>("NbGoodPv");
}

MakeZprimeMiniAodTreeMC::~MakeZprimeMiniAodTreeMC()
{
  delete rootFile_;
}

// ------------ method called once each job just before starting event loop  ------------
void MakeZprimeMiniAodTreeMC::beginJob()
{
  
  // go to *OUR* rootfile and book histograms
  rootFile_->cd();
  // Declare histograms
  mytree  = new TTree("tree","tr");
  //=============================================================
  //
  //           Create Branchs for Nb of event,run,lumi
  //
  //=============================================================
  mytree->Branch("event_runNo",  &Run,   "event_runNo/i");
  mytree->Branch("event_evtNo",  &Event, "event_evtNo/i");
  mytree->Branch("event_lumi",   &lumi,  "event_lumi/i");
  mytree->Branch("event_bunch",  &bunch, "event_bunch/i");
  //==================================================
  //
  //      Create Branchs for Electrons variables
  //
  //==================================================
  mytree->Branch("Ele_nbElectrons",&Ele_nbElectrons);
  mytree->Branch("Ele_isEcalDrivenSeed",&Ele_isEcalDrivenSeed);
  mytree->Branch("Ele_isPassConversionVeto",&Ele_isPassConversionVeto);
  mytree->Branch("Ele_charge",&Ele_charge);
  mytree->Branch("Ele_nbOfMissingHits",&Ele_nbOfMissingHits);
  mytree->Branch("Ele_nbVtx",&Ele_nbVtx);
  mytree->Branch("Ele_Et",&Ele_Et);
  mytree->Branch("Ele_EtFromCaloEn",&Ele_EtFromCaloEn);
  mytree->Branch("Ele_pt",&Ele_pt);
  mytree->Branch("Ele_thetaSC",&Ele_thetaSC);
  mytree->Branch("Ele_etaSC",&Ele_etaSC);
  mytree->Branch("Ele_phiSC",&Ele_phiSC);
  mytree->Branch("Ele_energySC",&Ele_energySC);
  mytree->Branch("Ele_preshowerEnergySC",&Ele_preshowerEnergySC);
  mytree->Branch("Ele_thetaTrack",&Ele_thetaTrack);
  mytree->Branch("Ele_etaTrack",&Ele_etaTrack);
  mytree->Branch("Ele_phiTrack",&Ele_phiTrack);
  mytree->Branch("Ele_hadronicOverEm",&Ele_hadronicOverEm);
  mytree->Branch("Ele_deltaEtaInSeedCluster",&Ele_deltaEtaInSeedCluster);
  mytree->Branch("Ele_deltaPhiInSeedCluster",&Ele_deltaPhiInSeedCluster);
  mytree->Branch("Ele_deltaEtaInSC",&Ele_deltaEtaInSC);
  mytree->Branch("Ele_deltaPhiInSC",&Ele_deltaPhiInSC);
  mytree->Branch("Ele_sigmaIetaIeta",&Ele_sigmaIetaIeta);
  mytree->Branch("Ele_e2x5Max",&Ele_e2x5Max);
  mytree->Branch("Ele_e1x5",&Ele_e1x5);
  mytree->Branch("Ele_frac15",&Ele_frac15);
  mytree->Branch("Ele_frac51",&Ele_frac51);
  mytree->Branch("Ele_e5x5",&Ele_e5x5);
  mytree->Branch("Ele3x3",&Ele_e3x3);
  mytree->Branch("Ele_e2x5MaxOver5x5",&Ele_e2x5MaxOver5x5);
  mytree->Branch("Ele_e1x5Over5x5",&Ele_e1x5Over5x5);
  mytree->Branch("Ele_sigmaIetaIetaFull5x5",&Ele_sigmaIetaIetaFull5x5);
  mytree->Branch("Ele_e2x5MaxFull5x5",&Ele_e2x5MaxFull5x5);
  mytree->Branch("Ele_e1x5Full5x5",&Ele_e1x5Full5x5);
  mytree->Branch("Ele_e5x5Full5x5",&Ele_e5x5Full5x5);
  mytree->Branch("Ele_e2x5MaxOver5x5Full5x5",&Ele_e2x5MaxOver5x5Full5x5);
  mytree->Branch("Ele_e1x5Over5x5Full5x5",&Ele_e1x5Over5x5Full5x5);
  mytree->Branch("Ele_e2x5Right",&Ele_e2x5Right);
  mytree->Branch("Ele_e2x5Left",&Ele_e2x5Left);
  mytree->Branch("Ele_e2x5Top",&Ele_e2x5Top);
  mytree->Branch("Ele_e2x5Bottom",&Ele_e2x5Bottom);
  mytree->Branch("Ele_eMax",&Ele_eMax);
  mytree->Branch("Ele_eRight",&Ele_eRight);
  mytree->Branch("Ele_eLeft",&Ele_eLeft);
  mytree->Branch("Ele_eTop",&Ele_eTop);
  mytree->Branch("Ele_eBottom",&Ele_eBottom);
  mytree->Branch("Ele_dxy",&Ele_dxy);
  mytree->Branch("Ele_dz",&Ele_dz);
  mytree->Branch("Ele_rhoIso",&Ele_rhoIso);
  mytree->Branch("Ele_fbrem",&Ele_fbrem);
  mytree->Branch("Ele_EoverP",&Ele_EoverP);
  mytree->Branch("Ele_Xposition",&Ele_Xposition);
  mytree->Branch("Ele_Yposition",&Ele_Yposition);
  mytree->Branch("Ele_EcalPlusHcald1iso",&Ele_EcalPlusHcald1iso);
  mytree->Branch("Ele_dr03EcalRecHitSumEt",&Ele_dr03EcalRecHitSumEt);
  mytree->Branch("Ele_dr03HcalDepth1TowerSumEt",&Ele_dr03HcalDepth1TowerSumEt);
  mytree->Branch("Ele_dr03HcalDepth1TowerSumEtBc",&Ele_dr03HcalDepth1TowerSumEtBc);
  mytree->Branch("Ele_hcalDepth1OverEcal",&Ele_hcalDepth1OverEcal);
  mytree->Branch("Ele_hcalDepth2OverEcal",&Ele_hcalDepth2OverEcal);
  mytree->Branch("Ele_dr03HcalDepth2TowerSumEt",&Ele_dr03HcalDepth2TowerSumEt);
  mytree->Branch("Ele_hcalDepth2TowerSumEtNoVeto",&Ele_hcalDepth2TowerSumEtNoVeto);
  mytree->Branch("Ele_hcalDepth1TowerSumEtNoVeto",&Ele_hcalDepth1TowerSumEtNoVeto);
  mytree->Branch("Ele_pfSumPhotonEt",&Ele_pfSumPhotonEt);
  mytree->Branch("Ele_pfSumChargedHadronPt",&Ele_pfSumChargedHadronPt);
  mytree->Branch("Ele_pfSumNeutralHadronEt",&Ele_pfSumNeutralHadronEt);
  mytree->Branch("Ele_pfSumPUPt",&Ele_pfSumPUPt);
  mytree->Branch("Ele_pfDeltaBeta",&Ele_pfDeltaBeta);
  mytree->Branch("Ele_rawId",&Ele_rawId);
  mytree->Branch("Ele_x",&Ele_x);
  mytree->Branch("Ele_y",&Ele_y);
  mytree->Branch("Ele_z",&Ele_z);
  mytree->Branch("Ele_zTrackPositionAtVtx",&Ele_zTrackPositionAtVtx);
  mytree->Branch("Ele_ieta",&Ele_ieta);
  mytree->Branch("Ele_phiWidth",&Ele_phiWidth);
  mytree->Branch("Ele_etaWidth",&Ele_etaWidth);
  mytree->Branch("Ele_dr03TkSumPt",&Ele_dr03TkSumPt);
  mytree->Branch("Ele_dr03TkSumPt_corrected",&Ele_dr03TkSumPt_corrected);
  mytree->Branch("Ele_nrSatCrys",&Ele_nrSatCrys);
  mytree->Branch("Ele_isPassHeepID",&Ele_isPassHeepID);
  //============================================================= 
  //
  //           Create Branchs for Muons match HLT variables
  //
  //=============================================================  
  mytree->Branch("HLT_nb", &HLT_nb);
  mytree->Branch("HLT_name", &HLT_name);
  mytree->Branch("HLT_isaccept", &HLT_isaccept);
  mytree->Branch("HLTObj_nbObj",&HLTObj_nbObj);
  mytree->Branch("HLTObj_pt",&HLTObj_pt);
  mytree->Branch("HLTObj_eta",&HLTObj_eta);
  mytree->Branch("HLTObj_phi",&HLTObj_phi);
  mytree->Branch("HLTObj_collection", &HLTObj_collection);
  //===================================================
  //
  //    Create vectors for gen jets variables
  //
  //===================================================
  mytree->Branch("iGenJet",&iGenJet);
  mytree->Branch("idGenJet",&idGenJet);
  mytree->Branch("statusGenJet",&statusGenJet);
  mytree->Branch("chargeGenJet",&chargeGenJet);
  mytree->Branch("ptGenJet",&ptGenJet);
  mytree->Branch("etaGenJet",&etaGenJet);
  mytree->Branch("phiGenJet",&phiGenJet);
  //===================================================
  //
  //    Create vectors for gen particles variables
  //
  //===================================================
  mytree->Branch("iGen",&iGen);
  mytree->Branch("idGen",&idGen);
  mytree->Branch("statusGen",&statusGen);
  mytree->Branch("ptGen",&ptGen);
  mytree->Branch("etaGen",&etaGen);
  mytree->Branch("phiGen",&phiGen);
  mytree->Branch("chargeGen",&chargeGen);
  mytree->Branch("EnergyGen",&EnergyGen);
  mytree->Branch("pxGen",&pxGen);
  mytree->Branch("pyGen",&pyGen);
  mytree->Branch("pzGen",&pzGen);
  //=============================================================  
  //
  //           Create Branchs for Muons variables 
  //
  //=============================================================
  mytree->Branch("Mu_nbMuon",&Mu_nbMuon);
  mytree->Branch("Mu_passOldMatchedStationsCut",&Mu_passOldMatchedStationsCut);
  mytree->Branch("Mu_passNewMatchedStationsCut",&Mu_passNewMatchedStationsCut);
  mytree->Branch("Mu_isTightMuon",&Mu_isTightMuon);
  mytree->Branch("Mu_isLooseMuon",&Mu_isLooseMuon);
  mytree->Branch("Mu_isGlobalMuon",&Mu_isGlobalMuon);
  mytree->Branch("Mu_isPF",&Mu_isPF);
  //mytree->Branch("Mu_isMuonsCleaned",&Mu_isMuonsCleaned);
  mytree->Branch("Mu_isTrackerMuon",&Mu_isTrackerMuon);
  mytree->Branch("Mu_et",&Mu_et);
  mytree->Branch("Mu_en",&Mu_en);
  mytree->Branch("Mu_pt",&Mu_pt);
  mytree->Branch("Mu_eta",&Mu_eta);
  mytree->Branch("Mu_phi",&Mu_phi);
  mytree->Branch("Mu_charge",&Mu_charge);
  mytree->Branch("Mu_ptTunePMuonBestTrack",&Mu_ptTunePMuonBestTrack);
  mytree->Branch("Mu_pxTunePMuonBestTrack",&Mu_pxTunePMuonBestTrack);
  mytree->Branch("Mu_pyTunePMuonBestTrack",&Mu_pyTunePMuonBestTrack);
  mytree->Branch("Mu_pzTunePMuonBestTrack",&Mu_pzTunePMuonBestTrack);
  mytree->Branch("Mu_pTunePMuonBestTrack",&Mu_pTunePMuonBestTrack);
  mytree->Branch("Mu_etaTunePMuonBestTrack",&Mu_etaTunePMuonBestTrack);
  mytree->Branch("Mu_phiTunePMuonBestTrack",&Mu_phiTunePMuonBestTrack);
  mytree->Branch("Mu_thetaTunePMuonBestTrack",&Mu_thetaTunePMuonBestTrack);
  mytree->Branch("Mu_chargeTunePMuonBestTrack",&Mu_chargeTunePMuonBestTrack);
  mytree->Branch("Mu_dPToverPTTunePMuonBestTrack",&Mu_dPToverPTTunePMuonBestTrack);
  mytree->Branch("Mu_absdxyTunePMuonBestTrack",&Mu_absdxyTunePMuonBestTrack);
  mytree->Branch("Mu_absdzTunePMuonBestTrack",&Mu_absdzTunePMuonBestTrack);
  mytree->Branch("Mu_ptInnerTrack",&Mu_ptInnerTrack);
  mytree->Branch("Mu_pxInnerTrack",&Mu_pxInnerTrack);
  mytree->Branch("Mu_pyInnerTrack",&Mu_pyInnerTrack);
  mytree->Branch("Mu_pzInnerTrack",&Mu_pzInnerTrack);
  mytree->Branch("Mu_pInnerTrack",&Mu_pInnerTrack);
  mytree->Branch("Mu_etaInnerTrack",&Mu_etaInnerTrack);
  mytree->Branch("Mu_phiInnerTrack",&Mu_phiInnerTrack);
  mytree->Branch("Mu_thetaInnerTrack",&Mu_thetaInnerTrack);
  mytree->Branch("Mu_chargeInnerTrack",&Mu_chargeInnerTrack);
  mytree->Branch("Mu_dPToverPTInnerTrack",&Mu_dPToverPTInnerTrack);
  mytree->Branch("Mu_absdxyInnerTrack",&Mu_absdxyInnerTrack);
  mytree->Branch("Mu_absdzInnerTrack",&Mu_absdzInnerTrack);
  mytree->Branch("Mu_normalizedChi2",&Mu_normalizedChi2);
  mytree->Branch("Mu_absdxy",&Mu_absdxy);
  mytree->Branch("Mu_absdz",&Mu_absdz);
  mytree->Branch("Mu_vtxMass",&Mu_vtxMass);
  mytree->Branch("Mu_vtxNormChi2",&Mu_vtxNormChi2);
  mytree->Branch("Mu_vtxMassLept",&Mu_vtxMassLept);
  mytree->Branch("Mu_numberOfMatchedStations",&Mu_numberOfMatchedStations);
  mytree->Branch("Mu_numberOfValidPixelHits",&Mu_numberOfValidPixelHits);
  mytree->Branch("Mu_numberOfValidMuonHits",&Mu_numberOfValidMuonHits);
  mytree->Branch("Mu_numberOftrackerLayersWithMeasurement",&Mu_numberOftrackerLayersWithMeasurement);
  mytree->Branch("Mu_emIso",&Mu_emIso);
  mytree->Branch("Mu_hadIso",&Mu_hadIso);
  mytree->Branch("Mu_trackiso",&Mu_trackiso);
  mytree->Branch("Mu_pfSumChargedHadronPt",&Mu_pfSumChargedHadronPt);
  mytree->Branch("Mu_pfSumNeutralHadronEt",&Mu_pfSumNeutralHadronEt);
  mytree->Branch("Mu_PFSumPhotonEt",&Mu_PFSumPhotonEt);
  mytree->Branch("Mu_pfSumPUPt",&Mu_pfSumPUPt);
  mytree->Branch("Mu_nbofpv",&Mu_nbofpv);
  mytree->Branch("Mu_patDeltaBeta",&Mu_patDeltaBeta);
  mytree->Branch("Mu_stationMask",&Mu_stationMask);
  mytree->Branch("Mu_numberOfMatchedRPCLayers",&Mu_numberOfMatchedRPCLayers);
  mytree->Branch("Mu_vtxMass30GeV",&Mu_vtxMass30GeV);
  mytree->Branch("Mu_vtxNormChi30GeV",&Mu_vtxNormChi30GeV);
  mytree->Branch("Mu_vtxMass30GeVLept",&Mu_vtxMass30GeVLept);
  //=============================================================
  //                   
  //           Create Branches for PF MET
  //
  //=============================================================
  mytree->Branch("GenMet_pt",&GenMet_pt,"GenMet_pt/D");
  //The default type1 corrected MET
  mytree->Branch("PFMet_et_cor",   &PFMet_et_cor,   "PFMet_et_cor/D");
  mytree->Branch("PFMet_pt_cor",   &PFMet_pt_cor,   "PFMet_pt_cor/D");
  mytree->Branch("PFMet_phi_cor",  &PFMet_phi_cor,  "PFMet_phi_cor/D");
  mytree->Branch("PFMet_en_cor",   &PFMet_en_cor,   "PFMet_en_cor/D");
  mytree->Branch("PFMet_px_cor",   &PFMet_px_cor,   "PFMet_px_cor/D");
  mytree->Branch("PFMet_py_cor",   &PFMet_py_cor,   "PFMet_py_cor/D");
  mytree->Branch("PFMet_pz_cor",   &PFMet_pz_cor,   "PFMet_pz_cor/D");
  mytree->Branch("PFMet_sumEt_cor",&PFMet_sumEt_cor,"PFMet_sumEt_cor/D");
  //The raw calo ETmiss  
  mytree->Branch("CaloMet_pt",   &CaloMet_pt,   "CaloMet_pt/D");
  mytree->Branch("CaloMet_phi",  &CaloMet_phi,  "CaloMet_phi/D");
  mytree->Branch("CaloMet_sumEt",&CaloMet_sumEt,"CaloMet_sumEt/D");
  mytree->Branch("PFMet_shiftedPt_JetEnUp",&PFMet_shiftedPt_JetEnUp,"PFMet_shiftedPt_JetEnUp/D");
  mytree->Branch("PFMet_shiftedPt_JetEnDown",&PFMet_shiftedPt_JetEnDown,"PFMet_shiftedPt_JetEnDown/D");
  //The Raw PF Met (un-corrected MET)
  /*mytree->Branch("PFMet_pt_uncor",   &PFMet_pt_uncor,   "PFMet_pt_uncor/D");
    mytree->Branch("PFMet_phi_uncor",  &PFMet_phi_uncor,  "PFMet_phi_uncor/D");
    mytree->Branch("PFMet_sumEt_uncor",&PFMet_sumEt_uncor,"PFMet_sumEt_uncor/D");*/
  //mytree->Branch("METSign",&METSign,"METSign/D");
  
  //Type I MET Uncertainties
  mytree->Branch("PFMet_METshiftedPtJetEnUp",  &pfMETshiftedPtJetEnUp,  "PFMet_METshiftedPtJetEnUp/D");
  mytree->Branch("PFMet_METshiftedPtJetEnDn",  &pfMETshiftedPtJetEnDn,  "PFMet_METshiftedPtJetEnDn/D");
  mytree->Branch("PFMet_METshiftedPtEleEnUp",  &pfMETshiftedPtEleEnUp,  "PFMet_METshiftedPtEleEnUp/D");
  mytree->Branch("PFMet_METshiftedPtEleEnDn",  &pfMETshiftedPtEleEnDn,  "PFMet_METshiftedPtEleEnDn/D");
  mytree->Branch("PFMet_METshiftedPtMuEnUp",   &pfMETshiftedPtMuEnUp,   "PFMet_METshiftedPtMuEnUp/D");
  mytree->Branch("PFMet_METshiftedPtMuEnDn",   &pfMETshiftedPtMuEnDn,   "PFMet_METshiftedPtMuEnDn/D");
  mytree->Branch("PFMet_METshiftedPtJetResUp", &pfMETshiftedPtJetResUp, "PFMet_METshiftedPtJetResUp/D");
  mytree->Branch("PFMet_METshiftedPtJetResDn", &pfMETshiftedPtJetResDn, "PFMet_METshiftedPtJetResDn/D");
  mytree->Branch("PFMet_METshiftedPtUnclEnUp", &pfMETshiftedPtUnclEnUp, "PFMet_METshiftedPtUnclEnUp/D");
  mytree->Branch("PFMet_METshiftedPtUnclEnDn", &pfMETshiftedPtUnclEnDn, "PFMet_METshiftedPtUnclEnDn/D");
  mytree->Branch("PFMet_METshiftedPtPhoEnUp",  &pfMETshiftedPtPhoEnUp,  "PFMet_METshiftedPtPhoEnUp/D");
  mytree->Branch("PFMet_METshiftedPtPhoEnDn",  &pfMETshiftedPtPhoEnDn,  "PFMet_METshiftedPtPhoEnDn/D");
  //mytree->Branch("PFMet_METshiftedPtJetEnUpSmear",  &pfMETshiftedPtJetEnUpSmear,  "PFMet_METshiftedPtJetEnUpSmear/D");
  //mytree->Branch("PFMet_METshiftedPtJetEnDnSmear",  &pfMETshiftedPtJetEnDnSmear,  "PFMet_METshiftedPtJetEnDnSmear/D");
  //mytree->Branch("PFMet_METUncertaintySize",        &pfMETUncertaintySize,  "PFMet_METUncertaintySize/D");
  //mytree->Branch("PFMet_METFullUncertaintySize",    &pfMETFullUncertaintySize,  "PFMet_METFullUncertaintySize/D");
  //=============================================================
  //
  // Create Branches for jets variables 
  //
  //============================================================= 
  mytree->Branch("jet_nb",&jet_nb);
  mytree->Branch("jet_charge",&jet_charge);
  mytree->Branch("jet_et",&jet_et);
  mytree->Branch("jet_pt",&jet_pt);
  mytree->Branch("jet_eta",&jet_eta);
  mytree->Branch("jet_phi",&jet_phi);
  mytree->Branch("jet_en",&jet_en);
  mytree->Branch("jet_theta",&jet_theta);
  mytree->Branch("jet_beta",&jet_beta);
  mytree->Branch("jet_pileup_mva_disc",&jet_pileup_mva_disc);
  //=============================================================
  //
  // Create Branches for Btagging variables 
  //
  //============================================================= 
  mytree->Branch("Nb_bDiscriminators",&Nb_bDiscriminators);
  mytree->Branch("jet_btag_pt",&jet_btag_pt);
  mytree->Branch("jet_btag_eta",&jet_btag_eta);
  mytree->Branch("jet_btag_phi",&jet_btag_phi);
  mytree->Branch("jet_btag_flavor",&jet_btag_flavor);
  mytree->Branch("jet_btag_pfCSVv2IVF_discriminator",&jet_btag_pfCSVv2IVF_discriminator);
  //====================================================
  //
  // Create vectors for Taus variables
  //
  //====================================================
  mytree->Branch("Nb_taus",&Nb_taus);
  mytree->Branch("Tau_pt",&Tau_pt);
  mytree->Branch("Tau_eta",&Tau_eta);
  mytree->Branch("Tau_phi",&Tau_phi);
  mytree->Branch("Tau_id",&Tau_id);
  mytree->Branch("Tau_LooseCombinedIsolationDeltaBetaCorr3Hits",&Tau_LooseCombinedIsolationDeltaBetaCorr3Hits);
  //===================================================
  //
  //    Create vectors for photons variables
  //
  //===================================================
  mytree->Branch("pfphoton_size",&pfphoton_size,"pfphoton_size/I");
  mytree->Branch("pfphoton_pt",&pfphoton_pt);
  mytree->Branch("pfphoton_eta",&pfphoton_eta);
  mytree->Branch("pfphoton_phi",&pfphoton_phi);
  mytree->Branch("pfphoton_theta",&pfphoton_theta);
  //mytree->Branch("pfphoton_PFchHad",&pfphoton_PFchHad);
  //mytree->Branch("pfphoton_PFneuHad",&pfphoton_PFneuHad);
  //mytree->Branch("pfphoton_PFphoton",&pfphoton_PFphoton);
  //mytree->Branch("pfphoton_PFPUchAllPart",&pfphoton_PFPUchAllPart);
  //mytree->Branch("pfphoton_pfphoton_PFX_rho",&pfphoton_pfphoton_PFX_rho);
  //=============================================================
  //                   
  //           Create Branchs for PileUp tree  
  //
  //=============================================================
  mytree->Branch("num_PU_vertices",&num_PU_vertices,"num_PU_vertices/I");
  mytree->Branch("PU_BunchCrossing",&PU_BunchCrossing,"PU_BunchCrossing/I");
  mytree->Branch("num_PU_gen_vertices",&num_PU_gen_vertices,"num_PU_gen_vertices/I");
  //=============================================================
  //                   
  //           Create Branch for Rho
  //
  //=============================================================
  mytree->Branch("Rho",&Rho,"Rho/F");
  //=============================================================
  //                   
  //           Create Branch for events reweighting
  //
  //=============================================================
  mytree->Branch("MC_weighting",&MC_weighting);
  //=============================================================
  //
  //           Create Branchs for Pimary Vertice variables
  //
  //=============================================================
  /*mytree->Branch("nbPv",&nbPv);
  mytree->Branch("Nbdof",&Nbdof);
  mytree->Branch("PositionRho",&PositionRho);
  mytree->Branch("PositionX",&PositionX);
  mytree->Branch("PositionY",&PositionY);
  mytree->Branch("PositionZ",&PositionZ);*/

}

void MakeZprimeMiniAodTreeMC::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //==============================================
  //=        Begin of the main program           =
  //============================================== 
  Run   = (unsigned int) iEvent.id().run();
  Event = (unsigned int) iEvent.id().event();
  lumi  = (unsigned int) iEvent.luminosityBlock();
  bunch = (unsigned int) iEvent.bunchCrossing();
  //===================== Handle For Primary Vertics ===========
  /*edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByToken(vtxToken_,pvHandle);
  const reco::VertexCollection &vertices = *pvHandle.product();
  bool GoodPv = PrimaryVertex(vertices);
  //cout << "At least one good primary vertex" << GoodPv << endl;
  if( GoodPv == 0 ) return; //the aim of that command is to select only events*/
  PatElectronTree(iEvent,iSetup);
  PatMuonTree(iEvent,iSetup);
  TriggerMatchingTree(iEvent,iSetup); 
  ComputeMuonMassVtx(iEvent,iSetup);
  ComputeMuonMassVtx30GeV(iEvent,iSetup);
  accessGenInfo(iEvent,iSetup);
  fillMET(iEvent);  
  JetsTree(iEvent,iSetup);
  GenJetTree(iEvent);
  ParticleFlowPhotonTree(iEvent,iSetup);
  EventsReWeighting(iEvent);
  fillPU(iEvent);
  fillRho(iEvent);
  BtaggingTree(iEvent);
  TauTree(iEvent);
  //==============================================
  //=        End of the main program             =
  //============================================== 
  mytree->Fill();
}
// ------------ method called once each job just after ending the event loop  ------------
void MakeZprimeMiniAodTreeMC::endJob() 
{
  // go to *OUR* root file and store histograms
  rootFile_->cd();
  mytree->Write();
  rootFile_->Close();
}
//define this as a plug-in
DEFINE_FWK_MODULE(MakeZprimeMiniAodTreeMC);
//=============================================================
//
//                        Start the Methods
//
//=============================================================
//=============================================================
//
// Method for finding good Primary Vertex
//
//=============================================================
bool MakeZprimeMiniAodTreeMC::PrimaryVertex(const reco::VertexCollection &vertices)
{
  int nbGoodPv = 0;
  bool result = false;
  for(reco::VertexCollection::const_iterator it=vertices.begin() ; it!=vertices.end() ; ++it)
    {
      if(it->ndof() > minNdof_ &&
	 ( (maxAbsZ_ <= 0.0) || fabs(it->position().z()) <= maxAbsZ_ ) &&
	 ( (maxd0_ <= 0.0) || fabs(it->position().rho()) <= maxd0_ ) ) nbGoodPv++;
    }
  if( nbGoodPv>=NbGoodPv_ ) result = true;
  return result;
}
//=============================================================
//
//                Method for Pat Electron Tree
//
//=============================================================
void MakeZprimeMiniAodTreeMC::PatElectronTree(const edm::Event& iEvent,const edm::EventSetup& es)
{ 
  unsigned int NbElectrons = -1 ;
  Ele_nbElectrons.clear();
  Ele_dr03TkSumPt.clear();
  Ele_dr03TkSumPt_corrected.clear();
  Ele_Et.clear();
  Ele_EtFromCaloEn.clear();
  Ele_rawId.clear();
  Ele_isEcalDrivenSeed.clear();
  Ele_isPassConversionVeto.clear();
  Ele_charge.clear();
  Ele_nbOfMissingHits.clear();
  Ele_nbVtx.clear();
  Ele_pt.clear();
  Ele_thetaSC.clear();
  Ele_etaSC.clear();
  Ele_phiSC.clear();
  Ele_energySC.clear();
  Ele_preshowerEnergySC.clear();
  Ele_thetaTrack.clear();
  Ele_etaTrack.clear();
  Ele_phiTrack.clear();
  Ele_hadronicOverEm.clear();
  Ele_deltaEtaInSeedCluster.clear();
  Ele_deltaPhiInSeedCluster.clear();      
  Ele_deltaEtaInSC.clear();
  Ele_deltaPhiInSC.clear();
  Ele_sigmaIetaIeta.clear();
  Ele_e2x5Max.clear();
  Ele_e1x5.clear();
  Ele_frac51.clear();
  Ele_frac15.clear();
  Ele_e5x5.clear();
  Ele_e3x3.clear();
  Ele_e2x5MaxOver5x5.clear();
  Ele_e1x5Over5x5.clear();
  Ele_sigmaIetaIetaFull5x5.clear();
  Ele_e2x5MaxFull5x5.clear();
  Ele_e1x5Full5x5.clear();
  Ele_e5x5Full5x5.clear();
  Ele_e2x5MaxOver5x5Full5x5.clear();
  Ele_e1x5Over5x5Full5x5.clear();
  Ele_e2x5Right.clear();
  Ele_e2x5Left.clear();
  Ele_e2x5Top.clear();
  Ele_e2x5Bottom.clear();
  Ele_eMax.clear();
  Ele_eRight.clear();
  Ele_eLeft.clear();
  Ele_eTop.clear();
  Ele_eBottom.clear();
  Ele_dxy.clear();
  Ele_dz.clear();
  Ele_rhoIso.clear();
  Ele_fbrem.clear();
  Ele_EoverP.clear();
  Ele_Xposition.clear();
  Ele_Yposition.clear();
  Ele_EcalPlusHcald1iso.clear();
  Ele_dr03EcalRecHitSumEt.clear();
  Ele_dr03HcalDepth1TowerSumEt.clear();
  Ele_dr03HcalDepth1TowerSumEtBc.clear();
  Ele_hcalDepth1OverEcal.clear();
  Ele_hcalDepth2OverEcal.clear();
  Ele_dr03HcalDepth2TowerSumEt.clear();
  Ele_hcalDepth2TowerSumEtNoVeto.clear();
  Ele_hcalDepth1TowerSumEtNoVeto.clear();                       
  Ele_pfSumPhotonEt.clear();
  Ele_pfSumChargedHadronPt.clear();
  Ele_pfSumNeutralHadronEt.clear();
  Ele_pfSumPUPt.clear();
  Ele_pfDeltaBeta.clear();
  Ele_x.clear();
  Ele_y.clear();
  Ele_z.clear();
  Ele_zTrackPositionAtVtx.clear();
  Ele_ieta.clear();
  Ele_phiWidth.clear();
  Ele_etaWidth.clear();
  Ele_nrSatCrys.clear();
  Ele_isPassHeepID.clear();
  // rho for isolation
  edm::Handle<double> rhoIso_h;
  iEvent.getByToken(rhoIsoInputTag_, rhoIso_h);
  double rhoIso = *(rhoIso_h.product());
  // primary vertex candidate collection
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  const reco::Vertex &PV = vertices->front();
  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);
  //------------ Rec Hits in EB ----------------
  edm::Handle<EcalRecHitCollection> ecalEB;
  iEvent.getByToken(tok_EB_,ecalEB);
  //------------ Rec Hits in EE ----------------  
  edm::Handle<EcalRecHitCollection> ecalEE;
  iEvent.getByToken(tok_EE_,ecalEE);
  std::auto_ptr<EcalClusterLazyTools> lazyTools_;
  lazyTools_ .reset(new EcalClusterLazyTools( iEvent, es, ecalRechitEBToken_, ecalRechitEEToken_ ));  
  edm::Handle<edm::View<pat::Electron> > eleHandle;
  iEvent.getByToken(eleToken_,eleHandle);
  for(auto& el : *eleHandle){
    if(el.caloEnergy() * sin(el.p4().theta()) < 20) continue;
    NbElectrons++;
    Ele_nbElectrons.push_back(NbElectrons); 
    Ele_Et.push_back(el.superCluster()->energy() * sin(el.p4().theta()));
    Ele_EtFromCaloEn.push_back(el.caloEnergy() * sin(el.p4().theta()));
    //================================================================
    //
    // Begin New piece of code by Sam
    //
    //================================================================
    //access new tracker isolation
    const float trkIso = el.userFloat("trkPtIso");
    //std::cout<<"trkIso(corr) ="<<trkIso<<std::endl;
    //std::cout<<"trkIso ="<<ele.dr03TkSumPt()<<std::endl;
    Ele_dr03TkSumPt_corrected.push_back(trkIso);
    Ele_dr03TkSumPt.push_back(el.dr03TkSumPt());
    //access # saturated crystals in the 5x5
    const float nrSatCrys = el.userInt("nrSatCrys");
    Ele_nrSatCrys.push_back(nrSatCrys);
    //access the HEEP ID pass / fail
    const bool heepID = el.userInt("heepElectronID_HEEPV70");
    Ele_isPassHeepID.push_back(heepID);
    //================================================================
    //
    // End New piece of code by Sam
    //
    //================================================================
    Ele_pt.push_back(el.pt()); 
    Ele_thetaSC.push_back(el.caloPosition().theta()); //theta SC
    Ele_etaSC.push_back(el.superCluster()->eta());    //eta SC
    Ele_phiSC.push_back(el.superCluster()->phi());    //phi SC
    Ele_phiWidth.push_back(el.superCluster()->phiWidth()); 
    Ele_etaWidth.push_back(el.superCluster()->etaWidth()); 
    Ele_energySC.push_back(el.superCluster()->energy()); //energy SC
    Ele_preshowerEnergySC.push_back(el.superCluster()->preshowerEnergy()); 
    Ele_thetaTrack.push_back(el.p4().theta()); //theta track
    Ele_etaTrack.push_back(el.p4().eta());     //eta track
    Ele_phiTrack.push_back(el.p4().phi());     //phi track
    Ele_x.push_back(el.p4().x());
    Ele_y.push_back(el.p4().y());
    Ele_z.push_back(el.p4().z());
    Ele_zTrackPositionAtVtx.push_back(el.TrackPositionAtVtx().Z());
    Ele_hadronicOverEm.push_back(el.hadronicOverEm());
    Ele_deltaEtaInSC.push_back(el.deltaEtaSuperClusterTrackAtVtx());
    Ele_deltaPhiInSC.push_back(el.deltaPhiSuperClusterTrackAtVtx());
    Ele_deltaEtaInSeedCluster.push_back(el.deltaEtaSeedClusterTrackAtVtx());
    Ele_deltaPhiInSeedCluster.push_back(el.deltaPhiSeedClusterTrackAtCalo());
    Ele_sigmaIetaIeta.push_back(el.sigmaIetaIeta());
    Ele_e2x5Max.push_back(el.e2x5Max());
    Ele_e1x5.push_back(el.e1x5());
    //Ele_e5x1.push_back(el.e5x1());
    //Ele_e3x3.push_back(el.e3x3());
    Ele_e5x5.push_back(el.e5x5());
    Ele_e2x5MaxOver5x5.push_back(el.e2x5Max()/el.e5x5());
    Ele_e1x5Over5x5.push_back(el.e1x5()/el.e5x5());
    Ele_sigmaIetaIetaFull5x5.push_back(el.full5x5_sigmaIetaIeta()); 
    Ele_e2x5MaxFull5x5.push_back(el.full5x5_e2x5Max());
    Ele_e1x5Full5x5.push_back(el.full5x5_e1x5());
    Ele_e5x5Full5x5.push_back(el.full5x5_e5x5());
    Ele_e2x5MaxOver5x5Full5x5.push_back(el.full5x5_e2x5Max()/el.full5x5_e5x5());
    Ele_e1x5Over5x5Full5x5.push_back(el.full5x5_e1x5()/el.full5x5_e5x5());
    //Ele_rawId.push_back(el.superCluster()->seed()->seed().id());
    EBDetId BarrelId = el.superCluster()->seed()->seed();
    //std::cout << "ieta =  " << BarrelId.ieta() << std::endl;
    //std::cout << "id =  " << BarrelId.rawId() << std::endl;
    Ele_rawId.push_back(BarrelId.rawId());
    Ele_ieta.push_back(BarrelId.ieta());
    Ele_e2x5Right.push_back(lazyTools_->e2x5Right(*(el.superCluster()->seed())));
    Ele_e2x5Left.push_back(lazyTools_->e2x5Left(*(el.superCluster()->seed())));
    Ele_e2x5Top.push_back(lazyTools_->e2x5Top(*(el.superCluster()->seed())));
    Ele_e2x5Bottom.push_back(lazyTools_->e2x5Bottom(*(el.superCluster()->seed())));
    Ele_eMax.push_back(lazyTools_->eMax(*(el.superCluster()->seed())));
    Ele_eRight.push_back(lazyTools_->eRight(*(el.superCluster()->seed())));
    Ele_eLeft.push_back(lazyTools_->eLeft(*(el.superCluster()->seed())));
    Ele_eTop.push_back(lazyTools_->eTop(*(el.superCluster()->seed())));
    Ele_eBottom.push_back(lazyTools_->eBottom(*(el.superCluster()->seed())));
    Ele_e3x3.push_back(lazyTools_->e3x3(*(el.superCluster()->seed())));
    Ele_frac51.push_back( lazyTools_->e5x1(*(el.superCluster()->seed()))/el.full5x5_e5x5() );
    Ele_frac15.push_back( lazyTools_->e1x5(*(el.superCluster()->seed()))/el.full5x5_e5x5() );
    Ele_nbVtx.push_back(vertices->size());
    if(vertices->size()>0){
      Ele_dxy.push_back(el.gsfTrack()->dxy(PV.position()));   
      Ele_dz.push_back(el.gsfTrack()->dz(PV.position())); 
    }
    else{
      Ele_dxy.push_back(el.gsfTrack()->dxy());   
      Ele_dz.push_back(el.gsfTrack()->dz());
    }
    Ele_isEcalDrivenSeed.push_back(el.ecalDrivenSeed());
    Ele_isPassConversionVeto.push_back(el.passConversionVeto());
    Ele_charge.push_back(el.gsfTrack()->charge());
    Ele_rhoIso.push_back(rhoIso);
    Ele_nbOfMissingHits.push_back(el.gsfTrack()->numberOfLostHits()); 
    Ele_fbrem.push_back(el.fbrem());
    Ele_EoverP.push_back(el.eSeedClusterOverP());
    Ele_Xposition.push_back(el.caloPosition().x());   
    Ele_Yposition.push_back(el.caloPosition().y()); 
    //------------- detector isolation -------------------------
    Ele_hcalDepth1OverEcal.push_back(el.hcalDepth1OverEcal());
    Ele_hcalDepth2OverEcal.push_back(el.hcalDepth2OverEcal());
    Ele_dr03HcalDepth2TowerSumEt.push_back(el.dr03HcalDepth2TowerSumEt());
    Ele_hcalDepth2TowerSumEtNoVeto.push_back(el.isolationVariables03().hcalDepth2TowerSumEt);// hcaldepht2 iso deposit with 
    // electron footprint removed
    Ele_hcalDepth1TowerSumEtNoVeto.push_back(el.isolationVariables03().hcalDepth1TowerSumEt);// hcaldepht1 iso deposit with 
    // electron footprint removed
    Ele_EcalPlusHcald1iso.push_back(el.dr03EcalRecHitSumEt() + el.dr03HcalDepth1TowerSumEt());
    Ele_dr03EcalRecHitSumEt.push_back(el.dr03EcalRecHitSumEt());
    Ele_dr03HcalDepth1TowerSumEt.push_back(el.dr03HcalDepth1TowerSumEt());
    Ele_dr03HcalDepth1TowerSumEtBc.push_back(el.dr03HcalDepth1TowerSumEtBc());
    //------------- PF isolation from pat::ele -------------------------
    Ele_pfSumPhotonEt.push_back(el.pfIsolationVariables().sumPhotonEt);
    Ele_pfSumChargedHadronPt.push_back(el.pfIsolationVariables().sumChargedHadronPt); 
    Ele_pfSumNeutralHadronEt.push_back(el.pfIsolationVariables().sumNeutralHadronEt);
    Ele_pfSumPUPt.push_back(el.pfIsolationVariables().sumPUPt);  
    // do deltaBeta
    double charged   = el.pfIsolationVariables().sumPhotonEt;
    double neutral   = el.pfIsolationVariables().sumNeutralHadronEt;
    double pileup    = el.pfIsolationVariables().sumPUPt;
    double deltaBeta = charged + std::max(0.0, neutral-0.5*pileup);
    Ele_pfDeltaBeta.push_back(deltaBeta);
  }    
}
//=============================================================
//
//            Method for Trigger Matching Tree
//
//=============================================================
void MakeZprimeMiniAodTreeMC::TriggerMatchingTree(const edm::Event& iEvent,const edm::EventSetup& es)
{ 
  int NbTriggers = 0;
  int NbTriggerObj = 0;
  HLT_nb.clear();
  HLT_name.clear();
  HLT_isaccept.clear();
  HLTObj_nbObj.clear();
  HLTObj_pt.clear();
  HLTObj_eta.clear();
  HLTObj_phi.clear();
  HLTObj_collection.clear();
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::TriggerObjectStandAloneCollection> trigobj_handle;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  std::string full_name;
  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  iEvent.getByToken(triggerObjects_, trigobj_handle);
  iEvent.getByToken(triggerPrescales_, triggerPrescales);
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  //std::cout << "\n === TRIGGER PATHS === " << std::endl;
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    if( names.triggerName(i) != "HLT_Mu27_v1" && names.triggerName(i) != "HLT_Mu27_v2"  &&
        names.triggerName(i) != "HLT_Mu27_v3" && names.triggerName(i) != "HLT_Mu27_v4"  &&
        names.triggerName(i) != "HLT_Mu27_v5" && names.triggerName(i) != "HLT_Mu27_v6"  &&
        names.triggerName(i) != "HLT_Mu27_v7" && names.triggerName(i) != "HLT_Mu27_v8"  &&
        names.triggerName(i) != "HLT_Mu27_v9" && names.triggerName(i) != "HLT_Mu27_v10" &&
        names.triggerName(i) != "HLT_TkMu27_v1" && names.triggerName(i) != "HLT_TkMu27_v2"  &&
        names.triggerName(i) != "HLT_TkMu27_v3" && names.triggerName(i) != "HLT_TkMu27_v4"  &&
        names.triggerName(i) != "HLT_TkMu27_v5" && names.triggerName(i) != "HLT_TkMu27_v6"  &&
        names.triggerName(i) != "HLT_TkMu27_v7" && names.triggerName(i) != "HLT_TkMu27_v8"  &&
        names.triggerName(i) != "HLT_TkMu27_v9" && names.triggerName(i) != "HLT_TkMu27_v10" &&
        names.triggerName(i) != "HLT_Mu50_v1" && names.triggerName(i) != "HLT_Mu50_v2"  &&
        names.triggerName(i) != "HLT_Mu50_v3" && names.triggerName(i) != "HLT_Mu50_v4"  &&
	names.triggerName(i) != "HLT_Mu50_v5" && names.triggerName(i) != "HLT_Mu50_v6"  &&
        names.triggerName(i) != "HLT_Mu50_v7" && names.triggerName(i) != "HLT_Mu50_v8"  &&
        names.triggerName(i) != "HLT_Mu50_v9" && names.triggerName(i) != "HLT_Mu50_v10" &&
        names.triggerName(i) != "HLT_TkMu50_v1" && names.triggerName(i) != "HLT_TkMu50_v2"  &&
        names.triggerName(i) != "HLT_TkMu50_v3" && names.triggerName(i) != "HLT_TkMu50_v4"  &&
        names.triggerName(i) != "HLT_TkMu50_v5" && names.triggerName(i) != "HLT_TkMu50_v6"  &&
        names.triggerName(i) != "HLT_TkMu50_v7" && names.triggerName(i) != "HLT_TkMu50_v8"  &&
        names.triggerName(i) != "HLT_TkMu50_v9" && names.triggerName(i) != "HLT_TkMu50_v10" && 
	names.triggerName(i) != "HLT_DoubleEle33_CaloIdL_MW_v1" && names.triggerName(i) != "HLT_DoubleEle33_CaloIdL_MW_v2"  &&
        names.triggerName(i) != "HLT_DoubleEle33_CaloIdL_MW_v3" && names.triggerName(i) != "HLT_DoubleEle33_CaloIdL_MW_v4"  &&
        names.triggerName(i) != "HLT_DoubleEle33_CaloIdL_MW_v5" && names.triggerName(i) != "HLT_DoubleEle33_CaloIdL_MW_v6"  &&
        names.triggerName(i) != "HLT_DoubleEle33_CaloIdL_MW_v7" && names.triggerName(i) != "HLT_DoubleEle33_CaloIdL_MW_v8"  &&
        names.triggerName(i) != "HLT_DoubleEle33_CaloIdL_MW_v9" && names.triggerName(i) != "HLT_DoubleEle33_CaloIdL_MW_v10" &&
	names.triggerName(i) != "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v1" && names.triggerName(i) != "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v2"  &&
        names.triggerName(i) != "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3" && names.triggerName(i) != "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v4"  &&
        names.triggerName(i) != "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v5" && names.triggerName(i) != "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v6"  &&
        names.triggerName(i) != "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v7" && names.triggerName(i) != "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v8"  &&
        names.triggerName(i) != "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v9" && names.triggerName(i) != "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v10" &&
	names.triggerName(i) != "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v1" && names.triggerName(i) != "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v2"  &&
        names.triggerName(i) != "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v3" && names.triggerName(i) != "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v4"  &&
        names.triggerName(i) != "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v5" && names.triggerName(i) != "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v6"  &&
        names.triggerName(i) != "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v7" && names.triggerName(i) != "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v8"  &&
        names.triggerName(i) != "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v9" && names.triggerName(i) != "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v10" &&
	names.triggerName(i) != "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v1" && names.triggerName(i) != "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v2"  &&
        names.triggerName(i) != "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v3" && names.triggerName(i) != "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v4"  &&
        names.triggerName(i) != "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v5" && names.triggerName(i) != "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v6"  &&
        names.triggerName(i) != "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v7" && names.triggerName(i) != "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v8"  &&
        names.triggerName(i) != "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v9" && names.triggerName(i) != "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v10" 
        ) continue;
    //std::string const& name = names.triggerName(i);
    //full_name = name;
    NbTriggers++;
    HLT_nb.push_back(NbTriggers);
    HLT_name.push_back(names.triggerName(i));   
    HLT_isaccept.push_back(triggerBits->accept(i));
    /*std::cout << "Trigger " << names.triggerName(i)  << 
      ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
      ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)") 
      << std::endl;*/
  }
  for (unsigned i = 0; i < trigobj_handle->size(); ++i) {
    pat::TriggerObjectStandAlone src = trigobj_handle->at(i);
    src.unpackPathNames(names);
    std::vector<std::string> const& pathnames = src.pathNames();
    for (unsigned j = 0; j < pathnames.size(); ++j) {
      //std::cout<<" full_name = "<< pathnames[j] << endl;
      //printf("pt = %f eta = %f phi = %f \n",src.pt(),src.eta(),src.phi());
      if( src.pt() < 20.0 ) continue;
      if( fabs(src.eta()) > 3.0 ) continue;
      if( pathnames[j] != "HLT_Mu27_v1" && pathnames[j] != "HLT_Mu27_v2"  &&
          pathnames[j] != "HLT_Mu27_v3" && pathnames[j] != "HLT_Mu27_v4"  &&
          pathnames[j] != "HLT_Mu27_v5" && pathnames[j] != "HLT_Mu27_v6"  &&
          pathnames[j] != "HLT_Mu27_v7" && pathnames[j] != "HLT_Mu27_v8"  &&
          pathnames[j] != "HLT_Mu27_v9" && pathnames[j] != "HLT_Mu27_v10" &&
          pathnames[j] != "HLT_TkMu27_v1" && pathnames[j] != "HLT_TkMu27_v2"  &&
          pathnames[j] != "HLT_TkMu27_v3" && pathnames[j] != "HLT_TkMu27_v4"  &&
          pathnames[j] != "HLT_TkMu27_v5" && pathnames[j] != "HLT_TkMu27_v6"  &&
          pathnames[j] != "HLT_TkMu27_v7" && pathnames[j] != "HLT_TkMu27_v8"  &&
          pathnames[j] != "HLT_TkMu27_v9" && pathnames[j] != "HLT_TkMu27_v10" &&
          pathnames[j] != "HLT_Mu50_v1" && pathnames[j] != "HLT_Mu50_v2"  &&
          pathnames[j] != "HLT_Mu50_v3" && pathnames[j] != "HLT_Mu50_v4"  &&
          pathnames[j] != "HLT_Mu50_v5" && pathnames[j] != "HLT_Mu50_v6"  &&
          pathnames[j] != "HLT_Mu50_v7" && pathnames[j] != "HLT_Mu50_v8"  &&
          pathnames[j] != "HLT_Mu50_v9" && pathnames[j] != "HLT_Mu50_v10" &&
          pathnames[j] != "HLT_TkMu50_v1" && pathnames[j] != "HLT_TkMu50_v2"  &&
          pathnames[j] != "HLT_TkMu50_v3" && pathnames[j] != "HLT_TkMu50_v4"  &&
          pathnames[j] != "HLT_TkMu50_v5" && pathnames[j] != "HLT_TkMu50_v6"  &&
          pathnames[j] != "HLT_TkMu50_v7" && pathnames[j] != "HLT_TkMu50_v8"  &&
          pathnames[j] != "HLT_TkMu50_v9" && pathnames[j] != "HLT_TkMu50_v10" && 
	  pathnames[j] != "HLT_DoubleEle33_CaloIdL_MW_v1" && pathnames[j] != "HLT_DoubleEle33_CaloIdL_MW_v2"  &&
          pathnames[j] != "HLT_DoubleEle33_CaloIdL_MW_v3" && pathnames[j] != "HLT_DoubleEle33_CaloIdL_MW_v4"  &&
          pathnames[j] != "HLT_DoubleEle33_CaloIdL_MW_v5" && pathnames[j] != "HLT_DoubleEle33_CaloIdL_MW_v6"  &&
          pathnames[j] != "HLT_DoubleEle33_CaloIdL_MW_v7" && pathnames[j] != "HLT_DoubleEle33_CaloIdL_MW_v8"  &&
          pathnames[j] != "HLT_DoubleEle33_CaloIdL_MW_v9" && pathnames[j] != "HLT_DoubleEle33_CaloIdL_MW_v10" &&
	  pathnames[j] != "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v1" && pathnames[j] != "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v2"  &&
          pathnames[j] != "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3" && pathnames[j] != "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v4"  &&
          pathnames[j] != "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v5" && pathnames[j] != "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v6"  &&
          pathnames[j] != "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v7" && pathnames[j] != "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v8"  &&
          pathnames[j] != "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v9" && pathnames[j] != "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v10" &&
	  pathnames[j] != "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v1" && pathnames[j] != "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v2"  &&
          pathnames[j] != "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v3" && pathnames[j] != "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v4"  &&
          pathnames[j] != "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v5" && pathnames[j] != "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v6"  &&
          pathnames[j] != "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v7" && pathnames[j] != "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v8"  &&
          pathnames[j] != "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v9" && pathnames[j] != "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v10" &&
	  pathnames[j] != "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v1" && pathnames[j] != "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v2"  &&
          pathnames[j] != "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v3" && pathnames[j] != "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v4"  &&
          pathnames[j] != "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v5" && pathnames[j] != "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v6"  &&
          pathnames[j] != "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v7" && pathnames[j] != "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v8"  &&
          pathnames[j] != "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v9" && pathnames[j] != "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v10"
          ) continue;
      HLTObj_nbObj.push_back(NbTriggerObj);
      HLTObj_pt.push_back(src.pt()); 
      HLTObj_eta.push_back(src.eta());
      HLTObj_phi.push_back(src.phi());
      HLTObj_collection.push_back(pathnames[j]);
    }
  }
}

//=============================================================
//
//     Method for Generated info the variables
//
//=============================================================
//Check recursively if any ancestor of particle is the given one
bool MakeZprimeMiniAodTreeMC::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
{
  //particle is already the ancestor
  if(ancestor == particle ) return true;
  //otherwise loop on mothers, if any and return true if the ancestor is found
  for(size_t i=0;i< particle->numberOfMothers();i++)
    {
      if(isAncestor(ancestor,particle->mother(i))) return true;
    }
  //if we did not return yet, then particle and ancestor are not relatives
  return false;
}

void MakeZprimeMiniAodTreeMC::accessGenInfo(const edm::Event& iEvent,const edm::EventSetup& es)
{ 
  int NbGenMuons  = 0;
  iGen.clear();
  idGen.clear();
  statusGen.clear();
  ptGen.clear();
  etaGen.clear();
  phiGen.clear();
  chargeGen.clear();
  EnergyGen.clear();
  pxGen.clear();
  pyGen.clear();
  pzGen.clear();
  // Pruned particles are the one containing "important" stuff
  Handle<edm::View<reco::GenParticle> > pruned;
  iEvent.getByToken(prunedGenToken_,pruned);
  if (!(pruned.isValid())) return; 
  // Packed particles are all the status 1, so usable to remake jets
  // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
  Handle<edm::View<pat::PackedGenParticle> > packed;
  iEvent.getByToken(packedGenToken_,packed);
  //let's try to find all status1 originating directly from a B meson decay 
  for(size_t i=0; i<pruned->size();i++){
    //if(abs((*pruned)[i].pdgId()) > 22 && abs((*pruned)[i].pdgId()) <24){
    if(abs((*pruned)[i].pdgId()) < BosonID_){
      const Candidate * Zprime = &(*pruned)[i];
      //std::cout << "PdgID: " << Zprime->pdgId() << " pt " << Zprime->pt() << " eta: " << Zprime->eta() << " phi: " << Zprime->phi() << std::endl;
      //std::cout << "  found daugthers: " << std::endl;
      for(size_t j=0; j<packed->size();j++){
	//get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection 
	const Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
	if(motherInPrunedCollection != nullptr && isAncestor( Zprime , motherInPrunedCollection)){
	  if( (*packed)[j].pt() < 20.0 ) continue;
	  if( fabs((*packed)[j].pdgId()) != ParticleID1_ && 
	      fabs((*packed)[j].pdgId()) != ParticleID2_ && 
	      fabs((*packed)[j].pdgId()) != ParticleID3_ ) continue;
	  if( (*packed)[j].status() > ParticleStatus_ ) continue;
	  NbGenMuons++;
	  iGen.push_back(NbGenMuons);
	  idGen.push_back((*packed)[j].pdgId());
	  statusGen.push_back((*packed)[j].status());
	  ptGen.push_back((*packed)[j].pt());
	  etaGen.push_back((*packed)[j].eta());
	  phiGen.push_back((*packed)[j].phi());
	  chargeGen.push_back((*packed)[j].charge());
	  EnergyGen.push_back((*packed)[j].energy());
	  pxGen.push_back((*packed)[j].px());
	  pyGen.push_back((*packed)[j].py());
	  pzGen.push_back((*packed)[j].pz());
	}
      }
    }
  }
}
//=============================================================
//
//                Method for Pat Muon Tree
//
//=============================================================
void MakeZprimeMiniAodTreeMC::PatMuonTree(const edm::Event& evt,const edm::EventSetup& es)
{  
  int NbMuons = 0;
  Mu_nbMuon.clear();
  Mu_isTightMuon.clear();
  Mu_isLooseMuon.clear();
  Mu_isGlobalMuon.clear();
  Mu_isPF.clear();
  Mu_isTrackerMuon.clear();
  Mu_isHighPtMuon.clear();
  Mu_en.clear();
  Mu_pt.clear();
  Mu_eta.clear();
  Mu_phi.clear();
  Mu_et.clear();
  Mu_charge.clear();
  Mu_numberOfMatchedStations.clear();
  Mu_trackiso.clear();
  Mu_pfSumChargedHadronPt.clear();
  Mu_pfSumNeutralHadronEt.clear();
  Mu_PFSumPhotonEt.clear();
  Mu_pfSumPUPt.clear();
  Mu_normalizedChi2.clear();
  Mu_numberOfValidPixelHits.clear();
  Mu_numberOftrackerLayersWithMeasurement.clear();
  Mu_numberOfValidMuonHits.clear();
  Mu_calEnergy.clear();
  Mu_emIso.clear();
  Mu_hadIso.clear();
  Mu_nbofpv.clear();
  Mu_ptTunePMuonBestTrack.clear();
  Mu_dPToverPTTunePMuonBestTrack.clear();
  Mu_pxTunePMuonBestTrack.clear();
  Mu_pyTunePMuonBestTrack.clear();
  Mu_pzTunePMuonBestTrack.clear();
  Mu_pTunePMuonBestTrack.clear();
  Mu_etaTunePMuonBestTrack.clear();
  Mu_phiTunePMuonBestTrack.clear();
  Mu_thetaTunePMuonBestTrack.clear();
  Mu_chargeTunePMuonBestTrack.clear();
  Mu_absdxyTunePMuonBestTrack.clear();
  Mu_absdzTunePMuonBestTrack.clear();  
  Mu_ptInnerTrack.clear();
  Mu_dPToverPTInnerTrack.clear();
  Mu_pxInnerTrack.clear();
  Mu_pyInnerTrack.clear();
  Mu_pzInnerTrack.clear();
  Mu_pInnerTrack.clear();
  Mu_etaInnerTrack.clear();
  Mu_phiInnerTrack.clear();
  Mu_thetaInnerTrack.clear();
  Mu_chargeInnerTrack.clear();
  Mu_absdxyInnerTrack.clear();
  Mu_absdzInnerTrack.clear();
  Mu_absdxy.clear();
  Mu_absdz.clear();
  Mu_patDeltaBeta.clear();
  Mu_passNewMatchedStationsCut.clear();
  Mu_passOldMatchedStationsCut.clear();
  Mu_stationMask.clear();
  Mu_numberOfMatchedRPCLayers.clear();
  //Mu_isMuonsCleaned.clear();
  // Get TransientTracks (for use in e.g. the vertex fit) for each of                   
  // the muon tracks, using e.g. the cocktail momentum.
  //edm::ESHandle<TransientTrackBuilder> ttkb;
  //es.get<TransientTrackRecord>().get("TransientTrackBuilder", ttkb);
  //std::vector<reco::TransientTrack> ttv;
  // primary vertex candidate collection
  edm::Handle<reco::VertexCollection> vertices;
  evt.getByToken(vtxToken_,vertices);
  const reco::Vertex &PV = vertices->front();
  // pat candidate collection
  edm::Handle<pat::MuonCollection> muons;
  evt.getByToken(muonToken_, muons);
  for (const pat::Muon &mu : *muons) {
    if( mu.pt() < 20.0 ) continue;
    //============= Parameters related to matched Gen info =====================
    if( !mu.innerTrack().isNonnull() ) continue;
    if( !mu.globalTrack().isNonnull() ) continue;
    //if( mu.isTrackerMuon()==false ) continue;
    NbMuons++;
    Mu_nbMuon.push_back(NbMuons);
    Mu_isLooseMuon.push_back(mu.isLooseMuon());
    Mu_isTightMuon.push_back(mu.isTightMuon(PV));
    Mu_isHighPtMuon.push_back(mu.isHighPtMuon(PV));
    Mu_isGlobalMuon.push_back(mu.isGlobalMuon());
    Mu_isPF.push_back(mu.isPFMuon());
    Mu_isTrackerMuon.push_back(mu.isTrackerMuon());
    Mu_nbofpv.push_back(vertices->size());
    //Mu_isMuonsCleaned.push_back(mu.userInt("muonsCleaned:oldPF"));
    //============== Parameters related to Kinematics =====================
    /*reco::TrackRef cktTrack;
    cktTrack = (muon::tevOptimized(mu, 200, 17., 40., 0.25)).first;
    Mu_ptcocktail.push_back(cktTrack->pt());
    Mu_dPToverPTcocktail.push_back(cktTrack->ptError()/cktTrack->pt());
    Mu_absdxy.push_back(fabs(cktTrack->dxy(PV.position())));
    Mu_absdz.push_back(fabs(cktTrack->dz(PV.position())));
    Mu_etaCocktail.push_back(cktTrack->eta());
    Mu_phiCocktail.push_back(cktTrack->phi());
    Mu_thetaCocktail.push_back(cktTrack->theta());
    Mu_chargeCocktail.push_back(cktTrack->charge());
    Mu_pxCocktail.push_back(cktTrack->px()); //px component of the track
    Mu_pyCocktail.push_back(cktTrack->py()); //py component of the track
    Mu_pzCocktail.push_back(cktTrack->pz()); //pz component of the track
    Mu_pCocktail.push_back(cktTrack->p()); //magnitude of momentum vector*/
    // part for TuneP Muon Best track
    const reco::TrackRef& tunePTrack  = mu.tunePMuonBestTrack();
    Mu_ptTunePMuonBestTrack.push_back(tunePTrack->pt());
    Mu_dPToverPTTunePMuonBestTrack.push_back(tunePTrack->ptError()/tunePTrack->pt());
    Mu_pxTunePMuonBestTrack.push_back(tunePTrack->px()); //px component of the track 
    Mu_pyTunePMuonBestTrack.push_back(tunePTrack->py()); //py component of the track 
    Mu_pzTunePMuonBestTrack.push_back(tunePTrack->pz()); //pz component of the track 
    Mu_pTunePMuonBestTrack.push_back(tunePTrack->p());   //magnitude of momentum vector 
    Mu_etaTunePMuonBestTrack.push_back(tunePTrack->eta());
    Mu_phiTunePMuonBestTrack.push_back(tunePTrack->phi());
    Mu_thetaTunePMuonBestTrack.push_back(tunePTrack->theta());
    Mu_chargeTunePMuonBestTrack.push_back(tunePTrack->charge());
    Mu_absdxyTunePMuonBestTrack.push_back(fabs(tunePTrack->dxy(PV.position())));
    Mu_absdzTunePMuonBestTrack.push_back(fabs(tunePTrack->dz(PV.position())));
    Mu_en.push_back(mu.energy());
    Mu_et.push_back(mu.et());
    Mu_pt.push_back(mu.pt());
    Mu_eta.push_back(mu.eta());
    Mu_phi.push_back(mu.phi());
    Mu_charge.push_back(mu.charge());
    Mu_ptInnerTrack.push_back(mu.innerTrack()->pt());    
    Mu_dPToverPTInnerTrack.push_back(mu.innerTrack()->ptError()/mu.innerTrack()->pt());
    Mu_pxInnerTrack.push_back(mu.innerTrack()->px()); 
    Mu_pyInnerTrack.push_back(mu.innerTrack()->py()); 
    Mu_pzInnerTrack.push_back(mu.innerTrack()->pz()); 
    Mu_pInnerTrack.push_back(mu.innerTrack()->p());   
    Mu_etaInnerTrack.push_back(mu.innerTrack()->eta());
    Mu_phiInnerTrack.push_back(mu.innerTrack()->phi());
    Mu_thetaInnerTrack.push_back(mu.innerTrack()->theta());
    Mu_chargeInnerTrack.push_back(mu.innerTrack()->charge());
    Mu_absdxyInnerTrack.push_back(fabs(mu.innerTrack()->dxy(PV.position())));
    Mu_absdzInnerTrack.push_back(fabs(mu.innerTrack()->dz(PV.position())));
    Mu_absdxy.push_back(fabs(mu.globalTrack()->dxy(PV.position())));
    Mu_absdz.push_back(fabs(mu.globalTrack()->dz(PV.position())));
    //====================== Parameters related to track quality =====================
    Mu_normalizedChi2.push_back(mu.globalTrack()->normalizedChi2());
    Mu_numberOfValidPixelHits.push_back(mu.globalTrack()->hitPattern().numberOfValidPixelHits());
    Mu_numberOfValidMuonHits.push_back(mu.globalTrack()->hitPattern().numberOfValidMuonHits());
    Mu_numberOftrackerLayersWithMeasurement.push_back(mu.globalTrack()->hitPattern().trackerLayersWithMeasurement());
    Mu_numberOfMatchedStations.push_back(mu.numberOfMatchedStations());
    //============= Parameters related to detector isolation =====================
    Mu_emIso.push_back(mu.isolationR03().emEt);
    Mu_hadIso.push_back(mu.isolationR03().hadEt);
    Mu_trackiso.push_back(mu.isolationR03().sumPt);
    //============= Parameters related to PF isolation =====================
    Mu_pfSumChargedHadronPt.push_back(mu.pfIsolationR03().sumChargedHadronPt);
    Mu_pfSumNeutralHadronEt.push_back(mu.pfIsolationR03().sumNeutralHadronEt);
    Mu_PFSumPhotonEt.push_back(mu.pfIsolationR03().sumPhotonEt);
    Mu_pfSumPUPt.push_back(mu.pfIsolationR03().sumPUPt);
    // do deltaBeta
    double charged   = mu.pfIsolationR03().sumChargedHadronPt;
    double neutral   = mu.pfIsolationR03().sumNeutralHadronEt;
    double pileup    = mu.pfIsolationR03().sumPUPt;
    double deltaBeta = (charged + std::max(0.0, neutral-0.5*pileup))/mu.pt();
    Mu_patDeltaBeta.push_back(deltaBeta);
    Mu_stationMask.push_back(mu.stationMask());
    Mu_numberOfMatchedRPCLayers.push_back(mu.numberOfMatchedRPCLayers());
    //New Proposed Muon Matched Stations Cut*:
    if(mu.numberOfMatchedStations()>1) {Mu_passOldMatchedStationsCut.push_back(true);}
    else {Mu_passOldMatchedStationsCut.push_back(false);}
    if (mu.numberOfMatchedStations()>1
	|| (mu.numberOfMatchedStations()==1 && !(mu.stationMask()==1 || mu.stationMask()==16))
	|| ((mu.numberOfMatchedStations()==1 && (mu.stationMask()==1 || mu.stationMask()==16)) && mu.numberOfMatchedRPCLayers()>2))  {Mu_passNewMatchedStationsCut.push_back(true);}
    else {Mu_passNewMatchedStationsCut.push_back(false);}
  }
}

void MakeZprimeMiniAodTreeMC::ComputeMuonMassVtx(const edm::Event& evt,const edm::EventSetup& es)
{ 
  int NbMuons1 = 0;
  int NbMuons2 = 0;
  int NbMuons3 = 0;
  Mu_vtxMass.clear();
  Mu_vtxNormChi2.clear();
  Mu_vtxMassLept.clear();
  float vtxNormChi1,DiMass1;
  float vtxNormChi2,DiMass2;
  float vtxNormChi3,DiMass3;
  // Get TransientTracks (for use in e.g. the vertex fit) for each of
  // the muon tracks, using e.g. the cocktail momentum.
  edm::ESHandle<TransientTrackBuilder> ttkb1;
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", ttkb1);
  std::vector<reco::TransientTrack> ttv1;
  edm::ESHandle<TransientTrackBuilder> ttkb2;
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", ttkb2);
  std::vector<reco::TransientTrack> ttv2;
  edm::ESHandle<TransientTrackBuilder> ttkb3;
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", ttkb3);
  std::vector<reco::TransientTrack> ttv3;
  //find the first high pt muon
  reco::TrackRef MuonBestTrack1;
  reco::TrackRef MuonBestTrack2;
  reco::TrackRef MuonBestTrack3;
  // primary vertex candidate collection
  edm::Handle<reco::VertexCollection> vertices;
  evt.getByToken(vtxToken_,vertices);
  const reco::Vertex &vertex = vertices->front();
  // pat candidate collection
  edm::Handle<pat::MuonCollection> muons;
  evt.getByToken(muonToken_, muons);
  for (const pat::Muon &mu : *muons) {
    if( mu.pt() < 5.0 ) continue;
    if( !mu.globalTrack().isNonnull() ) continue;
    if( mu.isTrackerMuon()==false ) continue;
    MuonBestTrack1 = mu.tunePMuonBestTrack();
    if( MuonBestTrack1->pt()>53.0 && (mu.isolationR03().sumPt/mu.innerTrack()->pt()<0.10) 
    	&& (MuonBestTrack1->ptError()/MuonBestTrack1->pt()<0.3) && 
        fabs(MuonBestTrack1->dxy(vertex.position()))<0.2 &&
        mu.globalTrack()->hitPattern().trackerLayersWithMeasurement()>5 && 
        mu.globalTrack()->hitPattern().numberOfValidPixelHits()>0 &&
        mu.globalTrack()->hitPattern().numberOfValidMuonHits()>0 &&
        (mu.numberOfMatchedStations()>1
	 || (mu.numberOfMatchedStations()==1 && !(mu.stationMask()==1 || mu.stationMask()==16))
	 || ((mu.numberOfMatchedStations()==1 && (mu.stationMask()==1 || mu.stationMask()==16)) && mu.numberOfMatchedRPCLayers()>2)) ){
      //cout << "PT 1mu= " << MuonBestTrack1->pt()  << " and charge =" <<  MuonBestTrack1->charge() << endl;
      NbMuons1++;
      if(NbMuons1>1) continue;
      //cout<<"==================================================================="<<endl;
      //printf ("pt1 = %f eta1 = %f phi1 = %f charge1 = %d \n",MuonBestTrack1->pt(),MuonBestTrack1->eta(),MuonBestTrack1->phi(),MuonBestTrack1->charge());
      ttv1.push_back(ttkb1->build(MuonBestTrack1));
      //find the second high pt muon
      for (const pat::Muon &mu2 : *muons) {
	if( mu2.pt() < 5.0 ) continue;
	if( !mu2.globalTrack().isNonnull() ) continue;
	if( mu2.isTrackerMuon()==false ) continue;
	MuonBestTrack2 = mu2.tunePMuonBestTrack();
	//cout << "PT 2mu= " << MuonBestTrack2->pt()  << " and charge =" <<  MuonBestTrack2->charge()<< endl;
        if(MuonBestTrack2->pt() == MuonBestTrack1->pt()) continue;
	if( MuonBestTrack2->pt()>53.0 && (mu2.isolationR03().sumPt/mu2.innerTrack()->pt()<0.10) 
	    && (MuonBestTrack2->ptError()/MuonBestTrack2->pt()<0.3) && 
	    fabs(MuonBestTrack2->dxy(vertex.position()))<0.2 &&
	    mu2.globalTrack()->hitPattern().trackerLayersWithMeasurement()>5 && 
	    mu2.globalTrack()->hitPattern().numberOfValidPixelHits()>0 &&
	    mu2.globalTrack()->hitPattern().numberOfValidMuonHits()>0 &&
	    (mu2.numberOfMatchedStations()>1
	     || (mu2.numberOfMatchedStations()==1 && !(mu2.stationMask()==1 || mu2.stationMask()==16))
	     || ((mu2.numberOfMatchedStations()==1 && (mu2.stationMask()==1 || mu2.stationMask()==16)) && mu2.numberOfMatchedRPCLayers()>2)) ){
	  //cout << "PT 2mu= " << MuonBestTrack2->pt()  << " and charge =" <<  MuonBestTrack2->charge() << endl;
	  NbMuons2++;
	  if(NbMuons2>1) continue;
	  ttv1.push_back(ttkb1->build(MuonBestTrack2));
	  //printf ("pt2 = %f eta2 = %f phi2 = %f charge2 = %d \n",MuonBestTrack2->pt(),MuonBestTrack2->eta(),MuonBestTrack2->phi(),MuonBestTrack2->charge());
	  KalmanVertexFitter kvf(true); // false means no smoothing which means no track re-fit
	  CachingVertex<5> vtx1 = kvf.vertex(ttv1);
	  vtxNormChi1 = vtx1.totalChiSquared()/vtx1.degreesOfFreedom();
	  InvariantMassFromVertex imfv1;
	  static const double muon_mass1 = 0.1056583;
	  Measurement1D mass1 = imfv1.invariantMass(vtx1, muon_mass1);
	  DiMass1 = mass1.value();
	  Mu_vtxMass.push_back(DiMass1);
	  Mu_vtxNormChi2.push_back(vtxNormChi1);
	  Mu_vtxMassLept.push_back(MuonBestTrack1->pt());
	  Mu_vtxMassLept.push_back(MuonBestTrack2->pt());
	  //cout<<"size(inside) = "<<ttv1.size()<<endl;
	  //printf ("normChi2 = %f mass = %f\n",vtxNormChi1,DiMass1);
	  /*cout << "PT 1mu= " << MuonBestTrack1->pt()  << " and charge =" <<  MuonBestTrack1->charge() << endl;
	    cout << "PT 2mu= " << MuonBestTrack2->pt()  << " and charge =" <<  MuonBestTrack2->charge() << endl;
	    cout << "normChi2 = mass = " << vtxNormChi1 << " " << DiMass1 << endl;*/
	  //find the third high pt muon
	  for (const pat::Muon &mu3 : *muons) {
	    if( mu3.pt() < 5.0 ) continue;
	    if( !mu3.globalTrack().isNonnull() ) continue;
	    if( mu3.isTrackerMuon()==false ) continue;
	    MuonBestTrack3 = mu3.tunePMuonBestTrack();
	    //cout << "PT 3mu= " << MuonBestTrack3->pt()  << " and charge =" <<  MuonBestTrack3->charge() << endl;
	    if(MuonBestTrack3->pt() == MuonBestTrack1->pt()) continue;	   
	    if(MuonBestTrack3->pt() == MuonBestTrack2->pt()) continue;	    
	    if( MuonBestTrack3->pt()>53.0 && (mu3.isolationR03().sumPt/mu3.innerTrack()->pt()<0.10) 
		&& (MuonBestTrack3->ptError()/MuonBestTrack3->pt()<0.3) && 
		fabs(MuonBestTrack3->dxy(vertex.position()))<0.2 &&
		mu3.globalTrack()->hitPattern().trackerLayersWithMeasurement()>5 && 
		mu3.globalTrack()->hitPattern().numberOfValidPixelHits()>0 &&
		mu3.globalTrack()->hitPattern().numberOfValidMuonHits()>0 &&
		(mu3.numberOfMatchedStations()>1
		 || (mu3.numberOfMatchedStations()==1 && !(mu3.stationMask()==1 || mu3.stationMask()==16))
		 || ((mu3.numberOfMatchedStations()==1 && (mu3.stationMask()==1 || mu3.stationMask()==16)) && mu3.numberOfMatchedRPCLayers()>2)) ){
	      //cout << "PT 3mu= " << MuonBestTrack3->pt()  << " and charge =" <<  MuonBestTrack3->charge() << endl;
	      NbMuons3++;
	      if(NbMuons3==0) continue;
	      //printf ("pt3 = %f eta3 = %f phi3 = %f charge3 = %d \n",MuonBestTrack3->pt(),MuonBestTrack3->eta(),MuonBestTrack3->phi(),MuonBestTrack3->charge());
	      //cout<<"NbMuons3 = "<<NbMuons3<<endl;
	      ttv2.push_back(ttkb2->build(MuonBestTrack1));
	      ttv2.push_back(ttkb2->build(MuonBestTrack3));
	      ttv3.push_back(ttkb3->build(MuonBestTrack2));
	      ttv3.push_back(ttkb3->build(MuonBestTrack3));
	      KalmanVertexFitter kvf(true); // false means no smoothing which means no track re-fit
	      CachingVertex<5> vtx2 = kvf.vertex(ttv2);
	      vtxNormChi2 = vtx2.totalChiSquared()/vtx2.degreesOfFreedom();
	      InvariantMassFromVertex imfv2;
	      static const double muon_mass2 = 0.1056583;
	      Measurement1D mass2 = imfv2.invariantMass(vtx2, muon_mass2);
	      DiMass2 = mass2.value();
	      Mu_vtxMass.push_back(DiMass2);
	      Mu_vtxNormChi2.push_back(vtxNormChi2);
	      Mu_vtxMassLept.push_back(MuonBestTrack1->pt());
	      Mu_vtxMassLept.push_back(MuonBestTrack3->pt());
	      /*cout << "PT 1mu= " << MuonBestTrack1->pt()  << " and charge =" <<  MuonBestTrack1->charge() << endl;
	      cout << "PT 3mu= " << MuonBestTrack3->pt()  << " and charge =" <<  MuonBestTrack3->charge() << endl;
	      cout << "normChi2 = %f mass = " << vtxNormChi2 << " " << DiMass2 << endl;*/
	      CachingVertex<5> vtx3 = kvf.vertex(ttv3);
	      vtxNormChi3 = vtx3.totalChiSquared()/vtx3.degreesOfFreedom();
	      InvariantMassFromVertex imfv3;
	      static const double muon_mass3 = 0.1056583;
	      Measurement1D mass3 = imfv3.invariantMass(vtx3, muon_mass3);
	      DiMass3 = mass3.value();
	      Mu_vtxMass.push_back(DiMass3);
	      Mu_vtxNormChi2.push_back(vtxNormChi3);
	      Mu_vtxMassLept.push_back(MuonBestTrack2->pt());
              Mu_vtxMassLept.push_back(MuonBestTrack3->pt());
	      /*cout << "PT 2mu= " << MuonBestTrack2->pt()  << " and charge =" <<  MuonBestTrack2->charge() << endl;
		cout <<  "PT 3mu= " << MuonBestTrack3->pt()  << " and charge =" <<  MuonBestTrack3->charge() << endl;
		cout << "normChi2 = %f mass = " << vtxNormChi3 << " " << DiMass3 << endl;*/
	      //cout<<"size(inside) = "<<ttv2.size()<<endl;
	      /*printf ("normChi2 = %f mass = %f\n",vtxNormChi2,DiMass2);
		printf ("normChi2 = %f mass = %f\n",vtxNormChi3,DiMass3);*/
	      //cout<<"==================================================================="<<endl;
	    }
	  }
	}
      }
    }
  }
}
void MakeZprimeMiniAodTreeMC::ComputeMuonMassVtx30GeV(const edm::Event& evt,const edm::EventSetup& es)
{ 
  int NbMuons1 = 0;
  int NbMuons2 = 0;
  int NbMuons3 = 0;
  Mu_vtxMass30GeV.clear();
  Mu_vtxNormChi30GeV.clear();
  Mu_vtxMass30GeVLept.clear();
  float vtxNormChi1,DiMass1;
  float vtxNormChi2,DiMass2;
  float vtxNormChi3,DiMass3;
  // Get TransientTracks (for use in e.g. the vertex fit) for each of
  // the muon tracks, using e.g. the cocktail momentum.
  edm::ESHandle<TransientTrackBuilder> ttkb1;
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", ttkb1);
  std::vector<reco::TransientTrack> ttv1;
  edm::ESHandle<TransientTrackBuilder> ttkb2;
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", ttkb2);
  std::vector<reco::TransientTrack> ttv2;
  edm::ESHandle<TransientTrackBuilder> ttkb3;
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", ttkb3);
  std::vector<reco::TransientTrack> ttv3;
  //find the first high pt muon
  reco::TrackRef MuonBestTrack1;
  reco::TrackRef MuonBestTrack2;
  reco::TrackRef MuonBestTrack3;
  // primary vertex candidate collection
  edm::Handle<reco::VertexCollection> vertices;
  evt.getByToken(vtxToken_,vertices);
  const reco::Vertex &vertex = vertices->front();
  // pat candidate collection
  edm::Handle<pat::MuonCollection> muons;
  evt.getByToken(muonToken_, muons);
  for (const pat::Muon &mu : *muons) {
    if( mu.pt() < 5.0 ) continue;
    if( !mu.globalTrack().isNonnull() ) continue;
    if( mu.isTrackerMuon()==false ) continue;
    MuonBestTrack1 = mu.tunePMuonBestTrack();
    if( MuonBestTrack1->pt()>30.0 && (mu.isolationR03().sumPt/mu.innerTrack()->pt()<0.10) 
    	&& (MuonBestTrack1->ptError()/MuonBestTrack1->pt()<0.3) && 
        fabs(MuonBestTrack1->dxy(vertex.position()))<0.2 &&
        mu.globalTrack()->hitPattern().trackerLayersWithMeasurement()>5 && 
        mu.globalTrack()->hitPattern().numberOfValidPixelHits()>0 &&
        mu.globalTrack()->hitPattern().numberOfValidMuonHits()>0 &&
        (mu.numberOfMatchedStations()>1
	 || (mu.numberOfMatchedStations()==1 && !(mu.stationMask()==1 || mu.stationMask()==16))
	 || ((mu.numberOfMatchedStations()==1 && (mu.stationMask()==1 || mu.stationMask()==16)) && mu.numberOfMatchedRPCLayers()>2)) ){
      NbMuons1++;
      if(NbMuons1>1) continue;
      ttv1.push_back(ttkb1->build(MuonBestTrack1));
      //find the second high pt muon
      for (const pat::Muon &mu2 : *muons) {
	if( mu2.pt() < 5.0 ) continue;
	if( !mu2.globalTrack().isNonnull() ) continue;
	if( mu2.isTrackerMuon()==false ) continue;
	MuonBestTrack2 = mu2.tunePMuonBestTrack();
        if(MuonBestTrack2->pt() == MuonBestTrack1->pt()) continue;
	if( MuonBestTrack2->pt()>30.0 && (mu2.isolationR03().sumPt/mu2.innerTrack()->pt()<0.10) 
	    && (MuonBestTrack2->ptError()/MuonBestTrack2->pt()<0.3) && 
	    fabs(MuonBestTrack2->dxy(vertex.position()))<0.2 &&
	    mu2.globalTrack()->hitPattern().trackerLayersWithMeasurement()>5 && 
	    mu2.globalTrack()->hitPattern().numberOfValidPixelHits()>0 &&
	    mu2.globalTrack()->hitPattern().numberOfValidMuonHits()>0 &&
	    (mu2.numberOfMatchedStations()>1
	     || (mu2.numberOfMatchedStations()==1 && !(mu2.stationMask()==1 || mu2.stationMask()==16))
	     || ((mu2.numberOfMatchedStations()==1 && (mu2.stationMask()==1 || mu2.stationMask()==16)) && mu2.numberOfMatchedRPCLayers()>2)) ){
	  NbMuons2++;
	  if(NbMuons2>1) continue;
	  ttv1.push_back(ttkb1->build(MuonBestTrack2));
	  KalmanVertexFitter kvf(true); // false means no smoothing which means no track re-fit
	  CachingVertex<5> vtx1 = kvf.vertex(ttv1);
	  vtxNormChi1 = vtx1.totalChiSquared()/vtx1.degreesOfFreedom();
	  InvariantMassFromVertex imfv1;
	  static const double muon_mass1 = 0.1056583;
	  Measurement1D mass1 = imfv1.invariantMass(vtx1, muon_mass1);
	  DiMass1 = mass1.value();
	  Mu_vtxMass30GeV.push_back(DiMass1);
	  Mu_vtxNormChi30GeV.push_back(vtxNormChi1);
	  Mu_vtxMass30GeVLept.push_back(MuonBestTrack1->pt());
	  Mu_vtxMass30GeVLept.push_back(MuonBestTrack2->pt());
	  //find the third high pt muon
	  for (const pat::Muon &mu3 : *muons) {
	    if( mu3.pt() < 5.0 ) continue;
	    if( !mu3.globalTrack().isNonnull() ) continue;
	    if( mu3.isTrackerMuon()==false ) continue;
	    MuonBestTrack3 = mu3.tunePMuonBestTrack();
	    if(MuonBestTrack3->pt() == MuonBestTrack1->pt()) continue;	   
	    if(MuonBestTrack3->pt() == MuonBestTrack2->pt()) continue;	    
	    if( MuonBestTrack3->pt()>30.0 && (mu3.isolationR03().sumPt/mu3.innerTrack()->pt()<0.10) 
		&& (MuonBestTrack3->ptError()/MuonBestTrack3->pt()<0.3) && 
		fabs(MuonBestTrack3->dxy(vertex.position()))<0.2 &&
		mu3.globalTrack()->hitPattern().trackerLayersWithMeasurement()>5 && 
		mu3.globalTrack()->hitPattern().numberOfValidPixelHits()>0 &&
		mu3.globalTrack()->hitPattern().numberOfValidMuonHits()>0 &&
		(mu3.numberOfMatchedStations()>1
		 || (mu3.numberOfMatchedStations()==1 && !(mu3.stationMask()==1 || mu3.stationMask()==16))
		 || ((mu3.numberOfMatchedStations()==1 && (mu3.stationMask()==1 || mu3.stationMask()==16)) && mu3.numberOfMatchedRPCLayers()>2)) ){
	      NbMuons3++;
	      if(NbMuons3==0) continue;
	      ttv2.push_back(ttkb2->build(MuonBestTrack1));
	      ttv2.push_back(ttkb2->build(MuonBestTrack3));
	      ttv3.push_back(ttkb3->build(MuonBestTrack2));
	      ttv3.push_back(ttkb3->build(MuonBestTrack3));
	      KalmanVertexFitter kvf(true); // false means no smoothing which means no track re-fit
	      CachingVertex<5> vtx2 = kvf.vertex(ttv2);
	      vtxNormChi2 = vtx2.totalChiSquared()/vtx2.degreesOfFreedom();
	      InvariantMassFromVertex imfv2;
	      static const double muon_mass2 = 0.1056583;
	      Measurement1D mass2 = imfv2.invariantMass(vtx2, muon_mass2);
	      DiMass2 = mass2.value();
	      Mu_vtxMass30GeV.push_back(DiMass2);
	      Mu_vtxNormChi30GeV.push_back(vtxNormChi2);
	      Mu_vtxMass30GeVLept.push_back(MuonBestTrack1->pt());
	      Mu_vtxMass30GeVLept.push_back(MuonBestTrack3->pt());
	      CachingVertex<5> vtx3 = kvf.vertex(ttv3);
	      vtxNormChi3 = vtx3.totalChiSquared()/vtx3.degreesOfFreedom();
	      InvariantMassFromVertex imfv3;
	      static const double muon_mass3 = 0.1056583;
	      Measurement1D mass3 = imfv3.invariantMass(vtx3, muon_mass3);
	      DiMass3 = mass3.value();
	      Mu_vtxMass30GeV.push_back(DiMass3);
	      Mu_vtxNormChi30GeV.push_back(vtxNormChi3);
	      Mu_vtxMass30GeVLept.push_back(MuonBestTrack2->pt());
              Mu_vtxMass30GeVLept.push_back(MuonBestTrack3->pt());
	    }
	  }
	}
      }
    }
  }
}
//=============================================================
//
//            Method for Pimary Vertices Tree
//
//=============================================================
void MakeZprimeMiniAodTreeMC::PrimaryVertexTree(const edm::Event& iEvent,const edm::EventSetup& es)
{
  int value3_ = 0;
  nbPv.clear();
  Nbdof.clear();
  PositionX.clear();
  PositionY.clear();
  PositionZ.clear();
  PositionRho.clear();
  edm::Handle<reco::VertexCollection> pvHandle; 
  iEvent.getByToken(vtxToken_,pvHandle);
  const reco::VertexCollection &vertices = *pvHandle.product();
  for(reco::VertexCollection::const_iterator it=vertices.begin() ; it!=vertices.end() ; ++it){
    value3_ ++;
    nbPv.push_back(value3_);
    Nbdof.push_back((*it).ndof());
    PositionX.push_back((*it).position().x());
    PositionY.push_back((*it).position().y());
    PositionZ.push_back((*it).position().z());
    PositionRho.push_back((*it).position().rho());
  }
}

//=============================================================
//
//                Method for PF MET Tree
//            mets = cms.InputTag("slimmedMETs"),
//=============================================================
void MakeZprimeMiniAodTreeMC::fillMET(const edm::Event& iEvent)
{
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
  const pat::MET &met = mets->front();
  //[1] https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETRun2Corrections
  //It is strongly recommended to use the Type-1 corrected MET, directly available from miniaod. 
  //In case it is needed to recompute the Type-1 MET, please follow this twiki: here you find how 
  //to get straight from miniAOD (slimmedMET and slimmedMETnoHF), how to recalculate correction 
  //and uncertainties with the latest and greatest JEC. 
  //You can find instructions for each analysis release. 
  //[2] https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD#Run2_Fall15_MiniAOD_v2_campaign
  //the slimmedMETsNoHF was dropped (the recommendation from the MET group is to use standard slimmedMETs in 76X)
  //The default type1 corrected MET
  PFMet_et_cor    = met.et();     
  PFMet_pt_cor    = met.pt();     
  PFMet_phi_cor   = met.phi();    
  PFMet_en_cor    = met.energy(); 
  PFMet_px_cor    = met.px();     
  PFMet_py_cor    = met.py();     
  PFMet_pz_cor    = met.pz();     
  PFMet_sumEt_cor = met.sumEt();
  //Type I MET Uncertainties
  pfMETshiftedPtJetEnUp      = met.shiftedPt(pat::MET::JetEnUp, pat::MET::Type1 ) ;
  pfMETshiftedPtJetEnDn      = met.shiftedPt(pat::MET::JetEnDown, pat::MET::Type1 ) ;
  pfMETshiftedPtEleEnUp      = met.shiftedPt(pat::MET::ElectronEnUp, pat::MET::Type1 ) ;
  pfMETshiftedPtEleEnDn      = met.shiftedPt(pat::MET::ElectronEnDown, pat::MET::Type1 ) ;
  pfMETshiftedPtMuEnUp       = met.shiftedPt(pat::MET::MuonEnUp, pat::MET::Type1 ) ;
  pfMETshiftedPtMuEnDn       = met.shiftedPt(pat::MET::MuonEnDown, pat::MET::Type1 ) ;
  pfMETshiftedPtJetResUp     = met.shiftedPt(pat::MET::JetResUp, pat::MET::Type1 ) ;
  pfMETshiftedPtJetResDn     = met.shiftedPt(pat::MET::JetResDown, pat::MET::Type1 ) ;
  pfMETshiftedPtUnclEnUp     = met.shiftedPt(pat::MET::UnclusteredEnUp, pat::MET::Type1 ) ;
  pfMETshiftedPtUnclEnDn     = met.shiftedPt(pat::MET::UnclusteredEnDown, pat::MET::Type1 ) ; 
  pfMETshiftedPtPhoEnUp      = met.shiftedPt(pat::MET::PhotonEnUp, pat::MET::Type1 ) ;
  pfMETshiftedPtPhoEnDn      = met.shiftedPt(pat::MET::PhotonEnDown, pat::MET::Type1 ) ; 
  //pfMETshiftedPtJetEnUpSmear = met.shiftedPhi(pat::MET::METUncertainty::JetResUpSmear);
  //pfMETshiftedPtJetEnDnSmear = met.shiftedPt(pat::MET::JetResDownSmear, pat::MET::Type1Smear ) ;
  //pfMETUncertaintySize       = met.shiftedPt(pat::MET::METUncertaintySize, pat::MET::Type1 ) ;
  //pfMETFullUncertaintySize   = met.shiftedPt(pat::MET::METFullUncertaintySize, pat::MET::Type1 ) ;

  //cout<<"jet smeared = "<<met.shiftedPhi(pat::MET::METUncertainty::JetResUpSmear)<<endl;
 
  /*
  printf("met_pt = %f pt:JEC up = %f pt:JEC dn = %f uncertainty1 = %f \n",
  met.pt(),met.shiftedPt(pat::MET::JetEnUp, pat::MET::Type1),met.shiftedPhi(pat::MET::JetEnDown, pat::MET::Type1),met.pt()-met.shiftedPt(pat::MET::JetEnUp, pat::MET::Type1));
  printf("met_pt = %f pt:Uncluster up = %f pt:Uncluster dn = %f uncertainty2 = %f \n",
  met.pt(),met.shiftedPt(pat::MET::UnclusteredEnUp, pat::MET::Type1),met.shiftedPhi(pat::MET::UnclusteredEnDown, pat::MET::Type1),met.pt()-met.shiftedPt(pat::MET::UnclusteredEnUp, pat::MET::Type1));
  */
  //std::cout<<"met uncer"<< pat::MET::METUncertainty(12) <<endl;
  //The Raw PF Met (un-corrected MET in CMSSW_7_4_12 or later, on 74X version 2 miniAODs)
  /*PFMet_pt_uncor    = met.uncorPt();     
  PFMet_phi_uncor   = met.uncorPhi();    
  PFMet_en_uncor    = met.uncorEnergy(); 
  PFMet_px_uncor    = met.uncorPx();     
  PFMet_py_uncor    = met.uncorPy();     
  PFMet_pz_uncor    = met.uncorPz();     
  PFMet_sumEt_uncor = met.uncorSumEt();
  //The Raw PF Met (un-corrected MET in CMSSW_7_4_11 or earlier, on 74X version 1 miniAODs)
  PFMet_pt_uncor    = met.uncorrectedPt();     
  PFMet_phi_uncor   = met.uncorrectedPhi();    
  PFMet_sumEt_uncor = met.uncorrectedSumEt();*/
  //The raw calo ETmiss  
  CaloMet_pt    = met.caloMETPt();     
  CaloMet_phi   = met.caloMETPhi();    
  CaloMet_sumEt = met.caloMETSumEt();
  PFMet_shiftedPt_JetEnUp   = met.shiftedPt(pat::MET::JetEnUp);
  PFMet_shiftedPt_JetEnDown = met.shiftedPt(pat::MET::JetEnDown);
  if (met.genMET() != NULL ) GenMet_pt = met.genMET()->pt();  
  /*edm::Handle<double> metsighandle;
  iEvent.getByToken(theMETSignificance_, metsighandle);
  METSign=*metsighandle;*/
}

//=============================================================
//
// Method for Gen Jet Tree
//
//=============================================================
void MakeZprimeMiniAodTreeMC::GenJetTree(const edm::Event& iEvent)
{ 
  int NbGenJets  = 0;
  iGenJet.clear();
  idGenJet.clear();
  statusGenJet.clear();
  ptGenJet.clear();
  etaGenJet.clear();
  phiGenJet.clear();
  chargeGenJet.clear();
  /**** GET GENJETS ****/
  edm::Handle<reco::GenJetCollection> h_genjets;
  iEvent.getByToken( EDMGenJetsToken_,h_genjets );
  if (!(h_genjets.isValid())) return; 
  std::vector<reco::GenJet> const &genjets = *h_genjets;
  for(size_t i=0; i<genjets.size();i++){
    if(genjets[i].pt()<20.0) continue;
    NbGenJets++;
    iGenJet.push_back(NbGenJets);
    idGenJet.push_back(genjets[i].pdgId());
    statusGenJet.push_back(genjets[i].status());
    ptGenJet.push_back(genjets[i].pt());
    etaGenJet.push_back(genjets[i].eta());
    phiGenJet.push_back(genjets[i].phi());
    chargeGenJet.push_back(genjets[i].charge());
  }
}
 
//=============================================================
//
//            Method for Jets Tree
//
//=============================================================
void MakeZprimeMiniAodTreeMC::JetsTree(const edm::Event& iEvent,const edm::EventSetup& es)
{
  jet_nb.clear();
  jet_charge.clear();
  jet_et.clear();
  jet_pt.clear();
  jet_en.clear();
  jet_eta.clear();
  jet_phi.clear();
  jet_theta.clear();
  jet_beta.clear();
  jet_pileup_mva_disc.clear();
  int jetnumber=0;
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);
  //int ijet = 0;
  for (const pat::Jet &j : *jets) {
    double in = 0, out = 0; 
    for (unsigned int id = 0, nd = j.numberOfDaughters(); id < nd; ++id) {
      const pat::PackedCandidate &dau = dynamic_cast<const pat::PackedCandidate &>(*j.daughter(id));
      if (dau.charge() == 0) continue;
      (fabs(dau.dz())<0.1 ? in : out) += dau.pt();
    }
    double sum = in + out;
    // Loose ID jet selection
    if( fabs(j.eta()) < 2.4 &&
	j.pt() > 20 &&
	j.neutralHadronEnergyFraction() < 0.99 &&
	j.neutralEmEnergyFraction() < 0.99 &&
	j.numberOfDaughters() > 1.0 &&
	j.chargedHadronEnergyFraction() > 0.0 &&
	j.chargedMultiplicity() > 0.0 &&
	j.chargedEmEnergyFraction() < 0.99 ) {
      //cout << "Jet passing the Tight ID with pT=" << j.pt() << " and eta=" << j.eta() << endl;
      jetnumber++;
      //printf("j = %d jet_pt = %f \n",jetnumber,j.pt());
      jet_nb.push_back(jetnumber);
      jet_charge.push_back(j.charge());
      jet_et.push_back(j.et());
      jet_pt.push_back(j.pt());
      jet_eta.push_back(j.eta());
      jet_phi.push_back(j.phi());
      jet_en.push_back(j.muonEnergy());
      jet_theta.push_back(j.theta());
      jet_beta.push_back(sum ? in/sum : 0);
      jet_pileup_mva_disc.push_back(j.userFloat("pileupJetId:fullDiscriminant"));
    }
    else continue;
  }   
}
//=============================================================
//
//            Method for Btagging Tree
//
//=============================================================
void MakeZprimeMiniAodTreeMC::BtaggingTree(const edm::Event& iEvent)
{
  int bDiscriminatorsNumber=0;
  Nb_bDiscriminators.clear();
  jet_btag_pt.clear();
  jet_btag_eta.clear();
  jet_btag_phi.clear();
  jet_btag_flavor.clear();
  jet_btag_pfCSVv2IVF_discriminator.clear();
  // define a jet handle
  edm::Handle<std::vector<pat::Jet> > jets;
  // get jets from the event
  iEvent.getByToken(jetToken_, jets);
  // loop over jets
  for( auto jet = jets->begin(); jet != jets->end(); ++jet )
    { 
      // fill discriminator variables
      int flavor = std::abs( jet->partonFlavour() );
      for( const std::string &bDiscr : bDiscriminators_ )
	{
	  //cout<<"Btag = "<<jet->bDiscriminator(bDiscr)<<endl;
          //if( flavor==0 ) continue; // skip jets with undefined flavor
	  //cout<<"Btag = "<<jet->bDiscriminator(bDiscr)<<endl;
          if( jet->pt()<20.0 || std::abs(jet->eta())>2.4 ) continue; // skip jets with low pT or outside the tracker acceptance
	  //cout<<"Btag = "<<jet->bDiscriminator(bDiscr)<<endl;
          bDiscriminatorsNumber++;
	  Nb_bDiscriminators.push_back(bDiscriminatorsNumber);
	  jet_btag_pt.push_back(jet->pt());
	  jet_btag_eta.push_back(jet->eta());
	  jet_btag_phi.push_back(jet->phi());
	  jet_btag_flavor.push_back(flavor);
	  jet_btag_pfCSVv2IVF_discriminator.push_back(jet->bDiscriminator(bDiscr));
	}
    }
}
//if( flavor==5 )      {jet_btag_pfCSVv2IVF_discriminator_b.push_back(j.bDiscriminator(bDiscr));} // b jet
//else if( flavor==4 ) {jet_btag_pfCSVv2IVF_discriminator_c.push_back(j.bDiscriminator(bDiscr));} // c jets
//else                 {jet_btag_pfCSVv2IVF_discriminator_udsg.push_back(j.bDiscriminator(bDiscr));}// light-flavor jet

//=============================================================
//
//            Method for particle Flow Tree
//
//=============================================================
void MakeZprimeMiniAodTreeMC::ParticleFlowPhotonTree(const edm::Event& iEvent,const edm::EventSetup& es)
{ 
  pfphoton_size = 0;
  pfphoton_pt.clear();
  pfphoton_eta.clear();
  pfphoton_phi.clear();
  pfphoton_theta.clear();
  //pfphoton_PFchHad.clear();
  //pfphoton_PFneuHad.clear();
  //pfphoton_PFphoton.clear();
  //pfphoton_PFPUchAllPart.clear();
  //pfphoton_PFX_rho.clear();
  edm::Handle<pat::PackedCandidateCollection> pfs;
  iEvent.getByToken(pfToken_, pfs);
  std::vector<const reco::Candidate *> leptons;
  edm::Handle<pat::PhotonCollection> photons;
  iEvent.getByToken(photonToken_, photons);
  for (const pat::Photon &pho : *photons) leptons.push_back(&pho);
  for (const reco::Candidate *lep : leptons) {
    // initialize sums
    double charged = 0, neutral = 0, pileup  = 0;
    // now get a list of the PF candidates used to build this lepton, so to exclude them
    std::vector<reco::CandidatePtr> footprint;
    for (unsigned int i = 0, n = lep->numberOfSourceCandidatePtrs(); i < n; ++i) {
      footprint.push_back(lep->sourceCandidatePtr(i));
    }
    // now loop on pf candidates
    //pfphoton_size=pfs->size();
    for (unsigned int i = 0, n = pfs->size(); i < n; ++i) {
      const pat::PackedCandidate &pf = (*pfs)[i];
      if (deltaR(pf,*lep) < 0.2) {
	// pfcandidate-based footprint removal
	if(std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfs,i)) != footprint.end()) continue;
	if( pf.pt() < 20.0 ) continue;
	if( fabs(pf.eta()) > 2.5 ) continue;
	pfphoton_pt.push_back(pf.pt());
	pfphoton_eta.push_back(pf.eta());
	pfphoton_phi.push_back(pf.phi());
	pfphoton_theta.push_back(pf.theta());
	pfphoton_size++;
	if (pf.charge() == 0) {
	  if (pf.pt() > 0.5) neutral += pf.pt();
	} else if (pf.fromPV() >= 2) {
	  charged += pf.pt();
	} else {
	  if (pf.pt() > 0.5) pileup += pf.pt();
	}
      }
    }
  }
}
//madgraph MC samples reweighing 
void MakeZprimeMiniAodTreeMC::EventsReWeighting(const edm::Event& evt){
  MC_weighting.clear();
  float EventWeight = 1.0;
  edm::Handle<GenEventInfoProduct> gen_ev_info;
  evt.getByToken(genInfoProductToken, gen_ev_info);
  if(!gen_ev_info.isValid()) return;
  EventWeight = gen_ev_info->weight();
  //std::cout<<"mc_weight = "<< gen_ev_info->weight() <<std::endl;
  float mc_weight = ( EventWeight > 0 ) ? 1 : -1;
  //std::cout<<"mc_weight = "<< mc_weight <<std::endl;
  MC_weighting.push_back(mc_weight);
}

void MakeZprimeMiniAodTreeMC::fillPU(const edm::Event& iEvent){
  edm::Handle<vector<PileupSummaryInfo> > PupInfo;
  iEvent.getByToken(PileupSrc_, PupInfo);
  if(!PupInfo.isValid()) return;
  for( vector<PileupSummaryInfo>::const_iterator cand = PupInfo->begin();cand != PupInfo->end(); ++ cand ) {
    num_PU_vertices = cand->getTrueNumInteractions();
    PU_BunchCrossing = cand->getBunchCrossing();
    num_PU_gen_vertices = cand->getPU_NumInteractions();
    //std::cout << " Pileup Information: bunchXing, nvtx: " << cand->getBunchCrossing() << " " << cand->getPU_NumInteractions() << std::endl;
  }
}
// rho for isolation
void MakeZprimeMiniAodTreeMC::fillRho(const edm::Event& evt){
  edm::Handle<double> rhoHandle;
  evt.getByToken(rhoToken,rhoHandle);
  if(!rhoHandle.isValid()) return;
  Rho = *rhoHandle;
}
//=============================================================
//
//            Method for Tau Tree
//
//=============================================================
void MakeZprimeMiniAodTreeMC::TauTree(const edm::Event& iEvent)
{
  int TausNumber=0;
  Nb_taus.clear();
  Tau_pt.clear();
  Tau_eta.clear();
  Tau_phi.clear();
  Tau_id.clear();
  Tau_LooseCombinedIsolationDeltaBetaCorr3Hits.clear();
  edm::Handle<pat::TauCollection> taus;
  iEvent.getByToken(tauToken_, taus);
  for (const pat::Tau &tau : *taus) {
    if(tau.pt() < 18.0) continue;
    if(fabs(tau.eta()) > 2.3) continue;
    /*printf("tau  with pt %4.1f, dxy signif %.1f, ID(byLooseCombinedIsolationDeltaBetaCorr3Hits) %.1f, lead candidate pt %.1f, pdgId %d \n",
      tau.pt(), tau.dxy_Sig(), tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"), tau.leadCand()->pt(), tau.leadCand()->pdgId());*/
    Nb_taus.push_back(TausNumber);
    Tau_pt.push_back(tau.pt());
    Tau_eta.push_back(tau.eta());
    Tau_phi.push_back(tau.phi());
    Tau_id.push_back(tau.pdgId());
    Tau_LooseCombinedIsolationDeltaBetaCorr3Hits.push_back(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
    TausNumber++;
  }
}
//====================================== the end of subrotines ======================================
/*
  [0] kgood=0,                   // channel ok, the energy and time measurement are reliable
  [1] kpoorreco,                 // the energy is available from the uncalibrechit, but approximate (bad shape, large chi2)
  [2] koutoftime,                // the energy is available from the uncalibrechit (sync reco), but the event is out of time
  [3] kfaultyhardware,           // the energy is available from the uncalibrechit, channel is faulty at some hardware level (e.g. noisy)
  [4] knoisy,                    // the channel is very noisy
  [5] kpoorcalib,                // the energy is available from the uncalibrechit, but the calibration of the channel is poor
  [6] ksaturated,                // saturated channel (recovery not tried)
  [7] kleadingedgerecovered,     // saturated channel: energy estimated from the leading edge before saturation
  [8] kNeighboursRecovered,      // saturated/isolated dead: energy estimated from neighbours
  [9] kTowerRecovered,           // channel in TT with no data link, info retrieved from Trigger Primitive
  [10] kDead,                     // channel is dead and any recovery fails
  [11] kKilled,                   // MC only flag: the channel is killed in the real detector
  [12] kTPSaturated,              // the channel is in a region with saturated TP
  [13] kL1SpikeFlag,              // the channel is in a region with TP with sFGVB = 0
  [14] kWeird,                    // the signal is believed to originate from an anomalous deposit (spike) 
  [15] kDiWeird,                  // the signal is anomalous, and neighbors another anomalous signal  
  [16] kHasSwitchToGain6,         // at least one data frame is in G6
  [17] kHasSwitchToGain1,         // at least one data frame is in G1
  [18] kUnknown                   // to ease the interface with functions returning flags.
*/
//====================================== The end of subrotines ======================================

