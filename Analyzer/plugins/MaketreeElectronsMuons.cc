//=========================================================================  
//      Make the root tree for Z' boson to Mu Mu or e e analysis          =
//                      (works only with AOD)                             =  
//                                                                        =
//                         CMSSW_7_4_5                                    =
//                                                                        =
//                 Written by Sherif Elgammal                             =
//                                                                        =
//                         30/03/2015                                     =
//=========================================================================
#include "ZprimeDiLeptons/Analyzer/interface/MaketreeElectronsMuons.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/GeometryVector/interface/Phi.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/View.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/Common/interface/RefCore.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "HLTrigger/HLTfilters/interface/HLTHighLevel.h"
using namespace edm;
using namespace std;
using namespace reco;

//====================== end a part for photons ==========================
#include <vector>
#include <set>
#include <stdio.h>
#include "TFile.h"
#include <math.h>
#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"

//========================================================================
MaketreeElectronsMuons::MaketreeElectronsMuons( const edm::ParameterSet& ps )
//========================================================================
{
  vertexSrc                      = ps.getParameter<edm::InputTag>("vertexCollection");
  maxAbsZ_                       = ps.getParameter<double>("maxAbsZ");
  maxd0_                         = ps.getParameter<double>("maxd0");
  minNdof_                       = ps.getParameter<int>("minndof");
  NbGoodPv_                      = ps.getParameter<int>("NbGoodPv");
  globalMuons_                   = ps.getParameter<edm::InputTag>("globalMuons");
  genParticlesColl_              = ps.getParameter<edm::InputTag>("genparticleCollection");
  token_globalMuons              = ps.getParameter<edm::InputTag>( "globalMuonTracks" );
  rhoIsoInputTag                 = ps.getParameter<edm::InputTag>("rhoIsoInputTag");
  reducedBarrelRecHitCollection_ = ps.getParameter<edm::InputTag>("reducedBarrelRecHitCollection");
  reducedEndcapRecHitCollection_ = ps.getParameter<edm::InputTag>("reducedEndcapRecHitCollection");
  GsfElectronCollection_         = ps.getParameter<std::string>("GsfElectronCollection");
  GsfElectronProducer_           = ps.getParameter<std::string>("GsfElectronProducer");
  // get names from module parameters, then derive slot numbers
  n_                  = 0;
  firstevent_         = true;  
  //produces<vector<std::string> >();
  // PG and FRC 06-07-11
  debug	=	ps.getUntrackedParameter<bool> ("debug", false);
  //pfMuons_                       = ps.getParameter<edm::InputTag>("pfMuons");
  //pfMuonToken_                   = ps.getParameter<edm::InputTag>("PFjetsCollection");
  //tokentracks_                   = ps.getParameter<edm::InputTag>("TrackCollectionTag");
  //theRecoLabel_                  = ps.getParameter<edm::InputTag>("inputTagMuonReco");
  //------------------------------------------------
  //       Get Parameters for HLT objects          -
  //------------------------------------------------
  //inputTag_                    = ps.getParameter<edm::InputTag> ("TriggerResultsTag");
  //HLTFilterName_               = ps.getUntrackedParameter<std::vector<std::string> >("HLTFilterName");
  triggerEvent                   = ps.getParameter<edm::InputTag>("triggerEvent");
  triggerFilter                  = ps.getParameter<std::string>("triggerFilter");
  //triggerFilter_asym             = ps.getParameter<std::vector<std::string> >("triggerFilterAsym");
  //HLTFilterName_                 = ps.getUntrackedParameter<std::string>("HLTFilterName");
  //===============================================================================================
  outputFile_ = ps.getParameter<std::string>("outputFile");
  rootFile_   = TFile::Open(outputFile_.c_str(),"RECREATE"); // open output file to store histograms
}

 
//========================================================================
MaketreeElectronsMuons::~MaketreeElectronsMuons()
//========================================================================    
{
  delete rootFile_;
}

//========================================================================
void MaketreeElectronsMuons::beginJob() {
  //======================================================================== 
  // go to *OUR* rootfile and book histograms
  rootFile_->cd();
  //========================================================
  mytree  = new TTree("tree","tr");
  //=============================================================
  //
  //           Create Branchs for Nb of event,run,lumi
  //
  //=============================================================
  TotalNbEvents = 0;
  NbEventsPass  = 0;
  Run       = 0;
  Event     = 0;
  lumi      = 0;
  bunch     = 0;
  mytree->Branch("NbEventsPassTrigger",&NbEventsPassTrigger);
  mytree->Branch("NbEventsPassTriggerandPVcond",&NbEventsPassTriggerandPVcond);
  mytree->Branch("event_runNo",  &Run,   "event_runNo/I");
  mytree->Branch("event_evtNo",  &Event, "event_evtNo/I");
  mytree->Branch("event_lumi",   &lumi,  "event_lumi/I");
  mytree->Branch("event_bunch",  &bunch, "event_bunch/I");
  //=============================================================
  //
  //           Create Branchs for Muons variables
  //
  //=============================================================
  mytree->Branch("Mu_nbMuon",&Mu_nbMuon);
  mytree->Branch("Mu_isGlobalMuon",&Mu_isGlobalMuon);
  mytree->Branch("Mu_isPF",&Mu_isPF);
  mytree->Branch("Mu_isGoodMuon",&Mu_isGoodMuon);
  mytree->Branch("Mu_isTrackerMuon",&Mu_isTrackerMuon);
  mytree->Branch("Mu_et",&Mu_et);
  mytree->Branch("Mu_en",&Mu_en);
  mytree->Branch("Mu_phi",&Mu_phi);
  mytree->Branch("Mu_eta",&Mu_eta);
  mytree->Branch("Mu_pt",&Mu_pt);
  // part for old TuneP track
  mytree->Branch("Mu_ptcocktail",&Mu_ptcocktail);
  mytree->Branch("Mu_etaCocktail",&Mu_etaCocktail);
  mytree->Branch("Mu_phiCocktail",&Mu_phiCocktail);
  mytree->Branch("Mu_thetaCocktail",&Mu_thetaCocktail);  
  mytree->Branch("Mu_pxCocktail",&Mu_pxCocktail);
  mytree->Branch("Mu_pyCocktail",&Mu_pyCocktail);
  mytree->Branch("Mu_pzCocktail",&Mu_pzCocktail);
  mytree->Branch("Mu_pCocktail",&Mu_pCocktail);
  mytree->Branch("Mu_dPToverPTcocktail",&Mu_dPToverPTcocktail);
  mytree->Branch("Mu_chargeCocktail",&Mu_chargeCocktail);
  mytree->Branch("Mu_absdxy",&Mu_absdxy);
  mytree->Branch("Mu_absdz",&Mu_absdz);
  mytree->Branch("Mu_normalizedChi2",&Mu_normalizedChi2);
  mytree->Branch("Mu_vtxMass",&Mu_vtxMass);
  mytree->Branch("Mu_vtxNormChi2",&Mu_vtxNormChi2);
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
  mytree->Branch("Mu_ptInnerTrack",&Mu_ptInnerTrack);  
  // part for old TuneP Muon Best Track
  mytree->Branch("Mu_ptTunePMuonBestTrack",&Mu_ptTunePMuonBestTrack);
  mytree->Branch("Mu_dPToverPTTunePMuonBestTrack",&Mu_dPToverPTTunePMuonBestTrack);
  mytree->Branch("Mu_pxTunePMuonBestTrack",&Mu_pxTunePMuonBestTrack);
  mytree->Branch("Mu_pyTunePMuonBestTrack",&Mu_pyTunePMuonBestTrack);
  mytree->Branch("Mu_pzTunePMuonBestTrack",&Mu_pzTunePMuonBestTrack);
  mytree->Branch("Mu_pTunePMuonBestTrack",&Mu_pTunePMuonBestTrack);
  mytree->Branch("Mu_etaTunePMuonBestTrack",&Mu_etaTunePMuonBestTrack);
  mytree->Branch("Mu_phiTunePMuonBestTrack",&Mu_phiTunePMuonBestTrack);
  mytree->Branch("Mu_thetaTunePMuonBestTrack",&Mu_thetaTunePMuonBestTrack);
  mytree->Branch("Mu_chargeTunePMuonBestTrack",&Mu_chargeTunePMuonBestTrack);
  mytree->Branch("Mu_absdxyTunePMuonBestTrack",&Mu_absdxyTunePMuonBestTrack);
  // part for DYT track
  mytree->Branch("Mu_ptDYTTrack",&Mu_ptDYTTrack);
  mytree->Branch("Mu_pxDYTTrack",&Mu_pxDYTTrack);
  mytree->Branch("Mu_pyDYTTrack",&Mu_pyDYTTrack);
  mytree->Branch("Mu_pzDYTTrack",&Mu_pzDYTTrack);
  mytree->Branch("Mu_pDYTTrack",&Mu_pDYTTrack);
  mytree->Branch("Mu_etaDYTTrack",&Mu_etaDYTTrack);
  mytree->Branch("Mu_phiDYTTrack",&Mu_phiDYTTrack);
  mytree->Branch("Mu_thetaDYTTrack",&Mu_thetaDYTTrack);
  mytree->Branch("Mu_chargeDYTTrack",&Mu_chargeDYTTrack);
  mytree->Branch("Mu_absdxyDYTTrack",&Mu_absdxyDYTTrack);
  mytree->Branch("Mu_dPToverPTDYTTrack",&Mu_dPToverPTDYTTrack);
  // part for Picky Track
  mytree->Branch("Mu_ptPickyTrack",&Mu_ptPickyTrack);
  mytree->Branch("Mu_pxPickyTrack",&Mu_pxPickyTrack);
  mytree->Branch("Mu_pyPickyTrack",&Mu_pyPickyTrack);
  mytree->Branch("Mu_pzPickyTrack",&Mu_pzPickyTrack);
  mytree->Branch("Mu_pPickyTrack",&Mu_pPickyTrack);
  mytree->Branch("Mu_etaPickyTrack",&Mu_etaPickyTrack);
  mytree->Branch("Mu_phiPickyTrack",&Mu_phiPickyTrack);
  mytree->Branch("Mu_thetaPickyTrack",&Mu_thetaPickyTrack);
  mytree->Branch("Mu_chargePickyTrack",&Mu_chargePickyTrack);
  mytree->Branch("Mu_absdxyPickyTrack",&Mu_absdxyPickyTrack);
  mytree->Branch("Mu_dPToverPTPickyTrack",&Mu_dPToverPTPickyTrack);
  // part for Muon Best Track
  mytree->Branch("Mu_ptMuonBestTrack",&Mu_ptMuonBestTrack);
  mytree->Branch("Mu_dPToverPTMuonBestTrack",&Mu_dPToverPTMuonBestTrack);
  mytree->Branch("Mu_pxMuonBestTrack",&Mu_pxMuonBestTrack);
  mytree->Branch("Mu_pyMuonBestTrack",&Mu_pyMuonBestTrack);
  mytree->Branch("Mu_pzMuonBestTrack",&Mu_pzMuonBestTrack);
  mytree->Branch("Mu_pMuonBestTrack",&Mu_pMuonBestTrack);
  mytree->Branch("Mu_etaMuonBestTrack",&Mu_etaMuonBestTrack);
  mytree->Branch("Mu_phiMuonBestTrack",&Mu_phiMuonBestTrack);
  mytree->Branch("Mu_thetaMuonBestTrack",&Mu_thetaMuonBestTrack);
  mytree->Branch("Mu_chargeMuonBestTrack",&Mu_chargeMuonBestTrack);
  mytree->Branch("Mu_absdxyMuonBestTrack",&Mu_absdxyMuonBestTrack);
  //=============================================================
  //
  //           Create Branchs for gen particles variables
  //
  //=============================================================
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
  mytree->Branch("NbOfDaughters",&NbOfDaughters);
  mytree->Branch("McZmass",&McZmass);
  //mytree->Branch("McZpt",&McZpt);
  //mytree->Branch("McZpx",&McZpx);
  //mytree->Branch("McZpy",&McZpy);
  //mytree->Branch("McZpz",&McZpz);
  //mytree->Branch("McZen",&McZen);
  //mytree->Branch("massGen",&massGen);
  //mytree->Branch("vxGen",&vxGen);
  //mytree->Branch("vyGen",&vyGen);
  //mytree->Branch("vzGen",&vzGen);
  //=============================================================
  //
  //           Create Branchs for Pimary Vertice variables
  //
  //=============================================================
  mytree->Branch("nbPv",&nbPv);
  mytree->Branch("Nbdof",&Nbdof);
  mytree->Branch("PositionRho",&PositionRho);
  mytree->Branch("PositionX",&PositionX);
  mytree->Branch("PositionY",&PositionY);
  mytree->Branch("PositionZ",&PositionZ);
  //=============================================================
  //
  //           Create Branchs for Muons match HLT variables
  //
  //=============================================================
  mytree->Branch("MuHLTMatch_nbMuonMatchHLT",&MuHLTMatch_nbMuonMatchHLT);
  mytree->Branch("MuHLTMatch_pt",&MuHLTMatch_pt);
  mytree->Branch("MuHLTMatch_eta",&MuHLTMatch_eta);
  mytree->Branch("MuHLTMatch_phi",&MuHLTMatch_phi);
  //=================================== For HEEP ELE =====================
  mytree->Branch("value",&value_,"value/I");
  mytree->Branch("iEle",&iEle);
  mytree->Branch("PreshowerEn1",&PreshowerEn1);
  mytree->Branch("EtaEle1",&EtaEle1);
  mytree->Branch("PhiEle1",&PhiEle1);
  mytree->Branch("EtaEleSC1",&EtaEleSC1);
  mytree->Branch("PhiEleSC1",&PhiEleSC1);
  mytree->Branch("EtEle1",&EtEle1);
  mytree->Branch("EtFromCaloEnEle1",&EtFromCaloEnEle1);
  mytree->Branch("EtCorrEle1",&EtCorrEle1);
  mytree->Branch("EnEleSC1",&EnEleSC1);
  mytree->Branch("EnCorrEleSC1",&EnCorrEleSC1);
  mytree->Branch("scEmaxEle1",&scEmaxEle1);
  mytree->Branch("scE25Ele1",&scE25Ele1);
  mytree->Branch("scE2x5RightEle1",&scE2x5RightEle1);
  mytree->Branch("scE2x5LeftEle1",&scE2x5LeftEle1);
  mytree->Branch("scE2x5TopEle1",&scE2x5TopEle1);
  mytree->Branch("scE2x5BottomEle1",&scE2x5BottomEle1);
  mytree->Branch("EleClassEle1",&EleClassEle1);
  mytree->Branch("ThetaEleSC1",&ThetaEleSC1);
  mytree->Branch("ThetaEle1",&ThetaEle1);
  mytree->Branch("chargeEle1",&chargeEle1);
  mytree->Branch("DeltaEtaEle1",&DeltaEtaEle1);
  mytree->Branch("DeltaPhiEle1",&DeltaPhiEle1);
  mytree->Branch("HoeEle1",&HoeEle1);
  mytree->Branch("EcalPlusHcal1IsoEle1",&EcalPlusHcal1IsoEle1);
  mytree->Branch("EcalPlusHcal1BcIsoEle1",&EcalPlusHcal1BcIsoEle1);
  mytree->Branch("SigmaIetaIetaEle1",&SigmaIetaIetaEle1);
  mytree->Branch("HcalDepth2IsoEle1",&HcalDepth2IsoEle1);
  mytree->Branch("TkSumPtIsoEle1",&TkSumPtIsoEle1);
  mytree->Branch("E2x5MaxOverE5x5Ele1",&E2x5MaxOverE5x5Ele1);
  mytree->Branch("E1x5OverE5x5Ele1",&E1x5OverE5x5Ele1);
  mytree->Branch("fbremEle1",&fbremEle1);
  mytree->Branch("EnSeedClusterOverPEle1",&EnSeedClusterOverPEle1);
  mytree->Branch("EcalIsoEle1",&EcalIsoEle1);
  mytree->Branch("Hcal1IsoEle1",&Hcal1IsoEle1);
  mytree->Branch("Hcal1BcIsoEle1",&Hcal1BcIsoEle1);
  mytree->Branch("xSCEle1",&xSCEle1);
  mytree->Branch("ySCEle1",&ySCEle1);
  mytree->Branch("EcalDrivenSeedEle1",&EcalDrivenSeedEle1);
  mytree->Branch("NbOfLostInnerHitsEle1",&NbOfLostInnerHitsEle1);
  mytree->Branch("E2x5MaxEle1",&E2x5MaxEle1);
  mytree->Branch("E1x5Ele1",&E1x5Ele1);
  mytree->Branch("ElectronFromConv",&ElectronFromConv);
  mytree->Branch("LxyMinFromConv",&LxyMinFromConv);
  mytree->Branch("ProbMinFromConv",&ProbMinFromConv);
  mytree->Branch("NHitsBeforeVtxMaxFromConv",&NHitsBeforeVtxMaxFromConv);
  mytree->Branch("Rho",&Rho);
  mytree->Branch("pxEle1",&pxEle1);
  mytree->Branch("pyEle1",&pyEle1);
  mytree->Branch("pzEle1",&pzEle1);
  mytree->Branch("gsftrackDxyEle1",&gsftrackDxyEle1);
  mytree->Branch("gsftrackDzEle1",&gsftrackDzEle1);
  mytree->Branch("gsftrackDxyVtxEle1",&gsftrackDxyVtxEle1);
  mytree->Branch("gsftrackDzVtxEle1",&gsftrackDzVtxEle1);
  mytree->Branch("gsftrackDxyErrorVtxEle1",&gsftrackDxyErrorVtxEle1);
  mytree->Branch("gsftrackDzErrorVtxEle1",&gsftrackDzErrorVtxEle1);
  mytree->Branch("scEnLeft",&scEnLeft);
  mytree->Branch("scEnRight",&scEnRight);
  mytree->Branch("scEnTop",&scEnTop);
  mytree->Branch("scEnBottom",&scEnBottom);
  mytree->Branch("scSeedTime",&scSeedTime);
}
  

//========================================================================================
void MaketreeElectronsMuons::analyze( const edm::Event& evt, const edm::EventSetup& es ) {
//======================================================================================
  using namespace edm; // needed for all fwk related classe
  //==============================================
  //=        Begin of the main program           =
  //============================================== 
  TotalNbEvents ++;
  NbEventsPassTrigger.clear();
  NbEventsPassTrigger.push_back(TotalNbEvents);
  // rho for isolation
  edm::Handle<double> rhoIso_h;
  evt.getByLabel(rhoIsoInputTag, rhoIso_h);
  double rhoIso = *(rhoIso_h.product());
  //std::cout<<"rhoIso = "<< rhoIso <<std::endl;
  //=========================== start for GSF Electron =============
  edm::Handle<reco::GsfElectronCollection> electrons;
  evt.getByLabel(GsfElectronProducer_, GsfElectronCollection_, electrons);
  const reco::GsfElectronCollection* gsfEle = electrons.product();
  //------------------- Start part for conversion rejection ---------------//
  edm::Handle<reco::BeamSpot> bsHandle;
  evt.getByLabel("offlineBeamSpot", bsHandle);
  const reco::BeamSpot &thebs = *bsHandle.product();
  edm::Handle<reco::ConversionCollection> hConversions;
  evt.getByLabel("allConversions", hConversions);
  //------------------- End part for conversion rejection ---------------//
  //===================== Handle For Muons =====================
  edm::Handle<reco::MuonCollection> muons;
  evt.getByLabel(globalMuons_, muons);
  const reco::MuonCollection* realMuons = muons.product();
  //===================== Handle For Primary Vertics ===========
  edm::Handle<reco::VertexCollection> pvHandle; 
  evt.getByLabel(vertexSrc,pvHandle);
  const reco::VertexCollection &vertices = *pvHandle.product();
  bool GoodPv  = PrimaryVertex(vertices);
  if( GoodPv == 0 ) return;    //the aim of that command is to select only events
                               // with one good reconstructed primary vertics
  NbEventsPass ++;
  NbEventsPassTriggerandPVcond.clear();
  NbEventsPassTriggerandPVcond.push_back(NbEventsPass);
  eventVertexPosition_ = GlobalPoint(0., 0., 0.);
  const reco::Vertex& thePrimaryEventVertex = (*pvHandle->begin());
  eventVertexPosition_ = GlobalPoint(thePrimaryEventVertex.x(), thePrimaryEventVertex.y(), thePrimaryEventVertex.z());
  Run   = evt.id().run();
  Event = evt.id().event();
  lumi  = evt.luminosityBlock();
  bunch = evt.bunchCrossing();
  IntialValues();
  GenParticleTree(evt);
  PrimaryVertexTree(vertices);
  MuonTree(evt,es,realMuons,thePrimaryEventVertex,vertices);  
  GsfEleTree(evt,es,gsfEle,hConversions,thebs,allowCkfMatch_,
             lxyMin_,probMin_,nHitsBeforeVtxMax_,
             LxyPhoConv, FitProbPhoConv,NHitsBeforeVtxPhoConv,
             rhoIso,vertices);
  //hlt matching
  TriggerMatchingTree(evt,realMuons);
  //==============================================
  //=        End of the main program             =
  //============================================== 
  mytree->Fill();
}

//========================================================================
void MaketreeElectronsMuons::endJob() {
//========================================================================
  // go to *OUR* root file and store histograms
  rootFile_->cd();

  mytree->Write();

  rootFile_->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(MaketreeElectronsMuons);
//=============================================================
//
//         Method for finding good Primary Vertex
//
//=============================================================
bool MaketreeElectronsMuons::PrimaryVertex(const reco::VertexCollection &vertices)
{
  int nbGoodPv = 0;
  bool result = false; 
  for(reco::VertexCollection::const_iterator it=vertices.begin() ; it!=vertices.end() ; ++it)
    {
      if(it->ndof() > minNdof_ && 
	 ( (maxAbsZ_ <= 0.0) || fabs(it->position().z()) <= maxAbsZ_ )  && 
	 ( (maxd0_ <= 0.0) || fabs(it->position().rho()) <= maxd0_ ) ) nbGoodPv++;
    }
  if( nbGoodPv>=NbGoodPv_ ) result = true;
  return result;
}
//=============================================================
//
//                Method for Reco Muon Tree
//
//=============================================================
void MaketreeElectronsMuons::MuonTree(const edm::Event& evt,const edm::EventSetup& es,
				      const reco::MuonCollection* muons,const reco::Vertex &vertex,
				      const reco::VertexCollection &verticess)
{  
  Mu_nbMuon.clear();
  Mu_isGlobalMuon.clear();
  Mu_isPF.clear();
  Mu_isGoodMuon.clear();
  Mu_isTrackerMuon.clear();
  Mu_etaCocktail.clear();
  Mu_ptTunePMuonBestTrack.clear();
  Mu_dPToverPTTunePMuonBestTrack.clear();
  Mu_phi.clear();
  Mu_eta.clear();
  Mu_pt.clear();
  Mu_et.clear();
  Mu_phiCocktail.clear();	
  Mu_en.clear();
  Mu_numberOfMatchedStations.clear();
  Mu_absdxy.clear();
  Mu_absdz.clear();
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
  Mu_chargeCocktail.clear();
  Mu_emIso.clear();
  Mu_hadIso.clear();
  Mu_thetaCocktail.clear();
  Mu_vtxMass.clear();
  Mu_vtxNormChi2.clear();
  Mu_ptcocktail.clear();
  Mu_dPToverPTcocktail.clear();
  Mu_nbofpv.clear();
  Mu_pxCocktail.clear();
  Mu_pyCocktail.clear();
  Mu_pzCocktail.clear();
  Mu_pCocktail.clear();
  Mu_pxTunePMuonBestTrack.clear();
  Mu_pyTunePMuonBestTrack.clear();
  Mu_pzTunePMuonBestTrack.clear();
  Mu_pTunePMuonBestTrack.clear();
  Mu_etaTunePMuonBestTrack.clear();
  Mu_ptInnerTrack.clear();
  Mu_phiTunePMuonBestTrack.clear();
  Mu_thetaTunePMuonBestTrack.clear();
  Mu_chargeTunePMuonBestTrack.clear();
  Mu_absdxyTunePMuonBestTrack.clear();
  Mu_ptDYTTrack.clear();
  Mu_pxDYTTrack.clear();
  Mu_pyDYTTrack.clear();
  Mu_pzDYTTrack.clear();
  Mu_pDYTTrack.clear(); 
  Mu_etaDYTTrack.clear();
  Mu_phiDYTTrack.clear();
  Mu_thetaDYTTrack.clear();
  Mu_chargeDYTTrack.clear();    
  Mu_absdxyDYTTrack.clear();
  Mu_dPToverPTDYTTrack.clear();
  Mu_ptPickyTrack.clear();
  Mu_pxPickyTrack.clear();
  Mu_pyPickyTrack.clear();
  Mu_pzPickyTrack.clear();
  Mu_pPickyTrack.clear(); 
  Mu_etaPickyTrack.clear();
  Mu_phiPickyTrack.clear();
  Mu_thetaPickyTrack.clear();
  Mu_chargePickyTrack.clear();    
  Mu_absdxyPickyTrack.clear();
  Mu_dPToverPTPickyTrack.clear();
  float vtxNormChi2,DiMass;
  Mu_ptMuonBestTrack.clear();
  Mu_dPToverPTMuonBestTrack.clear();
  Mu_pxMuonBestTrack.clear();
  Mu_pyMuonBestTrack.clear();
  Mu_pzMuonBestTrack.clear();
  Mu_pMuonBestTrack.clear(); 
  Mu_etaMuonBestTrack.clear();
  Mu_phiMuonBestTrack.clear();
  Mu_thetaMuonBestTrack.clear();
  Mu_chargeMuonBestTrack.clear();
  Mu_absdxyMuonBestTrack.clear();
  // Get TransientTracks (for use in e.g. the vertex fit) for each of
  // the muon tracks, using e.g. the cocktail momentum.
  edm::ESHandle<TransientTrackBuilder> ttkb;
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", ttkb);
  std::vector<reco::TransientTrack> ttv;
  for (reco::MuonCollection::const_iterator muonIt = muons->begin(); muonIt!=muons->end(); ++muonIt){
    //Varaibles needed for high pt Muons are taken from the twiki below
    //https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#HighPT_Muon
    reco::Muon recoMu = *muonIt;
    if (recoMu.innerTrack().isNull()) continue;
    if (recoMu.globalTrack().isNull()) continue;
    //if(recoMu.isTrackerMuon()==false) continue;
    NbMuons++;
    Mu_nbMuon.push_back(NbMuons);
    Mu_isGlobalMuon.push_back(recoMu.isGlobalMuon());
    Mu_isPF.push_back(recoMu.isPFMuon());
    Mu_isTrackerMuon.push_back(recoMu.isTrackerMuon()); 
    Mu_nbofpv.push_back(verticess.size());
    //============== Parameters related to Kinematics =====================
    reco::TrackRef cktTrack = (muon::tevOptimized(recoMu, 200, 17., 40., 0.25)).first; 
    const reco::TrackRef& tunePMBTrack = recoMu.tunePMuonBestTrack();
    Mu_en.push_back(recoMu.energy());
    Mu_et.push_back(recoMu.et());
    Mu_phi.push_back(recoMu.phi());
    Mu_eta.push_back(recoMu.eta());
    Mu_pt.push_back(recoMu.pt());
    Mu_ptInnerTrack.push_back(recoMu.innerTrack()->pt());
    Mu_ptcocktail.push_back(cktTrack->pt());
    Mu_dPToverPTcocktail.push_back(cktTrack->ptError()/cktTrack->pt());
    Mu_absdxy.push_back(fabs(cktTrack->dxy(vertex.position())));
    Mu_absdz.push_back(fabs(cktTrack->dz(vertex.position())));
    Mu_etaCocktail.push_back(cktTrack->eta());
    Mu_phiCocktail.push_back(cktTrack->phi());
    Mu_thetaCocktail.push_back(cktTrack->theta());
    Mu_chargeCocktail.push_back(cktTrack->charge());
    Mu_pxCocktail.push_back(cktTrack->px()); //px component of the track 
    Mu_pyCocktail.push_back(cktTrack->py()); //py component of the track 
    Mu_pzCocktail.push_back(cktTrack->pz()); //pz component of the track 
    Mu_pCocktail.push_back(cktTrack->p());   //magnitude of momentum vector
    // part for TuneP Muon Best track
    Mu_ptTunePMuonBestTrack.push_back(recoMu.tunePMuonBestTrack()->pt());
    Mu_pxTunePMuonBestTrack.push_back(recoMu.tunePMuonBestTrack()->px()); //px component of the track 
    Mu_pyTunePMuonBestTrack.push_back(recoMu.tunePMuonBestTrack()->py()); //py component of the track 
    Mu_pzTunePMuonBestTrack.push_back(recoMu.tunePMuonBestTrack()->pz()); //pz component of the track 
    Mu_pTunePMuonBestTrack.push_back(recoMu.tunePMuonBestTrack()->p());   //magnitude of momentum vector 
    Mu_etaTunePMuonBestTrack.push_back(recoMu.tunePMuonBestTrack()->eta());
    Mu_phiTunePMuonBestTrack.push_back(recoMu.tunePMuonBestTrack()->phi());
    Mu_thetaTunePMuonBestTrack.push_back(recoMu.tunePMuonBestTrack()->theta());
    Mu_chargeTunePMuonBestTrack.push_back(recoMu.tunePMuonBestTrack()->charge());
    Mu_dPToverPTTunePMuonBestTrack.push_back(recoMu.tunePMuonBestTrack()->ptError()/recoMu.tunePMuonBestTrack()->pt());
    Mu_absdxyTunePMuonBestTrack.push_back(fabs(recoMu.tunePMuonBestTrack()->dxy(vertex.position())));
    // part for Muon Best track
    Mu_ptMuonBestTrack.push_back(recoMu.muonBestTrack()->pt());
    Mu_pxMuonBestTrack.push_back(recoMu.muonBestTrack()->px()); //px component of the track 
    Mu_pyMuonBestTrack.push_back(recoMu.muonBestTrack()->py()); //py component of the track 
    Mu_pzMuonBestTrack.push_back(recoMu.muonBestTrack()->pz()); //pz component of the track 
    Mu_pMuonBestTrack.push_back(recoMu.muonBestTrack()->p());   //magnitude of momentum vector 
    Mu_etaMuonBestTrack.push_back(recoMu.muonBestTrack()->eta());
    Mu_phiMuonBestTrack.push_back(recoMu.muonBestTrack()->phi());
    Mu_thetaMuonBestTrack.push_back(recoMu.muonBestTrack()->theta());
    Mu_chargeMuonBestTrack.push_back(recoMu.muonBestTrack()->charge());
    Mu_dPToverPTMuonBestTrack.push_back(recoMu.muonBestTrack()->ptError()/recoMu.muonBestTrack()->pt());
    Mu_absdxyMuonBestTrack.push_back(fabs(recoMu.muonBestTrack()->dxy(vertex.position())));
    // part for DYT track
    /*
      Mu_ptDYTTrack.push_back(recoMu.dytTrack()->pt());
      Mu_pxDYTTrack.push_back(recoMu.dytTrack()->px());
      Mu_pyDYTTrack.push_back(recoMu.dytTrack()->py());
      Mu_pzDYTTrack.push_back(recoMu.dytTrack()->pz());
      Mu_pDYTTrack.push_back(recoMu.dytTrack()-> p());
      Mu_etaDYTTrack.push_back(recoMu.dytTrack()->eta());
      Mu_phiDYTTrack.push_back(recoMu.dytTrack()->phi());
      Mu_thetaDYTTrack.push_back(recoMu.dytTrack()->theta());
      Mu_chargeDYTTrack.push_back(recoMu.dytTrack()->charge());    
      Mu_absdxyDYTTrack.push_back(fabs(recoMu.dytTrack()->dxy(vertex.position())));
      Mu_dPToverPTDYTTrack.push_back(recoMu.dytTrack()->ptError()/recoMu.dytTrack()->pt());
    */
    // part for Picky track                          
    /*
      Mu_ptDYTTrack.push_back(recoMu.pickyTrack()->pt());
      Mu_pxDYTTrack.push_back(recoMu.pickyTrack()->px());
      Mu_pyDYTTrack.push_back(recoMu.pickyTrack()->py());
      Mu_pzDYTTrack.push_back(recoMu.pickyTrack()->pz());
      Mu_pDYTTrack.push_back(recoMu.pickyTrack()-> p());
      Mu_etaDYTTrack.push_back(recoMu.pickyTrack()->eta());
      Mu_phiDYTTrack.push_back(recoMu.pickyTrack()->phi());
      Mu_thetaDYTTrack.push_back(recoMu.pickyTrack()->theta());
      Mu_chargeDYTTrack.push_back(recoMu.pickyTrack()->charge());    
      Mu_absdxyDYTTrack.push_back(fabs(recoMu.pickyTrack()->dxy(vertex.position())));
      Mu_dPToverPTDYTTrack.push_back(recoMu.pickyTrack()->ptError()/recoMu.pickyTrack()->pt());
    */
    //====================== Parameters related to track quality =====================
    Mu_normalizedChi2.push_back(recoMu.globalTrack()->normalizedChi2());
    Mu_numberOfValidPixelHits.push_back(recoMu.globalTrack()->hitPattern().numberOfValidPixelHits());
    Mu_numberOfValidMuonHits.push_back(recoMu.globalTrack()->hitPattern().numberOfValidMuonHits());
    Mu_numberOftrackerLayersWithMeasurement.push_back(recoMu.globalTrack()->hitPattern().trackerLayersWithMeasurement());
    Mu_numberOfMatchedStations.push_back(recoMu.numberOfMatchedStations());
    //============= Parameters related to detector isolation =====================
    Mu_emIso.push_back(recoMu.isolationR03().emEt);
    Mu_hadIso.push_back(recoMu.isolationR03().hadEt);
    Mu_trackiso.push_back(recoMu.isolationR03().sumPt);
    //============= Parameters related to PF isolation =====================
    Mu_pfSumChargedHadronPt.push_back(recoMu.pfIsolationR03().sumChargedHadronPt);
    Mu_pfSumNeutralHadronEt.push_back(recoMu.pfIsolationR03().sumNeutralHadronEt);
    Mu_PFSumPhotonEt.push_back(recoMu.pfIsolationR03().sumPhotonEt);
    Mu_pfSumPUPt.push_back(recoMu.pfIsolationR03().sumPUPt);
    //------------------------------------------------------
    //   part to compute mass from vertex info             -   
    //   Given to me by Piotr                              -
    //------------------------------------------------------
    if( NbMuons > 0 && recoMu.tunePMuonBestTrack()->pt()>45.0 && 
	(recoMu.isolationR03().sumPt/recoMu.innerTrack()->pt()<0.10) && 
    	recoMu.tunePMuonBestTrack()->ptError()/recoMu.tunePMuonBestTrack()->pt()<0.3 && 
        fabs(recoMu.tunePMuonBestTrack()->dxy(vertex.position()))<0.2 &&
        recoMu.globalTrack()->hitPattern().trackerLayersWithMeasurement()>5 && 
        recoMu.globalTrack()->hitPattern().numberOfValidPixelHits()>0 &&
        recoMu.globalTrack()->hitPattern().numberOfValidMuonHits()>0 &&
        recoMu.numberOfMatchedStations()>1 ){
      //ttv.push_back(ttkb->build(cktTrack));
      ttv.push_back(ttkb->build(tunePMBTrack));
      //cout<<"size(outside) = "<<ttv.size()<<endl;
      if(ttv.size()>1){
	//if(ttv.size()==2){
	KalmanVertexFitter kvf(true); // false means no smoothing which means no track re-fit
	CachingVertex<5> vtx = kvf.vertex(ttv);
	vtxNormChi2 = vtx.totalChiSquared()/vtx.degreesOfFreedom();
	InvariantMassFromVertex imfv;
	static const double muon_mass = 0.1056583;
	Measurement1D mass = imfv.invariantMass(vtx, muon_mass);
	DiMass = mass.value();
	Mu_vtxMass.push_back(DiMass);
	Mu_vtxNormChi2.push_back(vtxNormChi2);
	nbTk++;
      }
    }
  }
}
//=============================================================
//
//            Method for Genrated Particles Tree
//
//=============================================================
void MaketreeElectronsMuons::GenParticleTree(const edm::Event& evt){
  edm::Handle<GenParticleCollection> genParticles;
  evt.getByLabel(genParticlesColl_, genParticles);
  //if (!(genParticles.isValid())) return;
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
  NbOfDaughters.clear();
  McZmass.clear();
  McZpt.clear();
  //massGen.clear();
  //vxGen.clear();
  //vyGen.clear();
  //vzGen.clear();
  //McZpx.clear();
  //McZpy.clear();
  //McZpz.clear();
  //McZen.clear();
  for(size_t i = 0; i < genParticles->size(); ++ i) {
    const GenParticle & p = (*genParticles)[i];
    if( fabs(p.pdgId()) > 13 ) continue; 
    if( p.status() > 3 ) continue; //consider only status 1,3
    //if( p.numberOfDaughters() < 2 ) continue;
    value2_++;
    iGen.push_back(value2_);
    math::XYZTLorentzVector Zboson;
    Zboson = p.p4();
    idGen.push_back(p.pdgId());           //int id = 
    statusGen.push_back(p.status());      //int st = 
    ptGen.push_back(p.pt());
    etaGen.push_back(p.eta());
    phiGen.push_back(p.phi());
    chargeGen.push_back(p.charge());
    EnergyGen.push_back(p.energy());
    pxGen.push_back(p.px());
    pyGen.push_back(p.py());
    pzGen.push_back(p.pz());
    NbOfDaughters.push_back(p.numberOfDaughters());
    McZmass.push_back(Zboson.M()); 
    //ROOT::Math::Boost CoMBoost(Zboson.BoostToCM());
    //const Candidate * mom = p.mother();
    //massGen.push_back(p.mass());
    //vxGen.push_back(p.vx());
    //vyGen.push_back(p.vy());
    //vzGen.push_back(p.vz());
    //McZpt.push_back(Zboson.Pt()); 
    //McZpx.push_back(Zboson.Px()); 
    //McZpy.push_back(Zboson.Py()); 
    //McZpz.push_back(Zboson.Pz());
    //McZen.push_back(Zboson.E());
  }
}
//=============================================================
//
//            Method for Pimary Vertices Tree
//
//=============================================================
void MaketreeElectronsMuons::PrimaryVertexTree(const reco::VertexCollection &vertices)
{
  nbPv.clear();
  Nbdof.clear();
  PositionX.clear();
  PositionY.clear();
  PositionZ.clear();
  PositionRho.clear();
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
//            Method for Trigger Matching Tree
//
//=============================================================
void MaketreeElectronsMuons::TriggerMatchingTree(const edm::Event& iEvent,const reco::MuonCollection* muons){
  MuHLTMatch_pt.clear();
  MuHLTMatch_eta.clear();
  MuHLTMatch_phi.clear();
  MuHLTMatch_nbMuonMatchHLT.clear();
  edm::Handle<trigger::TriggerEvent> handleTriggerEvent;
  iEvent.getByLabel(triggerEvent, handleTriggerEvent );
  const trigger::TriggerObjectCollection & toc(handleTriggerEvent->getObjects());
  size_t nMuHLT =0;
  for ( size_t ia = 0; ia < handleTriggerEvent->sizeFilters(); ++ ia) {
    const trigger::Keys & k = handleTriggerEvent->filterKeys(ia);
    for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) {
      std::string fullname = handleTriggerEvent->filterTag(ia).encode();
      //std::cout << " fullname == " << fullname << std::endl;
      std::string name;
      size_t p = fullname.find_first_of(':');
      if (p != std::string::npos) {
	name = fullname.substr(0, p);
      } else {
	name = fullname;
      }
      if (name == triggerFilter.c_str()) {
	//if (name == "HLTFilterName_") {
	nMuHLT++;
	MuHLTMatch_nbMuonMatchHLT.push_back(nMuHLT);
	MuHLTMatch_pt.push_back(toc[*ki].pt());
	MuHLTMatch_eta.push_back(toc[*ki].eta());
	MuHLTMatch_phi.push_back(toc[*ki].phi());
      } 
    }
  }
}

//=============================================================
//
//     Method for initializing values for the variables
//
//=============================================================
void MaketreeElectronsMuons::IntialValues()
{
  
  NbMuons   = 0;
  value2_   = 0;
  value3_   = 0;
  value_    = 0;
  nbTk      = 0;
}
//=============================================================
//
//                Method for Reco GSF Ele Tree
//
//=============================================================
void MaketreeElectronsMuons::GsfEleTree(const edm::Event& evt, const edm::EventSetup& es,const reco::GsfElectronCollection* GsfEle,
					const edm::Handle<reco::ConversionCollection> &convCol,const reco::BeamSpot &bspot,
					bool allowCkfMatch, float lxyMin, float probMin, unsigned int nHitsBeforeVtxMax,
					float LxyConv, float FitProbConv, unsigned int NHitsBeforeVtxConv,float RhoIsoValue,
					const reco::VertexCollection &vtxs)
{
  PreshowerEn1.clear();
  EtaEle1.clear();
  PhiEle1.clear();
  EtaEleSC1.clear();  
  PhiEleSC1.clear();
  EtEle1.clear();
  EtFromCaloEnEle1.clear();
  EtCorrEle1.clear();
  EnEleSC1.clear();
  EnCorrEleSC1.clear();
  scEmaxEle1.clear(); 
  scE25Ele1.clear();
  scE2x5RightEle1.clear();
  scE2x5LeftEle1.clear();
  scE2x5TopEle1.clear();
  scE2x5BottomEle1.clear();
  EleClassEle1.clear();
  ThetaEleSC1.clear();
  ThetaEle1.clear();
  chargeEle1.clear();
  DeltaEtaEle1.clear();
  DeltaPhiEle1.clear();
  HoeEle1.clear();
  //HoeBcEle1.clear();
  EcalPlusHcal1IsoEle1.clear();
  EcalPlusHcal1BcIsoEle1.clear();
  SigmaIetaIetaEle1.clear();
  HcalDepth2IsoEle1.clear();
  TkSumPtIsoEle1.clear();
  E2x5MaxOverE5x5Ele1.clear();
  E1x5OverE5x5Ele1.clear();
  fbremEle1.clear();
  EnSeedClusterOverPEle1.clear();
  EcalIsoEle1.clear();
  Hcal1IsoEle1.clear();
  Hcal1BcIsoEle1.clear();
  xSCEle1.clear();
  ySCEle1.clear();
  EcalDrivenSeedEle1.clear();
  NbOfLostInnerHitsEle1.clear();
  E2x5MaxEle1.clear();
  E1x5Ele1.clear();
  iEle.clear();
  //DistEle1.clear();
  //DcotEle1.clear();
  //RadiusEle1.clear();
  //gsftrackD0Ele1.clear();
  //gsftrackPhiEle1.clear();
  //vxEle1.clear();
  //vyEle1.clear();
  //vzEle1.clear();
  ElectronFromConv.clear();
  LxyMinFromConv.clear();
  ProbMinFromConv.clear();
  NHitsBeforeVtxMaxFromConv.clear();
  Rho.clear();
  pxEle1.clear();
  pyEle1.clear();
  pzEle1.clear();
  gsftrackDxyEle1.clear();
  gsftrackDzEle1.clear();
  gsftrackDxyVtxEle1.clear();
  gsftrackDzVtxEle1.clear();
  gsftrackDxyErrorVtxEle1.clear();
  gsftrackDzErrorVtxEle1.clear();
  scEnLeft.clear();
  scEnRight.clear();
  scEnTop.clear();
  scEnBottom.clear();
  scSeedTime.clear();
  //prepare electron cluster shapes extraction
  //std::auto_ptr<EcalClusterLazyTools> lazyTools_;
  //lazyTools_ .reset(new EcalClusterLazyTools( evt , es , reducedBarrelRecHitCollection_ , reducedEndcapRecHitCollection_ ));  
  for(reco::GsfElectronCollection::const_iterator iElectron = GsfEle->begin();  iElectron != GsfEle->end(); iElectron++ ) { 
    bool passconversionveto = hasMatchedConversion(*iElectron,convCol,bspot.position(),
                                                   allowCkfMatch,lxyMin,probMin,
                                                   nHitsBeforeVtxMax);
    for (ConversionCollection::const_iterator it = convCol->begin(); it!=convCol->end(); ++it) {
      allConversion(*it,bspot.position(),LxyConv,FitProbConv,NHitsBeforeVtxConv);
    }
    value_++;
    iEle.push_back(value_);
    //DetId xxx = (*iElectron).superCluster()->seed()->seed();
    //EBDetId elementId = GsfElectron::superCluster()->seed()->seed()->id();
    //std::cout << "ieta =  " << xxx.rawId() << std::endl;
    PreshowerEn1.push_back((*iElectron).superCluster()->preshowerEnergy());
    EtaEle1.push_back((*iElectron).p4().eta());
    PhiEle1.push_back((*iElectron).p4().phi());
    //EtaEleSC1.push_back((*iElectron).caloPosition().eta());
    //PhiEleSC1.push_back((*iElectron).caloPosition().phi());
    EtaEleSC1.push_back((*iElectron).superCluster()->eta());
    PhiEleSC1.push_back((*iElectron).superCluster()->phi());
    //EnEleSC1.push_back((*iElectron).caloEnergy());
    EnEleSC1.push_back((*iElectron).superCluster()->energy());
    EnCorrEleSC1.push_back((*iElectron).correctedEcalEnergy());
    EtFromCaloEnEle1.push_back((*iElectron).caloEnergy() * sin((*iElectron).p4().theta()));
    EtEle1.push_back((*iElectron).superCluster()->energy() * sin((*iElectron).p4().theta()));
    EtCorrEle1.push_back((*iElectron).correctedEcalEnergy() * sin((*iElectron).p4().theta()));
    /*
    scEmaxEle1.push_back(lazyTools_->eMax(*((*iElectron).superCluster()->seed()))) ;
    scEmaxEle1.push_back(lazyTools_->eMax(*((*iElectron).superCluster()->seed()))) ;
    scE25Ele1.push_back(lazyTools_->e5x5(*((*iElectron).superCluster()->seed())));
    scE2x5RightEle1.push_back(lazyTools_->e2x5Right(*((*iElectron).superCluster()->seed())));
    scE2x5LeftEle1.push_back(lazyTools_->e2x5Left(*((*iElectron).superCluster()->seed())));
    scE2x5TopEle1.push_back(lazyTools_->e2x5Top(*((*iElectron).superCluster()->seed())));
    scE2x5BottomEle1.push_back(lazyTools_->e2x5Bottom(*((*iElectron).superCluster()->seed())));
    scEnLeft.push_back(lazyTools_->eLeft(*(iElectron->superCluster()->seed()))) ;
    scEnRight.push_back(lazyTools_->eRight(*(iElectron->superCluster()->seed()))) ;
    scEnTop.push_back(lazyTools_->eTop(*(iElectron->superCluster()->seed()))) ;
    scEnBottom.push_back(lazyTools_->eBottom(*(iElectron->superCluster()->seed()))) ;
    scSeedTime.push_back(lazyTools_->BasicClusterSeedTime(*(iElectron->superCluster()->seed()))) ;
    //gsftrackDxyEle1.push_back((*iElectron).gsfTrack()->dxy(bspot.position()));
    //gsftrackDzEle1.push_back((*iElectron).gsfTrack()->dz(bspot.position()));
    */
    EleClassEle1.push_back((*iElectron).classification());
    ThetaEleSC1.push_back((*iElectron).caloPosition().theta());
    ThetaEle1.push_back((*iElectron).p4().theta());
    pxEle1.push_back((*iElectron).p4().px());
    pyEle1.push_back((*iElectron).p4().py());
    pzEle1.push_back((*iElectron).p4().pz());
    chargeEle1.push_back((*iElectron).gsfTrack()->charge());
    DeltaEtaEle1.push_back((*iElectron).deltaEtaSuperClusterTrackAtVtx());
    DeltaPhiEle1.push_back((*iElectron).deltaPhiSuperClusterTrackAtVtx());
    HoeEle1.push_back((*iElectron).hadronicOverEm());
    //HoeBcEle1.push_back((*iElectron).hcalOverEcalBc());
    EcalPlusHcal1IsoEle1.push_back((*iElectron).dr03EcalRecHitSumEt() + (*iElectron).dr03HcalDepth1TowerSumEt());
    EcalPlusHcal1BcIsoEle1.push_back((*iElectron).dr03EcalRecHitSumEt() + (*iElectron).dr03HcalDepth1TowerSumEtBc());
    SigmaIetaIetaEle1.push_back((*iElectron).sigmaIetaIeta());
    HcalDepth2IsoEle1.push_back((*iElectron).dr03HcalDepth2TowerSumEt());
    TkSumPtIsoEle1.push_back((*iElectron).dr03TkSumPt());
    E2x5MaxOverE5x5Ele1.push_back((*iElectron).e2x5Max()/(*iElectron).e5x5()); 
    E1x5OverE5x5Ele1.push_back((*iElectron).e1x5()/(*iElectron).e5x5()); 
    fbremEle1.push_back((*iElectron).fbrem());
    EnSeedClusterOverPEle1.push_back((*iElectron).eSeedClusterOverP());
    EcalIsoEle1.push_back((*iElectron).dr03EcalRecHitSumEt());
    Hcal1IsoEle1.push_back((*iElectron).dr03HcalDepth1TowerSumEt());
    Hcal1BcIsoEle1.push_back((*iElectron).dr03HcalDepth1TowerSumEtBc());
    xSCEle1.push_back((*iElectron).caloPosition().x());   
    ySCEle1.push_back((*iElectron).caloPosition().y());
    EcalDrivenSeedEle1.push_back((*iElectron).ecalDrivenSeed());
    //NbOfLostInnerHitsEle1.push_back((*iElectron).gsfTrack()->trackerExpectedHitsInner().numberOfLostHits());
    NbOfLostInnerHitsEle1.push_back((*iElectron).gsfTrack()->numberOfLostHits());
    E2x5MaxEle1.push_back((*iElectron).e2x5Max());
    E1x5Ele1.push_back((*iElectron).e1x5());
    ElectronFromConv.push_back(passconversionveto);
    ProbMinFromConv.push_back(FitProbConv);
    LxyMinFromConv.push_back(LxyConv);
    NHitsBeforeVtxMaxFromConv.push_back(NHitsBeforeVtxConv);
    Rho.push_back(RhoIsoValue);
    if (vtxs.size() >0){
      gsftrackDxyVtxEle1.push_back((*iElectron).gsfTrack()->dxy(vtxs.front().position()));   
      gsftrackDzVtxEle1.push_back((*iElectron).gsfTrack()->dz(vtxs.front().position()));
      gsftrackDxyErrorVtxEle1.push_back((*iElectron).gsfTrack()->dxyError());
      gsftrackDzErrorVtxEle1.push_back((*iElectron).gsfTrack()->dzError());
    }else {
      gsftrackDxyVtxEle1.push_back((*iElectron).gsfTrack()->dxy());
      gsftrackDzVtxEle1.push_back((*iElectron).gsfTrack()->dz());
      gsftrackDxyErrorVtxEle1.push_back((*iElectron).gsfTrack()->dxyError());
      gsftrackDzErrorVtxEle1.push_back((*iElectron).gsfTrack()->dzError());
    }
  }
}

bool MaketreeElectronsMuons::matchesConversion(const reco::GsfElectron &ele, const reco::Conversion &conv, bool allowCkfMatch)
{
  
  //check if a given GsfElectron matches a given conversion (no quality cuts applied)
  //matching is always attempted through the gsf track ref, and optionally attempted through the
  //closest ctf track ref
  const std::vector<edm::RefToBase<reco::Track> > &convTracks = conv.tracks();
  for (std::vector<edm::RefToBase<reco::Track> >::const_iterator it=convTracks.begin(); it!=convTracks.end(); ++it) {
    if ( ele.gsfTrack().isNonnull() && ele.gsfTrack().id()==it->id() && ele.gsfTrack().key()==it->key()) return true;
    else if ( allowCkfMatch && ele.closestCtfTrackRef().isNonnull() && ele.closestCtfTrackRef().id()==it->id() && ele.closestCtfTrackRef().key()==it->
	      key() ) return true;
  }
  return false;
}
bool MaketreeElectronsMuons::isGoodConversion(const reco::Conversion &conv, const math::XYZPoint &beamspot, float lxyMin, float probMin, unsigned int nHitsBeforeVtxMax)
{
  //Check if a given conversion candidate passes the conversion selection cuts
  const reco::Vertex &vtx = conv.conversionVertex();
  //vertex validity
  if (!vtx.isValid()) return false;
  //fit probability
  if (TMath::Prob( vtx.chi2(),  vtx.ndof() )<probMin) return false;
  //compute transverse decay length
  math::XYZVector mom(conv.refittedPairMomentum()); 
  double dbsx = vtx.x() - beamspot.x();
  double dbsy = vtx.y() - beamspot.y();
  double lxy = (mom.x()*dbsx + mom.y()*dbsy)/mom.rho();
  //transverse decay length  
  if ( lxy<lxyMin )
    return false;
  //loop through daughters to check nhitsbeforevtx
  for (std::vector<uint8_t>::const_iterator it = conv.nHitsBeforeVtx().begin(); it!=conv.nHitsBeforeVtx().end(); ++it) {
    if ( (*it)>nHitsBeforeVtxMax ) return false;
  }
  return true;
}


bool MaketreeElectronsMuons::hasMatchedConversion(const reco::GsfElectron &ele,
						  const edm::Handle<reco::ConversionCollection> &convCol,
						  const math::XYZPoint &beamspot, bool allowCkfMatch, float lxyMin, 
						  float probMin, unsigned int nHitsBeforeVtxMax)
{
  //check if a given electron candidate matches to at least one conversion candidate in the
  //collection which also passes the selection cuts, optionally match with the closestckf track in
  //in addition to just the gsf track (enabled in default arguments)
  for (ConversionCollection::const_iterator it = convCol->begin(); it!=convCol->end(); ++it) {
    if (!matchesConversion(ele, *it, allowCkfMatch)) continue;
    if (!isGoodConversion(*it,beamspot,lxyMin,probMin,nHitsBeforeVtxMax)) continue;
    return true;
  }
  return false;
}

void MaketreeElectronsMuons::allConversion(const reco::Conversion &conv, const math::XYZPoint &beamspot, 
					   float &LxyConv, float &FitProbConv, unsigned int &NHitsBeforeVtxConv)
{
  //Check if a given conversion candidate passes the conversion selection cuts
  const reco::Vertex &vtx = conv.conversionVertex();
  //vertex validity
  //fit probability
  FitProbConv = TMath::Prob( vtx.chi2(),  vtx.ndof() );
  //compute transverse decay length
  math::XYZVector mom(conv.refittedPairMomentum()); 
  double dbsx = vtx.x() - beamspot.x();
  double dbsy = vtx.y() - beamspot.y();
  double lxy = (mom.x()*dbsx + mom.y()*dbsy)/mom.rho();
  LxyConv = lxy;
  //loop through daughters to check nhitsbeforevtx
  for (std::vector<uint8_t>::const_iterator it = conv.nHitsBeforeVtx().begin(); it!=conv.nHitsBeforeVtx().end(); ++it) {
    NHitsBeforeVtxConv = (*it);
  }
}
