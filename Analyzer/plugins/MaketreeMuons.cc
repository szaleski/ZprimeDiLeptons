//=========================================================================  
//      Make the root tree for Z boson to Mu Mu analysis                  =  
//                                                                        =
//                         CMSSW_7_4_5                                    =
//                                                                        =
//       Written by Sherif Elgammal & Nicola De Filippis                  =
//                                                                        =
//                         02/08/2015                                     =
//=========================================================================
#include "ZprimeDiLeptons/Analyzer/interface/MaketreeMuons.h"
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
#include "DataFormats/JetReco/interface/PFJet.h"
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
#include <string>

//========================================================================
MaketreeMuons::MaketreeMuons( const edm::ParameterSet& ps )
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
  triggerEvent                   = ps.getParameter<edm::InputTag>("triggerEvent");
  triggerFilter                  = ps.getParameter<vector<string> >("triggerFilter");
  //inputTag_                    = ps.getParameter<edm::InputTag> ("TriggerResultsTag");
  thePFMETCollectionToken_       = ps.getParameter<edm::InputTag>("thePFMETCollectionToken");
  theMETSignificance_            = ps.getParameter<edm::InputTag>("METSignificance");
  thejetsTag_                    = ps.getParameter<edm::InputTag>("Jets");
  PileupSrc_                     = ps.getParameter<edm::InputTag>("PileupSrc");
  rhoIsoInputTag_                = ps.getParameter<edm::InputTag>("rhoIsoInputTag");
  genEventInfo_                  = ps.getParameter<edm::InputTag>("genEventInfo");
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
  //===============================================================================================
  outputFile_ = ps.getParameter<std::string>("outputFile");
  rootFile_   = TFile::Open(outputFile_.c_str(),"RECREATE"); // open output file to store histograms
}

 
//========================================================================
MaketreeMuons::~MaketreeMuons()
//========================================================================    
{
  delete rootFile_;
}

//========================================================================
void MaketreeMuons::beginJob() {
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
  mytree->Branch("Mu_vtxMassLept",&Mu_vtxMassLept);
  mytree->Branch("Mu_numberOfMatchedStations",&Mu_numberOfMatchedStations);
  mytree->Branch("Mu_numberOfValidPixelHits",&Mu_numberOfValidPixelHits);
  mytree->Branch("Mu_numberOfValidMuonHits",&Mu_numberOfValidMuonHits);
  mytree->Branch("Mu_numberOftrackerLayersWithMeasurement",&Mu_numberOftrackerLayersWithMeasurement);
  mytree->Branch("Mu_innerTK_numberOfValidPixelHits",&Mu_innerTK_numberOfValidPixelHits);
  mytree->Branch("Mu_innerTK_numberOfValidMuonHits",&Mu_innerTK_numberOfValidMuonHits);
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
  //           Create Branches for jets variables
  //
  mytree->Branch("jetnumber",&jetnumber,"jetnumber/I");
  mytree->Branch("jetcharge",&jetcharge);
  mytree->Branch("jetet",&jetet);
  mytree->Branch("jetpt",&jetpt);
  mytree->Branch("jeteta",&jeteta);
  mytree->Branch("jetphi",&jetphi);

  //=============================================================
  //                   
  //           Create Branches for Muons match HLT variables
  //
  //=============================================================
  mytree->Branch("MuHLTMatch_nbMuonMatchHLT",&MuHLTMatch_nbMuonMatchHLT);
  mytree->Branch("MuHLTMatch_pt",&MuHLTMatch_pt);
  mytree->Branch("MuHLTMatch_eta",&MuHLTMatch_eta);
  mytree->Branch("MuHLTMatch_phi",&MuHLTMatch_phi);
  mytree->Branch("MuHLTMatch_nbHLTMuonMatchReco",&MuHLTMatch_nbHLTMuonMatchReco);
  mytree->Branch("MuHLTMatch_Trigger_pt",&MuHLTMatch_Trigger_pt);
  mytree->Branch("MuHLTMatch_Trigger_eta",&MuHLTMatch_Trigger_eta);
  mytree->Branch("MuHLTMatch_Trigger_phi",&MuHLTMatch_Trigger_phi);
  //=============================================================
  //                   
  //           Create Branches for PF MET
  //
  //=============================================================
  mytree->Branch("PFMet_pt",&PFMet_pt,"PFMet_pt/F");
  mytree->Branch("PFMet_eta",&PFMet_eta,"PFMet_eta/F");
  mytree->Branch("PFMet_phi",&PFMet_phi,"PFMet_phi/F");
  mytree->Branch("PFMet_en",&PFMet_en,"PFMet_en/F");
  mytree->Branch("PFMet_px",&PFMet_px,"PFMet_px/F");
  mytree->Branch("PFMet_py",&PFMet_py,"PFMet_py/F");
  mytree->Branch("PFMet_pz",&PFMet_pz,"PFMet_pz/F");
  mytree->Branch("PFMet_sumEt",&PFMet_sumEt,"PFMet_sumEt/F");
  mytree->Branch("METSign",&METSign,"METSign/D");
  //=============================================================
  //                   
  //           Create Branchs for PileUp tree
  //
  //=============================================================
  mytree->Branch("num_PU_vertices",&num_PU_vertices,"num_PU_vertices/I");
  mytree->Branch("PU_BunchCrossing",&PU_BunchCrossing,"PU_BunchCrossing/I");
  //=============================================================
  //                   
  //           Create Branch for Rho
  //
  //=============================================================
  mytree->Branch("Rho",&Rho);
  //=============================================================
  //                   
  //           Create Branch for events reweighting
  //
  //=============================================================
  mytree->Branch("MC_weighting",&MC_weighting);
  
  
}
  

//========================================================================================
void MaketreeMuons::analyze( const edm::Event& evt, const edm::EventSetup& es ) {
//======================================================================================
  using namespace edm; // needed for all fwk related classe
  //==============================================
  //=        Begin of the main program           =
  //============================================== 

  TotalNbEvents ++;
  NbEventsPassTrigger.clear();
  NbEventsPassTrigger.push_back(TotalNbEvents);
  //===================== Handle For Muons =====================
  edm::Handle<reco::MuonCollection> muons;
  evt.getByLabel(globalMuons_, muons);
  const reco::MuonCollection* realMuons = muons.product();
  //cout << "Ciao" << endl;
  //===================== Handle For Primary Vertics ===========
  edm::Handle<reco::VertexCollection> pvHandle; 
  evt.getByLabel(vertexSrc,pvHandle);
  const reco::VertexCollection &vertices = *pvHandle.product();
  bool GoodPv  = PrimaryVertex(vertices);
  //cout << "At least one good primary vertex" << GoodPv << endl;
  if( GoodPv == 0 ) return;    //the aim of that command is to select only events
                               // with one good reconstructed primary vertics
  NbEventsPass ++;
  NbEventsPassTriggerandPVcond.clear();
  NbEventsPassTriggerandPVcond.push_back(NbEventsPass);
  eventVertexPosition_ = GlobalPoint(0., 0., 0.);
  const reco::Vertex& thePrimaryEventVertex = (*pvHandle->begin());
  eventVertexPosition_ = GlobalPoint(thePrimaryEventVertex.x(), thePrimaryEventVertex.y(), thePrimaryEventVertex.z());

  //===================== Handle For jets ===========
  edm::Handle<reco::PFJetCollection> jetsHandle;
  evt.getByLabel(thejetsTag_,jetsHandle);
  const reco::PFJetCollection &pfjets = *jetsHandle.product();
  
  // Run/Event/Lumi/bunch
  Run   = evt.id().run();
  Event = evt.id().event();
  lumi  = evt.luminosityBlock();
  bunch = evt.bunchCrossing();
  cout << "Run= " << Run << " Lumi= "<< lumi << " Event= " << Event << endl;
  IntialValues();
  GenParticleTree(evt);
  PrimaryVertexTree(vertices);
  JetsTree(pfjets);
  MuonTree(evt,es,realMuons,thePrimaryEventVertex,vertices);  
  ComputeMuonMassVtx(evt,es,realMuons,thePrimaryEventVertex);
  pfMETtree(evt);
  fillPU(evt);
  fillRho(evt);
  EventsReWeighting(evt);
  //hlt matching
  TriggerMatchingTree(evt,realMuons,thePrimaryEventVertex);
  //==============================================
  //=        End of the main program             =
  //============================================== 
  //mytree->SetWeight(w);
  mytree->Fill();
}

//========================================================================
void MaketreeMuons::endJob() {
//========================================================================
  // go to *OUR* root file and store histograms
  rootFile_->cd();
 
  mytree->Write();

  rootFile_->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(MaketreeMuons);

//=============================================================
//
//         Method for finding good Primary Vertex
//
//=============================================================
bool MaketreeMuons::PrimaryVertex(const reco::VertexCollection &vertices)
{
  int nbGoodPv = 0;
  bool result = false; 
  for(reco::VertexCollection::const_iterator it=vertices.begin() ; it!=vertices.end() ; ++it)
    {
      //cout << it->ndof() << " " << fabs(it->position().z()) << " " << fabs(it->position().rho()) << endl;
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
void MaketreeMuons::MuonTree(const edm::Event& evt,const edm::EventSetup& es,
			     const reco::MuonCollection* muons,const reco::Vertex &vertex,
			     const reco::VertexCollection &verticess)
{ 
  Mu_nbMuon.clear();
  Mu_isGlobalMuon.clear();
  Mu_isPF.clear();
  Mu_isTrackerMuon.clear();
  Mu_isGoodMuon.clear();
  Mu_etaCocktail.clear();
  Mu_ptTunePMuonBestTrack.clear();
  Mu_dPToverPTTunePMuonBestTrack.clear();
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
  Mu_innerTK_numberOfValidPixelHits.clear();
  Mu_innerTK_numberOfValidMuonHits.clear();
  Mu_calEnergy.clear();
  Mu_chargeCocktail.clear();
  Mu_emIso.clear();
  Mu_hadIso.clear();
  Mu_thetaCocktail.clear();
  //Mu_vtxMass.clear();
  //Mu_vtxNormChi2.clear();
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
  //float vtxNormChi2,DiMass;
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
  reco::TrackRef cktTrack;
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
    //reco::TrackRef cktTrack = (muon::tevOptimized(recoMu, 200, 17., 40., 0.25)).first; 
    //reco::TrackRef cktTrack = recoMu.innerTrack();
    // 2015B
    //if(fabs(recoMu.innerTrack()->eta())<1.2) cktTrack = (muon::tevOptimized(recoMu, 200, 17., 40., 0.25)).first;
    //else if(fabs(recoMu.innerTrack()->eta())>1.2) cktTrack = recoMu.innerTrack();
    // 2015C and D
    cktTrack = (muon::tevOptimized(recoMu, 200, 17., 40., 0.25)).first;
    Mu_en.push_back(recoMu.energy());
    Mu_et.push_back(recoMu.et());
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
    Mu_innerTK_numberOfValidPixelHits.push_back(recoMu.innerTrack()->hitPattern().numberOfValidPixelHits());
    Mu_innerTK_numberOfValidMuonHits.push_back(recoMu.innerTrack()->hitPattern().numberOfValidMuonHits());
    //============= Parameters related to detector isolation =====================
    Mu_emIso.push_back(recoMu.isolationR03().emEt);
    Mu_hadIso.push_back(recoMu.isolationR03().hadEt);
    Mu_trackiso.push_back(recoMu.isolationR03().sumPt);
    //============= Parameters related to PF isolation =====================
    Mu_pfSumChargedHadronPt.push_back(recoMu.pfIsolationR03().sumChargedHadronPt);
    Mu_pfSumNeutralHadronEt.push_back(recoMu.pfIsolationR03().sumNeutralHadronEt);
    Mu_PFSumPhotonEt.push_back(recoMu.pfIsolationR03().sumPhotonEt);
    Mu_pfSumPUPt.push_back(recoMu.pfIsolationR03().sumPUPt);
  }
}

bool MaketreeMuons::GenParticleRecoMatching(const edm::Event& evt,float eta2,float phi2){
  edm::Handle<GenParticleCollection> genParticles;
  evt.getByLabel(genParticlesColl_, genParticles);
  int NbMuons = 0;
  for(size_t i = 0; i < genParticles->size(); ++ i) {
    const GenParticle & p = (*genParticles)[i];
    //const GenParticle* gen_mom = static_cast<const GenParticle*> (p.mother());
    if( fabs(p.pdgId()) != 13 ) continue; //consider only mu
    if( p.status() != 1 ) continue;        //consider only status 1
    //h1_DaughterID_->Fill(gen_mom->pdgId());
    float DeltaR = delR(p.eta(),p.phi(),eta2,phi2);
    if(fabs(DeltaR)>0.15) continue;
    //printf ("pt(gen) = %f eta(gen) = %f phi(gen) = %f charge(gen) = %d \n",p.pt(),p.eta(),p.phi(),p.charge());
    NbMuons ++;
  }
  if(NbMuons > 0) {
    return true;         
  }
  else return false;
}

//=============================================================
//
//            Method for Genrated Particles Tree
//
//=============================================================
void MaketreeMuons::GenParticleTree(const edm::Event& evt){
  edm::Handle<GenParticleCollection> genParticles;
  evt.getByLabel(genParticlesColl_, genParticles);

   if (!(genParticles.isValid())) return;  

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
    if( fabs(p.pdgId()) != 13 ) continue; //consider only mu
    if( p.status() > 3 ) continue;        //consider only status 1,3
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
void MaketreeMuons::PrimaryVertexTree(const reco::VertexCollection &vertices)
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
//            Method for Jets Tree
//
//=============================================================
void MaketreeMuons::JetsTree(const reco::PFJetCollection &pfjets)
{
  jetcharge.clear();
  jetet.clear();
  jetpt.clear();
  jeteta.clear();
  jetphi.clear();
  jetnumber=0;

  // loose jetID selection

  for ( PFJetCollection::const_iterator i=pfjets.begin(); i!=pfjets.end(); i++) {    

    // loose ID jet selection
    if (fabs(i->eta()) < 2.4  ){
      if ( i->chargedHadronEnergyFraction() > 0. 
	   && i->chargedMultiplicity() > 0. 
	   && i->chargedEmEnergyFraction() < 0.99 ){
	//cout << "Jet passing the loose ID with pT=" << i->pt() << " and eta=" << i->eta() << endl; 	
	jetnumber++;
	jetcharge.push_back(i->charge());
	jetet.push_back(i->et());
	jetpt.push_back(i->pt());
	jeteta.push_back(i->eta());
	jetphi.push_back(i->phi());	
      }
    }
    else if (i->neutralHadronEnergyFraction() < 0.99 && 
	     i->neutralEmEnergyFraction() < 0.99 && 
	     i->getPFConstituents().size() > 1 && 
	     i->muonEnergyFraction()<0.8){
      //cout << "Jet passing the loose ID with pT=" << mIter->pt() << " and eta=" << mIter->eta() << endl; 
      jetnumber++;
      jetcharge.push_back(i->charge());
      jetet.push_back(i->et());
      jetpt.push_back(i->pt());
      jeteta.push_back(i->eta());
      jetphi.push_back(i->phi());
      }

  
    jetcharge.push_back(i->charge());
    jetet.push_back(i->et());
    jetpt.push_back(i->pt());
    jeteta.push_back(i->eta());
    jetphi.push_back(i->phi());
  }
}

//=============================================================
//
//            Method for Trigger Matching Tree
//
//=============================================================
bool MaketreeMuons::IsMuMatchedToHLTMu ( const reco::Muon &mu, std::vector<reco::Particle> HLTMu , std::vector<string> HLTMuNames, double DR, double DPtRel) {
  size_t dim =  HLTMu.size();
  size_t nPass=0;
  if (dim==0) return false;
  //cout << "Entering the function isMuMatchedToHLT" << endl;
  for (size_t k =0; k< dim; k++ ) {
    reco::TrackRef cktTrack = (muon::tevOptimized(mu, 200, 17., 40., 0.25)).first;
    //cout << "HLT matching= " << " Delta R= " << deltaR(HLTMu[k], mu) 
    //cout << "Delta R= " << delR(HLTMu[k].eta(),HLTMu[k].phi(),mu.eta(),mu.phi()) << endl;
    //cout  << " Delta pT= " << fabs(HLTMu[k].pt() - mu.innerTrack()->pt())/ HLTMu[k].pt() 
    /*
    cout  << " Delta pT= " << fabs(HLTMu[k].pt() - cktTrack->pt())/ HLTMu[k].pt() 
	  << " " << fabs(HLTMu[k].pt() - cktTrack->pt())/cktTrack->pt() 
	  << " " << fabs(HLTMu[k].pt() - mu.pt())/mu.pt() << " " << fabs(HLTMu[k].pt() - mu.pt())/HLTMu[k].pt() 
	  << endl;
     */
    //if (  (delR(HLTMu[k].eta(),HLTMu[k].phi(),cktTrack->eta(),cktTrack->phi())  < DR)   && (fabs(HLTMu[k].pt() - cktTrack->pt())/ cktTrack->pt()<DPtRel)){
    //cout << "HLT mu passing filter is= " << HLTMuNames[k].c_str() << " Delta R= " << delR(HLTMu[k].eta(),HLTMu[k].phi(),cktTrack->eta(),cktTrack->phi()) << " Delta pT= " << fabs(HLTMu[k].pt() - cktTrack->pt())/ cktTrack->pt() << endl;
    //if (  (delR(HLTMu[k].eta(),HLTMu[k].phi(),cktTrack->eta(),cktTrack->phi())  < DR)   && (fabs(HLTMu[k].pt() - cktTrack->pt())/ HLTMu[k].pt()<DPtRel)){
    // cout << "HLT mu passing filter is= " << HLTMuNames[k].c_str() << " Delta R= " << delR(HLTMu[k].eta(),HLTMu[k].phi(),cktTrack->eta(),cktTrack->phi()) << " Delta pT= " << fabs(HLTMu[k].pt() - cktTrack->pt())/ cktTrack->pt() << endl;
    //if (  (delR(HLTMu[k].eta(),HLTMu[k].phi(),mu.eta(),mu.phi())  < DR)   && (fabs(HLTMu[k].pt() - mu.pt())/ mu.pt()<DPtRel)){
    if (  (delR(HLTMu[k].eta(),HLTMu[k].phi(),mu.eta(),mu.phi())  < DR) ){

      //if (  (delR(HLTMu[k].eta(),HLTMu[k].phi(),cktTrack->eta(),cktTrack->phi())  < DR) ){ 
      //cout << "HLT mu passing filter is= " << HLTMuNames[k].c_str() << " Delta R= " << delR(HLTMu[k].eta(),HLTMu[k].phi(),mu.eta(),mu.phi()) << " Delta pT= " << fabs(HLTMu[k].pt() - mu.pt())/ mu.pt() << endl;
      nPass++ ;
      break;
    }
  }
  return (nPass>0);
}

float MaketreeMuons::delR(float eta1,float phi1,float eta2,float phi2){
  float mpi=3.14;
  float dp=std::abs(phi1-phi2);
  if (dp>mpi) dp-=float(2*mpi);
  return sqrt((eta1-eta2)*(eta1-eta2) + dp*dp);
}

void MaketreeMuons::TriggerMatchingTree(const edm::Event& iEvent,const reco::MuonCollection* muons,const reco::Vertex &vertex){
  double maxDeltaR_ = 0.2;
  double maxDPtRel_ = 1.0;
  //cout << "Start Trigger matching for muon" << endl;
  MuHLTMatch_nbHLTMuonMatchReco.clear();
  MuHLTMatch_Trigger_pt.clear();
  MuHLTMatch_Trigger_eta.clear();
  MuHLTMatch_Trigger_phi.clear();
  MuHLTMatch_pt.clear();
  MuHLTMatch_eta.clear();
  MuHLTMatch_phi.clear();
  MuHLTMatch_nbMuonMatchHLT.clear();
  std::vector<reco::Particle>  HLTMuMatched;
  std::vector<string> HLTMuMatchedNames;
  edm::Handle<trigger::TriggerEvent> handleTriggerEvent;
  iEvent.getByLabel(triggerEvent, handleTriggerEvent );
  const trigger::TriggerObjectCollection & toc(handleTriggerEvent->getObjects());
  size_t nMuHLT =0;
  size_t nbHLTmuons =0;
  for ( size_t ia = 0; ia < handleTriggerEvent->sizeFilters(); ++ ia) {
    const trigger::Keys & k = handleTriggerEvent->filterKeys(ia);
    for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) {
      std::string fullname = handleTriggerEvent->filterTag(ia).encode();
      // std::cout << " fullname == " << fullname << std::endl;
      std::string name;
      size_t p = fullname.find_first_of(':');
      if (p != std::string::npos) {
	name = fullname.substr(0, p);
      } else {
	name = fullname;
      }

      for (unsigned int l=0;l<triggerFilter.size();l++){
	if (name == triggerFilter.at(l).c_str()) {
	  nbHLTmuons++;
	  MuHLTMatch_nbHLTMuonMatchReco.push_back(nbHLTmuons);
	  HLTMuMatched.push_back(toc[*ki].particle());
	  HLTMuMatchedNames.push_back(name);
	  MuHLTMatch_Trigger_pt.push_back(toc[*ki].pt());
	  MuHLTMatch_Trigger_eta.push_back(toc[*ki].eta());
	  MuHLTMatch_Trigger_phi.push_back(toc[*ki].phi());
	  //cout << "Matching " << triggerFilter.c_str()  << endl;
	}
      }

    }
  }
  // start matching one muon passing high pt muon selection with a muon passing HLT
  for (reco::MuonCollection::const_iterator muonIt = muons->begin(); muonIt!=muons->end(); ++muonIt){
    reco::Muon mu = *muonIt;
    if (mu.innerTrack().isNull()) continue;
    if (mu.globalTrack().isNull()) continue;
    reco::TrackRef cktTrack = (muon::tevOptimized(mu, 200, 17., 40., 0.25)).first;
    if(cktTrack->pt()>53.0 && (mu.isolationR03().sumPt/mu.innerTrack()->pt()<0.10)
       && (cktTrack->ptError()/cktTrack->pt()<0.3) &&
       fabs(cktTrack->dxy(vertex.position()))<0.2 &&
       mu.globalTrack()->hitPattern().trackerLayersWithMeasurement()>5 &&
       mu.globalTrack()->hitPattern().numberOfValidPixelHits()>0 &&
       mu.globalTrack()->hitPattern().numberOfValidMuonHits()>0 &&
       mu.numberOfMatchedStations()>1){ 
      if (IsMuMatchedToHLTMu(mu,HLTMuMatched,HLTMuMatchedNames,maxDeltaR_,maxDPtRel_)==true){
        nMuHLT++;
        //if(nMuHLT > 1) continue;
	cout << "Muon HLT Matched with pT" <<  cktTrack->pt() << endl;	
	MuHLTMatch_nbMuonMatchHLT.push_back(nMuHLT);
	MuHLTMatch_pt.push_back(cktTrack->pt());
	MuHLTMatch_eta.push_back(cktTrack->eta());
	MuHLTMatch_phi.push_back(cktTrack->phi());
      }
    }
  }
}
    
void MaketreeMuons::ComputeMuonMassVtx(const edm::Event& evt,const edm::EventSetup& es,
				       const reco::MuonCollection* muons,const reco::Vertex &vertex)
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
  reco::TrackRef cktTrack1;
  reco::TrackRef cktTrack2;
  reco::TrackRef cktTrack3;
  for (reco::MuonCollection::const_iterator muonIt1 = muons->begin(); muonIt1!=muons->end(); ++muonIt1){
    reco::Muon recoMu1 = *muonIt1;
    if (recoMu1.innerTrack().isNull()) continue;
    if (recoMu1.globalTrack().isNull()) continue;
    if(recoMu1.isTrackerMuon()==false) continue;
    // for 2015B
    //if(fabs(recoMu1.innerTrack()->eta())<1.2) cktTrack1 = (muon::tevOptimized(recoMu1, 200, 17., 40., 0.25)).first;
    //else if(fabs(recoMu1.innerTrack()->eta())>1.2) cktTrack1 = recoMu1.innerTrack();
    //cktTrack1 = recoMu1.innerTrack();
    // for 2015C TuneP and pt 48-> 53
    cktTrack1 = (muon::tevOptimized(recoMu1, 200, 17., 40., 0.25)).first;    
    cout << "PT 1mu= " << cktTrack1->pt() << " and charge =" <<  cktTrack1->charge() << endl;
    if( cktTrack1->pt()>53.0 && (recoMu1.isolationR03().sumPt/recoMu1.innerTrack()->pt()<0.10) 
    	&& (cktTrack1->ptError()/cktTrack1->pt()<0.3) && 
        fabs(cktTrack1->dxy(vertex.position()))<0.2 &&
        recoMu1.globalTrack()->hitPattern().trackerLayersWithMeasurement()>5 && 
        recoMu1.globalTrack()->hitPattern().numberOfValidPixelHits()>0 &&
        recoMu1.globalTrack()->hitPattern().numberOfValidMuonHits()>0 &&
        recoMu1.numberOfMatchedStations()>1 ){
      //cout << "PT 1mu= " << cktTrack1->pt()  << " and charge =" <<  cktTrack1->charge() << endl;
      NbMuons1++;
      if(NbMuons1>1) continue;
      //cout<<"==================================================================="<<endl;
      //printf ("pt1 = %f eta1 = %f phi1 = %f charge1 = %d \n",cktTrack1->pt(),cktTrack1->eta(),cktTrack1->phi(),cktTrack1->charge());
      ttv1.push_back(ttkb1->build(cktTrack1));
      //find the second high pt muon
      for (reco::MuonCollection::const_iterator muonIt2 = muons->begin(); muonIt2!=muons->end(); ++muonIt2){
	reco::Muon recoMu2 = *muonIt2;
	if (recoMu2.innerTrack().isNull()) continue;
	if (recoMu2.globalTrack().isNull()) continue;
	if(recoMu2.isTrackerMuon()==false) continue;
	// for 2015B
	//if(fabs(recoMu2.innerTrack()->eta())<1.2) cktTrack2 = (muon::tevOptimized(recoMu2, 200, 17., 40., 0.25)).first;
	//else if(fabs(recoMu2.innerTrack()->eta())>1.2) cktTrack2 = recoMu2.innerTrack();
	//cktTrack2 = recoMu2.innerTrack();
	// for 2015B TuneP and pt 48 > 53
	cktTrack2 = (muon::tevOptimized(recoMu2, 200, 17., 40., 0.25)).first;
	cout << "PT 2mu= " << cktTrack2->pt()  << " and charge =" <<  cktTrack2->charge()<< endl;
        if(cktTrack2->pt() == cktTrack1->pt()) continue;
	
	if( cktTrack2->pt()>53.0 && (recoMu2.isolationR03().sumPt/recoMu2.innerTrack()->pt()<0.10) 
	    && (cktTrack2->ptError()/cktTrack2->pt()<0.3) && 
	    fabs(cktTrack2->dxy(vertex.position()))<0.2 &&
	    recoMu2.globalTrack()->hitPattern().trackerLayersWithMeasurement()>5 && 
	    recoMu2.globalTrack()->hitPattern().numberOfValidPixelHits()>0 &&
	    recoMu2.globalTrack()->hitPattern().numberOfValidMuonHits()>0 &&
	    recoMu2.numberOfMatchedStations()>1 ){
	  //cout << "PT 2mu= " << cktTrack2->pt()  << " and charge =" <<  cktTrack2->charge() << endl;
	  NbMuons2++;
	  if(NbMuons2>1) continue;
	  ttv1.push_back(ttkb1->build(cktTrack2));
	  //printf ("pt2 = %f eta2 = %f phi2 = %f charge2 = %d \n",cktTrack2->pt(),cktTrack2->eta(),cktTrack2->phi(),cktTrack2->charge());
	  KalmanVertexFitter kvf(true); // false means no smoothing which means no track re-fit
	  CachingVertex<5> vtx1 = kvf.vertex(ttv1);
	  vtxNormChi1 = vtx1.totalChiSquared()/vtx1.degreesOfFreedom();
	  InvariantMassFromVertex imfv1;
	  static const double muon_mass1 = 0.1056583;
	  Measurement1D mass1 = imfv1.invariantMass(vtx1, muon_mass1);

	  DiMass1 = mass1.value();
	  Mu_vtxMass.push_back(DiMass1);
	  Mu_vtxNormChi2.push_back(vtxNormChi1);
	  Mu_vtxMassLept.push_back(cktTrack1->pt());
	  Mu_vtxMassLept.push_back(cktTrack2->pt());
	  //cout<<"size(inside) = "<<ttv1.size()<<endl;
	  //printf ("normChi2 = %f mass = %f\n",vtxNormChi1,DiMass1);
	  cout << "PT 1mu= " << cktTrack1->pt()  << " and charge =" <<  cktTrack1->charge() << endl;
	  cout << "PT 2mu= " << cktTrack2->pt()  << " and charge =" <<  cktTrack2->charge() << endl;
	  cout << "normChi2 = mass = " << vtxNormChi1 << " " << DiMass1 << endl;
	  
	  //find the third high pt muon
	  for (reco::MuonCollection::const_iterator muonIt3 = muons->begin(); muonIt3!=muons->end(); ++muonIt3){
	    reco::Muon recoMu3 = *muonIt3;
	    if (recoMu3.innerTrack().isNull()) continue;
	    if (recoMu3.globalTrack().isNull()) continue;
	    if (recoMu3.isTrackerMuon()==false) continue;
	    //if(fabs(recoMu3.innerTrack()->eta())<1.2) cktTrack3 = (muon::tevOptimized(recoMu3, 200, 17., 40., 0.25)).first;
	    //else if(fabs(recoMu3.innerTrack()->eta())>1.2) cktTrack3 = recoMu3.innerTrack();
	    //cktTrack3 = recoMu3.innerTrack();
	    cktTrack3 = (muon::tevOptimized(recoMu3, 200, 17., 40., 0.25)).first;
	    cout << "PT 3mu= " << cktTrack3->pt()  << " and charge =" <<  cktTrack3->charge() << endl;
	    if(cktTrack3->pt() == cktTrack1->pt()) continue;	   
	    if(cktTrack3->pt() == cktTrack2->pt()) continue;	    

	    if( cktTrack3->pt()>53.0 && (recoMu3.isolationR03().sumPt/recoMu3.innerTrack()->pt()<0.10) 
		&& (cktTrack3->ptError()/cktTrack3->pt()<0.3) && 
		fabs(cktTrack3->dxy(vertex.position()))<0.2 &&
		recoMu3.globalTrack()->hitPattern().trackerLayersWithMeasurement()>5 && 
		recoMu3.globalTrack()->hitPattern().numberOfValidPixelHits()>0 &&
		recoMu3.globalTrack()->hitPattern().numberOfValidMuonHits()>0 &&
		recoMu3.numberOfMatchedStations()>1 ){
	      //cout << "PT 3mu= " << cktTrack3->pt()  << " and charge =" <<  cktTrack3->charge() << endl;
	      NbMuons3++;
	      if(NbMuons3==0) continue;
	      //printf ("pt3 = %f eta3 = %f phi3 = %f charge3 = %d \n",cktTrack3->pt(),cktTrack3->eta(),cktTrack3->phi(),cktTrack3->charge());
	      //cout<<"NbMuons3 = "<<NbMuons3<<endl;
	      ttv2.push_back(ttkb2->build(cktTrack1));
	      ttv2.push_back(ttkb2->build(cktTrack3));

	      ttv3.push_back(ttkb3->build(cktTrack2));
	      ttv3.push_back(ttkb3->build(cktTrack3));
	      KalmanVertexFitter kvf(true); // false means no smoothing which means no track re-fit

	      CachingVertex<5> vtx2 = kvf.vertex(ttv2);
	      vtxNormChi2 = vtx2.totalChiSquared()/vtx2.degreesOfFreedom();
	      InvariantMassFromVertex imfv2;
	      static const double muon_mass2 = 0.1056583;
	      Measurement1D mass2 = imfv2.invariantMass(vtx2, muon_mass2);
	      DiMass2 = mass2.value();
	      Mu_vtxMass.push_back(DiMass2);
	      Mu_vtxNormChi2.push_back(vtxNormChi2);
	      Mu_vtxMassLept.push_back(cktTrack1->pt());
	      Mu_vtxMassLept.push_back(cktTrack3->pt());

	      cout << "PT 1mu= " << cktTrack1->pt()  << " and charge =" <<  cktTrack1->charge() << endl;
	      cout << "PT 3mu= " << cktTrack3->pt()  << " and charge =" <<  cktTrack3->charge() << endl;
	      cout << "normChi2 = %f mass = " << vtxNormChi2 << " " << DiMass2 << endl;
	      CachingVertex<5> vtx3 = kvf.vertex(ttv3);
	      vtxNormChi3 = vtx3.totalChiSquared()/vtx3.degreesOfFreedom();
	      InvariantMassFromVertex imfv3;
	      static const double muon_mass3 = 0.1056583;
	      Measurement1D mass3 = imfv3.invariantMass(vtx3, muon_mass3);
	      DiMass3 = mass3.value();
	      Mu_vtxMass.push_back(DiMass3);
	      Mu_vtxNormChi2.push_back(vtxNormChi3);
	      Mu_vtxMassLept.push_back(cktTrack2->pt());
              Mu_vtxMassLept.push_back(cktTrack3->pt());
	      cout << "PT 2mu= " << cktTrack2->pt()  << " and charge =" <<  cktTrack2->charge() << endl;
	      cout <<  "PT 3mu= " << cktTrack3->pt()  << " and charge =" <<  cktTrack3->charge() << endl;
	      cout << "normChi2 = %f mass = " << vtxNormChi3 << " " << DiMass3 << endl;
	      

	      //cout<<"size(inside) = "<<ttv2.size()<<endl;
	      printf ("normChi2 = %f mass = %f\n",vtxNormChi2,DiMass2);
	      printf ("normChi2 = %f mass = %f\n",vtxNormChi3,DiMass3);
	      //cout<<"==================================================================="<<endl;
	    }
	  }
	}
      }
    }
  }
}
//=============================================================
//
//                Method for PF MET Tree
//
//=============================================================
void MaketreeMuons::pfMETtree(const edm::Event& evt)
{
  edm::Handle<reco::PFMETCollection> pfmetcoll;
  evt.getByLabel(thePFMETCollectionToken_, pfmetcoll);
  if(!pfmetcoll.isValid()) return;
  const PFMETCollection *pfmetcol = pfmetcoll.product();
  const PFMET *pfmet;
  pfmet = &(pfmetcol->front());
  PFMet_pt = pfmet->pt();
  PFMet_eta = pfmet->eta();
  PFMet_phi = pfmet->phi();
  PFMet_en = pfmet->energy();
  PFMet_px = pfmet->px();
  PFMet_py = pfmet->py();
  PFMet_pz = pfmet->pz();
  //scalar sum of transverse energy over all objects
  PFMet_sumEt = pfmet->sumEt();

  edm::Handle<double> metsighandle;
  evt.getByLabel(theMETSignificance_, metsighandle);
  METSign=*metsighandle;
}
 
void MaketreeMuons::fillPU(const edm::Event& iEvent){
  edm::Handle<vector<PileupSummaryInfo> > PupInfo;
  iEvent.getByLabel(PileupSrc_, PupInfo);
  if(!PupInfo.isValid()) return;
  for( vector<PileupSummaryInfo>::const_iterator cand = PupInfo->begin();cand != PupInfo->end(); ++ cand ) {
    num_PU_vertices = cand->getTrueNumInteractions();
    PU_BunchCrossing = cand->getBunchCrossing();
    //std::cout << " Pileup Information: bunchXing, nvtx: " << cand->getBunchCrossing() << " " << cand->getPU_NumInteractions() << std::endl;
    // num_PU_vertices=cand->getPU_NumInteractions();
  }
}

// rho for isolation
void MaketreeMuons::fillRho(const edm::Event& evt){
  Rho.clear();
  edm::Handle<double> rhoHandle;
  evt.getByLabel(rhoIsoInputTag_, rhoHandle);
  if(!rhoHandle.isValid()) return;
  float RhoIsoValue = *rhoHandle;
  //std::cout<<"rhoIso = "<< RhoIsoValue <<std::endl;
  Rho.push_back(RhoIsoValue);
}
//madgraph MC samples reweighing 

void MaketreeMuons::EventsReWeighting(const edm::Event& evt){
  MC_weighting.clear();
  float EventWeight = 1.0;
  edm::Handle<GenEventInfoProduct> gen_ev_info;
  evt.getByLabel(genEventInfo_, gen_ev_info);
  if(!gen_ev_info.isValid()) return;
  EventWeight = gen_ev_info->weight();
  //std::cout<<"mc_weight = "<< gen_ev_info->weight() <<std::endl;
  float mc_weight = ( EventWeight > 0 ) ? 1 : -1;
  //std::cout<<"mc_weight = "<< mc_weight <<std::endl;
  MC_weighting.push_back(mc_weight);
}
//=============================================================
//
//     Method for initializing values for the variables
//
//=============================================================
void MaketreeMuons::IntialValues()
{
  
  NbMuons   = 0;
  value2_   = 0;
  value3_   = 0;
  nbTk      = 0;
  
  
}
