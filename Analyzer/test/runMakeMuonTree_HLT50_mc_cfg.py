import FWCore.ParameterSet.Config as cms
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import *
#from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *

process = cms.Process("tree")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load('Configuration/Geometry/GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
process.load("RecoTracker.Configuration.RecoTracker_cff")
#process.load("RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilderWithoutRefit_cfi")
#process.load("TrackingTools.TrackRefitter.TracksToTrajectories_cff")
process.load('Configuration/StandardSequences/EndOfProcess_cff')


process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True) 
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)  

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
#'/store/data/Run2015B/SingleMuon/AOD/PromptReco-v1/000/251/638/00000/4E720749-F62A-E511-A3CB-02163E014166.root',
#'/store/data/Run2015B/SingleMuon/AOD/PromptReco-v1/000/251/562/00000/52C6E715-A12A-E511-8EC6-02163E012603.root'

#'/store/data/Run2015B/SingleMuon/AOD/PromptReco-v1/000/251/162/00000/C2C5E84D-4227-E511-8878-02163E01280D.root'

#'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/0033A97B-8707-E511-9D3B-008CFA1980B8.root'
#'/store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/00000/0803A0A4-53FD-E411-8198-002618943983.root'
#'file:009B2197-1E0B-E511-BE62-00266CFCCD94.root'
'file:0E580169-0F34-E511-B1B8-B8CA3A709648.root'
 )
)

##-------- Electron events of interest --------
process.HLTEle =cms.EDFilter("HLTHighLevel",
     TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
     # HLTPaths = cms.vstring("HLT_Mu45_eta2p1_v1","HLT_Mu50_v1"),
     HLTPaths = cms.vstring("HLT_Mu50_v*"),
     eventSetupPathsKey = cms.string(''),
     andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
     throw = cms.bool(False) # throw exception on unknown path names
 )

from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

process.load("RecoMET/METProducers.METSignificance_cfi")
process.load("RecoMET/METProducers.METSignificanceParams_cfi")
process.METSignificance.srcLeptons = cms.VInputTag(
       'gedGsfElectrons',
       'muons',
       'gedPhotons'
       )
process.METSignificance.srcPfJets            = cms.InputTag('ak4PFJets')
process.METSignificance.srcMet               = cms.InputTag('pfMet')
process.METSignificance.srcPFCandidates      = cms.InputTag('particleFlow')

# Global tag
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_74_V9', '')
# process.GlobalTag = GlobalTag(process.GlobalTag, '74X_dataRun2_Prompt_v0', '')

process.maketreeMuon = cms.EDAnalyzer("MaketreeMuons",
    outputFile               = cms.string('CMSSW745_Data2015_ZprimeMuMu_13TeV_tree.root'),
    genEventInfo             = cms.InputTag('generator'),
    rhoIsoInputTag           = cms.InputTag("kt6PFJetsForIsolation", "rho"),
    #rhoIsoInputTag          = cms.InputTag("kt6PFJetsCentral:rho"),
    PileupSrc                = cms.InputTag("addPileupInfo"),
    thePFMETCollectionToken  = cms.InputTag("pfMet"),
    METSignificance          = cms.InputTag("METSignificance","METSignificance"),                                  
    globalMuons              = cms.InputTag('muons'),
    globalMuonTracks         = cms.InputTag('globalMuons'),
    vertexCollection         = cms.InputTag('offlinePrimaryVertices'),
    genparticleCollection    = cms.InputTag("genParticles"),                                   
    TrackCollectionTag       = cms.InputTag("generalTracks"),
    Jets                     = cms.InputTag("ak4PFJets"),                                   
    # Trigger matching                                           
    triggerEvent          = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    # triggerFilter         = cms.vstring('hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered45e2p1Q','hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered50Q'),
    # triggerFilter         = cms.vstring('hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered45e2p1Q'),                                  
    #triggerFilter         = cms.vstring('hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered50Q'),
    triggerFilter         = cms.vstring('hltL3fL1sMu16orMu25L1f0L2f16QL3Filtered50Q'),
    #triggerFilter2        = cms.string('hltL3fL1sMu16L1f0L2f16QL3Filtered40Q'),
    maxAbsZ  = cms.double(24),	
    maxd0    = cms.double(2),
    minndof  = cms.int32(4),
    NbGoodPv = cms.int32(1)
)



#process.p = cms.Path(process.HLTEle*process.maketreeMuon)

process.p = cms.Path( process.METSignificance*process.maketreeMuon)

 
