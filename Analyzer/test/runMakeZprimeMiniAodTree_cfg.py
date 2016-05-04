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
    input = cms.untracked.int32(-1),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# reduce verbosity
#process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)  

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'/store/mc/RunIIFall15MiniAODv2/ZprimeToMuMu_M-5000_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/1A6ED956-2AB8-E511-9CC8-002590A3C95E.root',
'/store/mc/RunIIFall15MiniAODv2/ZprimeToMuMu_M-5000_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/48E4F234-2AB8-E511-A241-002590A83308.root',
'/store/mc/RunIIFall15MiniAODv2/ZprimeToMuMu_M-5000_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/A47DD552-2AB8-E511-8D40-0025905C96A4.root',
'/store/mc/RunIIFall15MiniAODv2/ZprimeToMuMu_M-5000_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/EAB39656-2AB8-E511-B1D2-002590A3C95E.root' 

#'/store/mc/RunIIFall15MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_6000_Inf/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/70000/2E4092CB-F2B8-E511-8312-34E6D7BEAF01.root',
#'/store/mc/RunIIFall15MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_6000_Inf/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/70000/508B4B0E-F3B8-E511-86AB-441EA1733A74.root',
#'/store/mc/RunIIFall15MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_6000_Inf/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/70000/B07C7308-F3B8-E511-AC37-D8D385AE8848.root'

#'/store/mc/RunIIFall15MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/1A60999E-B1B9-E511-BCF9-44A842CFD5FF.root',
#'/store/mc/RunIIFall15MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/42C56857-FCB8-E511-882F-001EC9AF0999.root',
#'/store/mc/RunIIFall15MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/A66FB321-F2B8-E511-90D6-02163E00F408.root',
#'/store/mc/RunIIFall15MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/C04A3D0B-02B9-E511-A240-1CC1DE192736.root' 

#'/store/data/Run2015D/SingleMuon/MINIAOD/16Dec2015-v1/10000/00006301-CAA8-E511-AD39-549F35AD8BC9.root'
#'/store/mc/RunIIFall15MiniAODv2/ZprimeToMuMu_M-5000_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/1A6ED956-2AB8-E511-9CC8-002590A3C95E.root'
    )
)
##-------- Electron events of interest --------
process.HLTEle =cms.EDFilter("HLTHighLevel",
     TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
     HLTPaths = cms.vstring("HLT_Mu50_*"),
     eventSetupPathsKey = cms.string(''),
     andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
     throw = cms.bool(False) # throw exception on unknown path names
 )

process.noscraping = cms.EDFilter("FilterOutScraping",
                applyfilter = cms.untracked.bool(True),
                debugOn = cms.untracked.bool(False),
                numtrack = cms.untracked.uint32(10),
                thresh = cms.untracked.double(0.25)
                )

# Global tag (data)
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '76X_dataRun2_v15', '')

# Global tag (MC)
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '76X_mcRun2_asymptotic_v12', '')

process.demo = cms.EDAnalyzer("MakeZprimeMiniAodTree",
    #outputFile = cms.string('CMSSW763_MC_DYtoEE_6000_inf_13TeV_pattuple.root'),
    outputFile = cms.string('CMSSW763_MC_DYtoMuMu_13TeV_pattuple200.root'), 
    scProducer = cms.InputTag("reducedEgamma:reducedSuperClusters"),
    vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons      = cms.InputTag("slimmedMuons"),
    electrons  = cms.InputTag("slimmedElectrons"),
    taus       = cms.InputTag("slimmedTaus"),
    photons    = cms.InputTag("slimmedPhotons"),
    jets       = cms.InputTag("slimmedJets"),
    fatjets    = cms.InputTag("slimmedJetsAK8"),
    mets       = cms.InputTag("slimmedMETs"),
    #mets       = cms.InputTag("slimmedMETsNoHF"),
    packed     = cms.InputTag("packedGenParticles"),
    pruned     = cms.InputTag("prunedGenParticles"),
    pfCands    = cms.InputTag("packedPFCandidates"),
    rhoIsoInputTag          = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp"),
    EBrecHitCollectionLabel = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    EErecHitCollectionLabel = cms.InputTag("reducedEgamma","reducedEERecHits"),
    ecalRechitEB            = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    ecalRechitEE            = cms.InputTag("reducedEgamma","reducedEERecHits"),
    #METSignificance        = cms.InputTag("METSignificance","METSignificance"),
    #generalTracksLabel      = cms.InputTag("generalTracks"),
    bits           = cms.InputTag("TriggerResults","","HLT"),
    prescales      = cms.InputTag("patTrigger"),
    objects        = cms.InputTag("selectedPatTrigger"),
    GenBosonID     = cms.int32(1000000),
    ParticleID     = cms.int32(13),
    ParticleStatus = cms.int32(25),
    maxAbsZ  = cms.double(24),
    maxd0    = cms.double(2),
    minndof  = cms.int32(4),
    NbGoodPv = cms.int32(1),
    #Analysis = cms.string('ZprimeToEE')
    Analysis = cms.string('ZprimeToMuMu')
)

process.p = cms.Path(process.demo)
#process.p = cms.Path(process.HLTEle * process.demo)
 
