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
#process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)  

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'/store/mc/RunIISpring16MiniAODv1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/20000/0017320C-7BFC-E511-9B2D-0CC47A4C8E34.root'
#'/store/mc/RunIISpring16MiniAODv1/TTTo2L2Nu_13TeV-powheg/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/00000/004183D9-4A07-E611-B87E-00266CF9BCC4.root'
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/02244373-7E03-E611-B581-003048F5B2B4.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/7ED6B748-F202-E611-8730-0CC47A4D7690.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/88C43FF0-6A03-E611-ABEA-002590E2DDC8.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/B8148C0A-CF03-E611-ACB1-0025905A608C.root'
'file:pickevents_278957_129_131269573.root'
#    'file:14A4CA95-F437-E611-B6DC-02163E01387F.root'
#'/store/data/Run2016G/SingleMuon/MINIAOD/PromptReco-v1/000/278/957/00000/EEA4821F-BE65-E611-8420-02163E0122D9.root'
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
##################################################################
#####################################################################
# Global tag (data)
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v11', '')

# Global tag (MC)
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_v3', '')

process.demo = cms.EDAnalyzer("MakeZprimeMiniAodTree",
    #outputFile = cms.string('Data.root'),
    outputFile = cms.string('CMSSW803_pickevents.root'), 
    #outputFile = cms.string('TTTo2L2Nu_13TeV-powheg.root'),
    #outputFile = cms.string('CMSSW803_MC_DYtoTauTau_13TeV_pattuple.root'), 
    scProducer = cms.InputTag("reducedEgamma:reducedSuperClusters"),
    vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons      = cms.InputTag("slimmedMuons"),
    electrons  = cms.InputTag("slimmedElectrons"),
    taus       = cms.InputTag("slimmedTaus"),
    photons    = cms.InputTag("slimmedPhotons"),
    #jets       = cms.InputTag("slimmedJetsCMSTopTagCHSPacked:SubJets"),
    jets       = cms.InputTag("slimmedJets"),
    #jets        = cms.InputTag("slimmedJetsPuppi"),
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
    JetSource               = cms.InputTag('slimmedGenJets'),
    PileupSrc               = cms.InputTag("slimmedAddPileupInfo"),
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
    bDiscriminators = cms.vstring(      # list of b-tag discriminators to access
        #'pfTrackCountingHighEffBJetTags',
        #'pfTtrackCountingHighPurBJetTags',
        #'pfJetProbabilityBJetTags',
        #'pfJetBProbabilityBJetTags',
        #'pfSimpleSecondaryVertexHighEffBJetTags',
        #'pfSimpleSecondaryVertexHighPurBJetTags',
        #'pfCombinedSecondaryVertexV2BJetTags',
        'pfCombinedInclusiveSecondaryVertexV2BJetTags'
        #'pfCombinedMVABJetTags'
    ),

    #Analysis = cms.string('ZprimeToEE')
    Analysis = cms.string('ZprimeToMuMu')
)

process.p = cms.Path(process.demo)



