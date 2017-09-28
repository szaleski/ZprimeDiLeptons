import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import *
#from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *
from PhysicsTools.SelectorUtils.tools.vid_id_tools import switchOnVIDElectronIdProducer
from PhysicsTools.SelectorUtils.tools.vid_id_tools import setupVIDElectronSelection
from PhysicsTools.SelectorUtils.tools.vid_id_tools import setupAllVIDIdsInModule
from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat

process = cms.Process("tree")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load('Configuration/Geometry/GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
process.load("RecoTracker.Configuration.RecoTracker_cff")
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#Track isolation correction
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("PhysicsTools.PatAlgos.slimming.primaryVertexAssociation_cfi")
process.load("RecoEgamma.ElectronIdentification.heepIdVarValueMapProducer_cfi")


#EGamma VID
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)




process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True) 
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# reduce verbosity
#process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(200000)  

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

'/store/mc/RunIISpring16MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_200_400/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/20000/0C242E56-BA3A-E611-9B2D-0242AC130004.root',



#'/store/data/Run2016B/DoubleEG/MINIAOD/03Feb2017_ver2-v2/50000/00054938-CEEA-E611-889E-0CC47A4D7650.root',
#'/store/data/Run2016B/DoubleEG/MINIAOD/03Feb2017_ver2-v2/50000/001BC16F-90EA-E611-87E2-F04DA27540BB.root',
#'/store/data/Run2016B/DoubleEG/MINIAOD/03Feb2017_ver2-v2/50000/001EB4EF-D3EA-E611-B94E-0CC47A4C8F26.root'

#'/store/data/Run2016G/SingleMuon/MINIAOD/23Sep2016-v1/90000/9695BF0B-8E97-E611-A1A4-00259044051C.root',
#'/store/data/Run2016G/SingleMuon/MINIAOD/23Sep2016-v1/90000/1C603D88-C597-E611-A0BE-008CFAF2931E.root'
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
#from HEEP.VID.tools import addHEEPV70ElesMiniAOD
#addHEEPV70ElesMiniAOD(process,useStdName=True)

from HEEP.VID.tools import addHEEPV70ElesMiniAOD
addHEEPV70ElesMiniAOD(process,useStdName=True)

#####################################################################
# Global tag (data)
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7', '')

# Global tag (MC)
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_v3', '')

process.demo = cms.EDAnalyzer("MakeZprimeMiniAodTreeHEEPData",
    outputFile = cms.string('Data.root'),
    eleTrkPtIsoLabel = cms.InputTag("heepIDVarValueMaps","eleTrkPtIso"),
    scProducer = cms.InputTag("reducedEgamma:reducedSuperClusters"),
    vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons      = cms.InputTag("slimmedMuons"),
    electrons  = cms.InputTag("slimmedElectrons"),
    eles       = cms.InputTag("slimmedElectrons"),
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
    ParticleID1    = cms.int32(13),
    ParticleID2    = cms.int32(11),
    ParticleID3    = cms.int32(15),
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

process.p = cms.Path( process.heepIDVarValueMaps * process.demo )



