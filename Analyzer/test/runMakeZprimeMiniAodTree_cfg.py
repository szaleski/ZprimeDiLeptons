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


'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/56465C92-C904-E611-82D3-90B11C08C1BA.root',
'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/5E6F359F-8603-E611-8D7E-02163E013C48.root',
'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/7495FE7B-EA03-E611-AB8C-001E675A690A.root',
'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/8634D093-F404-E611-9A09-20CF3027A5EB.root',
'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/86C9B04B-7B05-E611-A464-90B11C066D31.root',
'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/A8EBBB51-8703-E611-B004-02163E013E7B.root',
'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/AC897F36-8C04-E611-8815-02163E00F7B0.root',
'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/C03EB34E-7B05-E611-BB31-02163E00EABE.root'


#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/525F2B4D-9102-E611-9CE9-02163E01764A.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/9644B314-B402-E611-AA09-001517FB1990.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/A05BD472-9503-E611-B36E-02163E0176CE.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/C85FDED8-7B03-E611-9965-001E67DDD0B9.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/D460DAEA-D402-E611-B824-02163E00BC94.root'

#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/50000/0E375935-3F00-E611-B27D-B083FED73AA1.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/50000/3238F5AF-3F00-E611-9EBC-001E67A3FB9B.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/50000/402B923A-3F00-E611-B756-90B11C0BCE26.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/50000/44E6EA62-6EFF-E511-A666-000F530E4778.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/50000/46EB434D-A6FF-E511-9022-5065F381D2C1.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/50000/7C857A43-3F00-E611-8D3A-24BE05C48841.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/50000/928F2FBE-3F00-E611-8ED5-0CC47A0AD6C4.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/50000/DABCBA40-3F00-E611-8B9D-0025901D08D6.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/50000/F25403F2-09FF-E511-8579-002590D9D8A4.root'


#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/1872272D-3405-E611-BE67-002590D0AF86.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/18A3A23D-BD05-E611-BF49-1CC1DE192802.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/6425AA0F-8405-E611-A822-002590D60028.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/C22FCC41-8506-E611-ABA2-02163E017623.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/D8F4F223-3405-E611-90DF-0090FAA1AD04.root'


#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/4E45A151-D6FC-E511-9FE4-0090FAA58134.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/6CEA1236-D6FC-E511-8162-002590FD5122.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/9428B929-D6FC-E511-9E20-A0000420FE80.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/C0C8AE35-D6FC-E511-BE0C-0CC47A0AD63E.root'

#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/02244373-7E03-E611-B581-003048F5B2B4.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/7ED6B748-F202-E611-8730-0CC47A4D7690.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/88C43FF0-6A03-E611-ABEA-002590E2DDC8.root',
#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/B8148C0A-CF03-E611-ACB1-0025905A608C.root' 



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
#process.load('RecoMET.METFilters.metFilters_cff')
#process.Path_BunchSpacingproducer=cms.Path(process.bunchSpacingProducer)
#process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer * process.HBHENoiseFilter)
#process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer * process.HBHENoiseIsoFilter)
## process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)                                                                                                           
 
## process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)                                                                                     
 
#process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
## process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)  
## process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)                                                                                                       
 
#process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
## process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
#process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
## process.Flag_trackingFailureFilter = cms.Path(process.goodVertices + process.trackingFailureFilter)                                                                              
 
#process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
## process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)                                                                                                         
 
## process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)                                                                                                                     
 
## process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)                                                                           
 
## proces..Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)                                                                                                           
 
## and the sub-filters                                                                                                                                                              
 
# process.Flag_trkPOG_manystripclus53X = cms.Path(~manystripclus53X)                                                                                                                
 
# process.Flag_trkPOG_toomanystripclus53X = cms.Path(~toomanystripclus53X)                                                                                                          
 
# process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~logErrorTooManyClusters)
#####################################################################


# Global tag (data)
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '76X_dataRun2_v15', '')

# Global tag (MC)
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_v3', '')

process.demo = cms.EDAnalyzer("MakeZprimeMiniAodTree",
    #outputFile = cms.string('test.root'),
    outputFile = cms.string('CMSSW803_MC_DYtoMuMu2300to3500_13TeV_pattuple.root'), 
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

#process.p = cms.Path(process.Path_BunchSpacingproducer *
#                     process.Flag_HBHENoiseFilter *
#                     process.Flag_HBHENoiseIsoFilter *
#                     process.Flag_CSCTightHalo2015Filter *
#                     process.Flag_EcalDeadCellTriggerPrimitiveFilter *
#                     process.demo)
#process.p = cms.Path(process.HLTEle * process.demo)
#process.p = cms.Path(process.Path_BunchSpacingproducer * process.demo) 


#process.schedule = cms.Schedule( process.Path_BunchSpacingproducer,
#                                 process.Flag_HBHENoiseFilter,
#                                 process.Flag_HBHENoiseIsoFilter,
#                                 process.Flag_CSCTightHalo2015Filter,
#                                 process.Flag_EcalDeadCellTriggerPrimitiveFilter,
#				)

process.p = cms.Path(process.demo)



