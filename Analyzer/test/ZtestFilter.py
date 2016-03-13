
import FWCore.ParameterSet.Config as cms

process = cms.Process('ZFilter')


process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')

process.ztautauMCGenFilter = cms.EDFilter("HZZ4LeptonsMCGenFilter",
    genParticles                     = cms.InputTag("genParticles"),
    DebugHZZ4LeptonsMCGenFilter      = cms.bool(False),
    # 16= Z->tautau                 
    HZZ4LeptonsMCFilterLeptonFlavour = cms.int32(16),
    acceptance                       = cms.double(2.5)
)

process.myPath = cms.Path(process.ztautauMCGenFilter)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('CMSSW763_ZprimeMuMu_13TeV.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('test')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('myPath')
    )                               
 )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
'file:DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0033A97B-8707-E511-9D3B-008CFA1980B8.root'
                             )
                           )


#Endpath
process.o = cms.EndPath ( process.output )
