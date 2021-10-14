
import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")
#import os
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
# from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')
#process.GlobalTag.globaltag = cms.string('102X_dataRun2_v12')
#process.GlobalTag.globaltag = cms.string('106X_mc2017_realistic_v8')
process.GlobalTag.globaltag = cms.string('106X_upgrade2018_realistic_v15_L1v1')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)
#process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(3000))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(
#MiniAOD

## RE-REco
"/store/mc/RunIISummer19UL18MiniAOD/BdToJpsiKs_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/10000/00108360-B706-3147-BEF0-36D6B48EA850.root",
#"/store/mc/RunIISummer19UL18MiniAOD/BdToJpsiKs_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/10000/0765772C-471D-F443-A273-21DB12717D14.root",
))        

process.load("slimmedMuonsTriggerMatcher_cfi")
process.load("Psiks0Rootupler_cfi")
#process.rootuple.OnlyGen = cms.bool(True)
process.rootuple.isMC = cms.bool(True)
process.TFileService = cms.Service("TFileService",
#                                fileName = cms.string(FILE2),#'Rootuple_Bdtojpiks0_2018MC_MiniAOD.root'),
                                fileName = cms.string('RootupleMC_BdJPsiKs.root')
)


process.p = cms.Path(process.slimmedMuonsWithTriggerSequence*process.rootuple)
#process.p = cms.Path(process.rootuple)

