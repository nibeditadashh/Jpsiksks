import FWCore.ParameterSet.Config as cms
#import os
process = cms.Process("Rootuple")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
# from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')
process.GlobalTag.globaltag = cms.string('106X_dataRun2_v35')#28')
#process.GlobalTag.globaltag = cms.string('106X_dataRun2_v28')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)
#process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

#FILE1 = os.environ.get('FILE1')
#FILE2 = os.environ.get('FILE2')

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(

#'/store/data/Run2018C/Charmonium/MINIAOD/UL2018_MiniAODv2-v1/00000/02304131-2A08-594C-844E-EDB7070C2C3A.root'  
#'/store/data/Run2018C/Charmonium/MINIAOD/12Nov2019_UL2018-v1/00000/41EF51AA-8AF6-6448-84DE-6402EACC3547.root'
'/store/data/Run2018C/Charmonium/MINIAOD/UL2018_MiniAODv2-v1/00000/02304131-2A08-594C-844E-EDB7070C2C3A.root' 
#'/store/data/Run2016C/Charmonium/MINIAOD/17Jul2018-v1/20000/9C03CBE2-4B8B-E811-9299-0CC47AC17678.root',
#'/store/data/Run2016C/Charmonium/MINIAOD/
   )
)

process.load("slimmedMuonsTriggerMatcher_cfi")
process.load("Psiks0ks0Rootupler_cfi")
process.rootuple.isMC = cms.bool(False)
#process.rootuple.OnlyGen = cms.bool(False)
process.TFileService = cms.Service("TFileService",

      fileName = cms.string('Rootuple_Bdtojpiks0_2018_MiniAOD.root'),
)


#process.p = cms.Path(process.mySequence)
process.p = cms.Path(process.slimmedMuonsWithTriggerSequence*process.rootuple)

