from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'job_crab_data_Charmonium_finaljob_MINIAOD_CMSSW10218_18B_v1_priyanka'
config.General.workArea = 'crab_data_Charmonium_finaljob_MINIAOD_CMSSW10218_18B_v1_priyanka'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Psiks0Rootupler.py'
config.JobType.outputFiles = ['Rootuple_Bdtojpiks0_2018_MiniAOD.root']

config.Data.inputDataset = '/Charmonium/Run2018B-17Sep2018-v1/MINIAOD'
#config.Data.inputDataset = '/Charmonium/Run2018D-PromptReco-v2/MINIAOD'

config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'
config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 50

config.JobType.allowUndistributedCMSSW = True
# config.Data.ignoreLocality = True
# config.Site.whitelist = ['T2_CH_*', 'T2_UK_*', 'T2_IT_*', 'T2_US_*']

config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON_MuonPhys.txt'

#config.Data.runRange = '315252-316995' #Era A
config.Data.runRange = '317080-319310' #Era B
#config.Data.runRange = '319337-320065' #era C
#config.Data.runRange = '320673-325175' #EraD
                  

config.Data.outLFNDirBase = '/store/user/psadangi/'
config.Data.publication = False


config.Site.storageSite = 'T2_IN_TIFR'

