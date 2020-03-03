from CRABClient.UserUtilities import config
config = config()

config.General.requestName = "deltaRNtuplizer_stealth_t5"
config.General.workArea = "crab_workArea_deltaRNtuplizer_stealth_t5"
config.General.transferLogs = False

config.JobType.pluginName = "Analysis"
config.JobType.psetName = "deltaRAnalyzer_cfg.py"
config.JobType.pyCfgParams = ["eventProgenitor=gluino"]

config.Data.inputDataset = "/SMS-T5WgStealth_TuneCP2_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PUFall17Fast_94X_mc2017_realistic_v15-v1/MINIAODSIM"
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = "/store/group/lpcsusystealth/deltaRNtuples/t5"
config.Site.storageSite = "T3_US_FNALLPC"
