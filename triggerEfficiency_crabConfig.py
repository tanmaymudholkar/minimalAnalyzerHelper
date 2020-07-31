from CRABClient.UserUtilities import config
config = config()

# config.General.requestName = "triggerEfficiency"
# config.General.workArea = "crab_workArea_triggerEfficiency"
config.General.transferLogs = False

config.JobType.pluginName = "Analysis"
config.JobType.psetName = "triggerEfficiency_cfg.py"

config.Data.inputDataset = "/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = "/store/group/lpcsusystealth/triggerEfficiencyNTuples"
config.Site.storageSite = "T3_US_FNALLPC"
