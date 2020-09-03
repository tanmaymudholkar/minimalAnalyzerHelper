# Meant to be sourced

make

set -x
# ./makeTriggerEfficiencyPlots inputFile_triggerEfficiencies_noInvMassCut.txt /uscms/home/tmudholk/cmslpc_scratch/temp/triggerEfficiencies_noInvMassCut.root && xrdcp --verbose --force --path --streams 15 /uscms/home/tmudholk/cmslpc_scratch/temp/triggerEfficiencies_noInvMassCut.root ${EOSPREFIX}/store/group/lpcsusystealth/triggerEfficiencyAnalysis/triggerEfficiencies_noInvMassCut.root && rm -f /uscms/home/tmudholk/cmslpc_scratch/temp/triggerEfficiencies_noInvMassCut.root

# ./makeTriggerEfficiencyPlots inputFile_triggerEfficiencies_invMassCut60.txt /uscms/home/tmudholk/cmslpc_scratch/temp/triggerEfficiencies_invMassCut60.root && xrdcp --verbose --force --path --streams 15 /uscms/home/tmudholk/cmslpc_scratch/temp/triggerEfficiencies_invMassCut60.root ${EOSPREFIX}/store/group/lpcsusystealth/triggerEfficiencyAnalysis/triggerEfficiencies_invMassCut60.root && rm -f /uscms/home/tmudholk/cmslpc_scratch/temp/triggerEfficiencies_invMassCut60.root

# ./makeTriggerEfficiencyPlots inputFile_triggerEfficiencies_invMassCut90.txt /uscms/home/tmudholk/cmslpc_scratch/temp/triggerEfficiencies_invMassCut90.root && xrdcp --verbose --force --path --streams 15 /uscms/home/tmudholk/cmslpc_scratch/temp/triggerEfficiencies_invMassCut90.root ${EOSPREFIX}/store/group/lpcsusystealth/triggerEfficiencyAnalysis/triggerEfficiencies_invMassCut90.root && rm -f /uscms/home/tmudholk/cmslpc_scratch/temp/triggerEfficiencies_invMassCut90.root

# ./makeTriggerEfficiencyPlots inputFile_triggerEfficiencies_noVeto.txt /uscms/home/tmudholk/cmslpc_scratch/temp/triggerEfficiencies_noVeto.root && xrdcp --verbose --force --path --streams 15 /uscms/home/tmudholk/cmslpc_scratch/temp/triggerEfficiencies_noVeto.root ${EOSPREFIX}/store/group/lpcsusystealth/triggerEfficiencyAnalysis/triggerEfficiencies_noVeto.root && rm -f /uscms/home/tmudholk/cmslpc_scratch/temp/triggerEfficiencies_noVeto.root

./makeTriggerEfficiencyPlots inputFile_triggerEfficiencies_extraTriggers.txt /uscms/home/tmudholk/cmslpc_scratch/temp/triggerEfficiencies_extraTriggers.root && xrdcp --verbose --force --path --streams 15 /uscms/home/tmudholk/cmslpc_scratch/temp/triggerEfficiencies_extraTriggers.root ${EOSPREFIX}/store/group/lpcsusystealth/triggerEfficiencyAnalysis/triggerEfficiencies_extraTriggers.root && rm -f /uscms/home/tmudholk/cmslpc_scratch/temp/triggerEfficiencies_extraTriggers.root

set +x
