# Meant to be sourced

make

function spawn_makePlots {
    ./makeDeltaRPlots inputFile_deltaRNtuples_${1}.txt /uscms/home/tmudholk/cmslpc_scratch/temp/deltaRHistograms_${1}.root ${EOSPREFIX}/store/group/lpcsusystealth/MCGeneratedMasses/MCGeneratedMasses/MCGeneratedMasses_stealth_${1}Wg_savedObjects.root && xrdcp --verbose --force --path --streams 15 /uscms/home/tmudholk/cmslpc_scratch/temp/deltaRHistograms_${1}.root ${EOSPREFIX}/store/group/lpcsusystealth/deltaRNtuples/histograms/deltaRHistograms_${1}.root && rm -f /uscms/home/tmudholk/cmslpc_scratch/temp/deltaRHistograms_${1}.root
}

set -x

spawn_makePlots t5 > deltaRHistogramsProducer_t5_output.txt 2>&1 &
spawn_makePlots t6 > deltaRHistogramsProducer_t6_output.txt 2>&1 &

set +x
