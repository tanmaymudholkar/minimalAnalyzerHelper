#!/usr/bin/env python

import tmEOSUtils

targetsDict = {
    "inputFile_triggerEfficiencies_noInvMassCut.txt": "/store/group/lpcsusystealth/triggerEfficiencyNTuples/noInvMassCut",
    "inputFile_triggerEfficiencies_invMassCut60.txt": "/store/group/lpcsusystealth/triggerEfficiencyNTuples/invMassCut60",
    "inputFile_triggerEfficiencies_invMassCut90.txt": "/store/group/lpcsusystealth/triggerEfficiencyNTuples/invMassCut90",
    "inputFile_triggerEfficiencies_noVeto.txt": "/store/group/lpcsusystealth/triggerEfficiencyNTuples/noVeto",
}

for (outFilePath, inputDir) in targetsDict.items():
    print("Generating files for outFilePath: {oFP}, inputDir: {iD}".format(oFP=outFilePath, iD=inputDir))
    fileListGenerator = tmEOSUtils.generate_list_of_files_in_eos_path(eos_path=inputDir)
    outputFile = open(outFilePath, 'w')
    for fileName in fileListGenerator:
        outputFile.write(fileName + "\n")
    outputFile.close()
