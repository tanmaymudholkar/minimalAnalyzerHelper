#!/usr/bin/env python

import tmEOSUtils

fileListGenerator = tmEOSUtils.generate_list_of_files_in_eos_path(eos_path="/store/group/lpcsusystealth/triggerEfficiencyNTuples")
outputFile = open("inputFile_triggerEfficiencies.txt", 'w')
for fileName in fileListGenerator:
    outputFile.write(fileName + "\n")
outputFile.close()
