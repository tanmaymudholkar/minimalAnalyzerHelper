#!/usr/bin/env python

import tmEOSUtils

for tDesignation in range(5, 7):
    print("tDesignation: {tD}".format(tD=tDesignation))
    fileListGenerator = tmEOSUtils.generate_list_of_files_in_eos_path(eos_path="/store/group/lpcsusystealth/deltaRNtuples/t{tD}".format(tD=tDesignation))
    outputFile = open("inputFile_deltaRNtuples_t{tD}.txt".format(tD=tDesignation), 'w')
    for fileName in fileListGenerator:
        outputFile.write(fileName + "\n")
    outputFile.close()
