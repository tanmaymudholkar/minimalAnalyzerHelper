import sys

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')
options.register(name="inputFilesList",
                 default="none",
                 mult=VarParsing.multiplicity.singleton,
		 mytype=VarParsing.varType.string,
		 info="Newline-separated list of input files to run on.")
options.register(name="outputFilePath",
                 default="none",
                 mult=VarParsing.multiplicity.singleton,
		 mytype=VarParsing.varType.string,
		 info="Path to output file.")
options.parseArguments()

if (options.inputFilesList == "none"):
    sys.exit("ERROR: Argument \"inputFilesList\" needs to be set explicitly. Currently, it is set to \"{iFL}\"".format(iFL=options.inputFilesList))

process = cms.Process("MERGE")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

listOfInputFiles = []
inputFileNamesFileObject = open("{iFL}".format(iFL=options.inputFilesList), 'r')
for inputFileName in inputFileNamesFileObject:
    if (inputFileName[:5] != "file:" ):
        listOfInputFiles.append("root://cms-xrd-global.cern.ch/" + inputFileName.strip())
    else:
        listOfInputFiles.append(inputFileName.strip())
inputFileNamesFileObject.close()
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(*tuple(listOfInputFiles))
)

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("{oFP}".format(oFP=options.outputFilePath))
)

process.p = cms.Path()
process.e = cms.EndPath(process.out)
