import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

import sys

options = VarParsing ('analysis')
options.register(name="inputPath",
                 default="none",
                 mult=VarParsing.multiplicity.singleton,
		 mytype=VarParsing.varType.string,
		 info="Path to file containing list of MINIAOD sources.")
options.register(name="manualOutputPath",
                 default="none",
                 mult=VarParsing.multiplicity.singleton,
		 mytype=VarParsing.varType.string,
		 info="Path to output file.")
options.register(name="verbosity",
                 default=0,
		 mult=VarParsing.multiplicity.singleton,
		 mytype=VarParsing.varType.int,
		 info="Verbosity.")
options.parseArguments()

process = cms.Process("StealthTriggerEfficiency")
process.load("temp.StealthTriggerEfficiency.customLogger_cfi")
process.load("temp.StealthTriggerEfficiency.stealthTriggerEfficiency_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
process.stealthTriggerEfficiency.verbosity = options.verbosity

outputPath = "triggerEfficiency.root"
if not(options.manualOutputPath == "none"): outputPath = options.manualOutputPath
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputPath)
)

listOfInputFiles = []
if not(options.inputPath == "none"):
    inputFileNamesFileObject = open(options.inputPath, 'r')
    for inputFileName in inputFileNamesFileObject:
        if (inputFileName[:5] != "file:" ):
            listOfInputFiles.append("root://cms-xrd-global.cern.ch/" + inputFileName.strip())
        else:
            listOfInputFiles.append(inputFileName.strip())
    inputFileNamesFileObject.close()

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(*tuple(listOfInputFiles))
)

process.p = cms.Path(process.stealthTriggerEfficiency)
