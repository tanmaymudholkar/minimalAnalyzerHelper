import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

import sys

options = VarParsing ('analysis')
options.register(name="inputPath",
                 default="none",
                 mult=VarParsing.multiplicity.singleton,
		 mytype=VarParsing.varType.string,
		 info="Path to file containing list of MINIAOD sources.")
options.register(name="inputSingleFile",
                 default="none",
                 mult=VarParsing.multiplicity.singleton,
		 mytype=VarParsing.varType.string,
		 info="Single file to use as source.")
options.register(name="eventProgenitor",
                 default="event progenitor",
                 mult=VarParsing.multiplicity.singleton,
		 mytype=VarParsing.varType.string,
		 info="Event progenitor: \"squark\" or \"gluino\".")
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

if (not(options.inputPath == "none") and not(options.inputSingleFile == "none")):
    sys.exit("ERROR: options.inputPath and options.inputSingleFile cannot both be \"none\".")

process = cms.Process("GenLevelDeltaRAnalyzer")
process.load("temp.GenLevelDeltaRAnalyzer.customLogger_cfi")
# process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("temp.GenLevelDeltaRAnalyzer.genLevelDeltaRAnalyzer_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
if ((options.eventProgenitor == "squark") or (options.eventProgenitor == "gluino")):
    process.genLevelDeltaRAnalyzer.eventProgenitor = options.eventProgenitor
else:
    sys.exit("ERROR: options.eventProgenitor can be one of \"squark\" or \"gluino\".")
process.genLevelDeltaRAnalyzer.verbosity = options.verbosity

outputPath = "deltaRNtuple.root"
if not(options.manualOutputPath == "none"): outputPath = options.manualOutputPath
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputPath)
)

listOfInputFiles = []
if not(options.inputPath == "none"):
    with open(options.inputPath, 'r') as inputFileNamesFileObject:
        for inputFileName in inputFileNamesFileObject:
            if (inputFileName[:5] != "file:" ):
                listOfInputFiles.append("root://cms-xrd-global.cern.ch/" + inputFileName.strip())
            else:
                listOfInputFiles.append(inputFileName.strip())
elif not(options.inputSingleFile == "none"):
    listOfInputFiles.append("root://cms-xrd-global.cern.ch/" + options.inputSingleFile)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(*tuple(listOfInputFiles))
)

process.p = cms.Path(process.genLevelDeltaRAnalyzer)
