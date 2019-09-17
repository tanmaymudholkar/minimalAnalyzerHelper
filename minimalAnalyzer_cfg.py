import sys

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')
options.register(name="inputType",
                 default="none",
                 mult=VarParsing.multiplicity.singleton,
		 mytype=VarParsing.varType.string,
		 info="Input to run on. Currently supported: \"hgg\" or \"stealth\" or \"stealth2017_t5Wg\" or \"stealth2017_t6Wg\" or \"stealth_privateMC_fastsim\" or \"stealth_privateMC_fullsim\".")
options.register(name="verbosity",
                 default=0,
		 mult=VarParsing.multiplicity.singleton,
		 mytype=VarParsing.varType.int,
		 info="Verbosity.")
options.parseArguments()

if (options.inputType == "none"):
    sys.exit("inputType=" + options.inputType + " should be set")
# if not((options.inputType == "stealth") or (options.inputType == "hgg") or (options.inputType == "stealth2017_t5Wg") or (options.inputType == "stealth2017_t6Wg")):
#     sys.exit("inputType=" + options.inputType + " should be one of \"hgg\" or \"stealth\" or \"stealth2017_t5Wg\" or \"stealth2017_t6Wg\"") # Disabling for now, the list is getting too long :-)

process = cms.Process("MinimalAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("temp.MinimalMiniAODAnalyzer.minimalAnalyzer_cfi")
process.minimalAnalyzer.outputPath = ("output_{t}.root").format(t=options.inputType)
process.minimalAnalyzer.verbosity = options.verbosity
if (options.inputType == "hgg"):
    process.minimalAnalyzer.filterType = options.inputType
    process.minimalAnalyzer.selectJetsNearPhotons = True
else:
    process.minimalAnalyzer.filterType = "stealth"

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

listOfInputFiles = []
inputFileNamesFileObject = open(("inputFileList_{t}.txt").format(t=options.inputType), 'r')
for inputFileName in inputFileNamesFileObject:
    if (inputFileName[:5] != "file:" ):
        listOfInputFiles.append("root://cms-xrd-global.cern.ch/" + inputFileName.strip())
    else:
        listOfInputFiles.append(inputFileName.strip())
inputFileNamesFileObject.close()
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(*tuple(listOfInputFiles))
)

process.p = cms.Path(process.minimalAnalyzer)
