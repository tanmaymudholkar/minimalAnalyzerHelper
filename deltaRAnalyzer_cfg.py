import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')
# options.register(name="inputType",
#                  default="none",
#                  mult=VarParsing.multiplicity.singleton,
# 		 mytype=VarParsing.varType.string,
# 		 info="Input to run on. Currently supported: \"hgg\" or \"stealth\" or \"stealth2017_t5Wg\" or \"stealth2017_t6Wg\" or \"stealth_privateMC_fastsim\" or \"stealth_privateMC_fullsim\".")
options.register(name="inputPath",
                 default="none",
                 mult=VarParsing.multiplicity.singleton,
		 mytype=VarParsing.varType.string,
		 info="Path to file containing list of MINIAOD sources.")
options.register(name="outputPath",
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

process = cms.Process("GenLevelDeltaRAnalyzer")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("temp.GenLevelDeltaRAnalyzer.genLevelDeltaRAnalyzer_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
process.genLevelDeltaRAnalyzer.outputPath = options.outputPath

listOfInputFiles = []
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

process.p = cms.Path(process.genLevelDeltaRAnalyzer)
