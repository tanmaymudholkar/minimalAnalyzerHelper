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

options.register ('eventsToProcess',
                  '',
                  VarParsing.multiplicity.list,
                  VarParsing.varType.string,
                  "Events to process.")
# options.register(name="selectEventsFromFile",
#                  default="none",
#                  mult=VarParsing.multiplicity.singleton,
# 		 mytype=VarParsing.varType.string,
# 		 info="Only analyze events with IDs from this file.")

# options.register(name="eventProgenitor",
#                  default="event progenitor",
#                  mult=VarParsing.multiplicity.singleton,
# 		 mytype=VarParsing.varType.string,
# 		 info="Event progenitor: \"squark\" or \"gluino\".")
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
options.register(name="runParticleTreeDrawer",
                 default=False,
		 mult=VarParsing.multiplicity.singleton,
		 mytype=VarParsing.varType.bool,
		 info="Run particle tree drawer.")
options.register(name="useFullGenCollection",
                 default=False,
		 mult=VarParsing.multiplicity.singleton,
		 mytype=VarParsing.varType.bool,
		 info="GenParticles are obtained from the \"genParticles\" collection rather than the \"prunedGenParticles\" collection. For AOD inputs.")
options.parseArguments()

if (not(options.inputPath == "none") and not(options.inputSingleFile == "none")):
    sys.exit("ERROR: options.inputPath and options.inputSingleFile cannot both be \"none\".")

process = cms.Process("GenLevelDeltaRAnalyzer")
process.load("temp.GenLevelDeltaRAnalyzer.customLogger_cfi")
# process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("temp.GenLevelDeltaRAnalyzer.genLevelDeltaRAnalyzer_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
# if ((options.eventProgenitor == "squark") or (options.eventProgenitor == "gluino")):
#     process.genLevelDeltaRAnalyzer.eventProgenitor = options.eventProgenitor
# else:
#     sys.exit("ERROR: options.eventProgenitor can be one of \"squark\" or \"gluino\".")
process.genLevelDeltaRAnalyzer.verbosity = options.verbosity

if options.useFullGenCollection:
    process.genLevelDeltaRAnalyzer.prunedGenParticlesSrc = cms.InputTag("genParticles")

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
                # listOfInputFiles.append("root://cms-xrd-global.cern.ch/" + inputFileName.strip())
                listOfInputFiles.append("root://cmsxrootd.fnal.gov/" + inputFileName.strip())
            else:
                listOfInputFiles.append(inputFileName.strip())
elif not(options.inputSingleFile == "none"):
    # listOfInputFiles.append("root://cms-xrd-global.cern.ch/" + options.inputSingleFile)
    listOfInputFiles.append("root://cmsxrootd.fnal.gov/" + options.inputSingleFile)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(*tuple(listOfInputFiles))
)

# listOfEventsToProcess = []
# if not(options.selectEventsFromFile == "none"):
#     with open(options.selectEventsFromFile, 'r') as selected_events_file_handle:
#         for selected_events_line_raw in selected_events_file_handle:
#             selected_events_details = (selected_events_line_raw.strip()).split(",")
#             if not(len(selected_events_details) == 2):
#                 sys.exit("ERROR: line {l} in unexpected format.".format(l=selected_events_line_raw.strip()))
#             runID = int(selected_events_details[0])
#             eventID = int(selected_events_details[1])
#             print("Adding to list of events to process: {r}, {e}".format(r=runID, e=eventID))
#             listOfEventsToProcess.append(cms.untracked.EventID(runID, eventID))

# if (len(listOfEventsToProcess) > 0):
#     process.source.eventsToProcess = cms.untracked.VEventRange(*tuple(listOfEventsToProcess))
if options.eventsToProcess:
    process.source.eventsToProcess = cms.untracked.VEventRange(options.eventsToProcess)

if options.runParticleTreeDrawer:
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
    process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
                                       # src = cms.InputTag("prunedGenParticles"),
                                       src = process.genLevelDeltaRAnalyzer.prunedGenParticlesSrc,
                                       printP4 = cms.untracked.bool(True),
                                       printPtEtaPhi = cms.untracked.bool(True),
                                       printVertex = cms.untracked.bool(True),
                                       printStatus = cms.untracked.bool(True),
                                       printIndex = cms.untracked.bool(True)# ,
                                       # status = cms.untracked.vint32(1)
    )
    process.p = cms.Path(process.genLevelDeltaRAnalyzer * process.printTree)
else:
    process.p = cms.Path(process.genLevelDeltaRAnalyzer)
