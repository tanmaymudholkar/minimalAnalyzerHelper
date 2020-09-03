#!/usr/bin/env python

from __future__ import print_function, division

import ROOT
import sys
import tmGeneralUtils
sys.path.append("/uscms/home/tmudholk/private/stealth/STEALTH")
import stealthEnv

ROOT.gROOT.SetBatch(ROOT.kTRUE)

inputFilePathsDict = {
    # "noInvMassCut": "{eP}/store/group/lpcsusystealth/triggerEfficiencyAnalysis/triggerEfficiencies_noInvMassCut.root".format(eP=stealthEnv.EOSPrefix),
    # "invMassCut60": "{eP}/store/group/lpcsusystealth/triggerEfficiencyAnalysis/triggerEfficiencies_invMassCut60.root".format(eP=stealthEnv.EOSPrefix),
    # "invMassCut90": "{eP}/store/group/lpcsusystealth/triggerEfficiencyAnalysis/triggerEfficiencies_invMassCut90.root".format(eP=stealthEnv.EOSPrefix),
    # "noVeto": "{eP}/store/group/lpcsusystealth/triggerEfficiencyAnalysis/triggerEfficiencies_noVeto.root".format(eP=stealthEnv.EOSPrefix),
    "extraTriggers": "{eP}/store/group/lpcsusystealth/triggerEfficiencyAnalysis/triggerEfficiencies_extraTriggers.root".format(eP=stealthEnv.EOSPrefix)
}
print("inputFilePathsDict:")
tmGeneralUtils.prettyPrintDictionary(inputDict=inputFilePathsDict)

nPatterns = 8
targetEfficienciesDict = {"triggerEfficiency_patternIndex_{i}".format(i=patternIndex): "efficiency_pattern{i}".format(i=patternIndex) for patternIndex in range(nPatterns)}
print("targetEfficienciesDict:")
tmGeneralUtils.prettyPrintDictionary(inputDict=targetEfficienciesDict)

outputFolder = "{aR}/triggerEfficiencyPlots_minimal".format(aR=stealthEnv.analysisRoot)

for (inputLabel, inputSourcePath) in inputFilePathsDict.items():
    inputFile = ROOT.TFile.Open(inputSourcePath, "READ")
    if ((inputFile.IsOpen() == ROOT.kFALSE) or (inputFile.IsZombie())): sys.exit("ERROR: unable to open file with name {n}".format(n=inputFilePath))
    for (targetHistogramName, outputID) in targetEfficienciesDict.items():
        extractedEfficiency = ROOT.TEfficiency()
        extractedEfficiency.SetName(outputID)
        print("Fetching efficiency with ID {oID} for input label {iL}".format(oID=outputID, iL=inputLabel))
        inputFile.GetObject(targetHistogramName, extractedEfficiency)
        extractedEfficiency.SetName(outputID)
        if (extractedEfficiency):
            outputFileName = inputLabel + "_" + outputID
            outputCanvas = ROOT.TCanvas("c_" + outputFileName, "c_" + outputFileName)
            extractedEfficiency.Draw()
            ROOT.gPad.Update()
            extractedEfficiency.GetPaintedGraph().SetMinimum(0.7)
            extractedEfficiency.GetPaintedGraph().SetMaximum(1.05)
            ROOT.gPad.Update()
            outputCanvas.SaveAs("{oF}/{oFN}.pdf".format(oF=outputFolder, oFN=outputFileName))
        else:
            sys.exit("ERROR: Unable to find nonempty efficiency with ID {oID} for input label {iL}".format(oID=outputID, iL=inputLabel))
    inputFile.Close()
