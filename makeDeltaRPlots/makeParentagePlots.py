#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse
sys.path.append("/uscms/home/tmudholk/private/stealth/STEALTH")
import stealthEnv

inputArgumentsParser = argparse.ArgumentParser(description='Make MC parentage plots from output of GenLevelDeltaRAnalyzer module.')
inputArgumentsParser.add_argument('--inputPath', required=True, action='append', help='Path to sample.', type=str)
inputArgumentsParser.add_argument('--outputFolder', default="/uscms/home/tmudholk/nobackup/analysisAreas/MCParentage", help='Output directory in which to store plots.',type=str)
inputArgumentsParser.add_argument('--outputSuffix', required=True, help='Suffix to append to all plots.',type=str)
inputArguments = inputArgumentsParser.parse_args()

import ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

os.system("mkdir -p {o}".format(o=inputArguments.outputFolder))

parentage_dictionary = {
    0: "other",
    1: "LHE",
    2: "status 4 proton",
    3: "radiated off jet"
}

output_histogram_nmatches = ROOT.TH1F("nmatches", "Number of gen-level final state photons matched to two selected photons", 3, -0.5, 2.5)
output_histograms_parentage = {
    "leading": ROOT.TH1F("parentage_leading", "Parentage of leading reco photon", 4, -0.5, 3.5),
    "subLeading": ROOT.TH1F("parentage_subLeading", "Parentage of subleading reco photon", 4, -0.5, 3.5)
}

for output_label in output_histograms_parentage:
    for parentage_key in parentage_dictionary:
        output_histograms_parentage[output_label].GetXaxis().SetBinLabel(output_histograms_parentage[output_label].GetXaxis().FindFixBin(parentage_key), parentage_dictionary[parentage_key])

inputChain = ROOT.TChain("genLevelDeltaRAnalyzer/eventInfoTree")
inputChain.SetMaxTreeSize(100000000000) # 1 TB
for inputPath in inputArguments.inputPath:
    inputChain.Add(inputPath)
nEntries = inputChain.GetEntries()

for eventIndex in range(0, nEntries):
    treeStatus = inputChain.LoadTree(eventIndex)
    if (treeStatus < 0):
        break
    evtStatus = inputChain.GetEntry(eventIndex)
    if (evtStatus <= 0):
        continue
    # event_weight = (inputChain.b_MCXSecWeight)*(inputChain.genWeight)*(inputChain.b_evtPrefiringWeight)*(inputChain.b_evtphotonMCScaleFactor)*(inputChain.b_PUWeightNoSelection)
    # output_histogram_nmatches.Fill(inputChain.nMatchedFinalStatePhotons, event_weight)
    output_histogram_nmatches.Fill(inputChain.nMatchedFinalStatePhotons)
    for output_label in output_histograms_parentage:
        # output_histograms_parentage[output_label].Fill(getattr(inputChain, "parentage_{l}Photon".format(l=output_label)), event_weight)
        output_histograms_parentage[output_label].Fill(getattr(inputChain, "parentage_{l}Photon".format(l=output_label)))

for output_label in output_histograms_parentage:
    output_canvas = ROOT.TCanvas("output_{l}".format(l=output_label), "output_{l}".format(l=output_label), 1200, 1024)
    output_histograms_parentage[output_label].Draw()
    output_canvas.SaveAs("{o}/parentage_{l}Photon_{s}.pdf".format(o=inputArguments.outputFolder, l=output_label, s=inputArguments.outputSuffix))

output_canvas = ROOT.TCanvas("output_n_truth_matched", "output_n_truth_matched", 1024, 768)
output_histogram_nmatches.Draw()
output_canvas.SaveAs("{o}/n_truth_matched_photons_{s}.pdf".format(o=inputArguments.outputFolder, l=output_label, s=inputArguments.outputSuffix))

print("All Done!")
# ./makeDeltaRPlots/makeParentagePlots.py --inputPath deltaR_diphoton17.root --outputSuffix "diphoton17"
# ./makeDeltaRPlots/makeParentagePlots.py --inputPath deltaR_GJet17.root --outputSuffix "GJet17"
# ./makeDeltaRPlots/makeParentagePlots.py --inputPath deltaR_QCD17.root --outputSuffix "QCD17"
