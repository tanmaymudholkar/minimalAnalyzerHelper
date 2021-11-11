#!/usr/bin/env python

from __future__ import print_function, division

import re
import sys
sys.path.append("/uscms/home/tmudholk/private/stealth/STEALTH")
import stealthEnv

import ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)

EB_EE_TRANSITION_ETA = 1.479

inputFilePath = "{eP}/store/group/lpcsusystealth/deltaRNtuples/histograms/deltaRHistograms_t5.root".format(eP=stealthEnv.EOSPrefix)
outputFolder = "{aR}/deltaRPlots".format(aR=stealthEnv.analysisRoot)
inputFile = ROOT.TFile.Open(inputFilePath, "READ")
if ((inputFile.IsOpen() == ROOT.kFALSE) or (inputFile.IsZombie())): sys.exit("ERROR: unable to open file with name {n}".format(n=inputFilePath))

prefixesToExtract = [
    "eventInfo_photonPairDeltaR",
    "eventInfo_progenitor_eta",
    "eventInfo_neutralino_progenitor_child_eta",
    "eventInfo_neutralino_photon_mother_eta",
    "eventInfo_photon_eta",
    "deltaR_closestGenJet",
    "closestGenJet_fraction_hadronic",
    "deltaR_secondClosestGenJet",
    "secondClosestGenJet_fraction_hadronic"
]

parameterSpaceBinsToExtract = [
    (10, 73),
    (21, 5),
    (21, 7),
    (21, 9),
    (21, 17),
    (21, 25),
    (21, 33),
    (21, 41),
    (21, 49),
    (21, 57),
    (21, 65),
    (21, 73),
    (21, 81),
    (21, 89),
    (21, 97),
    (21, 105),
    (21, 113),
    (21, 121),
    (21, 129),
    (21, 137),
    (21, 145),
    (27, 73)
]

def postDraw(prefixToExtract):
    if ((re.match("eventInfo_progenitor_", prefixToExtract)) or
        (re.match("eventInfo_neutralino_", prefixToExtract)) or
        (re.match("eventInfo_photon_", prefixToExtract))):
        ROOT.gStyle.SetOptStat(1110)
        line_object = ROOT.TLine()
        line_object.SetLineColor(ROOT.kRed)
        line_object.DrawLine( EB_EE_TRANSITION_ETA, ROOT.gPad.GetUymin(),  EB_EE_TRANSITION_ETA, ROOT.gPad.GetUymax())
        line_object.DrawLine(-EB_EE_TRANSITION_ETA, ROOT.gPad.GetUymin(), -EB_EE_TRANSITION_ETA, ROOT.gPad.GetUymax())
        ROOT.gPad.Update()

for prefixToExtract in prefixesToExtract:
    runningCounter = 1
    for parameterSpaceBinToExtract in parameterSpaceBinsToExtract:
        extractedHistogram = ROOT.TH1F()
        histogramNameToExtract = "h_{pTE}_eventProgenitorBin_{ePB}_neutralinoBin_{nB}".format(pTE=prefixToExtract, ePB=parameterSpaceBinToExtract[0], nB=parameterSpaceBinToExtract[1])
        print("Fetching histogram: {h}".format(h=histogramNameToExtract))
        inputFile.GetObject(histogramNameToExtract, extractedHistogram)
        if (extractedHistogram):
            outputCanvas = ROOT.TCanvas("c_" + histogramNameToExtract, "c_" + histogramNameToExtract)
            extractedHistogram.Draw()
            postDraw(prefixToExtract)
            outputCanvas.SaveAs("{oF}/{fname}.pdf".format(oF=outputFolder, fname="{pTE}_{rC:04d}".format(pTE=prefixToExtract, rC=runningCounter)))
            runningCounter += 1
        else:
            sys.exit("ERROR: histogram initialized from path {p} is a nullptr.".format(p=histogramNameToExtract))

inputFile.Close()
