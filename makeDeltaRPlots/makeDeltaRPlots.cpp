#include <ios>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <map>
#include <vector>
#include <cassert>

#include "../../../STEALTH/eventSelection/include/MCTemplateReader.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TAxis.h"
#include "TH1F.h"

int main(int argc, char* argv[]) {
  assert(argc == 4);
  // argv[0]: name of executable
  // argv[1]: inputFilePath, file containing paths to deltaR ntuples
  // argv[2]: outputFilePath, path to output file
  // argv[3]: MCTemplatePath, path to MC template file
  gROOT->SetBatch();
  // std::streamsize original_precision = std::cout.precision();
  // (void)original_precision;

  std::string inputFilePath(argv[1]);
  std::cout << "Opening file containing paths to deltaR ntuples: " << inputFilePath << std::endl;
  std::ifstream inputFile_deltaRNtuples;
  inputFile_deltaRNtuples.open(inputFilePath);
  assert(inputFile_deltaRNtuples.is_open());
  std::string path_deltaRNtuples;
  TChain inputChain("genLevelDeltaRAnalyzer/deltaRTree");
  while (std::getline(inputFile_deltaRNtuples, path_deltaRNtuples)) {
    inputChain.Add(path_deltaRNtuples.c_str());
  }
  inputFile_deltaRNtuples.close();
  // Long64_t totalNEntries = inputChain.GetEntries();
  // std::cout << "Available nEntries = " << totalNEntries << std::endl;
  // assert(totalNEntries > 0);

  TFile *outputFile = TFile::Open(argv[2], "RECREATE");
  std::string MCTemplatePath(argv[3]);
  MCTemplateReader templateReader = MCTemplateReader(MCTemplatePath);

  // step1: create histograms
  std::map<int, std::map<int, TH1F> > histograms_deltaR_closestGenJet;
  std::map<int, std::map<int, TH1F> > histograms_deltaR_secondClosestGenJet;
  std::map<int, std::map<int, TH1F> > histograms_photonPT;
  std::map<int, std::map<int, TH1F> > histograms_closestGenJet_PT;
  std::map<int, std::map<int, TH1F> > histograms_closestGenJet_fraction_EM;
  std::map<int, std::map<int, TH1F> > histograms_closestGenJet_fraction_hadronic;
  std::map<int, std::map<int, TH1F> > histograms_closestGenJet_fraction_invisible;
  std::map<int, std::map<int, TH1F> > histograms_closestGenJet_fraction_aux;
  std::map<int, std::map<int, TH1F> > histograms_closestGenJet_totalFraction;
  std::map<int, std::map<int, TH1F> > histograms_secondClosestGenJet_PT;
  std::map<int, std::map<int, TH1F> > histograms_secondClosestGenJet_fraction_EM;
  std::map<int, std::map<int, TH1F> > histograms_secondClosestGenJet_fraction_hadronic;
  std::map<int, std::map<int, TH1F> > histograms_secondClosestGenJet_fraction_invisible;
  std::map<int, std::map<int, TH1F> > histograms_secondClosestGenJet_fraction_aux;
  std::map<int, std::map<int, TH1F> > histograms_secondClosestGenJet_totalFraction;
  for (int eventProgenitorBinIndex = 1; eventProgenitorBinIndex <= templateReader.nEventProgenitorMassBins; ++eventProgenitorBinIndex) {
    for (int neutralinoBinIndex = 1; neutralinoBinIndex <= templateReader.nNeutralinoMassBins; ++neutralinoBinIndex) {
      if (templateReader.isValidBin(eventProgenitorBinIndex, neutralinoBinIndex)) {
	std::string massBinID = "_eventProgenitorBin_" + std::to_string(eventProgenitorBinIndex) + "_neutralinoBin_" + std::to_string(neutralinoBinIndex);
	histograms_deltaR_closestGenJet[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_deltaR_closestGenJet" + massBinID).c_str(), ("h_deltaR_closestGenJet" + massBinID + ";deltaR, closest GenJet;truth-matched EB #gamma").c_str(), 200, -0.15, 3.85);
	histograms_deltaR_secondClosestGenJet[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_deltaR_secondClosestGenJet" + massBinID).c_str(), ("h_deltaR_secondClosestGenJet" + massBinID + ";deltaR, second-closest GenJet;truth-matched EB #gamma").c_str(), 200, -0.15, 3.85);
	histograms_photonPT[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_photonPT" + massBinID).c_str(), ("h_photonPT" + massBinID + ";photonPT;truth-matched EB #gamma").c_str(), 250, 0., 1500.);
	histograms_closestGenJet_PT[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_closestGenJet_PT" + massBinID).c_str(), ("h_closestGenJet_PT" + massBinID + ";PT, closest GenJet;truth-matched EB #gamma").c_str(), 250, 0., 1500.);
	histograms_closestGenJet_fraction_EM[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_closestGenJet_fraction_EM" + massBinID).c_str(), ("h_closestGenJet_fraction_EM" + massBinID + ";EM fraction, closest GenJet;truth-matched EB #gamma").c_str(), 103, -0.015, 1.015);
	histograms_closestGenJet_fraction_hadronic[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_closestGenJet_fraction_hadronic" + massBinID).c_str(), ("h_closestGenJet_fraction_hadronic" + massBinID + ";hadronic fraction, closest GenJet;truth-matched EB #gamma").c_str(), 103, -0.015, 1.015);
	histograms_closestGenJet_fraction_invisible[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_closestGenJet_fraction_invisible" + massBinID).c_str(), ("h_closestGenJet_fraction_invisible" + massBinID + ";invisible fraction, closest GenJet;truth-matched EB #gamma").c_str(), 103, -0.015, 1.015);
	histograms_closestGenJet_fraction_aux[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_closestGenJet_fraction_aux" + massBinID).c_str(), ("h_closestGenJet_fraction_aux" + massBinID + ";aux fraction, closest GenJet;truth-matched EB #gamma").c_str(), 103, -0.015, 1.015);
	histograms_closestGenJet_totalFraction[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_closestGenJet_totalFraction" + massBinID).c_str(), ("h_closestGenJet_totalFraction" + massBinID + ";total fraction, closestGenJet;truth-matched EB #gamma").c_str(), 103, -0.015, 1.015);
	histograms_secondClosestGenJet_PT[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_secondClosestGenJet_PT" + massBinID).c_str(), ("h_secondClosestGenJet_PT" + massBinID + ";PT, second-closest GenJet;truth-matched EB #gamma").c_str(), 250, 0., 1500.);
	histograms_secondClosestGenJet_fraction_EM[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_secondClosestGenJet_fraction_EM" + massBinID).c_str(), ("h_secondClosestGenJet_fraction_EM" + massBinID + ";EM fraction, second-closest GenJet;truth-matched EB #gamma").c_str(), 103, -0.015, 1.015);
	histograms_secondClosestGenJet_fraction_hadronic[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_secondClosestGenJet_fraction_hadronic" + massBinID).c_str(), ("h_secondClosestGenJet_fraction_hadronic" + massBinID + ";hadronic fraction, second-closest GenJet;truth-matched EB #gamma").c_str(), 103, -0.015, 1.015);
	histograms_secondClosestGenJet_fraction_invisible[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_secondClosestGenJet_fraction_invisible" + massBinID).c_str(), ("h_secondClosestGenJet_fraction_invisible" + massBinID + ";invisible fraction, second-closest GenJet;truth-matched EB #gamma").c_str(), 103, -0.015, 1.015);
	histograms_secondClosestGenJet_fraction_aux[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_secondClosestGenJet_fraction_aux" + massBinID).c_str(), ("h_secondClosestGenJet_fraction_aux" + massBinID + ";aux fraction, second-closest GenJet;truth-matched EB #gamma").c_str(), 103, -0.015, 1.015);
	histograms_secondClosestGenJet_totalFraction[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_secondClosestGenJet_totalFraction" + massBinID).c_str(), ("h_secondClosestGenJet_totalFraction" + massBinID + ";total fraction, second-closest GenJet;truth-matched EB #gamma").c_str(), 103, -0.015, 1.015);
      }
    }
  }

  // step 2: fill histograms
  TTreeReader inputChainReader(&inputChain);
  TTreeReaderValue<float> evt_eventProgenitorMass(inputChainReader, "evt_eventProgenitorMass");
  TTreeReaderValue<float> evt_neutralinoMass(inputChainReader, "evt_neutralinoMass");
  TTreeReaderValue<float> deltaR_closestGenJet(inputChainReader, "deltaR_closestGenJet");
  TTreeReaderValue<float> deltaR_secondClosestGenJet(inputChainReader, "deltaR_secondClosestGenJet");
  TTreeReaderValue<int>   photonMom_pdgId(inputChainReader, "photonMom_pdgId");
  TTreeReaderValue<float> photonPT(inputChainReader, "photonPT");
  TTreeReaderValue<float> closestGenJet_PT(inputChainReader, "closestGenJet_PT");
  TTreeReaderValue<float> closestGenJet_fraction_EM(inputChainReader, "closestGenJet_fraction_EM");
  TTreeReaderValue<float> closestGenJet_fraction_hadronic(inputChainReader, "closestGenJet_fraction_hadronic");
  TTreeReaderValue<float> closestGenJet_fraction_invisible(inputChainReader, "closestGenJet_fraction_invisible");
  TTreeReaderValue<float> closestGenJet_fraction_aux(inputChainReader, "closestGenJet_fraction_aux");
  TTreeReaderValue<float> closestGenJet_totalFraction(inputChainReader, "closestGenJet_totalFraction");
  TTreeReaderValue<float> secondClosestGenJet_PT(inputChainReader, "secondClosestGenJet_PT");
  TTreeReaderValue<float> secondClosestGenJet_fraction_EM(inputChainReader, "secondClosestGenJet_fraction_EM");
  TTreeReaderValue<float> secondClosestGenJet_fraction_hadronic(inputChainReader, "secondClosestGenJet_fraction_hadronic");
  TTreeReaderValue<float> secondClosestGenJet_fraction_invisible(inputChainReader, "secondClosestGenJet_fraction_invisible");
  TTreeReaderValue<float> secondClosestGenJet_fraction_aux(inputChainReader, "secondClosestGenJet_fraction_aux");
  TTreeReaderValue<float> secondClosestGenJet_totalFraction(inputChainReader, "secondClosestGenJet_totalFraction");
  Long64_t entryIndex = 0;
  TAxis eventProgenitorAxisReference = TAxis(templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass);
  TAxis neutralinoAxisReference = TAxis(templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);
  while(inputChainReader.Next()) {
    if (entryIndex%200000 == 0) std::cout << "Control at entryIndex = " << entryIndex << std::endl;
    ++entryIndex;

    int eventProgenitorBinIndex = eventProgenitorAxisReference.FindFixBin(*evt_eventProgenitorMass);
    int neutralinoBinIndex = neutralinoAxisReference.FindFixBin(*evt_neutralinoMass);
    if (!(templateReader.isValidBin(eventProgenitorBinIndex, neutralinoBinIndex))) continue;

    ((histograms_deltaR_closestGenJet.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(*deltaR_closestGenJet);
    ((histograms_deltaR_secondClosestGenJet.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(*deltaR_secondClosestGenJet);
    ((histograms_photonPT.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(*photonPT);
    ((histograms_closestGenJet_PT.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(*closestGenJet_PT);
    ((histograms_closestGenJet_fraction_EM.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(*closestGenJet_fraction_EM);
    ((histograms_closestGenJet_fraction_hadronic.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(*closestGenJet_fraction_hadronic);
    ((histograms_closestGenJet_fraction_invisible.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(*closestGenJet_fraction_invisible);
    ((histograms_closestGenJet_fraction_aux.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(*closestGenJet_fraction_aux);
    ((histograms_closestGenJet_totalFraction.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(*closestGenJet_totalFraction);
    ((histograms_secondClosestGenJet_PT.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(*secondClosestGenJet_PT);
    ((histograms_secondClosestGenJet_fraction_EM.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(*secondClosestGenJet_fraction_EM);
    ((histograms_secondClosestGenJet_fraction_hadronic.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(*secondClosestGenJet_fraction_hadronic);
    ((histograms_secondClosestGenJet_fraction_invisible.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(*secondClosestGenJet_fraction_invisible);
    ((histograms_secondClosestGenJet_fraction_aux.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(*secondClosestGenJet_fraction_aux);
    ((histograms_secondClosestGenJet_totalFraction.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(*secondClosestGenJet_totalFraction);
  }

  // step 3: save histograms to output file
  for (int eventProgenitorBinIndex = 1; eventProgenitorBinIndex <= templateReader.nEventProgenitorMassBins; ++eventProgenitorBinIndex) {
    for (int neutralinoBinIndex = 1; neutralinoBinIndex <= templateReader.nNeutralinoMassBins; ++neutralinoBinIndex) {
      if (templateReader.isValidBin(eventProgenitorBinIndex, neutralinoBinIndex)) {
	outputFile->WriteTObject(&((histograms_deltaR_closestGenJet.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
	outputFile->WriteTObject(&((histograms_deltaR_secondClosestGenJet.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
	outputFile->WriteTObject(&((histograms_photonPT.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
	outputFile->WriteTObject(&((histograms_closestGenJet_PT.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
	outputFile->WriteTObject(&((histograms_closestGenJet_fraction_EM.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
	outputFile->WriteTObject(&((histograms_closestGenJet_fraction_hadronic.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
	outputFile->WriteTObject(&((histograms_closestGenJet_fraction_invisible.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
	outputFile->WriteTObject(&((histograms_closestGenJet_fraction_aux.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
	outputFile->WriteTObject(&((histograms_closestGenJet_totalFraction.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
	outputFile->WriteTObject(&((histograms_secondClosestGenJet_PT.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
	outputFile->WriteTObject(&((histograms_secondClosestGenJet_fraction_EM.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
	outputFile->WriteTObject(&((histograms_secondClosestGenJet_fraction_hadronic.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
	outputFile->WriteTObject(&((histograms_secondClosestGenJet_fraction_invisible.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
	outputFile->WriteTObject(&((histograms_secondClosestGenJet_fraction_aux.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
	outputFile->WriteTObject(&((histograms_secondClosestGenJet_totalFraction.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
      }
    }
  }

  outputFile->Close();
  std::cout << "All done!" << std::endl;
  return 0;
}
