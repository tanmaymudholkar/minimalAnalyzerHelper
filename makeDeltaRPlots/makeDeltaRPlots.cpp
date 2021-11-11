#include <ios>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <map>
#include <vector>
#include <cassert>

#include "/uscms/home/tmudholk/private/stealth/STEALTH/eventSelection/include/MCTemplateReader.h"

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
  std::streamsize original_precision = std::cout.precision();
  // (void)original_precision;

  std::string inputFilePath(argv[1]);
  std::cout << "Opening file containing paths to deltaR ntuples: " << inputFilePath << std::endl;
  std::ifstream inputFile_deltaRNtuples;
  inputFile_deltaRNtuples.open(inputFilePath);
  assert(inputFile_deltaRNtuples.is_open());
  std::string path_deltaRNtuples;

  TChain inputChain_eventInfo("genLevelDeltaRAnalyzer/eventInfoTree");
  TChain inputChain_truePhoton("genLevelDeltaRAnalyzer/deltaRTree");

  while (std::getline(inputFile_deltaRNtuples, path_deltaRNtuples)) {
    inputChain_eventInfo.Add(path_deltaRNtuples.c_str());
    inputChain_truePhoton.Add(path_deltaRNtuples.c_str());
  }
  inputFile_deltaRNtuples.close();
  // Long64_t totalNEntries = inputChain_truePhoton.GetEntries();
  // std::cout << "Available nEntries = " << totalNEntries << std::endl;
  // assert(totalNEntries > 0);

  TFile *outputFile = TFile::Open(argv[2], "RECREATE");
  std::string MCTemplatePath(argv[3]);
  MCTemplateReader templateReader = MCTemplateReader(MCTemplatePath);

  // step1: create histograms
  std::map<int, std::map<int, TH1F> > histograms_eventInfo_photonPairDeltaR;
  std::map<int, std::map<int, TH1F> > histograms_eventInfo_progenitor_eta;
  std::map<int, std::map<int, TH1F> > histograms_eventInfo_neutralino_progenitor_child_eta;
  std::map<int, std::map<int, TH1F> > histograms_eventInfo_neutralino_photon_mother_eta;
  std::map<int, std::map<int, TH1F> > histograms_eventInfo_photon_eta;
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
        std::stringstream massDescriptionStringStream;
	massDescriptionStringStream << std::fixed << std::setprecision(1) << ", m_{#tilde{g}}: " << (templateReader.eventProgenitorMasses).at(eventProgenitorBinIndex) << "GeV, m_{#tilde{#chi}_{1}^{0}}: " << (templateReader.neutralinoMasses).at(neutralinoBinIndex) << "GeV" << std::setprecision(original_precision);
	std::string massDescriptionString = massDescriptionStringStream.str();
	histograms_eventInfo_photonPairDeltaR[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_eventInfo_photonPairDeltaR" + massBinID).c_str(), ("deltaR_photonPair" + massDescriptionString + ";deltaR, photon pair;events").c_str(), 200, -0.15, 3.85);
	histograms_eventInfo_progenitor_eta[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_eventInfo_progenitor_eta" + massBinID).c_str(), ("progenitor eta" + massDescriptionString + ";eta;events").c_str(), 300, -3.0, 3.0);
	histograms_eventInfo_neutralino_progenitor_child_eta[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_eventInfo_neutralino_progenitor_child_eta" + massBinID).c_str(), ("neutralino (progenitor child) eta" + massDescriptionString + ";eta;events").c_str(), 300, -3.0, 3.0);
	histograms_eventInfo_neutralino_photon_mother_eta[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_eventInfo_neutralino_photon_mother_eta" + massBinID).c_str(), ("neutralino (photon parent) eta" + massDescriptionString + ";eta;events").c_str(), 300, -3.0, 3.0);
	histograms_eventInfo_photon_eta[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_eventInfo_photon_eta" + massBinID).c_str(), ("photon eta" + massDescriptionString + ";photon eta;events").c_str(), 300, -3.0, 3.0);
	histograms_deltaR_closestGenJet[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_deltaR_closestGenJet" + massBinID).c_str(), ("deltaR_closestGenJet" + massDescriptionString + ";deltaR, closest GenJet;truth-matched EB #gamma").c_str(), 200, -0.15, 3.85);
	histograms_deltaR_secondClosestGenJet[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_deltaR_secondClosestGenJet" + massBinID).c_str(), ("deltaR_secondClosestGenJet" + massDescriptionString + ";deltaR, second-closest GenJet;truth-matched EB #gamma").c_str(), 200, -0.15, 3.85);
	histograms_photonPT[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_photonPT" + massBinID).c_str(), ("photonPT" + massDescriptionString + ";photonPT;truth-matched EB #gamma").c_str(), 250, 0., 1500.);
	histograms_closestGenJet_PT[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_closestGenJet_PT" + massBinID).c_str(), ("closestGenJet_PT" + massDescriptionString + ";PT, closest GenJet;truth-matched EB #gamma").c_str(), 250, 0., 1500.);
	histograms_closestGenJet_fraction_EM[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_closestGenJet_fraction_EM" + massBinID).c_str(), ("closestGenJet_fraction_EM" + massDescriptionString + ";EM fraction, closest GenJet;truth-matched EB #gamma").c_str(), 103, -0.015, 1.015);
	histograms_closestGenJet_fraction_hadronic[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_closestGenJet_fraction_hadronic" + massBinID).c_str(), ("closestGenJet_fraction_hadronic" + massDescriptionString + ";hadronic fraction, closest GenJet;truth-matched EB #gamma").c_str(), 103, -0.015, 1.015);
	histograms_closestGenJet_fraction_invisible[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_closestGenJet_fraction_invisible" + massBinID).c_str(), ("closestGenJet_fraction_invisible" + massDescriptionString + ";invisible fraction, closest GenJet;truth-matched EB #gamma").c_str(), 103, -0.015, 1.015);
	histograms_closestGenJet_fraction_aux[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_closestGenJet_fraction_aux" + massBinID).c_str(), ("closestGenJet_fraction_aux" + massDescriptionString + ";aux fraction, closest GenJet;truth-matched EB #gamma").c_str(), 103, -0.015, 1.015);
	histograms_closestGenJet_totalFraction[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_closestGenJet_totalFraction" + massBinID).c_str(), ("closestGenJet_totalFraction" + massDescriptionString + ";total fraction, closestGenJet;truth-matched EB #gamma").c_str(), 103, -0.015, 1.015);
	histograms_secondClosestGenJet_PT[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_secondClosestGenJet_PT" + massBinID).c_str(), ("secondClosestGenJet_PT" + massDescriptionString + ";PT, second-closest GenJet;truth-matched EB #gamma").c_str(), 250, 0., 1500.);
	histograms_secondClosestGenJet_fraction_EM[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_secondClosestGenJet_fraction_EM" + massBinID).c_str(), ("secondClosestGenJet_fraction_EM" + massDescriptionString + ";EM fraction, second-closest GenJet;truth-matched EB #gamma").c_str(), 103, -0.015, 1.015);
	histograms_secondClosestGenJet_fraction_hadronic[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_secondClosestGenJet_fraction_hadronic" + massBinID).c_str(), ("secondClosestGenJet_fraction_hadronic" + massDescriptionString + ";hadronic fraction, second-closest GenJet;truth-matched EB #gamma").c_str(), 103, -0.015, 1.015);
	histograms_secondClosestGenJet_fraction_invisible[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_secondClosestGenJet_fraction_invisible" + massBinID).c_str(), ("secondClosestGenJet_fraction_invisible" + massDescriptionString + ";invisible fraction, second-closest GenJet;truth-matched EB #gamma").c_str(), 103, -0.015, 1.015);
	histograms_secondClosestGenJet_fraction_aux[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_secondClosestGenJet_fraction_aux" + massBinID).c_str(), ("secondClosestGenJet_fraction_aux" + massDescriptionString + ";aux fraction, second-closest GenJet;truth-matched EB #gamma").c_str(), 103, -0.015, 1.015);
	histograms_secondClosestGenJet_totalFraction[eventProgenitorBinIndex][neutralinoBinIndex] = TH1F(("h_secondClosestGenJet_totalFraction" + massBinID).c_str(), ("secondClosestGenJet_totalFraction" + massDescriptionString + ";total fraction, second-closest GenJet;truth-matched EB #gamma").c_str(), 103, -0.015, 1.015);
      }
    }
  }

  // step 2: fill histograms
  Long64_t entryIndex = 0;
  TAxis eventProgenitorAxisReference = TAxis(templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass);
  TAxis neutralinoAxisReference = TAxis(templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);

  std::cout << "Beginning to fill event-info histograms..." << std::endl;
  TTreeReader inputChainReader_eventInfo(&inputChain_eventInfo);
  TTreeReaderValue<int> nKinematicStealthPhotons(inputChainReader_eventInfo, "nKinematicStealthPhotons");
  TTreeReaderValue<int> nEnergeticStealthPhotons(inputChainReader_eventInfo, "nEnergeticStealthPhotons");
  TTreeReaderArray<float> evt_eta_progenitor(inputChainReader_eventInfo, "eta_progenitor");
  TTreeReaderArray<float> evt_eta_neutralino_progenitor_child(inputChainReader_eventInfo, "eta_neutralino_progenitor_child");
  TTreeReaderArray<float> evt_eta_neutralino_photon_mother(inputChainReader_eventInfo, "eta_neutralino_photon_mother");
  TTreeReaderValue<float> evt_eta_photon_leading(inputChainReader_eventInfo, "eta_photon_leading");
  TTreeReaderValue<float> evt_eta_photon_subleading(inputChainReader_eventInfo, "eta_photon_subleading");
  TTreeReaderValue<float> evt_deltaR_photonPair(inputChainReader_eventInfo, "deltaR_photonPair");
  TTreeReaderValue<float> eventProgenitorMass(inputChainReader_eventInfo, "eventProgenitorMass");
  TTreeReaderValue<float> neutralinoMass(inputChainReader_eventInfo, "neutralinoMass");
  entryIndex = 0;
  while(inputChainReader_eventInfo.Next()) {
    if (entryIndex%200000 == 0) std::cout << "Control at entryIndex = " << entryIndex << std::endl;
    ++entryIndex;

    int eventProgenitorBinIndex = eventProgenitorAxisReference.FindFixBin(*eventProgenitorMass);
    int neutralinoBinIndex = neutralinoAxisReference.FindFixBin(*neutralinoMass);
    bool is_valid_bin = false;
    try {
      is_valid_bin = (templateReader.isValidBin(eventProgenitorBinIndex, neutralinoBinIndex));
    }
    catch (const std::out_of_range& exception_oor) {
      is_valid_bin = false;
    }
    if (!is_valid_bin) continue;

    if (*nEnergeticStealthPhotons == 2) {
      ((histograms_eventInfo_photon_eta.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(*evt_eta_photon_leading);
      ((histograms_eventInfo_photon_eta.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(*evt_eta_photon_subleading);
      // ((histograms_eventInfo_neutralino_eta.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(*evt_eta_neutralino_leading);
      // ((histograms_eventInfo_neutralino_eta.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(*evt_eta_neutralino_subleading);
      assert(!evt_eta_progenitor.IsEmpty());
      for (const float & eta : evt_eta_progenitor) {
	((histograms_eventInfo_progenitor_eta.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(eta);
      }
      assert(!evt_eta_neutralino_progenitor_child.IsEmpty());
      for (const float & eta : evt_eta_neutralino_progenitor_child) {
	((histograms_eventInfo_neutralino_progenitor_child_eta.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(eta);
      }
      assert(!evt_eta_neutralino_photon_mother.IsEmpty());
      for (const float & eta : evt_eta_neutralino_photon_mother) {
	((histograms_eventInfo_neutralino_photon_mother_eta.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(eta);
      }
    }
    if (*nKinematicStealthPhotons == 2) ((histograms_eventInfo_photonPairDeltaR.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)).Fill(*evt_deltaR_photonPair);
  }

  std::cout << "Beginning to fill photon histograms..." << std::endl;
  TTreeReader inputChainReader_truePhoton(&inputChain_truePhoton);
  TTreeReaderValue<float> evt_eventProgenitorMass(inputChainReader_truePhoton, "evt_eventProgenitorMass");
  TTreeReaderValue<float> evt_neutralinoMass(inputChainReader_truePhoton, "evt_neutralinoMass");
  TTreeReaderValue<float> deltaR_closestGenJet(inputChainReader_truePhoton, "deltaR_closestGenJet");
  TTreeReaderValue<float> deltaR_secondClosestGenJet(inputChainReader_truePhoton, "deltaR_secondClosestGenJet");
  TTreeReaderValue<int>   photonMom_pdgId(inputChainReader_truePhoton, "photonMom_pdgId");
  TTreeReaderValue<float> photonPT(inputChainReader_truePhoton, "photonPT");
  TTreeReaderValue<float> closestGenJet_PT(inputChainReader_truePhoton, "closestGenJet_PT");
  TTreeReaderValue<float> closestGenJet_fraction_EM(inputChainReader_truePhoton, "closestGenJet_fraction_EM");
  TTreeReaderValue<float> closestGenJet_fraction_hadronic(inputChainReader_truePhoton, "closestGenJet_fraction_hadronic");
  TTreeReaderValue<float> closestGenJet_fraction_invisible(inputChainReader_truePhoton, "closestGenJet_fraction_invisible");
  TTreeReaderValue<float> closestGenJet_fraction_aux(inputChainReader_truePhoton, "closestGenJet_fraction_aux");
  TTreeReaderValue<float> closestGenJet_totalFraction(inputChainReader_truePhoton, "closestGenJet_totalFraction");
  TTreeReaderValue<float> secondClosestGenJet_PT(inputChainReader_truePhoton, "secondClosestGenJet_PT");
  TTreeReaderValue<float> secondClosestGenJet_fraction_EM(inputChainReader_truePhoton, "secondClosestGenJet_fraction_EM");
  TTreeReaderValue<float> secondClosestGenJet_fraction_hadronic(inputChainReader_truePhoton, "secondClosestGenJet_fraction_hadronic");
  TTreeReaderValue<float> secondClosestGenJet_fraction_invisible(inputChainReader_truePhoton, "secondClosestGenJet_fraction_invisible");
  TTreeReaderValue<float> secondClosestGenJet_fraction_aux(inputChainReader_truePhoton, "secondClosestGenJet_fraction_aux");
  TTreeReaderValue<float> secondClosestGenJet_totalFraction(inputChainReader_truePhoton, "secondClosestGenJet_totalFraction");
  entryIndex = 0;
  while(inputChainReader_truePhoton.Next()) {
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
	// outputFile->WriteTObject(&((histograms_eventInfo_neutralino_eta.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
	outputFile->WriteTObject(&((histograms_eventInfo_progenitor_eta.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
	outputFile->WriteTObject(&((histograms_eventInfo_neutralino_progenitor_child_eta.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
	outputFile->WriteTObject(&((histograms_eventInfo_neutralino_photon_mother_eta.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
	outputFile->WriteTObject(&((histograms_eventInfo_photon_eta.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
	outputFile->WriteTObject(&((histograms_eventInfo_photonPairDeltaR.at(eventProgenitorBinIndex)).at(neutralinoBinIndex)));
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
