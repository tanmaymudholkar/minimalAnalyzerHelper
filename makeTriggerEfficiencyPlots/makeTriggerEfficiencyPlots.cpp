#include <ios>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <map>
#include <vector>
#include <cassert>

#include "../../CMSSW_10_2_10/src/temp/StealthTriggerEfficiency/interface/triggers.h"
#include "../../../STEALTH/eventSelection/include/constants.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TEfficiency.h"

int main(int argc, char* argv[]) {
  assert(argc == 3);
  // argv[0]: name of executable
  // argv[1]: inputFilePath, file containing paths to deltaR ntuples
  // argv[2]: outputFilePath, path to output file
  gROOT->SetBatch();

  std::string inputFilePath(argv[1]);
  std::cout << "Opening file containing paths to trigger efficiency ntuples: " << inputFilePath << std::endl;
  std::ifstream inputFile_triggerEfficiencyNtuples;
  inputFile_triggerEfficiencyNtuples.open(inputFilePath);
  assert(inputFile_triggerEfficiencyNtuples.is_open());
  std::string path_triggerEfficiencyNtuples;

  TChain inputChain_eventInfo("stealthTriggerEfficiency/eventInfoTree");

  while (std::getline(inputFile_triggerEfficiencyNtuples, path_triggerEfficiencyNtuples)) {
    inputChain_eventInfo.Add(path_triggerEfficiencyNtuples.c_str());
  }
  inputFile_triggerEfficiencyNtuples.close();
  // Long64_t totalNEntries = inputChain_truePhoton.GetEntries();
  // std::cout << "Available nEntries = " << totalNEntries << std::endl;
  // assert(totalNEntries > 0);

  TFile *outputFile = TFile::Open(argv[2], "RECREATE");

  // step 1: create efficiencies
  std::map<unsigned int, TEfficiency*> triggerEfficiencies;
  for (unsigned int patternIndex = 0; patternIndex < (triggerPatterns::patternsToSave).size(); ++patternIndex) {
    std::string efficiencyName = "triggerEfficiency_patternIndex_" + std::to_string(patternIndex);
    std::string efficiencyTitle = "Efficiency: " + (triggerPatterns::patternsToSave)[patternIndex];
    triggerEfficiencies[patternIndex] = new TEfficiency(efficiencyName.c_str(), (efficiencyTitle + ";leading photon pT;").c_str(), ((HLTEmulation::pTBinEdges).size()-1), &((HLTEmulation::pTBinEdges).at(0)));
    triggerEfficiencies[patternIndex]->SetName(efficiencyName.c_str());
    triggerEfficiencies[patternIndex]->SetTitle((efficiencyTitle + ";leading photon pT;").c_str());
  }

  // step 2: fill efficiencies
  Long64_t entryIndex = 0;
  std::cout << "Beginning to fill event-info histograms..." << std::endl;
  TTreeReader inputChainReader_eventInfo(&inputChain_eventInfo);
  TTreeReaderValue<float> eventRho(inputChainReader_eventInfo, "eventRho");
  TTreeReaderValue<float> pT_leadingPhoton(inputChainReader_eventInfo, "pT_leadingPhoton");
  TTreeReaderValue<float> eta_leadingPhoton(inputChainReader_eventInfo, "eta_leadingPhoton");
  TTreeReaderValue<float> pT_subLeadingPhoton(inputChainReader_eventInfo, "pT_subLeadingPhoton");
  TTreeReaderValue<float> eta_subLeadingPhoton(inputChainReader_eventInfo, "eta_subLeadingPhoton");
  TTreeReaderValue<bool> passesSelection(inputChainReader_eventInfo, "passesSelection");
  std::vector<TTreeReaderValue<bool> > triggerResults;
  for (unsigned int patternIndex = 0; patternIndex < (triggerPatterns::patternsToSave).size(); ++patternIndex) {
    std::string branchName = "passesTrigger_patternIndex_" + std::to_string(patternIndex);
    triggerResults.push_back(TTreeReaderValue<bool>(inputChainReader_eventInfo, branchName.c_str()));
  }

  entryIndex = 0;
  while(inputChainReader_eventInfo.Next()) {
    if (entryIndex%20000 == 0) std::cout << "Control at entryIndex = " << entryIndex << std::endl;
    ++entryIndex;
    if (!(*passesSelection)) continue;

    for (unsigned int patternIndex = 0; patternIndex < (triggerPatterns::patternsToSave).size(); ++patternIndex) {
      triggerEfficiencies.at(patternIndex)->Fill(*(triggerResults.at(patternIndex)), *pT_leadingPhoton);
    }
  }

  // step 3: save efficiencies to output file
  for (unsigned int patternIndex = 0; patternIndex < (triggerPatterns::patternsToSave).size(); ++patternIndex) {
    outputFile->WriteTObject(triggerEfficiencies.at(patternIndex));
  }

  outputFile->Close();
  std::cout << "All done!" << std::endl;
  return 0;
}
