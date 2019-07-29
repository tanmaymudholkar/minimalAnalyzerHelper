#include <ios>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <utility>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <map>

#include "TROOT.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TText.h"

// std::pair<double, double> getRatioFromTH1(TH1F* inputHist, const float& mediumLow, const float& mediumHigh, const float& fakeLow, const float& fakeHigh) {
//   int bin_mediumLow = inputHist->GetXaxis()->FindFixBin(mediumLow);
//   int bin_mediumHigh = inputHist->GetXaxis()->FindFixBin(mediumHigh);
//   int bin_fakeLow = inputHist->GetXaxis()->FindFixBin(fakeLow);
//   int bin_fakeHigh = inputHist->GetXaxis()->FindFixBin(fakeHigh);
//   double numeratorError = 0.;
//   double numerator = inputHist->IntegralAndError(1+bin_fakeLow, bin_fakeHigh, numeratorError);
//   double denominatorError = 0.;
//   double denominator = inputHist->IntegralAndError(bin_mediumLow, bin_mediumHigh, denominatorError);
//   if (denominator <= 0.) {
//     std::cout << "ERROR: unexpected denominator = " << denominator << std::endl;
//     std::exit(EXIT_FAILURE);
//   }
//   double ratio = numerator/denominator;
//   double ratioError = ratio*std::sqrt(std::pow(numeratorError/numerator, 2) + std::pow(denominatorError/denominator, 2));
//   std::pair<double, double> outPair = std::make_pair(ratio, ratioError);
//   return outPair;
// }

enum class IDEfficiencyType{NMinus1_nonTruthMatched=0, NMinus1, global, nIDEfficiencyTypes};

float getIDEfficiency(TH1F* inputHist, const float& cut) {
  int cutBin = inputHist->GetXaxis()->FindFixBin(cut);
  float binLowEdge = inputHist->GetXaxis()->GetBinLowEdge(cutBin);
  float cutBin_fraction_below = (cut-binLowEdge)/(inputHist->GetXaxis()->GetBinWidth(cutBin));
  float numerator = cutBin_fraction_below*inputHist->GetBinContent(cutBin) + inputHist->Integral(0, std::max(0, cutBin-1));
  float denominator = inputHist->Integral(0, 1+inputHist->GetXaxis()->GetNbins());
  return (numerator/denominator);
}

std::map<std::string, std::map<IDEfficiencyType, float> > getIDEfficienciesFromFile(const std::string& inputFileName, const std::map<std::string, float>& criteriaCuts, const std::string& inputType) {
  std::map<std::string, std::map<IDEfficiencyType, float> > IDEfficiencies;
  TFile *inputFile = TFile::Open(inputFileName.c_str(), "READ");
  for (auto&& criterionCutPair: criteriaCuts) {
    auto& criterionName = criterionCutPair.first;
    auto& cutValue = criterionCutPair.second;

    TH1F *h_criterionHist;
    inputFile->GetObject((criterionName+"_NMinus1").c_str(), h_criterionHist);
    IDEfficiencies[criterionName][IDEfficiencyType::NMinus1_nonTruthMatched] = getIDEfficiency(h_criterionHist, cutValue);

    TH1F *h_criterionHist_truthMatched;
    inputFile->GetObject((criterionName+"_NMinus1_TruthMatched").c_str(), h_criterionHist_truthMatched);
    IDEfficiencies[criterionName][IDEfficiencyType::NMinus1] = getIDEfficiency(h_criterionHist_truthMatched, cutValue);

    TH1F *h_criterionHist_global;
    inputFile->GetObject((criterionName+"_global_TruthMatched").c_str(), h_criterionHist_global);
    IDEfficiencies[criterionName][IDEfficiencyType::global] = getIDEfficiency(h_criterionHist_global, cutValue);

    std::string outputFileName = "plots/NMinus1/" + criterionName + "_" + inputType + "_TruthMatched.png";
    TLine *lines = new TLine();
    TCanvas *c = new TCanvas(("c_output_" + outputFileName).c_str(), "c_output");
    gPad->SetLogy();
    gStyle->SetOptStat(110010);
    h_criterionHist_truthMatched->Draw();
    lines->DrawLine(cutValue, h_criterionHist_truthMatched->GetMinimum(), cutValue, h_criterionHist_truthMatched->GetMaximum());
    c->SaveAs(outputFileName.c_str());

    outputFileName = "plots/global/" + criterionName + "_" + inputType + "_TruthMatched.png";
    c = new TCanvas(("c_output_" + outputFileName).c_str(), "c_output");
    gPad->SetLogy();
    gStyle->SetOptStat(110010);
    h_criterionHist_global->Draw();
    lines->DrawLine(cutValue, h_criterionHist_global->GetMinimum(), cutValue, h_criterionHist_global->GetMaximum());
    c->SaveAs(outputFileName.c_str());
  }
  TH1F *h_photonType;
  inputFile->GetObject("photonType", h_photonType);
  if (h_photonType) {
    // std::cout << "Overall efficiency: " << (h_photonType->GetBinContent(h_photonType->FindFixBin(2.0)))/(h_photonType->GetBinContent(h_photonType->FindFixBin(5.0))) << std::endl;
    IDEfficiencies["overall"][IDEfficiencyType::nIDEfficiencyTypes] = (h_photonType->GetBinContent(h_photonType->FindFixBin(2.0)))/(h_photonType->GetBinContent(h_photonType->FindFixBin(5.0)));
  }
  else {
    std::cout << "ERROR: histogram with name \"photonType\" not found!" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  inputFile->Close();
  return IDEfficiencies;
}

std::map<unsigned int, std::map<unsigned int, float> > getStepByStepEfficienciesFromFile(const std::string& inputFileName, const std::vector<std::string>& photonIDCriteria, const std::map<std::string, float>& criteriaCuts, const std::vector<std::vector<std::string> >& stepByStepSequences, const std::string& inputType) {
  std::map<unsigned int, std::map<unsigned int, float> > IDEfficiencies;
  TFile *inputFile = TFile::Open(inputFileName.c_str(), "READ");
  for (unsigned int sequenceIndex = 0; sequenceIndex < stepByStepSequences.size(); ++sequenceIndex) {
    std::string stepByStepPrefix = "sequence_";
    const std::vector<std::string>& sequence = stepByStepSequences.at(sequenceIndex);
    for (unsigned int stepIndex = 0; stepIndex < sequence.size(); ++stepIndex) {
      stepByStepPrefix += (sequence.at(stepIndex) + "_");
    }
    std::string outputFileNamePrefix = "plots/stepByStep/" + stepByStepPrefix;
    for (unsigned int stepIndex = 0; stepIndex < photonIDCriteria.size(); ++stepIndex) {
      unsigned int stepNumber = stepIndex + 1;
      const std::string& criterion = sequence.at(stepIndex);
      std::string histName = (stepByStepPrefix + "step" + std::to_string(stepNumber));
      TH1F* h_criterionHist;
      inputFile->GetObject(histName.c_str(), h_criterionHist);
      assert(std::string(h_criterionHist->GetXaxis()->GetTitle()) == criterion);
      IDEfficiencies[sequenceIndex][stepIndex] = getIDEfficiency(h_criterionHist, criteriaCuts.at(criterion));
      std::string outputFileName = outputFileNamePrefix + "step" + std::to_string(stepNumber) + "_" + inputType;
      TLine *lines = new TLine();
      TCanvas *c = new TCanvas(("c_output_" + outputFileName).c_str(), "c_output");
      outputFileName += ".png";
      gPad->SetLogy();
      gStyle->SetOptStat(110010);
      h_criterionHist->Draw();
      lines->DrawLine(criteriaCuts.at(criterion), h_criterionHist->GetMinimum(), criteriaCuts.at(criterion), h_criterionHist->GetMaximum());
      c->SaveAs(outputFileName.c_str());
    }
  }
  inputFile->Close();
  return IDEfficiencies;
}

std::map<std::string, std::map<std::string, float> > getCorrelationsFromFile(const std::string& inputFileName, const std::vector<std::string>& photonIDCriteria, const std::map<std::string, float>& criteriaCuts, const std::string& inputType) {
  std::map<std::string, std::map<std::string, float> > correlations;
  TFile *inputFile = TFile::Open(inputFileName.c_str(), "READ");
  for (unsigned int criterion1Index = 0; criterion1Index < (-1+photonIDCriteria.size()); ++criterion1Index) {
    const std::string& criterion1 = photonIDCriteria.at(criterion1Index);
    float cut_x = criteriaCuts.at(criterion1);
    for (unsigned int criterion2Index = (1+criterion1Index); criterion2Index < photonIDCriteria.size(); ++criterion2Index) {
      const std::string& criterion2 = photonIDCriteria.at(criterion2Index);
      float cut_y = criteriaCuts.at(criterion2);
      std::string global2DPlotName = criterion1 + "_" + criterion2 + "_global2D_TruthMatched";
      TH2F* h_global2DPlot;
      inputFile->GetObject(global2DPlotName.c_str(), h_global2DPlot);
      if (h_global2DPlot) {
        float N1 = 0.;
        float N2 = 0.;
        float N3 = 0.;
        float N4 = 0.;
        for (int xbincounter = 0; xbincounter <= 1+h_global2DPlot->GetXaxis()->GetNbins(); ++xbincounter) {
          float xBinCenter = h_global2DPlot->GetXaxis()->GetBinCenter(xbincounter);
          for (int ybincounter = 0; ybincounter <= 1+h_global2DPlot->GetYaxis()->GetNbins(); ++ybincounter) {
            float yBinCenter = h_global2DPlot->GetYaxis()->GetBinCenter(ybincounter);
            float binContent = h_global2DPlot->GetBinContent(xbincounter, ybincounter);
            if (xBinCenter < cut_x) {
              if (yBinCenter < cut_y) {
                N1 += binContent;
              }
              else {
                N4 += binContent;
              }
            }
            else {
              if (yBinCenter < cut_y) {
                N2 += binContent;
              }
              else {
                N3 += binContent;
              }
            }
          }
        }
        correlations[criterion1][criterion2] = ((N1*N3) - (N2*N4))/((N1*N3) + (N2*N4));
        TLine *lines = new TLine();
        TText *text = new TText();
        text->SetTextAlign(22);
        std::string outputFileName = "plots/global2D/global2D_" + criterion1 + "_" + criterion2 + "_" + inputType;
        TCanvas *c = new TCanvas(("c_output_" + outputFileName).c_str(), "c_output");
        outputFileName += ".png";
        gPad->SetLogz();
        gStyle->SetOptStat("ou");
        h_global2DPlot->Draw("colz");
        // h_mediumFakeCriteria->GetZaxis()->SetRangeUser(1, 1000);
        lines->DrawLine(criteriaCuts.at(criterion1), h_global2DPlot->GetYaxis()->GetXmin(), criteriaCuts.at(criterion1), h_global2DPlot->GetYaxis()->GetXmax());
        lines->DrawLine(h_global2DPlot->GetXaxis()->GetXmin(), criteriaCuts.at(criterion2), h_global2DPlot->GetXaxis()->GetXmax(), criteriaCuts.at(criterion2));
        text->DrawText(0.5*(h_global2DPlot->GetXaxis()->GetXmin() + criteriaCuts.at(criterion1)), 0.5*(h_global2DPlot->GetYaxis()->GetXmin() + criteriaCuts.at(criterion2)), "N1");
        text->DrawText(0.5*(criteriaCuts.at(criterion1) + h_global2DPlot->GetXaxis()->GetXmax()), 0.5*(h_global2DPlot->GetYaxis()->GetXmin() + criteriaCuts.at(criterion2)), "N2");
        text->DrawText(0.5*(criteriaCuts.at(criterion1) + h_global2DPlot->GetXaxis()->GetXmax()), 0.5*(criteriaCuts.at(criterion2) + h_global2DPlot->GetYaxis()->GetXmax()), "N3");
        text->DrawText(0.5*(h_global2DPlot->GetXaxis()->GetXmin() + criteriaCuts.at(criterion1)), 0.5*(criteriaCuts.at(criterion2) + h_global2DPlot->GetYaxis()->GetXmax()), "N4");
        c->SaveAs(outputFileName.c_str());
      }
      else {
        std::cout << "ERROR: plot with name \"" << global2DPlotName << "\" not found in input file " << inputFileName << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  }
  inputFile->Close();
  return correlations;
}

void getFakeToMediumRatioFromFile(const std::string& inputFileName, const std::string& histogramName, const float& chIso_cutMedium, const float& chIso_cutLoose, const float& sigmaietaieta_cutMedium, const float& sigmaietaieta_cutLoose, const bool& isTruthMatched, const std::string& inputType) {
  std::cout << "From file: " << inputFileName << std::endl;
  float nMedium = 0;
  float nFake = 0;
  TFile *inputFile = TFile::Open(inputFileName.c_str(), "READ");
  TH2F* h_mediumFakeCriteria;
  inputFile->GetObject(histogramName.c_str(), h_mediumFakeCriteria);
  if (h_mediumFakeCriteria) {
    // std::pair<double, double> ratioAndError = getRatioFromTH1(chIso, 1.141, 6.0);
    // std::cout << "Ratio = " << ratioAndError.first << " +/- " << ratioAndError.second << std::endl;
    for (int xbincounter = 0; xbincounter <= 1+h_mediumFakeCriteria->GetXaxis()->GetNbins(); ++xbincounter) {
      float bin_chIso = h_mediumFakeCriteria->GetXaxis()->GetBinCenter(xbincounter);
      for (int ybincounter = 0; ybincounter <= 1+h_mediumFakeCriteria->GetYaxis()->GetNbins(); ++ybincounter) {
        float bin_sigmaietaieta = h_mediumFakeCriteria->GetYaxis()->GetBinCenter(ybincounter);
        if ((bin_chIso < chIso_cutMedium) && (bin_sigmaietaieta < sigmaietaieta_cutMedium)) {
          nMedium += h_mediumFakeCriteria->GetBinContent(xbincounter, ybincounter);
        }
        else if ((bin_chIso < chIso_cutLoose) && (bin_sigmaietaieta < sigmaietaieta_cutLoose)) {
          nFake += h_mediumFakeCriteria->GetBinContent(xbincounter, ybincounter);
        }
      }
    }
  }
  else {
    std::cout << "ERROR: Unable to find histogram named: " << histogramName << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::cout << "Ratio loose/medium: " << nFake/nMedium << std::endl;
  TLine *lines = new TLine();
  TText *text = new TText();
  text->SetTextAlign(22);
  text->SetTextAngle(90);
  std::string outputFileName = "plots/NMinus2/mediumFakeCriteria_";
  if (isTruthMatched) outputFileName += "TruthMatched_";
  outputFileName += inputType;

  TCanvas *c = new TCanvas(("c_output_" + outputFileName).c_str(), "c_output");
  outputFileName += ".png";
  gPad->SetLogz();
  gStyle->SetOptStat(0);
  h_mediumFakeCriteria->Draw("colz");
  h_mediumFakeCriteria->GetZaxis()->SetRangeUser(1, 1000);
  lines->DrawLine(chIso_cutMedium, 0., chIso_cutMedium, 0.025);
  lines->DrawLine(chIso_cutLoose, 0., chIso_cutLoose, 0.025);
  lines->DrawLine(0., sigmaietaieta_cutMedium, 15., sigmaietaieta_cutMedium);
  lines->DrawLine(0., sigmaietaieta_cutLoose, 15., sigmaietaieta_cutLoose);
  text->SetTextColor(kBlue);
  text->DrawText(0.5*(0.+chIso_cutMedium), 0.5*(0.+sigmaietaieta_cutMedium), "medium");
  text->SetTextColor(kRed);
  text->DrawText(0.5*(chIso_cutMedium+chIso_cutLoose), 0.5*(0.+sigmaietaieta_cutMedium), "fake");
  text->DrawText(0.5*(chIso_cutMedium+chIso_cutLoose), 0.5*(sigmaietaieta_cutMedium+sigmaietaieta_cutLoose), "fake");
  text->DrawText(0.5*(0.+chIso_cutMedium), 0.5*(sigmaietaieta_cutMedium+sigmaietaieta_cutLoose), "fake");
  c->SaveAs(outputFileName.c_str());
  inputFile->Close();
}

int main(int argc, char* argv[]) {
  (void)argc;
  (void)argv;
  gROOT->SetBatch();
  std::streamsize original_precision = std::cout.precision();

  const float& chIso_cutMedium = 1.141;
  const float& chIso_cutLoose = 6.0;
  const float& sigmaietaieta_cutMedium = 0.01015;
  const float& sigmaietaieta_cutLoose = 0.02;

  std::cout << "Getting ratio of number of photons in fake range to number of photons in good range..." << std::endl;
  std::string prefix = "output_";
  std::string suffix = ".root";
  std::string inputFileName_hgg = prefix + "hgg" + suffix;
  std::cout << "Without truth matching: " << std::endl;
  getFakeToMediumRatioFromFile(inputFileName_hgg, "mediumFakeCriteria", chIso_cutMedium, chIso_cutLoose, sigmaietaieta_cutMedium, sigmaietaieta_cutLoose, false, "hgg");
  std::cout << "With truth matching: " << std::endl;
  getFakeToMediumRatioFromFile(inputFileName_hgg, "mediumFakeCriteria_TruthMatched", chIso_cutMedium, chIso_cutLoose, sigmaietaieta_cutMedium, sigmaietaieta_cutLoose, true, "hgg");
  std::string inputFileName_stealth = prefix + "stealth" + suffix;
  std::cout << "Without truth matching: " << std::endl;
  getFakeToMediumRatioFromFile(inputFileName_stealth, "mediumFakeCriteria", chIso_cutMedium, chIso_cutLoose, sigmaietaieta_cutMedium, sigmaietaieta_cutLoose, false, "stealth");
  std::cout << "With truth matching: " << std::endl;
  getFakeToMediumRatioFromFile(inputFileName_stealth, "mediumFakeCriteria_TruthMatched", chIso_cutMedium, chIso_cutLoose, sigmaietaieta_cutMedium, sigmaietaieta_cutLoose, true, "stealth");

  std::cout << "Now beginning to calculate ID efficiencies..." << std::endl;
  std::vector<std::string> photonIDCriteria = {"hOverE", "sigmaIEtaIEta", "chIso", "neutIso", "phoIso"};
  std::map<std::string, float> criteriaCuts = {
    {"hOverE", 0.02197},
    {"sigmaIEtaIEta", 0.01015},
    {"chIso", 1.141},
    {"neutIso", 1.0},
    {"phoIso", 1.0}
  };
  std::map<std::string, std::map<IDEfficiencyType, float> > efficiencies_hgg = getIDEfficienciesFromFile(inputFileName_hgg, criteriaCuts, "hgg");
  std::map<std::string, std::map<IDEfficiencyType, float> > efficiencies_stealth = getIDEfficienciesFromFile(inputFileName_stealth, criteriaCuts, "stealth");
  std::cout << "Efficiencies in hgg and stealth:" << std::endl;
  std::cout << "Criterion & (N-1) efficiency (hgg) & global efficiency (hgg) & (N-1) efficiency (stealth) & global efficiency (stealth)\\\\ \\hline \\hline" << std::endl;
  for (unsigned int criterionIndex = 0; criterionIndex < photonIDCriteria.size(); ++criterionIndex) {
    const std::string& criterion = photonIDCriteria.at(criterionIndex);
    std::cout << criterion << " & " << std::setprecision(3) << efficiencies_hgg[criterion][IDEfficiencyType::NMinus1] << " & " << efficiencies_hgg[criterion][IDEfficiencyType::global] << " & " << efficiencies_stealth[criterion][IDEfficiencyType::NMinus1] << " & " << efficiencies_stealth[criterion][IDEfficiencyType::global] << std::setprecision(original_precision) << "\\\\ \\hline" << std::endl;
  }
  std::cout << "Overall & \\multicolumn{2}{c||}{" << std::setprecision(3) << efficiencies_hgg["overall"][IDEfficiencyType::nIDEfficiencyTypes] << std::setprecision(original_precision) << "} & \\multicolumn{2}{c|}{" << std::setprecision(3) << efficiencies_stealth["overall"][IDEfficiencyType::nIDEfficiencyTypes] << std::setprecision(original_precision) << "}\\\\ \\hline" << std::endl;

  std::cout << "Now beginning to calculate step-by-step efficiencies..." << std::endl;
  std::vector<std::vector<std::string> > stepByStepSequences = {
    {"hOverE", "sigmaIEtaIEta", "neutIso", "phoIso", "chIso"},
    {"hOverE", "sigmaIEtaIEta", "neutIso", "chIso", "phoIso"},
    {"chIso", "hOverE", "sigmaIEtaIEta", "neutIso", "phoIso"}
  };
  std::map<unsigned int, std::map<unsigned int, float> > stepByStepEfficiencies_hgg = getStepByStepEfficienciesFromFile(inputFileName_hgg, photonIDCriteria, criteriaCuts, stepByStepSequences, "hgg");
  std::map<unsigned int, std::map<unsigned int, float> > stepByStepEfficiencies_stealth = getStepByStepEfficienciesFromFile(inputFileName_stealth, photonIDCriteria, criteriaCuts, stepByStepSequences, "stealth");
  for (unsigned int sequenceIndex = 0; sequenceIndex < stepByStepSequences.size(); ++sequenceIndex) {
    std::cout << "\\multicolumn{3}{|c|}{Sequence: ";
    const std::vector<std::string>& sequence = stepByStepSequences.at(sequenceIndex);
    for (unsigned int stepIndex = 0; stepIndex < sequence.size(); ++stepIndex) {
      std::cout << (sequence.at(stepIndex) + " $\rightarrow$ ");
    }
    std::cout << "}\\\\ \\hline \\hline" << std::endl;
    std::cout << "Step & efficiency(hgg) & efficiency(stealth) \\\\ \\hline \\hline" << std::endl;
    float overallEfficiency_hgg = 1.;
    float overallEfficiency_stealth = 1.;
    for (unsigned int stepIndex = 0; stepIndex < photonIDCriteria.size(); ++stepIndex) {
      unsigned int stepNumber = stepIndex + 1;
      const std::string& criterion = sequence.at(stepIndex);
      std::cout << stepNumber << " (" << criterion << ") & " << std::setprecision(3) << stepByStepEfficiencies_hgg[sequenceIndex][stepIndex] << " & " << stepByStepEfficiencies_stealth[sequenceIndex][stepIndex] << "\\\\ \\hline" << std::setprecision(original_precision) << std::endl;
      overallEfficiency_hgg *= stepByStepEfficiencies_hgg[sequenceIndex][stepIndex];
      overallEfficiency_stealth *= stepByStepEfficiencies_stealth[sequenceIndex][stepIndex];
    }
    std::cout << "\\hline" << std::endl;
    std::cout << "Overall & " << std::setprecision(3) << overallEfficiency_hgg << " & " << overallEfficiency_stealth << "\\\\ \\hline" << std::setprecision(original_precision) << std::endl;
  }
  
  std::cout << "Now beginning to calculate correlations..." << std::endl;
  // std::map<std::string, int> nHistBins = {
  //   {"hOverE", 500},
  //   {"sigmaIEtaIEta", 500},
  //   {"chIso", 1000},
  //   {"neutIso", 500},
  //   {"phoIso", 500}
  // };
  // std::map<std::string, float> lowerHistLimits = {
  //   {"hOverE", 0.},
  //   {"sigmaIEtaIEta", 0.},
  //   {"chIso", 0.},
  //   {"neutIso", 0.},
  //   {"phoIso", 0.}
  // };
  // std::map<std::string, float> upperHistLimits = {
  //   {"hOverE", 0.1},
  //   {"sigmaIEtaIEta", 0.025},
  //   {"chIso", 60.},
  //   {"neutIso", 2.},
  //   {"phoIso", 2.}
  // };
  std::map<std::string, std::map<std::string, float> > correlations_hgg = getCorrelationsFromFile(inputFileName_hgg, photonIDCriteria, criteriaCuts, "hgg");
  std::map<std::string, std::map<std::string, float> > correlations_stealth = getCorrelationsFromFile(inputFileName_stealth, photonIDCriteria, criteriaCuts, "stealth");
  std::cout << "Correlations in hgg and stealth:" << std::endl;
  for (unsigned int criterion1Index = 0; criterion1Index < (-1+photonIDCriteria.size()); ++criterion1Index) {
    const std::string& criterion1 = photonIDCriteria.at(criterion1Index);
    for (unsigned int criterion2Index = (1+criterion1Index); criterion2Index < photonIDCriteria.size(); ++criterion2Index) {
      const std::string& criterion2 = photonIDCriteria.at(criterion2Index);
      // Formatting LaTeX-style
      std::cout << criterion1 << " & " << criterion2 << " & " << std::setprecision(3) << correlations_hgg[criterion1][criterion2] << " & " << correlations_stealth[criterion1][criterion2] << std::setprecision(original_precision) << "\\\\ \\hline" << std::endl;
    }
  }

  std::cout << "All done!" << std::endl;
  return 0;
}

// TODO:
// Calculate Pearson coefficients as well
// Save correlation 2D plots in proper location
