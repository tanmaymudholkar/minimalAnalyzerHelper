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
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TText.h"
#include "TLegendEntry.h"
#include "TLegend.h"

enum class IDEfficiencyType{NMinus1_nonTruthMatched=0, NMinus1, global, nIDEfficiencyTypes};
enum class correlationType{customized=0, pearson, nCorrelationTypes};

float getIDEfficiency(TH1F* inputHist, const float& cut, const bool& invertCut) {
  int cutBin = inputHist->GetXaxis()->FindFixBin(cut);
  float binLowEdge = inputHist->GetXaxis()->GetBinLowEdge(cutBin);
  float cutBin_fraction_below = (cut-binLowEdge)/(inputHist->GetXaxis()->GetBinWidth(cutBin));
  float cutBin_fraction_above = 1.0-cutBin_fraction_below;
  float belowCut = cutBin_fraction_below*inputHist->GetBinContent(cutBin) + inputHist->Integral(0, std::max(0, cutBin-1));
  float aboveCut = cutBin_fraction_above*inputHist->GetBinContent(cutBin) + inputHist->Integral(std::min(1+cutBin, 1+inputHist->GetXaxis()->GetNbins()), 1+inputHist->GetXaxis()->GetNbins());
  float all = inputHist->Integral(0, 1+inputHist->GetXaxis()->GetNbins());
  if (invertCut) return (aboveCut/all);
  return (belowCut/all);
}

std::map<std::string, std::map<IDEfficiencyType, float> > getIDEfficienciesFromFile(const std::string& inputFileName, const std::map<std::string, float>& criteriaCuts, const std::string& inputType) {
  std::map<std::string, std::map<IDEfficiencyType, float> > IDEfficiencies;
  TFile *inputFile = TFile::Open(inputFileName.c_str(), "READ");
  for (auto&& criterionCutPair: criteriaCuts) {
    auto& criterionName = criterionCutPair.first;
    auto& cutValue = criterionCutPair.second;

    TH1F *h_criterionHist;
    inputFile->GetObject((criterionName+"_NMinus1").c_str(), h_criterionHist);
    h_criterionHist->StatOverflows(kTRUE);
    IDEfficiencies[criterionName][IDEfficiencyType::NMinus1_nonTruthMatched] = getIDEfficiency(h_criterionHist, cutValue, false);

    TH1F *h_criterionHist_truthMatched;
    inputFile->GetObject((criterionName+"_NMinus1_TruthMatched").c_str(), h_criterionHist_truthMatched);
    h_criterionHist_truthMatched->StatOverflows(kTRUE);
    IDEfficiencies[criterionName][IDEfficiencyType::NMinus1] = getIDEfficiency(h_criterionHist_truthMatched, cutValue, false);

    TH1F *h_criterionHist_global;
    inputFile->GetObject((criterionName+"_global_TruthMatched").c_str(), h_criterionHist_global);
    h_criterionHist_global->StatOverflows(kTRUE);
    IDEfficiencies[criterionName][IDEfficiencyType::global] = getIDEfficiency(h_criterionHist_global, cutValue, false);

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
  h_photonType->StatOverflows(kTRUE);
  if (h_photonType) {
    IDEfficiencies["overall"][IDEfficiencyType::nIDEfficiencyTypes] = (h_photonType->GetBinContent(h_photonType->FindFixBin(2.0)))/(h_photonType->GetBinContent(h_photonType->FindFixBin(5.0)));
  }
  else {
    std::cout << "ERROR: histogram with name \"photonType\" not found!" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  TH1F* h_phoET_TruthMatched;
  inputFile->GetObject("phoET_TruthMatched", h_phoET_TruthMatched);
  TH1F* h_phoET_passingID_TruthMatched;
  inputFile->GetObject("phoET_passingID_TruthMatched", h_phoET_passingID_TruthMatched);
  if (h_phoET_TruthMatched) {
    if (h_phoET_passingID_TruthMatched) {
      std::string outputFileName = "plots/efficiency/IDEfficiencies_overall_check_" + inputType;
      TCanvas *c = new TCanvas(("c_output_" + outputFileName).c_str(), "c_output");
      outputFileName += ".png";
      c->Divide(1, 2);
      c->cd(1);
      gPad->SetLogy();
      gStyle->SetOptStat(0);
      h_phoET_TruthMatched->SetLineColor(kBlue);
      h_phoET_TruthMatched->Draw();
      h_phoET_passingID_TruthMatched->SetLineColor(kRed);
      h_phoET_passingID_TruthMatched->Draw("same");
      c->cd(2);
      TEfficiency* net_ID_efficiency = new TEfficiency(*h_phoET_passingID_TruthMatched, *h_phoET_TruthMatched);
      gPad->SetLogy(0);
      net_ID_efficiency->Draw();
      gPad->Update();
      net_ID_efficiency->GetPaintedGraph()->SetMinimum(-0.2);
      net_ID_efficiency->GetPaintedGraph()->SetMaximum(1.2);
      gPad->Update();
      c->SaveAs(outputFileName.c_str());
    }
    else {
      std::cout << "ERROR: histogram with name \"phoET_passingID_TruthMatched\" not found!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  else {
    std::cout << "ERROR: histogram with name \"phoET_TruthMatched\" not found!" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  inputFile->Close();
  return IDEfficiencies;
}

std::map<std::string, std::map<IDEfficiencyType, float> > getHLTEfficienciesFromFile(const std::string& inputFileName, const std::map<std::string, float>& criteriaCuts, const std::vector<std::string>& criteriaToInvert, const std::string& inputType) {
  std::map<std::string, std::map<IDEfficiencyType, float> > IDEfficiencies;
  TFile *inputFile = TFile::Open(inputFileName.c_str(), "READ");
  for (auto&& criterionCutPair: criteriaCuts) {
    auto& criterionName = criterionCutPair.first;
    auto& cutValue = criterionCutPair.second;
    bool toInvert;
    if (std::find(criteriaToInvert.begin(), criteriaToInvert.end(), criterionName) == criteriaToInvert.end()) toInvert = false;
    else toInvert = true;

    TH1F *h_criterionHist;
    inputFile->GetObject((criterionName+"_NMinus1_HLT_TruthMatched_").c_str(), h_criterionHist);
    h_criterionHist->StatOverflows(kTRUE);
    IDEfficiencies[criterionName][IDEfficiencyType::NMinus1] = getIDEfficiency(h_criterionHist, cutValue, toInvert);

    TH1F *h_criterionHist_global;
    inputFile->GetObject((criterionName+"_global_HLT_TruthMatched_").c_str(), h_criterionHist_global);
    h_criterionHist_global->StatOverflows(kTRUE);
    IDEfficiencies[criterionName][IDEfficiencyType::global] = getIDEfficiency(h_criterionHist_global, cutValue, toInvert);

    std::string outputFileName = "plots/HLT/NMinus1/" + criterionName + "_" + inputType + "_TruthMatched.png";
    TLine *lines = new TLine();
    TCanvas *c = new TCanvas(("c_output_" + outputFileName).c_str(), "c_output");
    gPad->SetLogy();
    gStyle->SetOptStat(110010);
    h_criterionHist->Draw();
    lines->DrawLine(cutValue, h_criterionHist->GetMinimum(), cutValue, h_criterionHist->GetMaximum());
    c->SaveAs(outputFileName.c_str());

    outputFileName = "plots/HLT/global/" + criterionName + "_" + inputType + "_TruthMatched.png";
    c = new TCanvas(("c_output_" + outputFileName).c_str(), "c_output");
    gPad->SetLogy();
    gStyle->SetOptStat(110010);
    h_criterionHist_global->Draw();
    lines->DrawLine(cutValue, h_criterionHist_global->GetMinimum(), cutValue, h_criterionHist_global->GetMaximum());
    c->SaveAs(outputFileName.c_str());
  }

  std::string outputFileName = "plots/HLT/overallEfficiency";
  TCanvas *c = new TCanvas(("c_output_" + outputFileName).c_str(), "c_output");
  outputFileName += ".png";
  TEfficiency* globalEfficiency;
  inputFile->GetObject("HLTEfficiency_TruthMatched", globalEfficiency);
  globalEfficiency->Draw("colz");
  gPad->Update();
  // globalEfficiency->GetPaintedGraph()->SetMinimum(-0.2);
  // globalEfficiency->GetPaintedGraph()->SetMaximum(1.2);
  gPad->Update();
  c->SaveAs(outputFileName.c_str());
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
      h_criterionHist->StatOverflows(kTRUE);
      assert(std::string(h_criterionHist->GetXaxis()->GetTitle()) == criterion);
      IDEfficiencies[sequenceIndex][stepIndex] = getIDEfficiency(h_criterionHist, criteriaCuts.at(criterion), false);
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

std::map<correlationType, std::map<std::string, std::map<std::string, float> > > getCorrelationsFromFile(const std::string& inputFileName, const std::vector<std::string>& photonIDCriteria, const std::map<std::string, float>& criteriaCuts, const bool& isNMinus2, const std::string& inputType) {
  std::map<correlationType, std::map<std::string, std::map<std::string, float> > > correlations;
  TFile *inputFile = TFile::Open(inputFileName.c_str(), "READ");
  for (unsigned int criterion1Index = 0; criterion1Index < (-1+photonIDCriteria.size()); ++criterion1Index) {
    const std::string& criterion1 = photonIDCriteria.at(criterion1Index);
    float cut_x = criteriaCuts.at(criterion1);
    for (unsigned int criterion2Index = (1+criterion1Index); criterion2Index < photonIDCriteria.size(); ++criterion2Index) {
      const std::string& criterion2 = photonIDCriteria.at(criterion2Index);
      float cut_y = criteriaCuts.at(criterion2);
      std::string plot2DName = criterion1 + "_" + criterion2;
      if (isNMinus2) plot2DName += "_NMinus2";
      else plot2DName += "_global2D";
      plot2DName += "_TruthMatched";
      TH2F* h_global2DPlot;
      inputFile->GetObject(plot2DName.c_str(), h_global2DPlot);
      h_global2DPlot->StatOverflows(kTRUE);
      if (h_global2DPlot) {
        correlations[correlationType::pearson][criterion1][criterion2] = h_global2DPlot->GetCorrelationFactor();
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
        float correlationValue = ((N1*N3) - (N2*N4))/((N1*N3) + (N2*N4));
        correlations[correlationType::customized][criterion1][criterion2] = correlationValue;
        if (!(isNMinus2)) {
          TLine *lines = new TLine();
          TText *text = new TText();
          text->SetTextAlign(22);
          std::string outputFileName = "plots/global2D/global2D_" + criterion1 + "_" + criterion2 + "_" + inputType;
          TCanvas *c = new TCanvas(("c_output_" + outputFileName).c_str(), "c_output");
          outputFileName += ".png";
          gPad->SetLogz();
          gStyle->SetOptStat("ou");
          h_global2DPlot->Draw("colz");
          lines->DrawLine(criteriaCuts.at(criterion1), h_global2DPlot->GetYaxis()->GetXmin(), criteriaCuts.at(criterion1), h_global2DPlot->GetYaxis()->GetXmax());
          lines->DrawLine(h_global2DPlot->GetXaxis()->GetXmin(), criteriaCuts.at(criterion2), h_global2DPlot->GetXaxis()->GetXmax(), criteriaCuts.at(criterion2));
          text->DrawText(0.5*(h_global2DPlot->GetXaxis()->GetXmin() + criteriaCuts.at(criterion1)), 0.5*(h_global2DPlot->GetYaxis()->GetXmin() + criteriaCuts.at(criterion2)), "N1");
          text->DrawText(0.5*(criteriaCuts.at(criterion1) + h_global2DPlot->GetXaxis()->GetXmax()), 0.5*(h_global2DPlot->GetYaxis()->GetXmin() + criteriaCuts.at(criterion2)), "N2");
          text->DrawText(0.5*(criteriaCuts.at(criterion1) + h_global2DPlot->GetXaxis()->GetXmax()), 0.5*(criteriaCuts.at(criterion2) + h_global2DPlot->GetYaxis()->GetXmax()), "N3");
          text->DrawText(0.5*(h_global2DPlot->GetXaxis()->GetXmin() + criteriaCuts.at(criterion1)), 0.5*(criteriaCuts.at(criterion2) + h_global2DPlot->GetYaxis()->GetXmax()), "N4");
          c->SaveAs(outputFileName.c_str());
        }
      }
      else {
        std::cout << "ERROR: plot with name \"" << plot2DName << "\" not found in input file " << inputFileName << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  }
  inputFile->Close();
  return correlations;
}

void getFakeToMediumRatioFromFile(const std::string& inputFileName, const std::string& histogramName, const float& chIso_cutMedium, const float& chIso_cutLoose, const float& sigmaIEtaIEta_cutMedium, const float& sigmaIEtaIEta_cutLoose, const bool& isTruthMatched, const std::string& inputType) {
  std::cout << "From file: " << inputFileName << std::endl;
  float nMedium = 0;
  float nFake = 0;
  TFile *inputFile = TFile::Open(inputFileName.c_str(), "READ");
  TH2F* h_mediumFakeCriteria;
  inputFile->GetObject(histogramName.c_str(), h_mediumFakeCriteria);
  h_mediumFakeCriteria->StatOverflows(kTRUE);
  if (h_mediumFakeCriteria) {
    // std::pair<double, double> ratioAndError = getRatioFromTH1(chIso, 1.141, 6.0);
    // std::cout << "Ratio = " << ratioAndError.first << " +/- " << ratioAndError.second << std::endl;
    for (int xbincounter = 0; xbincounter <= 1+h_mediumFakeCriteria->GetXaxis()->GetNbins(); ++xbincounter) {
      float bin_sigmaIEtaIEta = h_mediumFakeCriteria->GetXaxis()->GetBinCenter(xbincounter);
      for (int ybincounter = 0; ybincounter <= 1+h_mediumFakeCriteria->GetYaxis()->GetNbins(); ++ybincounter) {
        float bin_chIso = h_mediumFakeCriteria->GetYaxis()->GetBinCenter(ybincounter);
        if ((bin_sigmaIEtaIEta < sigmaIEtaIEta_cutMedium) && (bin_chIso < chIso_cutMedium)) {
          nMedium += h_mediumFakeCriteria->GetBinContent(xbincounter, ybincounter);
        }
        else if ((bin_sigmaIEtaIEta < sigmaIEtaIEta_cutLoose) && (bin_chIso < chIso_cutLoose)) {
          nFake += h_mediumFakeCriteria->GetBinContent(xbincounter, ybincounter);
        }
      }
    }
  }
  else {
    std::cout << "ERROR: Unable to find histogram named: " << histogramName << std::endl;
    std::exit(EXIT_FAILURE);
  }
  float fakeToMediumRatio = nFake/nMedium;
  std::cout << "Ratio loose/medium: " << fakeToMediumRatio << std::endl;
  std::cout << "With this ratio, (leakage into control)/signal: " << fakeToMediumRatio*fakeToMediumRatio + 2*fakeToMediumRatio << std::endl;
  TLine *lines = new TLine();
  TText *text = new TText();
  text->SetTextAlign(22);
  std::string outputFileName = "plots/NMinus2/mediumFakeCriteria_";
  if (isTruthMatched) outputFileName += "TruthMatched_";
  outputFileName += inputType;

  TCanvas *c = new TCanvas(("c_output_" + outputFileName).c_str(), "c_output");
  outputFileName += ".png";
  gPad->SetLogz();
  gStyle->SetOptStat(0);
  h_mediumFakeCriteria->Draw("colz");
  h_mediumFakeCriteria->GetZaxis()->SetRangeUser(1, 1000);
  lines->DrawLine(sigmaIEtaIEta_cutMedium, h_mediumFakeCriteria->GetYaxis()->GetXmin(), sigmaIEtaIEta_cutMedium, h_mediumFakeCriteria->GetYaxis()->GetXmax());
  lines->DrawLine(sigmaIEtaIEta_cutLoose, h_mediumFakeCriteria->GetYaxis()->GetXmin(), sigmaIEtaIEta_cutLoose, h_mediumFakeCriteria->GetYaxis()->GetXmax());
  lines->DrawLine(h_mediumFakeCriteria->GetXaxis()->GetXmin(), chIso_cutMedium, h_mediumFakeCriteria->GetXaxis()->GetXmax(), chIso_cutMedium);
  lines->DrawLine(h_mediumFakeCriteria->GetXaxis()->GetXmin(), chIso_cutLoose, h_mediumFakeCriteria->GetXaxis()->GetXmax(), chIso_cutLoose);
  text->SetTextColor(kBlue);
  text->DrawText(0.5*(h_mediumFakeCriteria->GetXaxis()->GetXmin() + sigmaIEtaIEta_cutMedium), 0.5*(h_mediumFakeCriteria->GetYaxis()->GetXmin() + chIso_cutMedium), "medium");
  text->SetTextColor(kRed);
  text->DrawText(0.5*(sigmaIEtaIEta_cutMedium+sigmaIEtaIEta_cutLoose), 0.5*(h_mediumFakeCriteria->GetYaxis()->GetXmin() + chIso_cutMedium), "fake");
  text->DrawText(0.5*(sigmaIEtaIEta_cutMedium+sigmaIEtaIEta_cutLoose), 0.5*(chIso_cutMedium+chIso_cutLoose), "fake");
  text->DrawText(0.5*(h_mediumFakeCriteria->GetXaxis()->GetXmin()+sigmaIEtaIEta_cutMedium), 0.5*(chIso_cutMedium+chIso_cutLoose), "fake");
  c->SaveAs(outputFileName.c_str());
  inputFile->Close();
}

float getPredictedStepEfficiency(const float& correlationFactor, const float& globalEfficiency_X, const float& globalEfficiency_Y) {
  float sprime = (1+correlationFactor)/(1-correlationFactor);
  float gxprime = globalEfficiency_X/(1-globalEfficiency_X);
  float gyprime = globalEfficiency_Y/(1-globalEfficiency_Y);
  float sq_coeff = sprime*gyprime*(1+gxprime);
  assert(std::fabs(sq_coeff) > 0.00001);
  float linear_coeff = (gxprime*gyprime) + sprime*(gyprime-gxprime) - 1.0;
  float const_coeff = (-1.0)*(1.0+gxprime);
  float quadraticSolution = (((-1.0*linear_coeff) + std::sqrt((linear_coeff*linear_coeff)+(-4.0*sq_coeff*const_coeff)))/(2.0*sq_coeff));
  return (1.0/(1.0+quadraticSolution));
}

void saveEfficiencies(const std::map<std::string, std::string>& inputFiles, const std::vector<std::string>& photonIDCriteria, const std::vector<std::string>& efficiencyTypes, const std::map<std::string, int>& colors) {
  for (const auto& criterion: photonIDCriteria) {
    for (const auto& efficiencyType: efficiencyTypes) {
      std::string outputFileName = "plots/efficiency/IDEfficiencies_" + efficiencyType + "_" + criterion;
      TCanvas *c = new TCanvas(("c_output_" + outputFileName).c_str(), "c_output");
      outputFileName += ".png";
      gPad->SetLogy(0);
      TLegend *legend = new TLegend(0.6, 0.1, 0.9, 0.3);
      bool isFirstIteration = true;
      for (const auto& inputFilesElement: inputFiles) {
        const std::string& inputType = inputFilesElement.first;
        const std::string& inputFileName = inputFilesElement.second;
        std::cout << "Reading from: inputType = " << inputType << ", inputFileName = " << inputFileName << std::endl;
        TFile *inputFile = TFile::Open(inputFileName.c_str(), "READ");
        TEfficiency* efficiency;
        inputFile->GetObject((criterion + "_ETEfficiency_" + efficiencyType + "_TruthMatched").c_str(), efficiency);
        std::cout << "Name of histo: " << (criterion + "_ETEfficiency_" + efficiencyType + "_TruthMatched") << std::endl;
        if (isFirstIteration) {
          efficiency->Draw();
          gPad->Update();
          efficiency->GetPaintedGraph()->SetMinimum(-0.2);
          efficiency->GetPaintedGraph()->SetMaximum(1.2);
          gPad->Update();
          isFirstIteration = false;
        }
        else {
          efficiency->Draw("same");
          gPad->Update();
        }
        efficiency->SetLineColor(colors.at(inputType));
        efficiency->SetMarkerColor(colors.at(inputType));
        TLegendEntry* legendEntry = legend->AddEntry(efficiency, inputType.c_str());
        legendEntry->SetLineColor(colors.at(inputType));
        legendEntry->SetTextColor(colors.at(inputType));
        legendEntry->SetMarkerColor(colors.at(inputType));
        inputFile->Close();
      }
      legend->Draw();
      gPad->Update();
      c->SaveAs(outputFileName.c_str());
    }
  }
  std::string outputFileName = "plots/efficiency/IDEfficiencies_overall";
  TCanvas *c = new TCanvas(("c_output_" + outputFileName).c_str(), "c_output");
  outputFileName += ".png";
  gPad->SetLogy(0);
  TLegend *legend = new TLegend(0.6, 0.1, 0.9, 0.3);
  bool isFirstIteration = true;
  for (const auto& inputFilesElement: inputFiles) {
    const std::string& inputType = inputFilesElement.first;
    const std::string& inputFileName = inputFilesElement.second;
    TFile *inputFile = TFile::Open(inputFileName.c_str(), "READ");
    TEfficiency* efficiency;
    inputFile->GetObject("overallETEfficiency_TruthMatched", efficiency);
    if (isFirstIteration) {
      efficiency->Draw();
      gPad->Update();
      efficiency->GetPaintedGraph()->SetMinimum(-0.2);
      efficiency->GetPaintedGraph()->SetMaximum(1.2);
      gPad->Update();
      isFirstIteration = false;
    }
    else {
      efficiency->Draw("same");
      gPad->Update();
    }
    efficiency->SetLineColor(colors.at(inputType));
    efficiency->SetMarkerColor(colors.at(inputType));
    TLegendEntry* legendEntry = legend->AddEntry(efficiency, inputType.c_str());
    legendEntry->SetLineColor(colors.at(inputType));
    legendEntry->SetTextColor(colors.at(inputType));
    legendEntry->SetMarkerColor(colors.at(inputType));
    inputFile->Close();
  }
  legend->Draw();
  gPad->Update();
  c->SaveAs(outputFileName.c_str());
}

int main(int argc, char* argv[]) {
  (void)argc;
  (void)argv;
  gROOT->SetBatch();
  std::streamsize original_precision = std::cout.precision();

  const float& chIso_cutMedium = 1.141;
  const float& chIso_cutLoose = 6.0;
  const float& sigmaIEtaIEta_cutMedium = 0.01015;
  const float& sigmaIEtaIEta_cutLoose = 0.02;

  std::cout << "Getting ratio of number of photons in fake range to number of photons in good range..." << std::endl;
  std::string prefix = "output_";
  std::string suffix = ".root";
  std::vector<std::string> inputTypes = {"fastsim_lowMass", "fullsim_lowMass", "hgg"};
  assert(inputTypes.size() == 3);
  std::map<std::string, std::string> inputFileNames = {
    {"fastsim_lowMass", "output_stealth_privateMC_fastsim_lowMass.root"},
    {"fullsim_lowMass", "output_stealth_privateMC_fullsim_lowMass.root"},
    {"hgg", "output_hgg.root"}
  };
  assert(inputFileNames.size() == 3);
  std::map<std::string, int> colors = {
    {"fastsim_lowMass", kBlue},
    {"fullsim_lowMass", kRed},
    {"hgg", kGreen}
  };
  assert(colors.size() == 3);
  for (const auto& inputType: inputTypes) {
    std::cout << inputType << ": without truth matching: " << std::endl;
    getFakeToMediumRatioFromFile(inputFileNames[inputType], "mediumFakeCriteria", chIso_cutMedium, chIso_cutLoose, sigmaIEtaIEta_cutMedium, sigmaIEtaIEta_cutLoose, false, inputType);
    std::cout << inputType << ": with truth matching: " << std::endl;
    getFakeToMediumRatioFromFile(inputFileNames[inputType], "mediumFakeCriteria_TruthMatched", chIso_cutMedium, chIso_cutLoose, sigmaIEtaIEta_cutMedium, sigmaIEtaIEta_cutLoose, true, inputType);
  }

  std::cout << "Now beginning to calculate ID efficiencies..." << std::endl;
  std::vector<std::string> photonIDCriteria = {"hOverE", "sigmaIEtaIEta", "chIso", "neutIso", "phoIso"};
  std::map<std::string, unsigned int> photonIDCriterionIndices;
  for (unsigned int criterionIndex = 0; criterionIndex < photonIDCriteria.size(); ++criterionIndex) {
    const std::string& criterion = photonIDCriteria.at(criterionIndex);
    photonIDCriterionIndices[criterion] = criterionIndex;
  }
  std::map<std::string, float> criteriaCuts = {
    {"hOverE", 0.02197},
    {"sigmaIEtaIEta", 0.01015},
    {"chIso", 1.141},
    {"neutIso", 1.0},
    {"phoIso", 1.0}
  };

  std::map<std::string, std::map<std::string, std::map<IDEfficiencyType, float> > > efficiencies;
  for (const auto& inputType: inputTypes) {
    efficiencies[inputType] = getIDEfficienciesFromFile(inputFileNames.at(inputType), criteriaCuts, inputType);
  }

  std::cout << "(N-1) efficiencies in " << inputTypes.at(0) << ", " << inputTypes.at(1) << ", and " << inputTypes.at(2) << ":" << std::endl;
  std::cout << "\\begin{tabular}{|p{0.25\\textwidth}||p{0.2\\textwidth}|p{0.2\\textwidth}|>{\\columncolor[gray]{0.8}}p{0.2\\textwidth}||}" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "Criterion & efficiency (" << inputTypes.at(0) << ") & efficiency (" << inputTypes.at(1) << ") & efficiency (" << inputTypes.at(2) << ")\\\\ \\hline \\hline" << std::endl;
  for (unsigned int criterionIndex = 0; criterionIndex < photonIDCriteria.size(); ++criterionIndex) {
    const std::string& criterion = photonIDCriteria.at(criterionIndex);
    // std::cout << criterion << " & " << std::setprecision(3) << efficiencies[inputTypes.at(0)][criterion][IDEfficiencyType::NMinus1] << " & " << efficiencies[inputTypes.at(0)][criterion][IDEfficiencyType::global] << " & " << efficiencies[inputTypes.at(1)][criterion][IDEfficiencyType::NMinus1] << " & " << efficiencies[inputTypes.at(1)][criterion][IDEfficiencyType::global] << std::setprecision(original_precision) << "\\\\ \\hline" << std::endl;
    std::cout << criterion << " & " << std::setprecision(3) << efficiencies[inputTypes.at(0)][criterion][IDEfficiencyType::NMinus1] << " & " << efficiencies[inputTypes.at(1)][criterion][IDEfficiencyType::NMinus1] << " & " << efficiencies[inputTypes.at(2)][criterion][IDEfficiencyType::NMinus1] << std::setprecision(original_precision) << "\\\\ \\hline" << std::endl;
  }
  std::cout << "Overall & " << std::setprecision(3) << efficiencies[inputTypes.at(0)]["overall"][IDEfficiencyType::nIDEfficiencyTypes] << std::setprecision(original_precision) << " & " << std::setprecision(3) << efficiencies[inputTypes.at(1)]["overall"][IDEfficiencyType::nIDEfficiencyTypes] << std::setprecision(original_precision) << " & " << std::setprecision(3) << efficiencies[inputTypes.at(2)]["overall"][IDEfficiencyType::nIDEfficiencyTypes] << std::setprecision(original_precision) << "\\\\ \\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;

  std::cout << "Global efficiencies in " << inputTypes.at(0) << ", " << inputTypes.at(1) << ", and " << inputTypes.at(2) << ":" << std::endl;
  std::cout << "\\begin{tabular}{|p{0.25\\textwidth}||p{0.2\\textwidth}|p{0.2\\textwidth}|>{\\columncolor[gray]{0.8}}p{0.2\\textwidth}||}" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "Criterion & efficiency (" << inputTypes.at(0) << ") & efficiency (" << inputTypes.at(1) << ") & efficiency (" << inputTypes.at(2) << ")\\\\ \\hline \\hline" << std::endl;
  for (unsigned int criterionIndex = 0; criterionIndex < photonIDCriteria.size(); ++criterionIndex) {
    const std::string& criterion = photonIDCriteria.at(criterionIndex);
    // std::cout << criterion << " & " << std::setprecision(3) << efficiencies[inputTypes.at(0)][criterion][IDEfficiencyType::NMinus1] << " & " << efficiencies[inputTypes.at(0)][criterion][IDEfficiencyType::global] << " & " << efficiencies[inputTypes.at(1)][criterion][IDEfficiencyType::NMinus1] << " & " << efficiencies[inputTypes.at(1)][criterion][IDEfficiencyType::global] << std::setprecision(original_precision) << "\\\\ \\hline" << std::endl;
    std::cout << criterion << " & " << std::setprecision(3) << efficiencies[inputTypes.at(0)][criterion][IDEfficiencyType::global] << " & " << efficiencies[inputTypes.at(1)][criterion][IDEfficiencyType::global] << " & " << efficiencies[inputTypes.at(2)][criterion][IDEfficiencyType::global] << std::setprecision(original_precision) << "\\\\ \\hline" << std::endl;
  }
  std::cout << "Overall & " << std::setprecision(3) << efficiencies[inputTypes.at(0)]["overall"][IDEfficiencyType::nIDEfficiencyTypes] << std::setprecision(original_precision) << " & " << std::setprecision(3) << efficiencies[inputTypes.at(1)]["overall"][IDEfficiencyType::nIDEfficiencyTypes] << std::setprecision(original_precision) << " & " << std::setprecision(3) << efficiencies[inputTypes.at(2)]["overall"][IDEfficiencyType::nIDEfficiencyTypes] << std::setprecision(original_precision) << "\\\\ \\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;

  // std::vector<std::string> diphotonHLTCriteria = {"R9_leading", "hOverE_leading", "sigmaIEtaIEta_leading", "clusIso_leading", "trkIso_leading", "R9_subLeading", "hOverE_subLeading", "sigmaIEtaIEta_subLeading", "clusIso_subLeading", "trkIso_subLeading"};
  // std::map<std::string, float> hltCriteriaCuts = {
  //   {"R9_leading", 0.5},
  //   {"hOverE_leading", 0.12},
  //   {"sigmaIEtaIEta_leading", 0.015},
  //   {"clusIso_leading", 1.},
  //   {"trkIso_leading", 1.},
  //   {"R9_subLeading", 0.5},
  //   {"hOverE_subLeading", 0.12},
  //   {"sigmaIEtaIEta_subLeading", 0.015},
  //   {"clusIso_subLeading", 1.},
  //   {"trkIso_subLeading", 1.}
  // };
  // std::vector<std::string> HLTCriteriaToInvert = {"R9_leading", "R9_subLeading"};
  // std::map<std::string, std::map<std::string, std::map<IDEfficiencyType, float> > > hlt_efficiencies;
  // for (const auto& inputType: inputTypes) {
  //   hlt_efficiencies[inputType] = getHLTEfficienciesFromFile(inputFileNames.at(inputType), hltCriteriaCuts, HLTCriteriaToInvert, inputType);
  // }
  // std::cout << "HLT Efficiencies in " << inputTypes.at(0) << " and " << inputTypes.at(1) << ":" << std::endl;
  // std::cout << "\\begin{tabular}{|p{0.2\\textwidth}||p{0.135\\textwidth}|>{\\columncolor[gray]{0.8}}p{0.135\\textwidth}||p{0.135\\textwidth}|>{\\columncolor[gray]{0.8}}p{0.135\\textwidth}|}" << std::endl;
  // std::cout << "\\hline" << std::endl;
  // std::cout << "Criterion & (N-1) efficiency (" << inputTypes.at(0) << ") & global efficiency (" << inputTypes.at(0) << ") & (N-1) efficiency (" << inputTypes.at(1) << ") & global efficiency (" << inputTypes.at(1) << ")\\\\ \\hline \\hline" << std::endl;
  // for (unsigned int criterionIndex = 0; criterionIndex < diphotonHLTCriteria.size(); ++criterionIndex) {
  //   const std::string& criterion = diphotonHLTCriteria.at(criterionIndex);
  //   std::cout << criterion << " & " << std::setprecision(3) << hlt_efficiencies[inputTypes.at(0)][criterion][IDEfficiencyType::NMinus1] << " & " << hlt_efficiencies[inputTypes.at(0)][criterion][IDEfficiencyType::global] << " & " << hlt_efficiencies[inputTypes.at(1)][criterion][IDEfficiencyType::NMinus1] << " & " << hlt_efficiencies[inputTypes.at(1)][criterion][IDEfficiencyType::global] << std::setprecision(original_precision) << "\\\\ \\hline" << std::endl;
  // }
  // // std::cout << "Overall & \\multicolumn{2}{c||}{" << std::setprecision(3) << hlt_efficiencies[inputTypes.at(0)]["overall"][IDEfficiencyType::nIDEfficiencyTypes] << std::setprecision(original_precision) << "} & \\multicolumn{2}{c|}{" << std::setprecision(3) << hlt_efficiencies[inputTypes.at(1)]["overall"][IDEfficiencyType::nIDEfficiencyTypes] << std::setprecision(original_precision) << "}\\\\ \\hline" << std::endl;
  // std::cout << "\\end{tabular}" << std::endl;

  // std::cout << "Now beginning to calculate correlations..." << std::endl;
  // std::map<std::string, std::map<correlationType, std::map<std::string, std::map<std::string, float> > > > correlations;
  // std::map<std::string, std::map<correlationType, std::map<std::string, std::map<std::string, float> > > > correlations_NMinus2;
  // for (const auto& inputType: inputTypes) {
  //   correlations[inputType] = getCorrelationsFromFile(inputFileNames.at(inputType), photonIDCriteria, criteriaCuts, false, inputType);
  //   correlations_NMinus2[inputType] = getCorrelationsFromFile(inputFileNames.at(inputType), photonIDCriteria, criteriaCuts, true, inputType);
  // }

  // std::cout << "Correlations in " << inputTypes.at(0) << " and " << inputTypes.at(1) << " (customized):" << std::endl;
  // std::cout << "\\begin{tabular}{|p{0.2\\textwidth}|p{0.2\\textwidth}||p{0.15\\textwidth}|p{0.15\\textwidth}|}" << std::endl;
  // std::cout << "\\hline" << std::endl;
  // std::cout << "Criterion1 & Criterion2 & $s$ (" << inputTypes.at(0) << ") & $s$ (" << inputTypes.at(1) << ")\\\\ \\hline \\hline" << std::endl;
  // for (unsigned int criterion1Index = 0; criterion1Index < (-1+photonIDCriteria.size()); ++criterion1Index) {
  //   const std::string& criterion1 = photonIDCriteria.at(criterion1Index);
  //   for (unsigned int criterion2Index = (1+criterion1Index); criterion2Index < photonIDCriteria.size(); ++criterion2Index) {
  //     const std::string& criterion2 = photonIDCriteria.at(criterion2Index);
  //     // Formatting LaTeX-style
  //     std::cout << criterion1 << " & " << criterion2 << " & " << std::setprecision(3) << correlations[inputTypes.at(0)][correlationType::customized][criterion1][criterion2] << " & " << correlations[inputTypes.at(1)][correlationType::customized][criterion1][criterion2] << std::setprecision(original_precision) << "\\\\ \\hline" << std::endl;
  //   }
  // }
  // std::cout << "\\end{tabular}" << std::endl;

  // std::cout << "Correlations in " << inputTypes.at(0) << " and " << inputTypes.at(1) << " (Pearson):" << std::endl;
  // std::cout << "\\hline" << std::endl;
  // std::cout << "Criterion1 & Criterion2 & $s$ (" << inputTypes.at(0) << ") & $s$ (" << inputTypes.at(1) << ")\\\\ \\hline \\hline" << std::endl;
  // for (unsigned int criterion1Index = 0; criterion1Index < (-1+photonIDCriteria.size()); ++criterion1Index) {
  //   const std::string& criterion1 = photonIDCriteria.at(criterion1Index);
  //   for (unsigned int criterion2Index = (1+criterion1Index); criterion2Index < photonIDCriteria.size(); ++criterion2Index) {
  //     const std::string& criterion2 = photonIDCriteria.at(criterion2Index);
  //     // Formatting LaTeX-style
  //     std::cout << criterion1 << " & " << criterion2 << " & " << std::setprecision(3) << correlations[inputTypes.at(0)][correlationType::pearson][criterion1][criterion2] << " & " << correlations[inputTypes.at(1)][correlationType::pearson][criterion1][criterion2] << std::setprecision(original_precision) << "\\\\ \\hline" << std::endl;
  //   }
  // }

  // std::cout << "Now beginning to calculate step-by-step efficiencies..." << std::endl;
  // std::vector<std::vector<std::string> > stepByStepSequences = {
  //   {"chIso", "hOverE", "sigmaIEtaIEta", "neutIso", "phoIso"},
  //   {"hOverE", "chIso", "sigmaIEtaIEta", "neutIso", "phoIso"},
  //   {"hOverE", "sigmaIEtaIEta", "chIso", "neutIso", "phoIso"},
  //   {"hOverE", "sigmaIEtaIEta", "neutIso", "chIso", "phoIso"},
  //   {"hOverE", "sigmaIEtaIEta", "neutIso", "phoIso", "chIso"},
  //   {"chIso", "phoIso", "hOverE", "sigmaIEtaIEta", "neutIso"},
  //   {"phoIso", "chIso", "hOverE", "sigmaIEtaIEta", "neutIso"}
  // };
  // std::map<std::string, std::map<unsigned int, std::map<unsigned int, float> > > stepByStepEfficiencies;
  // for (const auto& inputType: inputTypes) {
  //   stepByStepEfficiencies[inputType] = getStepByStepEfficienciesFromFile(inputFileNames.at(inputType), photonIDCriteria, criteriaCuts, stepByStepSequences, inputType);
  // }
  // for (unsigned int sequenceIndex = 0; sequenceIndex < stepByStepSequences.size(); ++sequenceIndex) {
  //   std::cout << "\\begin{tabular}{|p{0.275\\textwidth}||p{0.13\\textwidth}|>{\\columncolor[gray]{0.8}}p{0.13\\textwidth}||p{0.13\\textwidth}|>{\\columncolor[gray]{0.8}}p{0.13\\textwidth}|}" << std::endl;
  //   std::cout << "\\hline" << std::endl;
  //   std::cout << "\\multicolumn{5}{|c|}{Sequence: ";
  //   const std::vector<std::string>& sequence = stepByStepSequences.at(sequenceIndex);
  //   for (unsigned int stepIndex = 0; stepIndex < sequence.size(); ++stepIndex) {
  //     std::cout << sequence.at(stepIndex);
  //     if (stepIndex < (sequence.size()-1)) std::cout << " $\\rightarrow$ ";
  //   }
  //   std::cout << "}\\\\ \\hline \\hline" << std::endl;
  //   std::cout << "Step & observed efficiency (" << inputTypes.at(0) << ") & predicted efficiency (" << inputTypes.at(0) << ") & observed efficiency (" << inputTypes.at(1) << ") & predicted efficiency (" << inputTypes.at(1) << ") \\\\ \\hline \\hline" << std::endl;
  //   std::map<std::string, float> overallEfficiencies;
  //   std::map<std::string, float> predictedOverallEfficiencies;
  //   for (const auto& inputType: inputTypes) {
  //     overallEfficiencies[inputType] = 1.;
  //     predictedOverallEfficiencies[inputType] = 1.;
  //   }
  //   std::string headCriterion = sequence.at(0);
  //   unsigned int headCriterionIndex = photonIDCriterionIndices.at(headCriterion);
  //   for (unsigned int stepIndex = 0; stepIndex < photonIDCriteria.size(); ++stepIndex) {
  //     unsigned int stepNumber = stepIndex + 1;
  //     const std::string& criterion = sequence.at(stepIndex);
  //     for (const auto& inputType: inputTypes) {
  // 	overallEfficiencies[inputType] *= stepByStepEfficiencies[inputType][sequenceIndex][stepIndex];
  //     }
  //     unsigned int currentCriterionIndex = photonIDCriterionIndices.at(criterion);
  //     std::map<std::string, float> predictedEfficiencies;
  //     for (const auto& inputType: inputTypes) {
  // 	predictedEfficiencies[inputType] = -1.;
  //     }
  //     if (stepIndex == 0) {
  // 	for (const auto& inputType: inputTypes) {
  // 	  predictedEfficiencies[inputType] = efficiencies[inputType][criterion][IDEfficiencyType::global];
  // 	}
  //     }
  //     else if (stepIndex <= 3) {
  //       unsigned int criterionIndex1;
  //       unsigned int criterionIndex2;
  //       if (currentCriterionIndex < headCriterionIndex) {
  //         criterionIndex1 = currentCriterionIndex;
  //         criterionIndex2 = headCriterionIndex;
  //       }
  //       else {
  //         criterionIndex1 = headCriterionIndex;
  //         criterionIndex2 = currentCriterionIndex;
  //       }
  //       if (stepIndex == 1) {
  // 	  for (const auto& inputType: inputTypes) {
  // 	    predictedEfficiencies[inputType] = getPredictedStepEfficiency(correlations[inputType][correlationType::customized][photonIDCriteria.at(criterionIndex1)][photonIDCriteria.at(criterionIndex2)], efficiencies[inputType][headCriterion][IDEfficiencyType::global], efficiencies[inputType][criterion][IDEfficiencyType::global]);
  // 	  }
  //       }
  //       else if (stepIndex == 3) {
  // 	  for (const auto& inputType: inputTypes) {
  // 	    predictedEfficiencies[inputType] = getPredictedStepEfficiency(correlations_NMinus2[inputType][correlationType::customized][photonIDCriteria.at(criterionIndex1)][photonIDCriteria.at(criterionIndex2)], efficiencies[inputType][headCriterion][IDEfficiencyType::global], efficiencies[inputType][criterion][IDEfficiencyType::global]);
  // 	  }
  //       }
  //       else if (stepIndex == 2) {
  // 	  for (const auto& inputType: inputTypes) {
  // 	    predictedEfficiencies[inputType] = 0.5*(getPredictedStepEfficiency(correlations_NMinus2[inputType][correlationType::customized][photonIDCriteria.at(criterionIndex1)][photonIDCriteria.at(criterionIndex2)], efficiencies[inputType][headCriterion][IDEfficiencyType::global], efficiencies[inputType][criterion][IDEfficiencyType::global]) + getPredictedStepEfficiency(correlations[inputType][correlationType::customized][photonIDCriteria.at(criterionIndex1)][photonIDCriteria.at(criterionIndex2)], efficiencies[inputType][headCriterion][IDEfficiencyType::global], efficiencies[inputType][criterion][IDEfficiencyType::global]));
  // 	  }
  //       }
  //     }
  //     else {
  // 	for (const auto& inputType: inputTypes) {
  // 	  predictedEfficiencies[inputType] = efficiencies[inputType][criterion][IDEfficiencyType::NMinus1];
  // 	}
  //     }
  //     for (const auto& inputType: inputTypes) {
  // 	predictedOverallEfficiencies[inputType] *= predictedEfficiencies.at(inputType);
  //     }
  //     std::cout << stepNumber << " (" << criterion << ") & " << std::setprecision(3) << stepByStepEfficiencies[inputTypes.at(0)][sequenceIndex][stepIndex] << " & " << predictedEfficiencies.at(inputTypes.at(0)) << " & " << stepByStepEfficiencies[inputTypes.at(1)][sequenceIndex][stepIndex] << " & " << predictedEfficiencies.at(inputTypes.at(1)) << "\\\\ \\hline" << std::setprecision(original_precision) << std::endl;
  //   }
  //   std::cout << "\\hline" << std::endl;
  //   std::cout << "Overall & " << std::setprecision(3) << overallEfficiencies.at(inputTypes.at(0)) << " & " << predictedOverallEfficiencies.at(inputTypes.at(0)) << " & " << overallEfficiencies.at(inputTypes.at(1)) << " & " << predictedOverallEfficiencies.at(inputTypes.at(1)) << "\\\\ \\hline" << std::setprecision(original_precision) << std::endl;
  //   std::cout << "\\end{tabular}" << std::endl;
  //   std::cout << std::endl << std::endl;
  // }

  // std::cout << "Now saving ID efficiencies..." << std::endl;
  // std::vector<std::string> efficiencyTypes = {"global", "NMinus1"};
  // saveEfficiencies(inputFiles, photonIDCriteria, efficiencyTypes, colors);

  std::cout << "All done!" << std::endl;
  return 0;
}
