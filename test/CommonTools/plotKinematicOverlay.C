#include <iostream>
#include <iomanip>
#include "TH1.h"
#include "TROOT.h"
#include "THStack.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <vector>
#include <string>

using namespace std;

vector<string> years = {
    "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016preVFP_Muon20.root",
    "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016postVFP_Muon20.root",
    "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2017_Muon20.root",
    "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2018_Muon20.root"
};

vector<string> yearLabels = {"2016preVFP", "2016postVFP", "2017", "2018"};
vector<int> colors = {kRed, kBlue, kGreen, kMagenta};

void plotKinematicOverlay() {
    gROOT->SetBatch(kTRUE);

    vector<string> variables = {"lepPt", "lepEta", "lepPhi"};
    vector<string> xLabels = {"Muon p_{T} [GeV]", "Muon #eta", "Muon #phi"};

    TFile *yearFiles[years.size()];
    TTree *yearTrees[years.size()];

    vector<TH1D *> yearHists[years.size()];

    // Prepare histograms for each variable with custom binning
    for (size_t i = 0; i < variables.size(); i++) {
        for (size_t j = 0; j < years.size(); j++) {
            yearFiles[j] = TFile::Open(years[j].c_str());
            yearTrees[j] = (TTree *)yearFiles[j]->Get("signalTree");

            if (variables[i] == "lepPt") {
                Double_t ptBins[] = {20, 40, 60, 80, 100, 120, 140, 160, 200, 250, 300, 500,700, 900, 1200, 1500};
                TH1D *yearHist = new TH1D((yearLabels[j] + "_" + variables[i]).c_str(), "", 15, ptBins);
                yearHists[j].push_back(yearHist);
            } else if (variables[i] == "lepEta") {
                TH1D *yearHist = new TH1D((yearLabels[j] + "_" + variables[i]).c_str(), "", 90, -3, 3);
                yearHists[j].push_back(yearHist);
            } else if (variables[i] == "lepPhi") {
                TH1D *yearHist = new TH1D((yearLabels[j] + "_" + variables[i]).c_str(), "", 48, 0, 3.2);
                yearHists[j].push_back(yearHist);
            }
        }
    }

    // Read data and fill histograms
    float lepPt, lepEta, lepPhi;
    for (size_t j = 0; j < years.size(); j++) {
        yearTrees[j]->SetBranchAddress("lepPt", &lepPt);
        yearTrees[j]->SetBranchAddress("lepEta", &lepEta);
        yearTrees[j]->SetBranchAddress("lepPhi", &lepPhi);

        Long64_t nEntries = yearTrees[j]->GetEntries();

        for (Long64_t i = 0; i < nEntries; i++) {
            yearTrees[j]->GetEntry(i);
            yearHists[j][0]->Fill(lepPt);
            yearHists[j][1]->Fill(lepEta);
            yearHists[j][2]->Fill(lepPhi);
        }
    }
     for (size_t j = 0; j < years.size(); j++) {
        for (size_t i = 0; i < variables.size(); i++) {
            yearHists[j][i]->Scale(1.0 / yearHists[j][i]->Integral());
        }
    }

    // Create overlay plots
    for (size_t i = 0; i < variables.size(); i++) {
        TCanvas *c = new TCanvas(("canvas_" + variables[i]).c_str(), (variables[i] + " Canvas").c_str(), 800, 600);
        c->SetLogy();  // Set Y-axis to log scale

        TLegend *legend = new TLegend(0.65, 0.67, 0.85, 0.87);
        legend->SetTextSize(0.03);
        legend->SetBorderSize(0);
	        legend->AddEntry((TObject*)0, "#bf{Unit normalised}", "");

        double maxHistHeight = 0;
        for (size_t j = 0; j < years.size(); j++) {
            yearHists[j][i]->SetLineColor(colors[j]);
            yearHists[j][i]->SetLineWidth(1);
            yearHists[j][i]->Draw(j == 0 ? "hist" : "hist same");
            legend->AddEntry(yearHists[j][i], yearLabels[j].c_str(), "l");

            // Track maximum height for Y-axis adjustment
            double currentMax = yearHists[j][i]->GetMaximum();
            if (currentMax > maxHistHeight) maxHistHeight = currentMax;
        }
	

        // Adjust Y-axis maximum
        yearHists[0][i]->SetMaximum(maxHistHeight * 3.5);

        // Add x-axis title
        yearHists[0][i]->GetXaxis()->SetTitle(xLabels[i].c_str());

        // Remove stat box
        gStyle->SetOptStat(0);

        // Draw the legend
        legend->Draw();

        // Save the canvas
        c->SaveAs(("muon_"+variables[i] + "_data.pdf").c_str());
    }
}
