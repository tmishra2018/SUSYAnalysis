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

const int nBG = 6; // Number of background files
vector<string> bgFiles = {
    "/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_DYJetsToLL_2018.root",
    "/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_TTJets_2018.root",
    "/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_WGToLNuG_2018.root",
    "/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_WJetsToLNu_2018.root",
    "/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_WW_2018.root",
    "/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_WZ_2018.root"
};

vector<float> crossSections = {6424.0, 734.6, 551.008, 61526.7, 75.95, 27.59}; // Cross sections in pb
vector<string> bgLegend = {"DYJetsToLL", "TTJets", "WGToLNuG", "WJetsToLNu", "WW", "WZ"};

string signalFile = "/uscms/home/tmishra/nobackup/signal_trees/resTree_T5WG_2018.root";
float signalCrossSection = 0.001057; // Cross section in pb
vector<int> colors = {kPink + 1, kTeal + 9, kYellow, kGray, kOrange, kCyan};

void plotKinematicStack() {
    gROOT->SetBatch(kTRUE);

    vector<string> variables = {"lepPt", "lepEta", "lepPhi"};
    vector<string> xLabels = {"Muon p_{T} [GeV]", "Muon #eta", "Muon #phi"};

    TFile *bgFilesT[nBG];
    TTree *bgTrees[nBG];
    TFile *signalFileT = TFile::Open(signalFile.c_str());
    TTree *signalTree = (TTree *)signalFileT->Get("mgTree");

    vector<TH1D *> bgHists[nBG];
    vector<TH1D *> signalHists;
    vector<THStack *> stacks;

    // Prepare histograms for each variable with custom binning
    for (size_t i = 0; i < variables.size(); i++) {
        THStack *stack = new THStack((variables[i] + "_stack").c_str(), "");
        stacks.push_back(stack);

        // Prepare signal histogram with weights
        if (variables[i] == "lepPt") {
	    Double_t ptBins[] = {20, 40, 60, 80, 100, 120, 140, 160, 200, 250, 300, 500,700, 900, 1200, 1500};
            TH1D *sigHist = new TH1D(("signal_" + variables[i]).c_str(), "", 15, ptBins);
            signalHists.push_back(sigHist);
        } else if (variables[i] == "lepEta") {
            TH1D *sigHist = new TH1D(("signal_" + variables[i]).c_str(), "", 90, -3, 3);
            signalHists.push_back(sigHist);
        } else if (variables[i] == "lepPhi") {
            TH1D *sigHist = new TH1D(("signal_" + variables[i]).c_str(), "", 48, 0, 3.2);
            signalHists.push_back(sigHist);
        }

        // Load background files and histograms
        for (int j = 0; j < nBG; j++) {
            bgFilesT[j] = TFile::Open(bgFiles[j].c_str());
            bgTrees[j] = (TTree *)bgFilesT[j]->Get("signalTree");

            if (variables[i] == "lepPt") {
	    	Double_t ptBins[] = {20, 40, 60, 80, 100, 120, 140, 160, 200, 250, 300, 500,700, 900, 1200, 1500};
                TH1D *bgHist = new TH1D(("bg_" + to_string(j) + "_" + variables[i]).c_str(), "", 15, ptBins);
                bgHists[j].push_back(bgHist);
            } else if (variables[i] == "lepEta") {
                TH1D *bgHist = new TH1D(("bg_" + to_string(j) + "_" + variables[i]).c_str(), "", 90, -3, 3);
                bgHists[j].push_back(bgHist);
            } else if (variables[i] == "lepPhi") {
                TH1D *bgHist = new TH1D(("bg_" + to_string(j) + "_" + variables[i]).c_str(), "", 48, 0, 3.2);
                bgHists[j].push_back(bgHist);
            }
        }
    }

    // Read data and fill histograms with weights
    float lepPt, lepEta, lepPhi;
    for (int j = 0; j < nBG; j++) {
        bgTrees[j]->SetBranchAddress("lepPt", &lepPt);
        bgTrees[j]->SetBranchAddress("lepEta", &lepEta);
        bgTrees[j]->SetBranchAddress("lepPhi", &lepPhi);

        Long64_t nEntries = bgTrees[j]->GetEntries();
        float XS_weight = 59.816229 * 1000 * crossSections[j] / nEntries;

        for (Long64_t i = 0; i < nEntries; i++) {
            bgTrees[j]->GetEntry(i);
            bgHists[j][0]->Fill(lepPt, XS_weight);
            bgHists[j][1]->Fill(lepEta, XS_weight);
            bgHists[j][2]->Fill(lepPhi, XS_weight);
        }
    }

    // Signal data
    signalTree->SetBranchAddress("lepPt", &lepPt);
    signalTree->SetBranchAddress("lepEta", &lepEta);
    signalTree->SetBranchAddress("lepPhi", &lepPhi);

    Long64_t nSignalEntries = signalTree->GetEntries();
    float signalXS_weight = 59.816229 * 1000 * signalCrossSection / nSignalEntries;

    for (Long64_t i = 0; i < nSignalEntries; i++) {
        signalTree->GetEntry(i);
        signalHists[0]->Fill(lepPt, signalXS_weight * 10000); // Scale signal by 1000
        signalHists[1]->Fill(lepEta, signalXS_weight * 10000);
        signalHists[2]->Fill(lepPhi, signalXS_weight * 10000);
    }

    // Create stacked plots
    for (size_t i = 0; i < variables.size(); i++) {
        vector<pair<TH1D *, double>> integralOrder;
        for (int j = 0; j < nBG; j++) {
            integralOrder.push_back(make_pair(bgHists[j][i], bgHists[j][i]->Integral()));
        }

        // Sort histograms by integral value
        sort(integralOrder.begin(), integralOrder.end(), [](const pair<TH1D *, double> &a, const pair<TH1D *, double> &b) {
            return a.second > b.second;
        });

        TCanvas *c = new TCanvas(("canvas_" + variables[i]).c_str(), (variables[i] + " Canvas").c_str(), 800, 600);
        c->SetLogy();  // Set Y-axis to log scale

        TLegend *legend = new TLegend(0.65, 0.67, 0.85, 0.87);
	legend->SetTextSize(0.01); 
	legend->SetBorderSize(0);


        double maxStackHeight = 0;
        for (size_t j = integralOrder.size(); j > 0; j--) {
            integralOrder[j - 1].first->SetFillColor(colors[j - 1]);
            integralOrder[j - 1].first->SetLineColor(kBlack);
            stacks[i]->Add(integralOrder[j - 1].first);
		legend->AddEntry(integralOrder[integralOrder.size() - j].first, bgLegend[integralOrder.size() - j].c_str(), "f");

            // Track maximum height for Y-axis adjustment
            double currentMax = integralOrder[j - 1].first->GetMaximum();
            if (currentMax > maxStackHeight) maxStackHeight = currentMax;
        }

        // Signal histogram
        signalHists[i]->SetLineColor(kRed);
        signalHists[i]->SetLineStyle(2);
	legend->AddEntry(signalHists[i], "T5WG (scaled x10^{3})", "l");
	stacks[i]->SetMinimum(10000);

        // Draw the stack and signal
        stacks[i]->Draw("hist");
        stacks[i]->GetXaxis()->SetTitle(xLabels[i].c_str());
        stacks[i]->GetYaxis()->SetTitle("Events");
        stacks[i]->SetMaximum(maxStackHeight * 70.0); // Set Y-axis maximum to 2.5 times the stack maximum

        signalHists[i]->Draw("hist same");

        // Draw the legend
        legend->SetTextSize(0.03);
        legend->Draw();

        // Save the canvas
        c->SaveAs(("muon_"+variables[i] + "_SM-processes.pdf").c_str());
    }
}

