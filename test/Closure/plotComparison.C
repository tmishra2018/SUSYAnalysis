#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>

void plotComparison() {
    // Define years and colors for plotting
    TString years[4] = {"2016postVFP", "2016preVFP", "2017", "2018"};
    int colors[4] = {kRed, kBlue, kGreen, kBlack}; 

    // File paths for eγ (electron + photon) channel
    TString files_eg[4] = {
        "/eos/uscms/store/user/tmishra/elefakepho/Closure/Ele_fake_pho_Direct-Simulation_eg_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/elefakepho/Closure/Ele_fake_pho_Direct-Simulation_eg_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/elefakepho/Closure/Ele_fake_pho_Direct-Simulation_eg_2017.root",
        "/eos/uscms/store/user/tmishra/elefakepho/Closure/Ele_fake_pho_Direct-Simulation_eg_2018.root"
    };

    // File paths for μγ (muon + photon) channel
    TString files_mg[4] = {
        "/eos/uscms/store/user/tmishra/elefakepho/Closure/Ele_fake_pho_Direct-Simulation_mg_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/elefakepho/Closure/Ele_fake_pho_Direct-Simulation_mg_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/elefakepho/Closure/Ele_fake_pho_Direct-Simulation_mg_2017.root",
        "/eos/uscms/store/user/tmishra/elefakepho/Closure/Ele_fake_pho_Direct-Simulation_mg_2018.root"
    };

    // Histogram names and corresponding X-axis labels
    TString histNames[3] = {"p_MET", "p_Mt", "p_HT"};
    TString xAxisLabels[3] = {"MET", "MT", "HT"};

    // Loop over eγ and μγ channels
    for (int ch = 0; ch < 2; ch++) {
        TString channel = (ch == 0) ? "eg" : "mg";
        TString* files = (ch == 0) ? files_eg : files_mg;

        for (int v = 0; v < 3; v++) {
            TCanvas *c = new TCanvas(Form("c_%s_%s", channel.Data(), histNames[v].Data()), "", 800, 600);
            c->SetMargin(0.15, 0.05, 0.15, 0.1); // Increased margins for titles

            c->SetTitle("electron-fake-photon events from simulation"); // Set canvas title
            TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);

            bool firstHist = true;
            double maxY = 0; // Track maximum Y value to set range properly
            std::vector<TH1D*> histograms;

            for (int i = 0; i < 4; i++) {
                TFile *file = TFile::Open(files[i], "READ");
                if (!file || file->IsZombie()) {
                    std::cerr << "Error opening file: " << files[i] << std::endl;
                    continue;
                }

                TH1D *hist = (TH1D*) file->Get(histNames[v]);
                if (!hist) {
                    std::cerr << "Histogram " << histNames[v] << " not found in file: " << files[i] << std::endl;
                    file->Close();
                    continue;
                }

                // Print debug info
                std::cout << "Histogram " << histNames[v] << " from " << years[i] 
                          << " Integral: " << hist->Integral() << std::endl;

                // Scale to unity
                if (hist->Integral() > 0) {
                    hist->Scale(1.0 / hist->Integral());
                }

                histograms.push_back(hist);

                // Find the max Y value for setting axis limits
                if (hist->GetMaximum() > maxY) {
                    maxY = hist->GetMaximum();
                }

                // Style settings
                hist->SetLineColor(colors[i]);
                hist->SetLineWidth(3);
                hist->GetXaxis()->SetTitle(xAxisLabels[v]); // Set X-axis title
                hist->GetXaxis()->SetRangeUser(0, 500); 
                hist->GetXaxis()->SetLimits(0, 500);
                hist->GetXaxis()->SetTitleSize(0.07);
                hist->GetXaxis()->SetTitleOffset(1.4); // Adjust offset

                hist->GetYaxis()->SetTitle("Normalized Events");
                hist->GetYaxis()->SetTitleSize(0.07);
                hist->GetYaxis()->SetTitleOffset(1.2); // Y-axis title offset

                // Draw histograms
                if (firstHist) {
                    hist->Draw("hist");
                    firstHist = false;
                } else {
                    hist->Draw("hist same");
                }

                // Add to legend
                legend->AddEntry(hist, years[i], "l");
            }

            // Adjust Y-axis max for visibility
            for (auto& hist : histograms) {
                hist->SetMaximum(1.1 * maxY); // Increase Y max slightly for visibility
            }

            // Update canvas and draw legend
            c->Update();
            gPad->Modified();
            gPad->Update();
            legend->Draw();
            c->SaveAs(Form("%s_Comparison_%s-channel.png", xAxisLabels[v].Data(), channel.Data()));

            // Clean up memory
            delete c;
        }
    }
}

