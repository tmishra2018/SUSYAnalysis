#include <iostream>
#include <TFile.h>
#include <TH1D.h>

void printLatexTable() {
    // List of files for each year
    const char* files_egamma[] = {
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2017.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2018.root"
    };
    
    const char* files_mgamma[] = {
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016preVFP_Muon20.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016postVFP_Muon20.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2017_Muon20.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2018_Muon20.root"
    };

    // Array to store corresponding years
    const char* years[] = {"2016preVFP", "2016postVFP", "2017", "2018"};

    // Start printing LaTeX tables
    for (int year_idx = 0; year_idx < 4; ++year_idx) {
        std::cout << "\\begin{table*}[tbp]" << std::endl;
        std::cout << "  \\begin{center}" << std::endl;
        std::cout << "    \\caption{Event selection and the number of events after successive cuts for " << years[year_idx] << "." << std::endl;
        std::cout << "    \\label{tab:cutflow_" << years[year_idx] << "}" << std::endl;
        std::cout << "    \\begin{tabular}{|l|l|l|}" << std::endl;
        std::cout << "      \\hline" << std::endl;
        std::cout << "      Cut & \\Pe\\Pgg~Channel (egamma) & \\Pgm\\Pgg~Channel (mgamma) \\\\" << std::endl;
        std::cout << "      \\hline" << std::endl;
        std::cout << "      \\hline" << std::endl;

        // Open ROOT files for egamma and mgamma channels
        TFile *file_egamma = TFile::Open(files_egamma[year_idx]);
        TFile *file_mgamma = TFile::Open(files_mgamma[year_idx]);

        // Read the histograms
        TH1D* hist_egamma = (TH1D*)file_egamma->Get("p_eventcount");
        TH1D* hist_mgamma = (TH1D*)file_mgamma->Get("p_eventcount");

        // Get bin contents
        double egamma_bin1 = hist_egamma->GetBinContent(1);
        double mgamma_bin1 = hist_mgamma->GetBinContent(1);

        std::cout << "      Total & " << egamma_bin1 << " & " << mgamma_bin1 << " \\\\" << std::endl;

        // Print HLT
        double egamma_bin2 = hist_egamma->GetBinContent(2);
        double mgamma_bin2 = hist_mgamma->GetBinContent(2);
        std::cout << "      HLT & " << egamma_bin2 << " (" << (egamma_bin2 / egamma_bin1) * 100 << "\\%) & " 
                  << mgamma_bin2 << " (" << (mgamma_bin2 / mgamma_bin1) * 100 << "\\%) \\\\" << std::endl;

        // Print \(\geq 1\) good \(\gamma\)
        double egamma_bin3 = hist_egamma->GetBinContent(3);
        double mgamma_bin3 = hist_mgamma->GetBinContent(3);
        std::cout << "      $\\geq 1$ good $\\gamma$ & " << egamma_bin3 << " (" << (egamma_bin3 / egamma_bin1) * 100 << "\\%) & " 
                  << mgamma_bin3 << " (" << (mgamma_bin3 / mgamma_bin1) * 100 << "\\%) \\\\" << std::endl;

        // Print \(\geq 1\) good \(\ell\)
        double egamma_bin4 = hist_egamma->GetBinContent(4);
        double mgamma_bin4 = hist_mgamma->GetBinContent(4);
        std::cout << "      $\\geq 1$ good $\\ell$ & " << egamma_bin4 << " (" << (egamma_bin4 / egamma_bin1) * 100 << "\\%) & " 
                  << mgamma_bin4 << " (" << (mgamma_bin4 / mgamma_bin1) * 100 << "\\%) \\\\" << std::endl;

        // Print \(\Delta R(\gamma, l) > 0.8\)
        double egamma_bin5 = hist_egamma->GetBinContent(5);
        double mgamma_bin5 = hist_mgamma->GetBinContent(5);
        std::cout << "      $\\Delta R(\\gamma, l) > 0.8$ & " << egamma_bin5 << " (" << (egamma_bin5 / egamma_bin1) * 100 << "\\%) & " 
                  << mgamma_bin5 << " (" << (mgamma_bin5 / mgamma_bin1) * 100 << "\\%) \\\\" << std::endl;

        // Print Z veto
        double egamma_bin6 = hist_egamma->GetBinContent(6);
        double mgamma_bin6 = hist_mgamma->GetBinContent(6);
        std::cout << "      Z veto & " << egamma_bin6 << " (" << (egamma_bin6 / egamma_bin1) * 100 << "\\%) & " 
                  << mgamma_bin6 << " (" << (mgamma_bin6 / mgamma_bin1) * 100 << "\\%) \\\\" << std::endl;

        // Print MET filters
        double egamma_bin7 = hist_egamma->GetBinContent(7);
        double mgamma_bin7 = hist_mgamma->GetBinContent(7);
        std::cout << "      MET filters & " << egamma_bin7 << " (" << (egamma_bin7 / egamma_bin1) * 100 << "\\%) & " 
                  << mgamma_bin7 << " (" << (mgamma_bin7 / mgamma_bin1) * 100 << "\\%) \\\\" << std::endl;

        std::cout << "      \\hline" << std::endl;
        std::cout << "    \\end{tabular}" << std::endl;
        std::cout << "  \\end{center}" << std::endl;
        std::cout << "\\end{table*}" << std::endl;

        // Close the files
        file_egamma->Close();
        file_mgamma->Close();
    }
}


