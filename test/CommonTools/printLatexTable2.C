#include <iostream>
#include <TFile.h>
#include <TH1D.h>

void printLatexTable2() {
    const char* years[] = {"2016preVFP", "2016postVFP", "2017", "2018"};

    const char* files_egamma[] = {
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2017.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2018.root"};

    const char* signal_egamma[][4] = {
        {"/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1500-1_2016preVFP.root",
         "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1500-1_2016postVFP.root",
         "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1500-1_2017.root",
         "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1500-1_2018.root"},
        {"/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-800_2016preVFP.root",
         "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-800_2016postVFP.root",
         "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-800_2017.root",
         "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-800_2018.root"},
        {"/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-1600_2016preVFP.root",
         "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-1600_2016postVFP.root",
         "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-1600_2017.root",
         "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-1600_2018.root"}};

    const char* files_mgamma[] = {
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016preVFP_Muon20.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016postVFP_Muon20.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2017_Muon20.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2018_Muon20_GT36.root"};

    for (int year_idx = 0; year_idx < 4; ++year_idx) {
        for (int channel = 0; channel < 2; ++channel) {
            const char* channel_name = channel ? "mgamma" : "egamma";
            const char* data_file = channel ? files_mgamma[year_idx] : files_egamma[year_idx];
            const char** signal_files = channel ? signal_egamma[0] : signal_egamma[0];

            TFile *file_data = TFile::Open(data_file);
            TH1D* hist_data = (TH1D*)file_data->Get("p_eventcount");

            std::cout << "\\begin{table*}[tbp]" << std::endl;
            std::cout << "  \\begin{center}" << std::endl;
            std::cout << "    \\caption{Event selection and the number of events after successive cuts for " << years[year_idx] << " (" << channel_name << ").}" << std::endl;
            std::cout << "    \\label{tab:cutflow_" << years[year_idx] << "_" << channel_name << "}" << std::endl;
            std::cout << "    \\begin{tabular}{|l|l|l|l|l|}" << std::endl;
            std::cout << "      \\hline" << std::endl;
            std::cout << "      Cut & Data & T5Wg(1500, 1) & T5Wg(1800, 800) & T5Wg(1800, 1600) \\\\" << std::endl;
            std::cout << "      \\hline" << std::endl;

            for (int bin = 1; bin <= 7; ++bin) {
                double data_bin = hist_data->GetBinContent(bin);
                std::cout << "      Cut " << bin << " & " << data_bin;
                for (int sig = 0; sig < 3; ++sig) {
                    TFile *file_sig = TFile::Open(signal_files[year_idx]);
                    TH1D* hist_sig = (TH1D*)file_sig->Get("p_eventcount");
                    double sig_bin = hist_sig->GetBinContent(bin);
                    std::cout << " & " << sig_bin << " (" << (sig_bin / data_bin) * 100 << "\\%)";
                    file_sig->Close();
                }
                std::cout << " \\\\" << std::endl;
            }

            std::cout << "      \\hline" << std::endl;
            std::cout << "    \\end{tabular}" << std::endl;
            std::cout << "  \\end{center}" << std::endl;
            std::cout << "\\end{table*}" << std::endl;

            file_data->Close();
        }
    }
}
