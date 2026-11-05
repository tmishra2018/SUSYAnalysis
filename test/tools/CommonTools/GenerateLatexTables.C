#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include <iomanip>

void printLatexTable(const char* dataFiles[], const char* signal1Files[], const char* signal2Files[], const char* signal3Files[], const char* years[], const char* channel) {
    const char* cutNames[] = {
        "Total",
        "HLT",
        "$\\geq 1$ good $\\gamma$",
        "$\\geq 1$ good $\\ell$",
        "$\\Delta R(\\gamma, l) > 0.8$",
 //       "MET filters", "Z veto" // mg channel
        "Z veto", "MET filters" // eg channel
    };

    std::cout << "\\begin{table*}[tbp]" << std::endl;
    std::cout << "  \\begin{center}" << std::endl;
    std::cout << "    \\caption{Event selection and the number of events after successive cuts for " << channel << " channel.}" << std::endl;
    std::cout << "    \\label{tab:cutflow_" << channel << "}" << std::endl;
    std::cout << "    \\resizebox{\\textwidth}{!}{" << std::endl; 
    std::cout << "    \\begin{tabular}{|l|l|l|l|l|l|}" << std::endl;
    std::cout << "      \\hline" << std::endl;
    std::cout << "      Year & Cut & Data & T5Wg(1500, 1) & T5Wg(1800, 800) & T5Wg(1800, 1600) \\\\" << std::endl;
    std::cout << "      \\hline" << std::endl;

    for (int year_idx = 0; year_idx < 4; ++year_idx) {
        TFile *file_data = TFile::Open(dataFiles[year_idx]);
        TFile *file_signal1 = TFile::Open(signal1Files[year_idx]);
        TFile *file_signal2 = TFile::Open(signal2Files[year_idx]);
        TFile *file_signal3 = TFile::Open(signal3Files[year_idx]);

        TH1D* hist_data = (TH1D*)file_data->Get("p_eventcount");
        TH1D* hist_signal1 = (TH1D*)file_signal1->Get("p_eventcount");
        TH1D* hist_signal2 = (TH1D*)file_signal2->Get("p_eventcount");
        TH1D* hist_signal3 = (TH1D*)file_signal3->Get("p_eventcount");

        double total_data = hist_data->GetBinContent(1);
        double total_signal1 = hist_signal1->GetBinContent(1);
        double total_signal2 = hist_signal2->GetBinContent(1);
        double total_signal3 = hist_signal3->GetBinContent(1);

        for (int bin = 1; bin <= hist_data->GetNbinsX() && bin <= 7; ++bin) {
		 int bin_data = static_cast<int>(hist_data->GetBinContent(bin));
            int bin_signal1 = static_cast<int>(hist_signal1->GetBinContent(bin));
            int bin_signal2 = static_cast<int>(hist_signal2->GetBinContent(bin));
            int bin_signal3 = static_cast<int>(hist_signal3->GetBinContent(bin));

	    if (bin == 1) {
                std::cout << "      " << years[year_idx] << " & " << cutNames[bin - 1] << " & " << bin_data << " (" << std::fixed << std::setprecision(2) << (bin_data / total_data) * 100 << "\\%) & "
                          << bin_signal1 << " (" << std::fixed << std::setprecision(2) << (bin_signal1 / total_signal1) * 100 << "\\%) & "
                          << bin_signal2 << " (" << std::fixed << std::setprecision(2) << (bin_signal2 / total_signal2) * 100 << "\\%) & "
                          << bin_signal3 << " (" << std::fixed << std::setprecision(2) << (bin_signal3 / total_signal3) * 100 << "\\%) \\\\" << std::endl;
            } else {
                std::cout << "      & " << cutNames[bin - 1] << " & " << bin_data << " (" << std::fixed << std::setprecision(2) << (bin_data / total_data) * 100 << "\\%) & "
                          << bin_signal1 << " (" << std::fixed << std::setprecision(2) << (bin_signal1 / total_signal1) * 100 << "\\%) & "
                          << bin_signal2 << " (" << std::fixed << std::setprecision(2) << (bin_signal2 / total_signal2) * 100 << "\\%) & "
                          << bin_signal3 << " (" << std::fixed << std::setprecision(2) << (bin_signal3 / total_signal3) * 100 << "\\%) \\\\" << std::endl;
            }
        }
        std::cout << "      \\hline" << std::endl;

        file_data->Close();
        file_signal1->Close();
        file_signal2->Close();
        file_signal3->Close();
    }

    std::cout << "    \\end{tabular}" << std::endl;
    std::cout << "    \\}" << std::endl; 
    std::cout << "  \\end{center}" << std::endl;
    std::cout << "\\end{table*}" << std::endl;
}


void GenerateLatexTables() {
    const char* years[] = {"2016preVFP", "2016postVFP", "2017", "2018"};

    const char* files_egamma[] = {
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2017.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2018.root"};

    const char* files_mgamma[] = {
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016postVFP.root",
	"/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2017.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2018.root"
    };

    const char* signal_1_egamma[] = {
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1500-1_2016preVFP.root",
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1500-1_2016postVFP.root",
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1500-1_2017.root",
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1500-1_2017.root"};

    const char* signal_2_egamma[] = {
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-800_2016preVFP.root",
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-800_2016postVFP.root",
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-800_2017.root",
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-800_2017.root"};

    const char* signal_3_egamma[] = {
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-1600_2016preVFP.root",
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-1600_2016postVFP.root",
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-1600_2017.root",
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-1600_2017.root"};

    const char* signal_1_mgamma[] = {
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1500-1_2016preVFP.root",
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1500-1_2016postVFP.root",
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1500-1_2017.root",
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1500-1_2017.root"};

    const char* signal_2_mgamma[] = {
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1800-800_2016preVFP.root",
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1800-800_2016postVFP.root",
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1800-800_2017.root",
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1800-800_2017.root"};

    const char* signal_3_mgamma[] = {
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1800-1600_2016preVFP.root",
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1800-1600_2016postVFP.root",
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1800-1600_2017.root",
          "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1800-1600_2017.root"};

 //   printLatexTable(files_egamma, signal_1_egamma, signal_2_egamma, signal_3_egamma, years, "egamma");
    printLatexTable(files_mgamma, signal_1_mgamma, signal_2_mgamma, signal_3_mgamma, years, "mgamma");
}
