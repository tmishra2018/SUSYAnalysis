#include <iostream>
#include <TFile.h>
#include <TH1D.h>
bool channel = 1;

void printLatexTable() {
    const char* files_egamma[] = {
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2017.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2018.root"};

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
    
    const char* files_mgamma[] = {
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016preVFP_Muon20.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016postVFP_Muon20.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2017_Muon20.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2018_Muon20.root"};
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



    const char* years[] = {"2016preVFP", "2016postVFP", "2017", "2018"};

    for (int year_idx = 0; year_idx < 4; ++year_idx) {
        std::cout << "\\begin{table*}[tbp]" << std::endl;
        std::cout << "  \\begin{center}" << std::endl;
        std::cout << "    \\caption{Event selection and the number of events after successive cuts for " << years[year_idx] << "." << std::endl;
        std::cout << "    \\label{tab:cutflow_" << years[year_idx] << "}" << std::endl;
        std::cout << "    \\begin{tabular}{|l|l|l|}" << std::endl;
        std::cout << "      \\hline" << std::endl;
        std::cout << "      Cut & Data & T5Wg(1500, 1)  & T5Wg (1800, 1600) & T5Wg (1800, 800) \\\\" << std::endl;
        std::cout << "      \\hline" << std::endl;
        std::cout << "      \\hline" << std::endl;

        TFile *file_egamma = TFile::Open(files_egamma[year_idx]);
        TFile *file_mgamma = TFile::Open(files_mgamma[year_idx]);

        TH1D* hist_egamma = (TH1D*)file_egamma->Get("p_eventcount");
        TH1D* hist_mgamma = (TH1D*)file_mgamma->Get("p_eventcount");

        double egamma_bin1 = hist_egamma->GetBinContent(1);
        double mgamma_bin1 = hist_mgamma->GetBinContent(1);

        std::cout << "      Total & " << egamma_bin1 << " & " << mgamma_bin1 << " \\\\" << std::endl;

        double egamma_bin2 = hist_egamma->GetBinContent(2);
        double mgamma_bin2 = hist_mgamma->GetBinContent(2);
        std::cout << "      HLT & " << egamma_bin2 << " (" << (egamma_bin2 / egamma_bin1) * 100 << "\\%) & " 
                  << mgamma_bin2 << " (" << (mgamma_bin2 / mgamma_bin1) * 100 << "\\%) \\\\" << std::endl;

        double egamma_bin3 = hist_egamma->GetBinContent(3);
        double mgamma_bin3 = hist_mgamma->GetBinContent(3);
        std::cout << "      $\\geq 1$ good $\\gamma$ & " << egamma_bin3 << " (" << (egamma_bin3 / egamma_bin1) * 100 << "\\%) & " 
                  << mgamma_bin3 << " (" << (mgamma_bin3 / mgamma_bin1) * 100 << "\\%) \\\\" << std::endl;

        double egamma_bin4 = hist_egamma->GetBinContent(4);
        double mgamma_bin4 = hist_mgamma->GetBinContent(4);
        std::cout << "      $\\geq 1$ good $\\ell$ & " << egamma_bin4 << " (" << (egamma_bin4 / egamma_bin1) * 100 << "\\%) & " 
                  << mgamma_bin4 << " (" << (mgamma_bin4 / mgamma_bin1) * 100 << "\\%) \\\\" << std::endl;

        double egamma_bin5 = hist_egamma->GetBinContent(5);
        double mgamma_bin5 = hist_mgamma->GetBinContent(5);
        std::cout << "      $\\Delta R(\\gamma, l) > 0.8$ & " << egamma_bin5 << " (" << (egamma_bin5 / egamma_bin1) * 100 << "\\%) & " 
                  << mgamma_bin5 << " (" << (mgamma_bin5 / mgamma_bin1) * 100 << "\\%) \\\\" << std::endl;

        double egamma_bin6 = hist_egamma->GetBinContent(6);
        double mgamma_bin6 = hist_mgamma->GetBinContent(6);
        std::cout << "      Z veto & " << egamma_bin6 << " (" << (egamma_bin6 / egamma_bin1) * 100 << "\\%) & " 
                  << mgamma_bin6 << " (" << (mgamma_bin6 / mgamma_bin1) * 100 << "\\%) \\\\" << std::endl;

        double egamma_bin7 = hist_egamma->GetBinContent(7);
        double mgamma_bin7 = hist_mgamma->GetBinContent(7);
        std::cout << "      MET filters & " << egamma_bin7 << " (" << (egamma_bin7 / egamma_bin1) * 100 << "\\%) & " 
                  << mgamma_bin7 << " (" << (mgamma_bin7 / mgamma_bin1) * 100 << "\\%) \\\\" << std::endl;

        std::cout << "      \\hline" << std::endl;
        std::cout << "    \\end{tabular}" << std::endl;
        std::cout << "  \\end{center}" << std::endl;
        std::cout << "\\end{table*}" << std::endl;

        file_egamma->Close();
        file_mgamma->Close();
    }
}


