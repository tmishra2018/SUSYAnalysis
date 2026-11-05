#include <iostream>
#include <iomanip>
#include <TFile.h>
#include <TH1D.h>

// root -l 
// .L printLatexTable.C+
// printLatexTable(1); // eγ
// printLatexTable(2); // μγ


double pct(double num, double den) {
    if (den <= 0) return 0.0;
    return 100.0 * num / den;
}

void printLatexTable(int channel = 1) {

    const char* years[] = {"2016preVFP", "2016postVFP", "2017", "2018"};

    /* =======================
       File lists by channel
       ======================= */

    const char* data_egamma[] = {
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2017.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2018.root"
    };

    const char* data_mgamma[] = {
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2017.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2018.root"
    };

    const char* sig1500_eg[] = {
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1500-1_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1500-1_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1500-1_2017.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1500-1_2018.root"
    };

    const char* sig1500_mg[] = {
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1500-1_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1500-1_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1500-1_2017.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1500-1_2018.root"
    };

    const char* sig1800_800_eg[] = {
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-800_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-800_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-800_2017.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-800_2018.root"
    };

    const char* sig1800_800_mg[] = {
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1800-800_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1800-800_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1800-800_2017.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1800-800_2018.root"
    };

    const char* sig1800_1600_eg[] = {
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-1600_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-1600_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-1600_2017.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_eg_T5Wg_1800-1600_2018.root"
    };

    const char* sig1800_1600_mg[] = {
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1800-1600_2016preVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1800-1600_2016postVFP.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1800-1600_2017.root",
        "/eos/uscms/store/user/tmishra/eg_mg_treesData/sigtree_mg_T5Wg_1800-1600_2018.root"
    };

    /* =======================
       Select channel
       ======================= */

    const char** dataFiles;
    const char** sig1500;
    const char** sig1800_800;
    const char** sig1800_1600;

    const char* channelName;
    const char* channelLabel;

    if (channel == 1) {
        dataFiles       = data_egamma;
        sig1500         = sig1500_eg;
        sig1800_800     = sig1800_800_eg;
        sig1800_1600    = sig1800_1600_eg;
        channelName     = "$e\\gamma$";
        channelLabel    = "egamma";
    } else {
        dataFiles       = data_mgamma;
        sig1500         = sig1500_mg;
        sig1800_800     = sig1800_800_mg;
        sig1800_1600    = sig1800_1600_mg;
        channelName     = "$\\mu\\gamma$";
        channelLabel    = "mgamma";
    }

    /* =======================
       LaTeX header
       ======================= */

    std::cout << "\\begin{table*}[tbp]\n"
              << "  \\begin{center}\n"
              << "    \\caption{Event selection and the number of events after successive cuts for "
              << channelName << " channel.}\n"
              << "    \\label{tab:cutflow_" << channelLabel << "}\n"
              << "\\resizebox{\\textwidth}{!}{\n"
              << "    \\begin{tabular}{|l|l|l|l|l|l|}\n"
              << "      \\hline\n"
              << "      Year & Cut & Data & T5Wg(1500, 1) & T5Wg(1800, 800) & T5Wg(1800, 1600) \\\\\n"
              << "      \\hline\n";

    std::cout << std::fixed << std::setprecision(2);

    const char* cuts[] = {
        "Total",
        "HLT",
        "$\\geq 1$ good $\\gamma$",
        "$\\geq 1$ good $\\ell$",
        "$\\geq 1$ good $\\ell$ and $\\geq 1$ good $\\gamma$",
        "$\\Delta R(\\gamma, l) > 0.8$",
        "Z veto",
        "MET filters"
    };

    for (int y = 0; y < 4; ++y) {

        TFile *fD = TFile::Open(dataFiles[y]);
        TFile *f1 = TFile::Open(sig1500[y]);
        TFile *f2 = TFile::Open(sig1800_800[y]);
        TFile *f3 = TFile::Open(sig1800_1600[y]);

        TH1D *hD = (TH1D*)fD->Get("p_eventcount");
        TH1D *h1 = (TH1D*)f1->Get("p_eventcount");
        TH1D *h2 = (TH1D*)f2->Get("p_eventcount");
        TH1D *h3 = (TH1D*)f3->Get("p_eventcount");

        double D0 = hD->GetBinContent(1);
        double S10 = h1->GetBinContent(1);
        double S20 = h2->GetBinContent(1);
        double S30 = h3->GetBinContent(1);

        for (int i = 0; i < 8; ++i) {
            int bin = i + 1;

            std::cout << (i == 0 ? years[y] : "")
                      << " & " << cuts[i]
                      << " & " << hD->GetBinContent(bin)  << " (" << pct(hD->GetBinContent(bin),  D0)  << "\\%)"
                      << " & " << h1->GetBinContent(bin)  << " (" << pct(h1->GetBinContent(bin), S10) << "\\%)"
                      << " & " << h2->GetBinContent(bin)  << " (" << pct(h2->GetBinContent(bin), S20) << "\\%)"
                      << " & " << h3->GetBinContent(bin)  << " (" << pct(h3->GetBinContent(bin), S30) << "\\%) \\\\\n";
        }

        std::cout << "      \\hline\n";

        fD->Close(); f1->Close(); f2->Close(); f3->Close();
    }

    std::cout << "    \\end{tabular}\n"
              << "}\n"
              << "  \\end{center}\n"
              << "\\end{table*}\n";
}
