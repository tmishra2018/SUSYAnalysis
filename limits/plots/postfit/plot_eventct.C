#include <sstream>
#include <iomanip>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <TROOT.h>
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TF3.h"
#include "TH1D.h"
#include "TH1F.h"
#include "THStack.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "Math/QuantFuncMathCore.h"

#include "../../../../cmsstyle/src/cmsstyle.C"
#include "../../../include/analysis_commoncode.h"

void getSignalBinLabels(int i, TString &channel, TString &ptmiss,
                        TString &ht, TString &ptgamma)
{
    if (i < 1 || i > 36) {
        channel = "ERR"; ptmiss = "ERR"; ht = "ERR"; ptgamma = "ERR";
        return;
    }
    channel = (i <= 18) ? "$\\mu\\gamma$" : "$e\\gamma$";

    int idx  = i - 1;          // 0-based index
    int cell = idx % 9;         // position within a 9-bin ptgamma block

    // MET bin (rows of 3)
    if      (cell <= 2) ptmiss = "120--200";
    else if (cell <= 5) ptmiss = "200--400";
    else                ptmiss = "$>400$";

    // HT bin (columns within each MET row)
    if      (cell % 3 == 0) ht = "$<100$";
    else if (cell % 3 == 1) ht = "100--400";
    else                     ht = "$>400$";

    // pT(gamma) bin (two 9-bin blocks per channel)
    if ( (i >= 1 && i <= 9) || (i >= 19 && i <= 27) )
        ptgamma = "$<200$";
    else
        ptgamma = "$>200$";
}


std::string to_string_with_precision(double value, int precision = 1)
{
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;
    std::string result = out.str();
    // Strip trailing ".0" so e.g. "3.0" becomes "3"
    if (result.find(".0") != std::string::npos)
        result = result.substr(0, result.find(".0"));
    return result;
}


void getGarwoodInterval(int N, double &L, double &U, double CL = 0.6827)
{
    double alpha = 1.0 - CL;
    L = (N == 0) ? 0.0 : ROOT::Math::gamma_quantile(alpha / 2.0,       N,   1.0);
    U =                  ROOT::Math::gamma_quantile(1.0 - alpha / 2.0, N+1, 1.0);
}


bool total16 = false;
void readNormPostFit(const TString &postfitDir, const TString &proc, const TString &yearTag,
                      int nbins, std::vector<double> &val, std::vector<double> &err,
                      bool usePostfit = true)
{
    val.assign(nbins, 0.0);
    err.assign(nbins, 0.0);
    TString fname = Form("%s/NormPostFit_%s_%s.txt", postfitDir.Data(), proc.Data(), yearTag.Data());
    std::ifstream in(fname.Data());
    if (!in.is_open()) {
        std::cerr << "ERROR: cannot open " << fname << " -- did you run run_postfit_part1.sh "
                  << "for this year? All bins for this process will read back as 0." << std::endl;
        return;
    }
    double pre_y, pre_e, post_y, post_e;
    int i = 0;
    while (i < nbins && (in >> pre_y >> pre_e >> post_y >> post_e)) {
        val[i] = usePostfit ? post_y : pre_y;
        err[i] = usePostfit ? post_e : pre_e;
        i++;
    }
    if (i != nbins) {
        std::cerr << "WARNING: " << fname << " has " << i << " bin(s), expected " << nbins
                  << " -- remaining bins left at 0." << std::endl;
    }
}


void plot_eventct(int NBIN, TString postfitDir = "", TString yearTag = "2018")
{
    // ----- ROOT style setup -----
    gROOT->SetBatch(true);          // run without displaying graphics
    gStyle->SetOptStat(0);          // suppress statistics box on histograms
    cmsstyle::setCMSStyle();        // CMS style (https://github.com/cms-cat/cmsstyle), replaces setTDRStyle()
    SetSignalConfig();              // set RunYear, preVFP etc. from config
    gStyle->SetErrorX(0);           // no horizontal error bars on data points
    gSystem->Load("../../../lib/libAnaClasses.so");  // load analysis shared library

    std::string whichVFP;
    if      (yearTag == "2016preVFP")  { RunYear = 2016; preVFP = 1; whichVFP = "preVFP";  }
    else if (yearTag == "2016postVFP") { RunYear = 2016; preVFP = 0; whichVFP = "postVFP"; }
    else if (yearTag == "2017")        { RunYear = 2017; preVFP = -1; whichVFP = "";       }
    else if (yearTag == "2018")        { RunYear = 2018; preVFP = -1; whichVFP = "";       }
    else if (yearTag == "Run2")        { RunYear = 678;  preVFP = -1; whichVFP = "";       }
    else {
        std::cerr << "ERROR: unrecognized yearTag '" << yearTag << "' -- expected one of "
                  << "2016preVFP/2016postVFP/2017/2018/Run2. Aborting." << std::endl;
        return;
    }

    // total16 flag: merge both 2016 eras into a single combined dataset
    if (RunYear == 2016 && total16 == true) whichVFP = "";

    cout << "Year: " << RunYear << "\t" << whichVFP << "\t(yearTag=" << yearTag << ")" << endl;

    const int NBINS = NBIN * 2;

    std::vector<double> qcd_val, qcd_err, ele_val, ele_err, jet_val, jet_err,
                         vg_val,  vg_err,  rare_val, rare_err, tot_val, tot_err;
    readNormPostFit(postfitDir, "qcdfakelep",      yearTag, NBINS, qcd_val,  qcd_err);
    readNormPostFit(postfitDir, "elefakepho",      yearTag, NBINS, ele_val,  ele_err);
    readNormPostFit(postfitDir, "jetfakepho",      yearTag, NBINS, jet_val,  jet_err);
    readNormPostFit(postfitDir, "VGamma",          yearTag, NBINS, vg_val,   vg_err);
    readNormPostFit(postfitDir, "ttGam_VVGam",     yearTag, NBINS, rare_val, rare_err);
    readNormPostFit(postfitDir, "total_background",yearTag, NBINS, tot_val,  tot_err);

    // ----- Data (observed) signal-region event count files -----
    TFile *mgfile_sig  = TFile::Open(Form("/uscms_data/d3/tmishra/Output/signalTree_mg_signal_%d%s.root",      RunYear, whichVFP.c_str()));
    TFile *egfile_sig  = TFile::Open(Form("/uscms_data/d3/tmishra/Output/signalTree_egamma_signal_%d%s.root",  RunYear, whichVFP.c_str()));

    TH1D *eg_sig = (TH1D*)egfile_sig->Get("p_eventcount");
    TH1D *mg_sig = (TH1D*)mgfile_sig->Get("p_eventcount");

    TH1D *h_sig = new TH1D("h_sig", "", NBIN*2, 1, NBIN*2+1);
    h_sig->SetBinErrorOption(TH1::kPoisson);
    for (int ibin = 1; ibin <= NBIN; ibin++) {
        h_sig->SetBinContent(ibin,       mg_sig->GetBinContent(ibin));  // mg channel
        h_sig->SetBinContent(ibin+NBIN,  eg_sig->GetBinContent(ibin));  // eg channel
    }

    TFile *file_t5   = TFile::Open(Form("/uscms/home/tmishra/nobackup/signal_trees/signalTree_T5WG_%d%s.root",    RunYear, whichVFP.c_str()));
    TH1D  *p_t5wg_signal_raw   = (TH1D*)file_t5->Get("p_t5wg_1700_1000_nom");
    TH1D  *p_t5wg2_signal_raw  = (TH1D*)file_t5->Get("p_t5wg_1900_100_nom");   // additional T5Wg mass point

    TFile *file_tchi = TFile::Open(Form("/uscms/home/tmishra/nobackup/signal_trees/signalTree_TChiWG_%d%s.root",  RunYear, whichVFP.c_str()));
    TH1D  *p_tchiwg_signal_raw = (TH1D*)file_tchi->Get("p_tchiwg_800_nom");

    TH1D *p_t5wg_signal   = new TH1D("p_t5wg_signal",   "", NBIN*2, 1, NBIN*2+1);
    TH1D *p_t5wg2_signal  = new TH1D("p_t5wg2_signal",  "", NBIN*2, 1, NBIN*2+1);
    TH1D *p_tchiwg_signal = new TH1D("p_tchiwg_signal",  "", NBIN*2, 1, NBIN*2+1);
    for (int ibin = 1; ibin <= NBIN*2; ibin++) {
        p_t5wg_signal->SetBinContent(ibin,   p_t5wg_signal_raw->GetBinContent(ibin));
        p_t5wg_signal->SetBinError(ibin,     p_t5wg_signal_raw->GetBinError(ibin));
        p_t5wg2_signal->SetBinContent(ibin,  p_t5wg2_signal_raw->GetBinContent(ibin));
        p_t5wg2_signal->SetBinError(ibin,    p_t5wg2_signal_raw->GetBinError(ibin));
        p_tchiwg_signal->SetBinContent(ibin, p_tchiwg_signal_raw->GetBinContent(ibin));
        p_tchiwg_signal->SetBinError(ibin,   p_tchiwg_signal_raw->GetBinError(ibin));
    }

    TH1D          *h_bkg           = new TH1D("h_bkg",           "", NBIN*2, 1, NBIN*2+1);
    TH1D          *h_bkg_elefakepho= new TH1D("h_bkg_elefakepho","", NBIN*2, 1, NBIN*2+1);
    TH1D          *h_bkg_jetfakepho= new TH1D("h_bkg_jetfakepho","", NBIN*2, 1, NBIN*2+1);
    TH1D          *h_bkg_qcdfakepho= new TH1D("h_bkg_qcdfakepho","", NBIN*2, 1, NBIN*2+1);
    TH1D          *h_bkg_VGamma    = new TH1D("h_bkg_VGamma",    "", NBIN*2, 1, NBIN*2+1);
    // Individual (non-cumulative) per-component histograms, for the real THStack built via
    // cmsstyle::buildTHStack() below - h_bkg/h_bkg_elefakepho/h_bkg_jetfakepho/h_bkg_VGamma
    // above stay as-is (still used for the tex table and unchanged elsewhere), these are
    // additional, separate histograms just for the stacked plot.
    TH1D          *h_bkg_rare_indiv= new TH1D("h_bkg_rare_indiv","", NBIN*2, 1, NBIN*2+1);
    TH1D          *h_bkg_vg_indiv  = new TH1D("h_bkg_vg_indiv",  "", NBIN*2, 1, NBIN*2+1);
    TH1D          *h_bkg_jet_indiv = new TH1D("h_bkg_jet_indiv", "", NBIN*2, 1, NBIN*2+1);
    TH1D          *h_bkg_ele_indiv = new TH1D("h_bkg_ele_indiv", "", NBIN*2, 1, NBIN*2+1);
    TH1D          *h_bkg_qcd_indiv = new TH1D("h_bkg_qcd_indiv", "", NBIN*2, 1, NBIN*2+1);

    h_bkg_rare_indiv->SetDirectory(0);
    h_bkg_vg_indiv->SetDirectory(0);
    h_bkg_jet_indiv->SetDirectory(0);
    h_bkg_ele_indiv->SetDirectory(0);
    h_bkg_qcd_indiv->SetDirectory(0);
    TGraphErrors  *error_bkg       = new TGraphErrors(NBIN*2);
    TGraphErrors  *error_ratio     = new TGraphErrors(NBIN*2);
    TGraphAsymmErrors *ratio       = new TGraphAsymmErrors();
    TGraphAsymmErrors *g_data      = new TGraphAsymmErrors(NBIN*2);


    ofstream tex;
    tex.open(Form("logs/SignalRegionTable_%d%s.tex", RunYear, whichVFP.c_str()));
    tex << "\\begin{table}[htbp]\n";
    tex << "\\centering\n";
    tex << "\\caption{Observed number of events and predicted SM backgrounds in the search bins for the "
        << RunYear << " " << whichVFP
        << " data-taking period. "
        << "The predicted yields are shown with their systematic uncertainties. "
        << "The total SM background uncertainty is obtained by adding the individual background uncertainties in quadrature.}\n";
    tex << "\\label{tab:sr_yields_" << RunYear << whichVFP << "}\n";
    tex << "\\resizebox{\\textwidth}{!}{%\n";
    tex << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|}\n";
    tex << "\\hline\n";
    tex << "Bin & channel & $p_T^{miss}$ & $H_T$ & $p_T^{\\gamma}$ & "
        << "$e\\to\\gamma$ fakes & jet$\\to\\gamma$ fakes & jet$\\to\\ell$ fakes & "
        << "$W\\gamma$/$Z\\gamma$ & Electroweak & Total bkg & Obs \\\\\n";
    tex << "\\hline\n";

    double bkgContent_elefakepho(0), bkgContent_jetfakepho(0), bkgContent_jetfakelep(0);
    double bkgContent_VGamma(0), bkgContent_rare(0), bkgContent_Total(0), dataEvents(0);
    double bkgContent_elefakepho_mg(0), bkgContent_jetfakepho_mg(0), bkgContent_jetfakelep_mg(0);
    double bkgContent_VGamma_mg(0), bkgContent_rare_mg(0), bkgContent_Total_mg(0), dataEvents_mg(0);

    // Plain-text version of the same per-bin table (all backgrounds + data)
    ofstream txt;
    txt.open(Form("logs/SignalRegionTable_%d%s.txt", RunYear, whichVFP.c_str()));
    txt << "# " << yearTag << " postfit yields per search bin (value +- postfit unc.)\n";
    txt << Form("%-4s %-4s %-9s %-9s %-6s %18s %18s %18s %18s %18s %20s %6s %8s\n",
                "Bin", "Ch", "pTmiss", "HT", "pTg",
                "e->g misid", "j->g misid", "misid lepton", "Wg/Zg", "ttg/VVg",
                "Total bkg", "Obs", "Obs/Pred");
    // Channel subtotals: central values only (postfit uncertainties are correlated
    // across bins, so a quadrature sum here would not be correct)
    auto writeTotalRows = [&](const char *label, double ele, double jet, double qcd,
                              double vg, double rare, double tot, double obs) {
        tex << "\\hline\n"
            << "\\multicolumn{5}{|c|}{" << label << "} & "
            << Form("%.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.0f", ele, jet, qcd, vg, rare, tot, obs)
            << " \\\\\n";
        TString plain = label;
        plain.ReplaceAll("$\\mu\\gamma$", "mug").ReplaceAll("$e\\gamma$", "eg");
        txt << Form("%-35s %18.2f %18.2f %18.2f %18.2f %18.2f %20.2f %6.0f %8.3f\n",
                    plain.Data(), ele, jet, qcd, vg, rare, tot, obs, tot > 0 ? obs / tot : 0.0);
    };

    for (int ibin = 1; ibin <= 2*NBIN; ibin++) {
        int idx = ibin - 1;  // 0-based index into the NormPostFit vectors

        double v_qcd = qcd_val[idx],  e_qcd = qcd_err[idx];
        double v_ele = ele_val[idx],  e_ele = ele_err[idx];
        double v_jet = jet_val[idx],  e_jet = jet_err[idx];
        double v_vg  = vg_val[idx],   e_vg  = vg_err[idx];
        double v_rare= rare_val[idx], e_rare= rare_err[idx];
        double bkgContent = tot_val[idx];
        double totalerror = tot_err[idx];

        double bkgqcdfakelep = v_qcd;
        double bkgelefakepho = bkgqcdfakelep + v_ele;
        double bkgjetfakepho = bkgelefakepho  + v_jet;
        double bkgVGamma     = bkgjetfakepho  + v_vg;

        double nObs = h_sig->GetBinContent(ibin);

        TString channel, ptmiss, ht, ptgamma;
        getSignalBinLabels(ibin, channel, ptmiss, ht, ptgamma);
        tex << ibin << " & " << channel << " & " << ptmiss << " & " << ht << " & " << ptgamma << " & "
            << Form("%.2f $\\pm$ %.2f", v_ele,      e_ele)  << " & "
            << Form("%.2f $\\pm$ %.2f", v_jet,      e_jet)  << " & "
            << Form("%.2f $\\pm$ %.2f", v_qcd,      e_qcd)  << " & "
            << Form("%.2f $\\pm$ %.2f", v_vg,       e_vg)   << " & "
            << Form("%.2f $\\pm$ %.2f", v_rare,     e_rare) << " & "
            << Form("%.2f $\\pm$ %.2f", bkgContent, totalerror) << " & "
            << Form("%.0f", nObs)
            << " \\\\\n";

        // Strip LaTeX markup for the plain-text table
        TString ch_t = (ibin <= NBIN) ? "mug" : "eg";
        TString met_t = ptmiss, ht_t = ht, ptg_t = ptgamma;
        for (TString *s : {&met_t, &ht_t, &ptg_t}) { s->ReplaceAll("$", ""); s->ReplaceAll("--", "-"); }
        txt << Form("%-4d %-4s %-9s %-9s %-6s %8.2f +- %6.2f %8.2f +- %6.2f %8.2f +- %6.2f %8.2f +- %6.2f %8.2f +- %6.2f %10.2f +- %6.2f %6.0f %8.3f\n",
                    ibin, ch_t.Data(), met_t.Data(), ht_t.Data(), ptg_t.Data(),
                    v_ele, e_ele, v_jet, e_jet, v_qcd, e_qcd, v_vg, e_vg, v_rare, e_rare,
                    bkgContent, totalerror, nObs, bkgContent > 0 ? nObs / bkgContent : 0.0);

        // ----- Fill drawing histograms -----
        h_bkg->SetBinContent(ibin, bkgContent);
        h_bkg->SetBinError(ibin, totalerror);
        h_bkg_qcdfakepho->SetBinContent(ibin, bkgqcdfakelep);
        h_bkg_elefakepho->SetBinContent(ibin, bkgelefakepho);
        h_bkg_jetfakepho->SetBinContent(ibin, bkgjetfakepho);
        h_bkg_VGamma->SetBinContent(ibin, bkgVGamma);
        // Individual (non-cumulative) fills for the real THStack:
        h_bkg_rare_indiv->SetBinContent(ibin, v_rare);
        h_bkg_vg_indiv->SetBinContent(ibin,   v_vg);
        h_bkg_jet_indiv->SetBinContent(ibin,  v_jet);
        h_bkg_ele_indiv->SetBinContent(ibin,  v_ele);
        h_bkg_qcd_indiv->SetBinContent(ibin,  v_qcd);
        error_bkg->SetPoint(ibin-1, ibin + 0.5, bkgContent);
        error_bkg->SetPointError(ibin-1, 0.5, totalerror);
        error_ratio->SetPoint(ibin-1, ibin + 0.5, 1.0);
        error_ratio->SetPointError(ibin-1, 0.5, bkgContent > 0 ? totalerror / bkgContent : 0.0);

        // ----- Accumulate per-channel yield totals -----
        if (ibin <= NBIN) {
            dataEvents_mg            += h_sig->GetBinContent(ibin);
            bkgContent_elefakepho_mg += v_ele;
            bkgContent_jetfakepho_mg += v_jet;
            bkgContent_jetfakelep_mg += v_qcd;
            bkgContent_VGamma_mg     += v_vg;
            bkgContent_rare_mg       += v_rare;
            bkgContent_Total_mg      += bkgContent;
        } else {
            dataEvents            += h_sig->GetBinContent(ibin);
            bkgContent_elefakepho += v_ele;
            bkgContent_jetfakepho += v_jet;
            bkgContent_jetfakelep += v_qcd;
            bkgContent_VGamma     += v_vg;
            bkgContent_rare       += v_rare;
            bkgContent_Total      += bkgContent;
        }

        if (ibin == NBIN) {  // mu+gamma subtotal, then divider before e+gamma half
            writeTotalRows("Total $\\mu\\gamma$", bkgContent_elefakepho_mg, bkgContent_jetfakepho_mg,
                           bkgContent_jetfakelep_mg, bkgContent_VGamma_mg, bkgContent_rare_mg,
                           bkgContent_Total_mg, dataEvents_mg);
            tex << "\\hline\n";
        }
    }

    writeTotalRows("Total $e\\gamma$", bkgContent_elefakepho, bkgContent_jetfakepho,
                   bkgContent_jetfakelep, bkgContent_VGamma, bkgContent_rare,
                   bkgContent_Total, dataEvents);
    writeTotalRows("Total", bkgContent_elefakepho + bkgContent_elefakepho_mg,
                   bkgContent_jetfakepho + bkgContent_jetfakepho_mg,
                   bkgContent_jetfakelep + bkgContent_jetfakelep_mg,
                   bkgContent_VGamma + bkgContent_VGamma_mg, bkgContent_rare + bkgContent_rare_mg,
                   bkgContent_Total + bkgContent_Total_mg, dataEvents + dataEvents_mg);

    tex << "\\hline\n\\end{tabular}}\n\\end{table}\n";
    tex.close();
    txt.close();

    {
        std::string srYearTag = std::to_string(RunYear) + whichVFP;
        std::ofstream srout(Form("SR_counts_SUS-24-011_%s.txt", srYearTag.c_str()));
        srout << "\"\"\"\n"
              << "This file contains the observed data, background and uncertainties "
              << "for each signal region. Its purpose is to be imported when computing "
              << "the simplified likelihood\n"
              << "\"\"\"\n\n"
              << "counts = {}\n"
              << "counts[\"cms_sus_24_011\"] = {}\n\n"
              << "#NEED TO ADD THE BIN AS A VALUE\n";

        for (int ibin = 1; ibin <= 2*NBIN; ibin++) {
            int idx = ibin - 1;
            double nb_bin   = tot_val[idx];
            double totalErr = tot_err[idx];
            int    nobs     = (int)h_sig->GetBinContent(ibin);

            srout << "counts[\"cms_sus_24_011\"][" << ibin << "] = {"
                  << "\"n_obs\": " << nobs << ", "
                  << "\"n_bkg\": " << Form("%.4f", nb_bin) << ", "
                  << "\"uncertainty\": " << Form("%.4f", totalErr) << "}\n";
        }
        srout.close();
    }

    cout << "\nPer-bin background breakdown (postfit):\n";
    for (int ibin = 1; ibin <= 2*NBIN; ibin++) {
        int idx = ibin - 1;
        cout << std::fixed << std::setprecision(2)
             << ibin << "\t"
             << ele_val[idx]  << "\t"
             << jet_val[idx]  << "\t"
             << qcd_val[idx]  << "\t"
             << vg_val[idx]   << "\t"
             << rare_val[idx] << "\n";
    }
    cout << "\n--- mu+gamma channel ---\n";
    cout << "e->gamma fake  : " << bkgContent_elefakepho_mg << "\n";
    cout << "jet->gamma fake: " << bkgContent_jetfakepho_mg << "\n";
    cout << "fake lepton    : " << bkgContent_jetfakelep_mg << "\n";
    cout << "V+gamma        : " << bkgContent_VGamma_mg << "\n";
    cout << "Electroweak    : " << bkgContent_rare_mg << "\n";
    cout << "Total SM bkg   : " << bkgContent_Total_mg << "  |  Data: " << dataEvents_mg << "\n";
    if (dataEvents_mg > 0)
        cout << "Discrepancy    : " << 100*(bkgContent_Total_mg-dataEvents_mg)/dataEvents_mg << " %\n";
    cout << "\n--- e+gamma channel ---\n";
    cout << "e->gamma fake  : " << bkgContent_elefakepho << "\n";
    cout << "jet->gamma fake: " << bkgContent_jetfakepho << "\n";
    cout << "fake lepton    : " << bkgContent_jetfakelep << "\n";
    cout << "V+gamma        : " << bkgContent_VGamma << "\n";
    cout << "Electroweak    : " << bkgContent_rare << "\n";
    cout << "Total SM bkg   : " << bkgContent_Total << "  |  Data: " << dataEvents << "\n";
    if (dataEvents > 0)
        cout << "Discrepancy    : " << 100*(bkgContent_Total-dataEvents)/dataEvents << " %\n";

    cmsstyle::TCmsCanvas *can = cmsstyle::cmsCanvas("can", 1, 2*NBIN+1, 0.05, 1000000,
                                               "", "", false, 1200, 800);
    can->SetRightMargin(0.20);
    gStyle->SetPadLeftMargin(0.10);
    can->SetLeftMargin(0.10);
    can->SetRightMargin(0.04);
    can->SetTopMargin(0.08);
    can->SetBottomMargin(0.1);
    can->cd();

    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
    setTopPad(pad1);
    // BUG FIX: pad1 never had an explicit left/right margin, so it inherited
    // whatever the current gStyle default was -- too tight to fit the Y-axis
    // title at the offset used below, which clipped it out of the canvas.
    pad1->SetLeftMargin(0.12);
    pad1->SetRightMargin(0.04);
    pad1->Draw();
    pad1->cd();
    gPad->SetLogy();  // logarithmic y-axis for main plot

    // ----- Axis ranges -----
    h_sig->SetMarkerStyle(20);
    h_sig->GetYaxis()->SetRangeUser(0.05, 1000000);
    h_sig->GetXaxis()->SetRangeUser(1, 2*NBIN+1);
    h_sig->GetYaxis()->SetTitle("Events / bin");
    h_sig->GetYaxis()->SetTitleOffset(1.4);
    h_bkg->GetYaxis()->SetRangeUser(0.05, 1000000);
    h_bkg->GetXaxis()->SetRangeUser(1, 2*NBIN+1);
    h_bkg->GetYaxis()->SetTitle("Events / bin");
    h_bkg->GetYaxis()->SetTitleOffset(1.4);

    // ----- Background stack colours -----

    const int col_ewk = TColor::GetColor("#b9ac70");
    const int col_jetfake = TColor::GetColor("#92dadd");
    const int col_elefake = TColor::GetColor("#bd1f01");
    const int col_qcdfake = TColor::GetColor("#3f90da");
    const int col_vgamma  = TColor::GetColor("#ffa90e ");

    h_bkg->SetLineColor(kBlack);
    h_bkg->Draw("AXIS");

    THStack *hs_bkg = cmsstyle::buildTHStack(std::vector<TH1*>({h_bkg_qcd_indiv, h_bkg_jet_indiv, h_bkg_ele_indiv, h_bkg_vg_indiv,  h_bkg_rare_indiv}),
        {col_qcdfake, col_jetfake, col_elefake, col_vgamma, col_ewk}, "STACK");

    cmsstyle::cmsObjectDraw(hs_bkg, "same");

    const int col_TChiWG    = TColor::GetColor("#832db6");
    const int col_T5Wg  = TColor::GetColor("#e76300"); 
    const int col_T5Wg2 = TColor::GetColor("#FFD700");   // T5Wg (1900,100): gold, dash-dotted

    p_tchiwg_signal->SetLineColor(col_TChiWG);      // brown/olive
    p_tchiwg_signal->SetLineWidth(3);
    p_tchiwg_signal->SetLineStyle(2);
    p_tchiwg_signal->SetFillStyle(0);
    p_tchiwg_signal->Draw("HIST same");
    
    p_t5wg_signal->SetLineColor(col_T5Wg);         // blue
    p_t5wg_signal->SetLineWidth(3);
    p_t5wg_signal->SetFillStyle(0);
    p_t5wg_signal->Draw("HIST same");

    p_t5wg2_signal->SetLineColor(col_T5Wg2);
    p_t5wg2_signal->SetLineWidth(3);
    p_t5wg2_signal->SetLineStyle(4);   // dash-dot, so it is distinct from the solid T5Wg (1700,1000)
    p_t5wg2_signal->SetFillStyle(0);
    p_t5wg2_signal->Draw("HIST same");

    for (int ibin = 1; ibin <= 2*NBIN; ibin++) {
        int N = TMath::Nint(h_sig->GetBinContent(ibin));  // round to nearest integer count
        double L, U;
        getGarwoodInterval(N, L, U);
        double x = h_sig->GetBinCenter(ibin);
        g_data->SetPoint(ibin-1, x, N);
        g_data->SetPointEXlow(ibin-1,  0);
        g_data->SetPointEXhigh(ibin-1, 0);
        g_data->SetPointEYlow(ibin-1,  N - L);
        g_data->SetPointEYhigh(ibin-1, U - N);
    }
    g_data->SetMarkerStyle(20);
    g_data->SetMarkerSize(1.0);
    g_data->SetLineWidth(1);
    g_data->SetLineColor(kBlack);

    // Data points drawn last so they appear on top
    cmsstyle::cmsObjectDraw(g_data, "P");

    // ----- Legend (https://github.com/cms-cat/cmsstyle: cmsLeg/addToLegend) -----
    TLegend *leg = cmsstyle::cmsLeg(0.15, 0.82, 0.92, 0.9);
    leg->SetNColumns(7);
    h_bkg_rare_indiv->SetMarkerSize(0);
    h_bkg_vg_indiv->SetMarkerSize(0);
    h_bkg_jet_indiv->SetMarkerSize(0);
    h_bkg_ele_indiv->SetMarkerSize(0);
    h_bkg_qcd_indiv->SetMarkerSize(0);
    error_bkg->SetMarkerSize(0);
    error_bkg->SetLineWidth(0);
    cmsstyle::addToLegend(leg, {
        {g_data,           {"Data",                        "ep"}},
        {h_bkg_rare_indiv, {"t#bar{t}#gamma / WW#gamma / WZ#gamma", "f"}},
        {h_bkg_vg_indiv,   {"W#gamma / Z#gamma",            "f"}},
        {h_bkg_ele_indiv,  {"e->#gamma misid.",                "f"}},
        {h_bkg_jet_indiv,  {"j->#gamma misid.",                "f"}},
        {h_bkg_qcd_indiv,  {"Misid. lepton",                 "f"}},
        {error_bkg,        {"Unc.",                          "f"}}

        //{p_t5wg_signal,    {"T5Wg (1700, 1000)",             "l"}},
        //{p_t5wg2_signal,   {"T5Wg (1900, 100)",              "l"}},
        //{p_tchiwg_signal,  {"TChiWg (800)",                  "l"}}
    });
    leg->Draw("same");

    // ----- Total uncertainty band -----
    error_bkg->SetFillColor(12);
    error_bkg->SetFillStyle(3345);
    error_bkg->Draw("E2 same");

    // ----- Dividing lines between kinematic regions -----
    TLine *line_mg_pt1  = new TLine(10,0,10,150000);  line_mg_pt1->SetLineStyle(2);
    TLine *line_mg_pt2  = new TLine(19,0,19,150000);  line_mg_pt2->SetLineStyle(2); line_mg_pt2->SetLineWidth(4);
    TLine *line_eg_pt1  = new TLine(28,0,28,150000);  line_eg_pt1->SetLineStyle(2);
    TLine *line_eg_pt2  = new TLine(37,0,37,150000);  line_eg_pt2->SetLineStyle(2);
    TLine *line_mg_met1 = new TLine(4,0,4,10000);      line_mg_met1->SetLineStyle(3);
    TLine *line_mg_met2 = new TLine(7,0,7,10000);      line_mg_met2->SetLineStyle(3);
    TLine *line_mg_met3 = new TLine(13,0,13,10000);    line_mg_met3->SetLineStyle(3);
    TLine *line_mg_met4 = new TLine(16,0,16,10000);    line_mg_met4->SetLineStyle(3);
    TLine *line_eg_met1 = new TLine(22,0,22,10000);    line_eg_met1->SetLineStyle(3);
    TLine *line_eg_met2 = new TLine(25,0,25,10000);    line_eg_met2->SetLineStyle(3);
    TLine *line_eg_met3 = new TLine(31,0,31,10000);    line_eg_met3->SetLineStyle(3);
    TLine *line_eg_met4 = new TLine(34,0,34,10000);    line_eg_met4->SetLineStyle(3);
    line_mg_pt1->Draw("same");  line_mg_pt2->Draw("same");
    line_eg_pt1->Draw("same");  line_eg_pt2->Draw("same");
    line_mg_met1->Draw("same"); line_mg_met2->Draw("same");
    line_mg_met3->Draw("same"); line_mg_met4->Draw("same");
    line_eg_met1->Draw("same"); line_eg_met2->Draw("same");
    line_eg_met3->Draw("same"); line_eg_met4->Draw("same");

    TLatex *latex = new TLatex();
    latex->SetTextSize(0.04);
    latex->DrawLatex(2.5,  50000, "#mu#gamma, p_{T}^{#gamma} < 200 GeV");
    latex->DrawLatex(11.5, 50000, "#mu#gamma, p_{T}^{#gamma} > 200 GeV");
    latex->DrawLatex(20.5, 50000, "e#gamma, p_{T}^{#gamma} < 200 GeV");
    latex->DrawLatex(29.5, 50000, "e#gamma, p_{T}^{#gamma} > 200 GeV");
    latex->DrawLatex(4,  15000, "p^{miss}_{T} (GeV)");
    latex->DrawLatex(13,  15000, "p^{miss}_{T} (GeV)");
    latex->DrawLatex(22,  15000, "p^{miss}_{T} (GeV)");
    latex->DrawLatex(31,  15000, "p^{miss}_{T} (GeV)");
    // pTmiss category labels (lowest category is 120-200 GeV; ARC Fig. 4 comment)
    latex->SetTextSize(0.033);
    for (double x0 : {1.0, 10.0, 19.0, 28.0}) {
        latex->DrawLatex(x0 + 0.62, 5000, "120-200");
        latex->DrawLatex(x0 + 3.25, 5000, "200-400");
        latex->DrawLatex(x0 + 6.4,  5000, "> 400");
    }
    latex->SetTextSize(0.04);

    if      (RunYear == 2016 && preVFP == 1)  CMS_lumi(pad1, 1, 1, 11);  // 2016 preVFP
    else if (RunYear == 2016 && preVFP == 0)  CMS_lumi(pad1, 2, 1, 11);  // 2016 postVFP
    else if (RunYear == 2017)                  CMS_lumi(pad1, 3, 1, 11);  // 2017
    else if (RunYear == 2018)                  CMS_lumi(pad1, 4, 1, 11);  // 2018
    else if (RunYear == 678)                   CMS_lumi(pad1, 6, 1, 11);  // full Run 2 (v2 addition)

    if (TH1 *frame1 = (TH1*)pad1->GetPrimitive("h_bkg")) {
        frame1->GetYaxis()->SetTitle("Events / bin");
        frame1->GetYaxis()->SetTitleFont(43);
        frame1->GetYaxis()->SetTitleSize(26);
        frame1->GetYaxis()->SetTitleOffset(1.6);
        frame1->GetYaxis()->SetLabelFont(43);
        frame1->GetYaxis()->SetLabelSize(20);
    }
    pad1->Modified();
    pad1->RedrawAxis();  // redraw axes on top of histogram fills (v2 addition)
    pad1->Update();
    can->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.35);
    // Horizontal labels need less vertical room than rotated ones did.
    pad2->SetBottomMargin(0.32);
    // BUG FIX: match pad1's left/right margins so the ratio pad's axis lines
    // up vertically with the main pad, and the Y-title has room to render.
    pad2->SetLeftMargin(0.12);
    pad2->SetRightMargin(0.04);
    pad2->Draw();
    pad2->cd();
    TH1F *ratioframe = new TH1F("frame", ";Search bin; Data / Pred.", 36, 1, 37);
    ratioframe->GetXaxis()->SetRangeUser(1, 37);
    ratioframe->GetYaxis()->SetTitleOffset(1.4);
    ratioframe->GetYaxis()->SetRangeUser(0, 2.0);
    std::ostringstream ratiotitle;
    for (unsigned i = 1; i <= 36; i++) {
        ratiotitle.str("");
        ratiotitle << i;
        ratioframe->GetXaxis()->SetBinLabel(i, ratiotitle.str().c_str());
    }
    ratioframe->Draw();
 
    // Flat reference line at ratio = 1
    TLine *flatratio = new TLine(1, 1, 37, 1);
 
    for (int ibin = 1; ibin <= 2*NBIN; ibin++) {
        int    N      = TMath::Nint(h_sig->GetBinContent(ibin));  // round to nearest integer count
        double bkgval = h_bkg->GetBinContent(ibin);
        // Protect against zero background
        if (bkgval <= 0) continue;
        ratio->SetPoint(ibin-1, ibin + 0.5, (double)N / bkgval);
        double L, U;
        getGarwoodInterval(N, L, U);
        ratio->SetPointEYlow( ibin-1, (N - L) / bkgval);
        ratio->SetPointEYhigh(ibin-1, (U - N) / bkgval);
    }
    ratio->SetMarkerStyle(20);
    ratio->SetLineColor(kBlack);
    ratio->GetYaxis()->SetTitleOffset(1.4);
    ratio->GetYaxis()->SetRangeUser(0, 2.0);
    ratio->Draw("P same");
 
    // Background uncertainty band on ratio
    error_ratio->SetFillColor(12);
    error_ratio->SetFillStyle(3345);
    error_ratio->Draw("E2 same");
 
    flatratio->Draw();

    if (TH1 *frame2 = (TH1*)pad2->GetPrimitive("frame")) {
        frame2->GetXaxis()->SetTitle("Search bin");
        frame2->GetYaxis()->SetTitle("Data / Pred.");
        frame2->GetXaxis()->SetTitleFont(43);
        frame2->GetYaxis()->SetTitleFont(43);
        frame2->GetXaxis()->SetTitleSize(26);
        frame2->GetYaxis()->SetTitleSize(26);
        frame2->GetXaxis()->SetTitleOffset(1.2);
        frame2->GetYaxis()->SetTitleOffset(1.6);
        frame2->GetXaxis()->SetLabelFont(43);
        frame2->GetYaxis()->SetLabelFont(43);
        frame2->GetXaxis()->SetLabelSize(13);
        frame2->GetXaxis()->LabelsOption("h");
        frame2->GetYaxis()->SetLabelSize(20);
        frame2->GetYaxis()->SetNdivisions(505);
    }
    pad2->Modified();
    pad2->Update();
 
    can->Update();
    can->SaveAs(Form("/uscms_data/d3/tmishra/Output/SRbins_T5WG_%d%s_postfit_obs.png", RunYear, whichVFP.c_str()));
    can->SaveAs(Form("/uscms_data/d3/tmishra/Output/SRbins_T5WG_%d%s_postfit_obs.pdf", RunYear, whichVFP.c_str()));
 
    TH1D *test_h = new TH1D("test_h", "", 18, 0, 18);
    test_h->SetBinErrorOption(TH1::kPoisson);
    for (unsigned ibin = 1; ibin <= 18; ibin++) test_h->SetBinContent(ibin, ibin-1);
}
