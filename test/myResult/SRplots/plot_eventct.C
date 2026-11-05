// ----- Standard library includes -----
#include <sstream>
#include <iomanip>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

// ----- ROOT includes -----
#include <TROOT.h>
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TF3.h"
#include "TH1D.h"
#include "TH1F.h"
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

// ----- Analysis framework includes -----
// Provides: setTDRStyle(), SetSignalConfig(), CMS_lumi(), setTopPad(),
//           and global variables RunYear, preVFP
#include "../../include/analysis_commoncode.h"


// =============================================================================
// Helper: Build human-readable labels for each of the 36 signal bins.
//   Bins 1-18  -> muon+gamma (mg) channel
//   Bins 19-36 -> electron+gamma (eg) channel
//   Within each 18-bin block:
//     Bins 1-9  : pT(gamma) < 200 GeV
//     Bins 10-18: pT(gamma) > 200 GeV
//     Within each 9-bin block, the MET x HT grid is:
//       cell 0,1,2 -> MET [120,200] ; HT [<100, 100-400, >400]
//       cell 3,4,5 -> MET [200,400] ; HT [<100, 100-400, >400]
//       cell 6,7,8 -> MET [>400]    ; HT [<100, 100-400, >400]
// (Only used for the LaTeX SR table — from v1.)
// =============================================================================
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


// =============================================================================
// Helper: Format a double to a string with given decimal precision,
//         stripping trailing ".0" for clean LaTeX output.
// =============================================================================
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


// =============================================================================
// Helper: CMS-recommended Garwood (Poisson) confidence interval for an
// observed integer count, as described in
//   https://twiki.cern.ch/twiki/bin/view/CMS/PoissonErrorBars
// L, U are returned as bounds on the TRUE mean given N observed (i.e. the
// usual Garwood/Clopper-Pearson-type interval built from the gamma
// quantile function). The caller subtracts N to get the (asymmetric)
// bin error sizes: errLow = N - L, errHigh = U - N.
//
// This single helper is used for BOTH the upper-panel data points and the
// ratio-panel data points below, so the two panels are guaranteed to use
// an identical interval definition (fixes the prefit/postfit inconsistency
// noted in review: error bars on the same observed points should not
// differ between panels/plots).
//
// Reference sanity-check values (CL = 0.6827, i.e. 1-sigma equivalent):
//   N=0 -> [0.00, 1.84]   N=1 -> [0.17, 2.30]   N=2 -> [0.71, 3.30]
// =============================================================================
void getGarwoodInterval(int N, double &L, double &U, double CL = 0.6827)
{
    double alpha = 1.0 - CL;
    L = (N == 0) ? 0.0 : ROOT::Math::gamma_quantile(alpha / 2.0,       N,   1.0);
    U =                  ROOT::Math::gamma_quantile(1.0 - alpha / 2.0, N+1, 1.0);
}


// =============================================================================
// Global flag (from v1):
//   total16 = true  -> treat 2016preVFP + 2016postVFP as a single 2016 dataset
//                       (whichVFP forced to "")
//   total16 = false -> use preVFP / postVFP split as set by analysis_commoncode
// =============================================================================
bool total16 = false;


// =============================================================================
// Main function
//   NBIN : number of signal bins per lepton channel (typically 18)
//          Total histogram range will be [0, 2*NBIN]:
//            bins 1..NBIN       -> muon+gamma (mg)
//            bins NBIN+1..2*NBIN -> electron+gamma (eg)
// =============================================================================
void plot_eventct(int NBIN)
{
    // ----- ROOT style setup -----
    gROOT->SetBatch(true);          // run without displaying graphics
    gStyle->SetOptStat(0);          // suppress statistics box on histograms
    setTDRStyle();                  // CMS TDR style
    SetSignalConfig();              // set RunYear, preVFP etc. from config
    gStyle->SetErrorX(0);           // no horizontal error bars on data points
    gSystem->Load("../../lib/libAnaClasses.so");  // load analysis shared library

    // -------------------------------------------------------------------------
    // Determine the VFP suffix for 2016 data file naming.
    // For other years the suffix is empty.
    // -------------------------------------------------------------------------
    std::string whichVFP;
    if      (RunYear == 2016 && preVFP == 1) whichVFP = "preVFP";
    else if (RunYear == 2016 && preVFP == 0) whichVFP = "postVFP";
    else                                     whichVFP = "";

    // total16 flag: merge both 2016 eras into a single combined dataset
    if (RunYear == 2016 && total16 == true) whichVFP = "";

    cout << "Year: " << RunYear << "\t" << whichVFP << endl;

    // =========================================================================
    // Open background input ROOT files
    //   mg* -> muon+gamma stream
    //   eg* -> electron+gamma stream
    //   Background categories:
    //     ele  : e -> gamma fake background
    //     jet  : jet -> gamma fake background
    //     qcd  : QCD fake lepton background
    //     VG   : V+gamma (Wgamma / Zgamma) background
    //     rare : rare electroweak backgrounds (ttgamma, WWgamma, WZgamma, ...)
    // =========================================================================
    TFile *mgfile_ele  = TFile::Open(Form("/uscms_data/d3/tmishra/Output/signalTree_mg_eleBkg_%d%s.root",    RunYear, whichVFP.c_str()));
    TFile *mgfile_jet  = TFile::Open(Form("/uscms_data/d3/tmishra/Output/signalTree_mg_jetbkg_%d%s.root",    RunYear, whichVFP.c_str()));
    TFile *mgfile_qcd  = TFile::Open(Form("/uscms_data/d3/tmishra/Output/signalTree_mg_qcd_%d%s.root",       RunYear, whichVFP.c_str()));
    TFile *mgfile_VG   = TFile::Open(Form("/uscms_data/d3/tmishra/Output/signalTree_mg_VGBkg_%d%s.root",     RunYear, whichVFP.c_str()));
    TFile *mgfile_rare = TFile::Open(Form("/uscms_data/d3/tmishra/Output/signalTree_mg_rareBkg_%d%s.root",   RunYear, whichVFP.c_str()));
    TFile *egfile_ele  = TFile::Open(Form("/uscms_data/d3/tmishra/Output/signalTree_egamma_eleBkg_%d%s.root",  RunYear, whichVFP.c_str()));
    TFile *egfile_jet  = TFile::Open(Form("/uscms_data/d3/tmishra/Output/signalTree_egamma_jetbkg_%d%s.root",  RunYear, whichVFP.c_str()));
    TFile *egfile_qcd  = TFile::Open(Form("/uscms_data/d3/tmishra/Output/signalTree_egamma_qcd_%d%s.root",     RunYear, whichVFP.c_str()));
    TFile *egfile_VG   = TFile::Open(Form("/uscms_data/d3/tmishra/Output/signalTree_egamma_VGBkg_%d%s.root",   RunYear, whichVFP.c_str()));
    TFile *egfile_rare = TFile::Open(Form("/uscms_data/d3/tmishra/Output/signalTree_egamma_rareBkg_%d%s.root", RunYear, whichVFP.c_str()));

    // ----- Data (observed) signal-region event count files -----
    TFile *mgfile_sig  = TFile::Open(Form("/uscms_data/d3/tmishra/Output/signalTree_mg_signal_%d%s.root",      RunYear, whichVFP.c_str()));
    TFile *egfile_sig  = TFile::Open(Form("/uscms_data/d3/tmishra/Output/signalTree_egamma_signal_%d%s.root",  RunYear, whichVFP.c_str()));

    // =========================================================================
    // Read data (observed) event count histograms and fill combined histogram
    //   h_sig: bins 1..NBIN filled from mg, bins NBIN+1..2*NBIN from eg
    // =========================================================================
    TH1D *eg_sig = (TH1D*)egfile_sig->Get("p_eventcount");
    TH1D *mg_sig = (TH1D*)mgfile_sig->Get("p_eventcount");

    TH1D *h_sig = new TH1D("h_sig", "", NBIN*2, 1, NBIN*2+1);
    // NOTE: kPoisson here only affects h_sig's own (unused-for-drawing)
    // GetBinError*(); the actual displayed error bars on the data points
    // (both panels) are computed explicitly via getGarwoodInterval() below,
    // to guarantee consistency between the upper panel and the ratio panel.
    h_sig->SetBinErrorOption(TH1::kPoisson);
    for (int ibin = 1; ibin <= NBIN; ibin++) {
        h_sig->SetBinContent(ibin,       mg_sig->GetBinContent(ibin));  // mg channel
        h_sig->SetBinContent(ibin+NBIN,  eg_sig->GetBinContent(ibin));  // eg channel
    }

    // =========================================================================
    // Open SUSY signal model files (v2 addition)
    //   T5Wg  : gluino pair -> wino + photon, benchmark (m_gl=1700, m_LSP=1000 GeV)
    //   TChiWg: chargino-neutralino -> W + photon, benchmark (m=800 GeV)
    // =========================================================================
    TFile *file_t5   = TFile::Open(Form("/uscms/home/tmishra/nobackup/signal_trees/signalTree_T5WG_%d%s.root",    RunYear, whichVFP.c_str()));
    TH1D  *p_t5wg_signal_raw   = (TH1D*)file_t5->Get("p_t5wg_1700_1000_nom");

    TFile *file_tchi = TFile::Open(Form("/uscms/home/tmishra/nobackup/signal_trees/signalTree_TChiWG_%d%s.root",  RunYear, whichVFP.c_str()));
    TH1D  *p_tchiwg_signal_raw = (TH1D*)file_tchi->Get("p_tchiwg_800_nom");

    // BUG FIX: signal histograms from file are binned as (NBIN*2, 0, NBIN*2)
    // so their bin i has centre i-0.5.  All background histograms in this
    // script were fixed to (NBIN*2, 1, NBIN*2+1) (centre at i+0.5), shifting
    // the axis by one unit.  Copy signal bin contents into new histograms with
    // the corrected range so they align with the background stack on the pad.
    TH1D *p_t5wg_signal   = new TH1D("p_t5wg_signal",   "", NBIN*2, 1, NBIN*2+1);
    TH1D *p_tchiwg_signal = new TH1D("p_tchiwg_signal",  "", NBIN*2, 1, NBIN*2+1);
    for (int ibin = 1; ibin <= NBIN*2; ibin++) {
        p_t5wg_signal->SetBinContent(ibin,   p_t5wg_signal_raw->GetBinContent(ibin));
        p_t5wg_signal->SetBinError(ibin,     p_t5wg_signal_raw->GetBinError(ibin));
        p_tchiwg_signal->SetBinContent(ibin, p_tchiwg_signal_raw->GetBinContent(ibin));
        p_tchiwg_signal->SetBinError(ibin,   p_tchiwg_signal_raw->GetBinError(ibin));
    }

    // =========================================================================
    // Create combined (mg+eg) background and error histograms / graphs
    //   h_bkg*       : stacked background components for display
    //   error_bkg    : TGraphErrors carrying total (stat+syst) error band on bkg
    //   error_ratio  : same error band, normalised to bkg, for ratio pad
    //   ratio        : TGraphAsymmErrors for observed/bkg with Poisson data errors
    // =========================================================================
    TH1D          *h_bkg           = new TH1D("h_bkg",           "", NBIN*2, 1, NBIN*2+1);
    TH1D          *h_bkg_elefakepho= new TH1D("h_bkg_elefakepho","", NBIN*2, 1, NBIN*2+1);
    TH1D          *h_bkg_jetfakepho= new TH1D("h_bkg_jetfakepho","", NBIN*2, 1, NBIN*2+1);
    TH1D          *h_bkg_qcdfakepho= new TH1D("h_bkg_qcdfakepho","", NBIN*2, 1, NBIN*2+1);
    TH1D          *h_bkg_VGamma    = new TH1D("h_bkg_VGamma",    "", NBIN*2, 1, NBIN*2+1);
    TGraphErrors  *error_bkg       = new TGraphErrors(NBIN*2);
    TGraphErrors  *error_ratio     = new TGraphErrors(NBIN*2);
    TGraphAsymmErrors *ratio       = new TGraphAsymmErrors();
    // g_data: observed data points in the UPPER panel, drawn with explicit
    // Garwood/Poisson asymmetric error bars (CMS recommendation), built via
    // getGarwoodInterval() so they match the ratio-panel data errors exactly.
    TGraphAsymmErrors *g_data      = new TGraphAsymmErrors(NBIN*2);

    // =========================================================================
    // Declare combined (mg+eg) systematic uncertainty histograms.
    // Each background category has a set of systematic variations.
    // The naming convention is:
    //   h_<bkg>_syserr_<source>
    // where <bkg>  in {elefakepho, jetfakepho, qcdfakelep, VGamma, rare}
    //       <source> in {jes, jer, esf, scale, eleshape, jetshape, qcdshape,
    //                    xs, lumi, isr, pu, qcdmuF, qcdmuR, pdf}
    //
    // For data-driven backgrounds (ele/jet/qcd fake) the dominant uncertainty
    // is the shape uncertainty (eleshape / jetshape / qcdshape / scale).
    // For MC backgrounds (VGamma, rare) full sets of JEC/JER/PDF/QCD-scale
    // uncertainties are included.
    // =========================================================================

    // --- e->gamma fake ---
    TH1D *h_elefakepho_norm            = new TH1D("h_elefakepho_norm",            "h_elefakepho_norm", NBIN*2, 1, NBIN*2+1);
    TH1D *h_elefakepho_controlsample   = new TH1D("h_elefakepho_controlsample",   "h_elefakepho_controlsample", NBIN*2, 1, NBIN*2+1);
    TH1D *h_elefakepho_transferfactor  = new TH1D("h_elefakepho_transferfactor",  "h_elefakepho_transferfactor", NBIN*2, 1, NBIN*2+1);
    TH1D *h_elefakepho_syserr_jes      = new TH1D("h_elefakepho_syserr_jes",      "h_elefakepho_syserr_jes", NBIN*2, 1, NBIN*2+1);
    TH1D *h_elefakepho_syserr_jer      = new TH1D("h_elefakepho_syserr_jer",      "h_elefakepho_syserr_jer", NBIN*2, 1, NBIN*2+1);
    TH1D *h_elefakepho_syserr_esf      = new TH1D("h_elefakepho_syserr_esf",      "h_elefakepho_syserr_esf", NBIN*2, 1, NBIN*2+1);
    TH1D *h_elefakepho_syserr_scale    = new TH1D("h_elefakepho_syserr_scale",    "h_elefakepho_syserr_scale", NBIN*2, 1, NBIN*2+1);
    TH1D *h_elefakepho_syserr_eleshape = new TH1D("h_elefakepho_syserr_e_to_pho_syst", "h_elefakepho_syserr_eleshape", NBIN*2, 1, NBIN*2+1);
    TH1D *h_elefakepho_syserr_jetshape = new TH1D("h_elefakepho_syserr_j_to_pho_syst", "h_elefakepho_syserr_jetshape", NBIN*2, 1, NBIN*2+1);
    TH1D *h_elefakepho_syserr_qcdshape = new TH1D("h_elefakepho_syserr_fakelep_shape", "h_elefakepho_syserr_qcdshape", NBIN*2, 1, NBIN*2+1);
    TH1D *h_elefakepho_syserr_xs       = new TH1D("h_elefakepho_syserr_xs",       "h_elefakepho_syserr_xs", NBIN*2, 1, NBIN*2+1);
    TH1D *h_elefakepho_syserr_lumi     = new TH1D("h_elefakepho_syserr_lumi",     "h_elefakepho_syserr_lumi", NBIN*2, 1, NBIN*2+1);
    TH1D *h_elefakepho_syserr_isr      = new TH1D("h_elefakepho_syserr_isr",      "h_elefakepho_syserr_isr", NBIN*2, 1, NBIN*2+1);

    // --- jet->gamma fake ---
    TH1D *h_jetfakepho_norm            = new TH1D("h_jetfakepho_norm",            "h_jetfakepho_norm", NBIN*2, 1, NBIN*2+1);
    TH1D *h_jetfakepho_controlsample   = new TH1D("h_jetfakepho_controlsample",   "h_jetfakepho_controlsample", NBIN*2, 1, NBIN*2+1);
    TH1D *h_jetfakepho_transferfactor  = new TH1D("h_jetfakepho_transferfactor",  "h_jetfakepho_transferfactor", NBIN*2, 1, NBIN*2+1);
    TH1D *h_jetfakepho_syserr_jes      = new TH1D("h_jetfakepho_syserr_jes",      "h_jetfakepho_syserr_jes", NBIN*2, 1, NBIN*2+1);
    TH1D *h_jetfakepho_syserr_jer      = new TH1D("h_jetfakepho_syserr_jer",      "h_jetfakepho_syserr_jer", NBIN*2, 1, NBIN*2+1);
    TH1D *h_jetfakepho_syserr_esf      = new TH1D("h_jetfakepho_syserr_esf",      "h_jetfakepho_syserr_esf", NBIN*2, 1, NBIN*2+1);
    TH1D *h_jetfakepho_syserr_scale    = new TH1D("h_jetfakepho_syserr_scale",    "h_jetfakepho_syserr_scale", NBIN*2, 1, NBIN*2+1);
    TH1D *h_jetfakepho_syserr_eleshape = new TH1D("h_jetfakepho_syserr_e_to_pho_syst", "h_jetfakepho_syserr_eleshape", NBIN*2, 1, NBIN*2+1);
    TH1D *h_jetfakepho_syserr_jetshape = new TH1D("h_jetfakepho_syserr_j_to_pho_syst", "h_jetfakepho_syserr_jetshape", NBIN*2, 1, NBIN*2+1);
    TH1D *h_jetfakepho_syserr_qcdshape = new TH1D("h_jetfakepho_syserr_fakelep_shape", "h_jetfakepho_syserr_qcdshape", NBIN*2, 1, NBIN*2+1);
    TH1D *h_jetfakepho_syserr_xs       = new TH1D("h_jetfakepho_syserr_xs",       "h_jetfakepho_syserr_xs", NBIN*2, 1, NBIN*2+1);
    TH1D *h_jetfakepho_syserr_lumi     = new TH1D("h_jetfakepho_syserr_lumi",     "h_jetfakepho_syserr_lumi", NBIN*2, 1, NBIN*2+1);
    TH1D *h_jetfakepho_syserr_isr      = new TH1D("h_jetfakepho_syserr_isr",      "h_jetfakepho_syserr_isr", NBIN*2, 1, NBIN*2+1);

    // --- QCD fake lepton ---
    TH1D *h_qcdfakelep_norm            = new TH1D("h_qcdfakelep_norm",            "h_qcdfakelep_norm", NBIN*2, 1, NBIN*2+1);
    TH1D *h_qcdfakelep_controlsample   = new TH1D("h_qcdfakelep_controlsample",   "h_qcdfakelep_controlsample", NBIN*2, 1, NBIN*2+1);
    TH1D *h_qcdfakelep_transferfactor  = new TH1D("h_qcdfakelep_transferfactor",  "h_qcdfakelep_transferfactor", NBIN*2, 1, NBIN*2+1);
    TH1D *h_qcdfakelep_syserr_jes      = new TH1D("h_qcdfakelep_syserr_jes",      "h_qcdfakelep_syserr_jes", NBIN*2, 1, NBIN*2+1);
    TH1D *h_qcdfakelep_syserr_jer      = new TH1D("h_qcdfakelep_syserr_jer",      "h_qcdfakelep_syserr_jer", NBIN*2, 1, NBIN*2+1);
    TH1D *h_qcdfakelep_syserr_esf      = new TH1D("h_qcdfakelep_syserr_esf",      "h_qcdfakelep_syserr_esf", NBIN*2, 1, NBIN*2+1);
    TH1D *h_qcdfakelep_syserr_scale    = new TH1D("h_qcdfakelep_syserr_scale",    "h_qcdfakelep_syserr_scale", NBIN*2, 1, NBIN*2+1);
    TH1D *h_qcdfakelep_syserr_eleshape = new TH1D("h_qcdfakelep_syserr_e_to_pho_syst", "h_qcdfakelep_syserr_eleshape", NBIN*2, 1, NBIN*2+1);
    TH1D *h_qcdfakelep_syserr_jetshape = new TH1D("h_qcdfakelep_syserr_j_to_pho_syst", "h_qcdfakelep_syserr_jetshape", NBIN*2, 1, NBIN*2+1);
    TH1D *h_qcdfakelep_syserr_qcdshape = new TH1D("h_qcdfakelep_syserr_fakelep_shape", "h_qcdfakelep_syserr_qcdshape", NBIN*2, 1, NBIN*2+1);
    TH1D *h_qcdfakelep_syserr_xs       = new TH1D("h_qcdfakelep_syserr_xs",       "h_qcdfakelep_syserr_xs", NBIN*2, 1, NBIN*2+1);
    TH1D *h_qcdfakelep_syserr_lumi     = new TH1D("h_qcdfakelep_syserr_lumi",     "h_qcdfakelep_syserr_lumi", NBIN*2, 1, NBIN*2+1);
    TH1D *h_qcdfakelep_syserr_isr      = new TH1D("h_qcdfakelep_syserr_isr",      "h_qcdfakelep_syserr_isr", NBIN*2, 1, NBIN*2+1);

    // --- V+gamma (Wgamma / Zgamma) MC background ---
    // Additional uncertainties: pileup (pu), QCD scale (qcdmuF), PDF
    TH1D *h_VGamma_norm            = new TH1D("h_VGamma_norm",            "h_VGamma_norm", NBIN*2, 1, NBIN*2+1);
    TH1D *h_VGamma_syserr_jes      = new TH1D("h_VGamma_syserr_jes",      "h_VGamma_syserr_jes", NBIN*2, 1, NBIN*2+1);
    TH1D *h_VGamma_syserr_jer      = new TH1D("h_VGamma_syserr_jer",      "h_VGamma_syserr_jer", NBIN*2, 1, NBIN*2+1);
    TH1D *h_VGamma_syserr_esf      = new TH1D("h_VGamma_syserr_esf",      "h_VGamma_syserr_esf", NBIN*2, 1, NBIN*2+1);
    TH1D *h_VGamma_syserr_scale    = new TH1D("h_VGamma_syserr_scale",    "h_VGamma_syserr_scale", NBIN*2, 1, NBIN*2+1);
    TH1D *h_VGamma_syserr_eleshape = new TH1D("h_VGamma_syserr_e_to_pho_syst", "h_VGamma_syserr_eleshape", NBIN*2, 1, NBIN*2+1);
    TH1D *h_VGamma_syserr_jetshape = new TH1D("h_VGamma_syserr_j_to_pho_syst", "h_VGamma_syserr_jetshape", NBIN*2, 1, NBIN*2+1);
    TH1D *h_VGamma_syserr_qcdshape = new TH1D("h_VGamma_syserr_fakelep_shape", "h_VGamma_syserr_qcdshape", NBIN*2, 1, NBIN*2+1);
    TH1D *h_VGamma_syserr_xs       = new TH1D("h_VGamma_syserr_xs",       "h_VGamma_syserr_xs", NBIN*2, 1, NBIN*2+1);
    TH1D *h_VGamma_syserr_lumi     = new TH1D("h_VGamma_syserr_lumi",     "h_VGamma_syserr_lumi", NBIN*2, 1, NBIN*2+1);
    TH1D *h_VGamma_syserr_isr      = new TH1D("h_VGamma_syserr_isr",      "h_VGamma_syserr_isr", NBIN*2, 1, NBIN*2+1);
    TH1D *h_VGamma_syserr_pu       = new TH1D("h_VGamma_syserr_pu",       "h_VGamma_syserr_pu", NBIN*2, 1, NBIN*2+1);   // pileup (v1 only)
    TH1D *h_VGamma_syserr_qcdmuF   = new TH1D("h_VGamma_syserr_qcdmuF",   "h_VGamma_syserr_qcdmuF", NBIN*2, 1, NBIN*2+1);
    TH1D *h_VGamma_syserr_pdf      = new TH1D("h_VGamma_syserr_pdf",      "h_VGamma_syserr_pdf", NBIN*2, 1, NBIN*2+1);

    // --- Rare electroweak backgrounds (ttgamma, WWgamma, WZgamma, ...) ---
    // Additional uncertainties: pileup (pu), QCD scale (qcdmuR, qcdmuF), PDF
    TH1D *h_rare_norm            = new TH1D("h_rare_norm",            "h_rare_norm", NBIN*2, 1, NBIN*2+1);
    TH1D *h_rare_syserr_jes      = new TH1D("h_rare_syserr_jes",      "h_rare_syserr_jes", NBIN*2, 1, NBIN*2+1);
    TH1D *h_rare_syserr_jer      = new TH1D("h_rare_syserr_jer",      "h_rare_syserr_jer", NBIN*2, 1, NBIN*2+1);
    TH1D *h_rare_syserr_esf      = new TH1D("h_rare_syserr_esf",      "h_rare_syserr_esf", NBIN*2, 1, NBIN*2+1);
    TH1D *h_rare_syserr_scale    = new TH1D("h_rare_syserr_scale",    "h_rare_syserr_scale", NBIN*2, 1, NBIN*2+1);
    TH1D *h_rare_syserr_eleshape = new TH1D("h_rare_syserr_e_to_pho_syst", "h_rare_syserr_eleshape", NBIN*2, 1, NBIN*2+1);
    TH1D *h_rare_syserr_jetshape = new TH1D("h_rare_syserr_j_to_pho_syst", "h_rare_syserr_jetshape", NBIN*2, 1, NBIN*2+1);
    TH1D *h_rare_syserr_qcdshape = new TH1D("h_rare_syserr_fakelep_shape", "h_rare_syserr_qcdshape", NBIN*2, 1, NBIN*2+1);
    TH1D *h_rare_syserr_xs       = new TH1D("h_rare_syserr_xs",       "h_rare_syserr_xs", NBIN*2, 1, NBIN*2+1);
    TH1D *h_rare_syserr_lumi     = new TH1D("h_rare_syserr_lumi",     "h_rare_syserr_lumi", NBIN*2, 1, NBIN*2+1);
    TH1D *h_rare_syserr_isr      = new TH1D("h_rare_syserr_isr",      "h_rare_syserr_isr", NBIN*2, 1, NBIN*2+1);
    TH1D *h_rare_syserr_pu       = new TH1D("h_rare_syserr_pu",       "h_rare_syserr_pu", NBIN*2, 1, NBIN*2+1);
    TH1D *h_rare_syserr_qcdmuR   = new TH1D("h_rare_syserr_qcdmuR",   "h_rare_syserr_qcdmuR", NBIN*2, 1, NBIN*2+1);
    TH1D *h_rare_syserr_qcdmuF   = new TH1D("h_rare_syserr_qcdmuF",   "h_rare_syserr_qcdmuF", NBIN*2, 1, NBIN*2+1);
    TH1D *h_rare_syserr_pdf      = new TH1D("h_rare_syserr_pdf",      "h_rare_syserr_pdf", NBIN*2, 1, NBIN*2+1);

    // =========================================================================
    // Read per-channel systematic histograms from input files
    //
    // eg channel histograms (from egfile_*)
    // =========================================================================

    // e->gamma fake, eg channel
    TH1D *eg_elefakepho_norm            = (TH1D*)egfile_ele->Get("eg_elefakepho_norm");
    TH1D *eg_elefakepho_controlsample   = (TH1D*)egfile_ele->Get("eg_elefakepho_controlsample");
    TH1D *eg_elefakepho_transferfactor  = (TH1D*)egfile_ele->Get("eg_elefakepho_transferfactor");
    TH1D *eg_elefakepho_syserr_jes      = (TH1D*)egfile_ele->Get("eg_elefakepho_syserr_jes");
    TH1D *eg_elefakepho_syserr_jer      = (TH1D*)egfile_ele->Get("eg_elefakepho_syserr_jer");
    TH1D *eg_elefakepho_syserr_esf      = (TH1D*)egfile_ele->Get("eg_elefakepho_syserr_esf");
    TH1D *eg_elefakepho_syserr_scale    = (TH1D*)egfile_ele->Get("eg_elefakepho_syserr_scale");
    TH1D *eg_elefakepho_syserr_eleshape = (TH1D*)egfile_ele->Get("eg_elefakepho_syserr_e_to_pho");
    TH1D *eg_elefakepho_syserr_jetshape = (TH1D*)egfile_ele->Get("eg_elefakepho_syserr_j_to_pho");
    TH1D *eg_elefakepho_syserr_qcdshape = (TH1D*)egfile_ele->Get("eg_elefakepho_syserr_j_to_lep");
    TH1D *eg_elefakepho_syserr_xs       = (TH1D*)egfile_ele->Get("eg_elefakepho_syserr_xs");
    TH1D *eg_elefakepho_syserr_lumi     = (TH1D*)egfile_ele->Get("eg_elefakepho_syserr_lumi");
    TH1D *eg_elefakepho_syserr_isr      = (TH1D*)egfile_ele->Get("eg_elefakepho_syserr_isr");
    TH1D *eg_VGamma_syserr_pu           = (TH1D*)egfile_VG->Get("eg_VGamma_syserr_pu");   // VGamma PU from VG file

    // jet->gamma fake, eg channel
    TH1D *eg_jetfakepho_norm            = (TH1D*)egfile_jet->Get("eg_jetfakepho_norm");
    TH1D *eg_jetfakepho_controlsample   = (TH1D*)egfile_jet->Get("eg_jetfakepho_controlsample");
    TH1D *eg_jetfakepho_transferfactor  = (TH1D*)egfile_jet->Get("eg_jetfakepho_transferfactor");
    TH1D *eg_jetfakepho_syserr_jes      = (TH1D*)egfile_jet->Get("eg_jetfakepho_syserr_jes");
    TH1D *eg_jetfakepho_syserr_jer      = (TH1D*)egfile_jet->Get("eg_jetfakepho_syserr_jer");
    TH1D *eg_jetfakepho_syserr_esf      = (TH1D*)egfile_jet->Get("eg_jetfakepho_syserr_esf");
    TH1D *eg_jetfakepho_syserr_scale    = (TH1D*)egfile_jet->Get("eg_jetfakepho_syserr_scale");
    TH1D *eg_jetfakepho_syserr_eleshape = (TH1D*)egfile_jet->Get("eg_jetfakepho_syserr_e_to_pho");
    TH1D *eg_jetfakepho_syserr_jetshape = (TH1D*)egfile_jet->Get("eg_jetfakepho_syserr_j_to_pho");
    TH1D *eg_jetfakepho_syserr_qcdshape = (TH1D*)egfile_jet->Get("eg_jetfakepho_syserr_j_to_lep");
    TH1D *eg_jetfakepho_syserr_xs       = (TH1D*)egfile_jet->Get("eg_jetfakepho_syserr_xs");
    TH1D *eg_jetfakepho_syserr_lumi     = (TH1D*)egfile_jet->Get("eg_jetfakepho_syserr_lumi");
    TH1D *eg_jetfakepho_syserr_isr      = (TH1D*)egfile_jet->Get("eg_jetfakepho_syserr_isr");
    TH1D *eg_rare_syserr_pu             = (TH1D*)egfile_rare->Get("eg_rare_syserr_pu");   // rare PU from rare file

    // QCD fake lepton, eg channel
    TH1D *eg_qcdfakelep_norm            = (TH1D*)egfile_qcd->Get("eg_qcdfakelep_norm");
    TH1D *eg_qcdfakelep_controlsample   = (TH1D*)egfile_qcd->Get("eg_qcdfakelep_controlsample");
    TH1D *eg_qcdfakelep_transferfactor  = (TH1D*)egfile_qcd->Get("eg_qcdfakelep_transferfactor");
    TH1D *eg_qcdfakelep_syserr_jes      = (TH1D*)egfile_qcd->Get("eg_qcdfakelep_syserr_jes");
    TH1D *eg_qcdfakelep_syserr_jer      = (TH1D*)egfile_qcd->Get("eg_qcdfakelep_syserr_jer");
    TH1D *eg_qcdfakelep_syserr_esf      = (TH1D*)egfile_qcd->Get("eg_qcdfakelep_syserr_esf");
    TH1D *eg_qcdfakelep_syserr_scale    = (TH1D*)egfile_qcd->Get("eg_qcdfakelep_syserr_scale");
    TH1D *eg_qcdfakelep_syserr_eleshape = (TH1D*)egfile_qcd->Get("eg_qcdfakelep_syserr_e_to_pho");
    TH1D *eg_qcdfakelep_syserr_jetshape = (TH1D*)egfile_qcd->Get("eg_qcdfakelep_syserr_j_to_pho");
    TH1D *eg_qcdfakelep_syserr_qcdshape = (TH1D*)egfile_qcd->Get("eg_qcdfakelep_syserr_j_to_lep");
    TH1D *eg_qcdfakelep_syserr_xs       = (TH1D*)egfile_qcd->Get("eg_qcdfakelep_syserr_xs");
    TH1D *eg_qcdfakelep_syserr_lumi     = (TH1D*)egfile_qcd->Get("eg_qcdfakelep_syserr_lumi");
    TH1D *eg_qcdfakelep_syserr_isr      = (TH1D*)egfile_qcd->Get("eg_qcdfakelep_syserr_isr");

    // V+gamma, eg channel
    TH1D *eg_VGamma_norm            = (TH1D*)egfile_VG->Get("eg_VGamma_norm");
    TH1D *eg_VGamma_syserr_jes      = (TH1D*)egfile_VG->Get("eg_VGamma_syserr_jes");
    TH1D *eg_VGamma_syserr_jer      = (TH1D*)egfile_VG->Get("eg_VGamma_syserr_jer");
    TH1D *eg_VGamma_syserr_esf      = (TH1D*)egfile_VG->Get("eg_VGamma_syserr_esf");
    TH1D *eg_VGamma_syserr_scale    = (TH1D*)egfile_VG->Get("eg_VGamma_syserr_scale");
    TH1D *eg_VGamma_syserr_eleshape = (TH1D*)egfile_VG->Get("eg_VGamma_syserr_eleshape");
    TH1D *eg_VGamma_syserr_jetshape = (TH1D*)egfile_VG->Get("eg_VGamma_syserr_jetshape");
    TH1D *eg_VGamma_syserr_qcdshape = (TH1D*)egfile_VG->Get("eg_VGamma_syserr_qcdshape");
    TH1D *eg_VGamma_syserr_xs       = (TH1D*)egfile_VG->Get("eg_VGamma_syserr_xs");
    TH1D *eg_VGamma_syserr_lumi     = (TH1D*)egfile_VG->Get("eg_VGamma_syserr_lumi");
    TH1D *eg_VGamma_syserr_isr      = (TH1D*)egfile_VG->Get("eg_VGamma_syserr_isr");
    TH1D *eg_VGamma_syserr_pdf      = (TH1D*)egfile_VG->Get("eg_VGamma_syserr_pdf");
    TH1D *eg_VGamma_syserr_qcdmuF   = (TH1D*)egfile_VG->Get("eg_VGamma_syserr_qcdmuF");

    // Rare EWK, eg channel
    TH1D *eg_rare_norm            = (TH1D*)egfile_rare->Get("eg_rare_norm");
    TH1D *eg_rare_syserr_jes      = (TH1D*)egfile_rare->Get("eg_rare_syserr_jes");
    TH1D *eg_rare_syserr_jer      = (TH1D*)egfile_rare->Get("eg_rare_syserr_jer");
    TH1D *eg_rare_syserr_esf      = (TH1D*)egfile_rare->Get("eg_rare_syserr_esf");
    TH1D *eg_rare_syserr_scale    = (TH1D*)egfile_rare->Get("eg_rare_syserr_scale");
    TH1D *eg_rare_syserr_eleshape = (TH1D*)egfile_rare->Get("eg_rare_syserr_eleshape");
    TH1D *eg_rare_syserr_jetshape = (TH1D*)egfile_rare->Get("eg_rare_syserr_jetshape");
    TH1D *eg_rare_syserr_qcdshape = (TH1D*)egfile_rare->Get("eg_rare_syserr_qcdshape");
    TH1D *eg_rare_syserr_xs       = (TH1D*)egfile_rare->Get("eg_rare_syserr_xs");
    TH1D *eg_rare_syserr_lumi     = (TH1D*)egfile_rare->Get("eg_rare_syserr_lumi");
    TH1D *eg_rare_syserr_isr      = (TH1D*)egfile_rare->Get("eg_rare_syserr_isr");
    TH1D *eg_rare_syserr_pdf      = (TH1D*)egfile_rare->Get("eg_rare_syserr_pdf");
    TH1D *eg_rare_syserr_qcdmuF   = (TH1D*)egfile_rare->Get("eg_rare_syserr_qcdmuF");
    TH1D *eg_rare_syserr_qcdmuR   = (TH1D*)egfile_rare->Get("eg_rare_syserr_qcdmuR");

    // =========================================================================
    // mg channel histograms (from mgfile_*)
    // =========================================================================

    // e->gamma fake, mg channel
    TH1D *mg_elefakepho_norm            = (TH1D*)mgfile_ele->Get("mg_elefakepho_norm");
    TH1D *mg_elefakepho_controlsample   = (TH1D*)mgfile_ele->Get("mg_elefakepho_controlsample");
    TH1D *mg_elefakepho_transferfactor  = (TH1D*)mgfile_ele->Get("mg_elefakepho_transferfactor");
    TH1D *mg_elefakepho_syserr_jes      = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_jes");
    TH1D *mg_elefakepho_syserr_jer      = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_jer");
    TH1D *mg_elefakepho_syserr_esf      = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_esf");
    TH1D *mg_elefakepho_syserr_scale    = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_scale");
    TH1D *mg_elefakepho_syserr_eleshape = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_e_to_pho");
    TH1D *mg_elefakepho_syserr_jetshape = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_j_to_pho");
    TH1D *mg_elefakepho_syserr_qcdshape = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_j_to_lep");
    TH1D *mg_elefakepho_syserr_xs       = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_xs");
    TH1D *mg_elefakepho_syserr_lumi     = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_lumi");
    TH1D *mg_elefakepho_syserr_isr      = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_isr");

    // jet->gamma fake, mg channel
    TH1D *mg_jetfakepho_norm            = (TH1D*)mgfile_jet->Get("mg_jetfakepho_norm");
    TH1D *mg_jetfakepho_controlsample   = (TH1D*)mgfile_jet->Get("mg_jetfakepho_controlsample");
    TH1D *mg_jetfakepho_transferfactor  = (TH1D*)mgfile_jet->Get("mg_jetfakepho_transferfactor");
    TH1D *mg_jetfakepho_syserr_jes      = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_jes");
    TH1D *mg_jetfakepho_syserr_jer      = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_jer");
    TH1D *mg_jetfakepho_syserr_esf      = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_esf");
    TH1D *mg_jetfakepho_syserr_scale    = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_scale");
    TH1D *mg_jetfakepho_syserr_eleshape = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_e_to_pho");
    TH1D *mg_jetfakepho_syserr_jetshape = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_j_to_pho");
    TH1D *mg_jetfakepho_syserr_qcdshape = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_j_to_lep");
    TH1D *mg_jetfakepho_syserr_xs       = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_xs");
    TH1D *mg_jetfakepho_syserr_lumi     = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_lumi");
    TH1D *mg_jetfakepho_syserr_isr      = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_isr");

    // QCD fake lepton, mg channel
    TH1D *mg_qcdfakelep_norm            = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_norm");
    TH1D *mg_qcdfakelep_controlsample   = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_controlsample");
    TH1D *mg_qcdfakelep_transferfactor  = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_transferfactor");
    TH1D *mg_qcdfakelep_syserr_jes      = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_jes");
    TH1D *mg_qcdfakelep_syserr_jer      = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_jer");
    TH1D *mg_qcdfakelep_syserr_esf      = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_esf");
    TH1D *mg_qcdfakelep_syserr_scale    = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_scale");
    TH1D *mg_qcdfakelep_syserr_eleshape = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_e_to_pho");
    TH1D *mg_qcdfakelep_syserr_jetshape = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_j_to_pho");
    TH1D *mg_qcdfakelep_syserr_qcdshape = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_j_to_lep");
    TH1D *mg_qcdfakelep_syserr_xs       = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_xs");
    TH1D *mg_qcdfakelep_syserr_lumi     = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_lumi");
    TH1D *mg_qcdfakelep_syserr_isr      = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_isr");

    // V+gamma, mg channel
    TH1D *mg_VGamma_norm            = (TH1D*)mgfile_VG->Get("mg_VGamma_norm");
    TH1D *mg_VGamma_syserr_jes      = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_jes");
    TH1D *mg_VGamma_syserr_jer      = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_jer");
    TH1D *mg_VGamma_syserr_esf      = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_esf");
    TH1D *mg_VGamma_syserr_scale    = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_scale");
    TH1D *mg_VGamma_syserr_eleshape = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_eleshape");
    TH1D *mg_VGamma_syserr_jetshape = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_jetshape");
    TH1D *mg_VGamma_syserr_qcdshape = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_qcdshape");
    TH1D *mg_VGamma_syserr_xs       = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_xs");
    TH1D *mg_VGamma_syserr_lumi     = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_lumi");
    TH1D *mg_VGamma_syserr_isr      = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_isr");
    TH1D *mg_VGamma_syserr_pu       = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_pu");
    TH1D *mg_VGamma_syserr_pdf      = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_pdf");
    TH1D *mg_VGamma_syserr_qcdmuF   = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_qcdmuF");

    // Rare EWK, mg channel
    TH1D *mg_rare_norm            = (TH1D*)mgfile_rare->Get("mg_rare_norm");
    TH1D *mg_rare_syserr_jes      = (TH1D*)mgfile_rare->Get("mg_rare_syserr_jes");
    TH1D *mg_rare_syserr_jer      = (TH1D*)mgfile_rare->Get("mg_rare_syserr_jer");
    TH1D *mg_rare_syserr_esf      = (TH1D*)mgfile_rare->Get("mg_rare_syserr_esf");
    TH1D *mg_rare_syserr_scale    = (TH1D*)mgfile_rare->Get("mg_rare_syserr_scale");
    TH1D *mg_rare_syserr_eleshape = (TH1D*)mgfile_rare->Get("mg_rare_syserr_eleshape");
    TH1D *mg_rare_syserr_jetshape = (TH1D*)mgfile_rare->Get("mg_rare_syserr_jetshape");
    TH1D *mg_rare_syserr_qcdshape = (TH1D*)mgfile_rare->Get("mg_rare_syserr_qcdshape");
    TH1D *mg_rare_syserr_xs       = (TH1D*)mgfile_rare->Get("mg_rare_syserr_xs");
    TH1D *mg_rare_syserr_lumi     = (TH1D*)mgfile_rare->Get("mg_rare_syserr_lumi");
    TH1D *mg_rare_syserr_isr      = (TH1D*)mgfile_rare->Get("mg_rare_syserr_isr");
    TH1D *mg_rare_syserr_pu       = (TH1D*)mgfile_rare->Get("mg_rare_syserr_pu");
    TH1D *mg_rare_syserr_pdf      = (TH1D*)mgfile_rare->Get("mg_rare_syserr_pdf");
    TH1D *mg_rare_syserr_qcdmuF   = (TH1D*)mgfile_rare->Get("mg_rare_syserr_qcdmuF");
    TH1D *mg_rare_syserr_qcdmuR   = (TH1D*)mgfile_rare->Get("mg_rare_syserr_qcdmuR");

    // =========================================================================
    // Fill combined (mg+eg) histograms:
    //   bins 1..NBIN       <- mg channel
    //   bins NBIN+1..2*NBIN <- eg channel
    // Also propagate statistical (bin) errors for the normalisation histograms
    // and transfer factors (used later in the table).
    // =========================================================================
    double bkgContent(0), bkgError(0), statError(0), totalerror(0);

    for (int ibin = 1; ibin <= NBIN; ibin++) {

        // ---- mg bins ----
        h_elefakepho_norm->SetBinContent(ibin,            mg_elefakepho_norm->GetBinContent(ibin));
        h_elefakepho_controlsample->SetBinContent(ibin,   mg_elefakepho_controlsample->GetBinContent(ibin));
        h_elefakepho_transferfactor->SetBinContent(ibin,  mg_elefakepho_transferfactor->GetBinContent(ibin));
        h_elefakepho_syserr_jes->SetBinContent(ibin,      mg_elefakepho_syserr_jes->GetBinContent(ibin));
        h_elefakepho_syserr_jer->SetBinContent(ibin,      mg_elefakepho_syserr_jer->GetBinContent(ibin));
        h_elefakepho_syserr_esf->SetBinContent(ibin,      mg_elefakepho_syserr_esf->GetBinContent(ibin));
        h_elefakepho_syserr_scale->SetBinContent(ibin,    mg_elefakepho_syserr_scale->GetBinContent(ibin));
        h_elefakepho_syserr_eleshape->SetBinContent(ibin, mg_elefakepho_syserr_eleshape->GetBinContent(ibin));
        h_elefakepho_syserr_jetshape->SetBinContent(ibin, mg_elefakepho_syserr_jetshape->GetBinContent(ibin));
        h_elefakepho_syserr_qcdshape->SetBinContent(ibin, mg_elefakepho_syserr_qcdshape->GetBinContent(ibin));
        h_elefakepho_syserr_xs->SetBinContent(ibin,       mg_elefakepho_syserr_xs->GetBinContent(ibin));
        h_elefakepho_syserr_lumi->SetBinContent(ibin,     mg_elefakepho_syserr_lumi->GetBinContent(ibin));
        h_elefakepho_syserr_isr->SetBinContent(ibin,      mg_elefakepho_syserr_isr->GetBinContent(ibin));

        h_jetfakepho_norm->SetBinContent(ibin,            mg_jetfakepho_norm->GetBinContent(ibin));
        h_jetfakepho_controlsample->SetBinContent(ibin,   mg_jetfakepho_controlsample->GetBinContent(ibin));
        h_jetfakepho_transferfactor->SetBinContent(ibin,  mg_jetfakepho_transferfactor->GetBinContent(ibin));
        h_jetfakepho_syserr_jes->SetBinContent(ibin,      mg_jetfakepho_syserr_jes->GetBinContent(ibin));
        h_jetfakepho_syserr_jer->SetBinContent(ibin,      mg_jetfakepho_syserr_jer->GetBinContent(ibin));
        h_jetfakepho_syserr_esf->SetBinContent(ibin,      mg_jetfakepho_syserr_esf->GetBinContent(ibin));
        h_jetfakepho_syserr_scale->SetBinContent(ibin,    mg_jetfakepho_syserr_scale->GetBinContent(ibin));
        h_jetfakepho_syserr_eleshape->SetBinContent(ibin, mg_jetfakepho_syserr_eleshape->GetBinContent(ibin));
        h_jetfakepho_syserr_jetshape->SetBinContent(ibin, mg_jetfakepho_syserr_jetshape->GetBinContent(ibin));
        h_jetfakepho_syserr_qcdshape->SetBinContent(ibin, mg_jetfakepho_syserr_qcdshape->GetBinContent(ibin));
        h_jetfakepho_syserr_xs->SetBinContent(ibin,       mg_jetfakepho_syserr_xs->GetBinContent(ibin));
        h_jetfakepho_syserr_lumi->SetBinContent(ibin,     mg_jetfakepho_syserr_lumi->GetBinContent(ibin));
        h_jetfakepho_syserr_isr->SetBinContent(ibin,      mg_jetfakepho_syserr_isr->GetBinContent(ibin));

        h_qcdfakelep_norm->SetBinContent(ibin,            mg_qcdfakelep_norm->GetBinContent(ibin));
        h_qcdfakelep_controlsample->SetBinContent(ibin,   mg_qcdfakelep_controlsample->GetBinContent(ibin));
        h_qcdfakelep_transferfactor->SetBinContent(ibin,  mg_qcdfakelep_transferfactor->GetBinContent(ibin));
        h_qcdfakelep_syserr_jes->SetBinContent(ibin,      mg_qcdfakelep_syserr_jes->GetBinContent(ibin));
        h_qcdfakelep_syserr_jer->SetBinContent(ibin,      mg_qcdfakelep_syserr_jer->GetBinContent(ibin));
        h_qcdfakelep_syserr_esf->SetBinContent(ibin,      mg_qcdfakelep_syserr_esf->GetBinContent(ibin));
        h_qcdfakelep_syserr_scale->SetBinContent(ibin,    mg_qcdfakelep_syserr_scale->GetBinContent(ibin));
        h_qcdfakelep_syserr_eleshape->SetBinContent(ibin, mg_qcdfakelep_syserr_eleshape->GetBinContent(ibin));
        h_qcdfakelep_syserr_jetshape->SetBinContent(ibin, mg_qcdfakelep_syserr_jetshape->GetBinContent(ibin));
        h_qcdfakelep_syserr_qcdshape->SetBinContent(ibin, mg_qcdfakelep_syserr_qcdshape->GetBinContent(ibin));
        h_qcdfakelep_syserr_xs->SetBinContent(ibin,       mg_qcdfakelep_syserr_xs->GetBinContent(ibin));
        h_qcdfakelep_syserr_lumi->SetBinContent(ibin,     mg_qcdfakelep_syserr_lumi->GetBinContent(ibin));
        h_qcdfakelep_syserr_isr->SetBinContent(ibin,      mg_qcdfakelep_syserr_isr->GetBinContent(ibin));

        h_VGamma_norm->SetBinContent(ibin,            mg_VGamma_norm->GetBinContent(ibin));
        h_VGamma_syserr_jes->SetBinContent(ibin,      mg_VGamma_syserr_jes->GetBinContent(ibin));
        h_VGamma_syserr_jer->SetBinContent(ibin,      mg_VGamma_syserr_jer->GetBinContent(ibin));
        h_VGamma_syserr_esf->SetBinContent(ibin,      mg_VGamma_syserr_esf->GetBinContent(ibin));
        h_VGamma_syserr_scale->SetBinContent(ibin,    mg_VGamma_syserr_scale->GetBinContent(ibin));
        h_VGamma_syserr_eleshape->SetBinContent(ibin, mg_VGamma_syserr_eleshape->GetBinContent(ibin));
        h_VGamma_syserr_jetshape->SetBinContent(ibin, mg_VGamma_syserr_jetshape->GetBinContent(ibin));
        h_VGamma_syserr_qcdshape->SetBinContent(ibin, mg_VGamma_syserr_qcdshape->GetBinContent(ibin));
        h_VGamma_syserr_xs->SetBinContent(ibin,       mg_VGamma_syserr_xs->GetBinContent(ibin));
        h_VGamma_syserr_lumi->SetBinContent(ibin,     mg_VGamma_syserr_lumi->GetBinContent(ibin));
        h_VGamma_syserr_isr->SetBinContent(ibin,      mg_VGamma_syserr_isr->GetBinContent(ibin));
        h_VGamma_syserr_pu->SetBinContent(ibin,       mg_VGamma_syserr_pu->GetBinContent(ibin));
        h_VGamma_syserr_qcdmuF->SetBinContent(ibin,   mg_VGamma_syserr_qcdmuF->GetBinContent(ibin));
        h_VGamma_syserr_pdf->SetBinContent(ibin,      mg_VGamma_syserr_pdf->GetBinContent(ibin));

        h_rare_norm->SetBinContent(ibin,            mg_rare_norm->GetBinContent(ibin));
        h_rare_syserr_jes->SetBinContent(ibin,      mg_rare_syserr_jes->GetBinContent(ibin));
        h_rare_syserr_jer->SetBinContent(ibin,      mg_rare_syserr_jer->GetBinContent(ibin));
        h_rare_syserr_esf->SetBinContent(ibin,      mg_rare_syserr_esf->GetBinContent(ibin));
        h_rare_syserr_scale->SetBinContent(ibin,    mg_rare_syserr_scale->GetBinContent(ibin));
        h_rare_syserr_eleshape->SetBinContent(ibin, mg_rare_syserr_eleshape->GetBinContent(ibin));
        h_rare_syserr_jetshape->SetBinContent(ibin, mg_rare_syserr_jetshape->GetBinContent(ibin));
        h_rare_syserr_qcdshape->SetBinContent(ibin, mg_rare_syserr_qcdshape->GetBinContent(ibin));
        h_rare_syserr_xs->SetBinContent(ibin,       mg_rare_syserr_xs->GetBinContent(ibin));
        h_rare_syserr_lumi->SetBinContent(ibin,     mg_rare_syserr_lumi->GetBinContent(ibin));
        h_rare_syserr_isr->SetBinContent(ibin,      mg_rare_syserr_isr->GetBinContent(ibin));
        h_rare_syserr_pu->SetBinContent(ibin,       mg_rare_syserr_pu->GetBinContent(ibin));
        h_rare_syserr_qcdmuR->SetBinContent(ibin,   mg_rare_syserr_qcdmuR->GetBinContent(ibin));
        h_rare_syserr_qcdmuF->SetBinContent(ibin,   mg_rare_syserr_qcdmuF->GetBinContent(ibin));
        h_rare_syserr_pdf->SetBinContent(ibin,      mg_rare_syserr_pdf->GetBinContent(ibin));

        // ---- eg bins (offset by NBIN) ----
        h_elefakepho_norm->SetBinContent(ibin+NBIN,            eg_elefakepho_norm->GetBinContent(ibin));
        h_elefakepho_controlsample->SetBinContent(ibin+NBIN,   eg_elefakepho_controlsample->GetBinContent(ibin));
        h_elefakepho_transferfactor->SetBinContent(ibin+NBIN,  eg_elefakepho_transferfactor->GetBinContent(ibin));
        h_elefakepho_syserr_jes->SetBinContent(ibin+NBIN,      eg_elefakepho_syserr_jes->GetBinContent(ibin));
        h_elefakepho_syserr_jer->SetBinContent(ibin+NBIN,      eg_elefakepho_syserr_jer->GetBinContent(ibin));
        h_elefakepho_syserr_esf->SetBinContent(ibin+NBIN,      eg_elefakepho_syserr_esf->GetBinContent(ibin));
        h_elefakepho_syserr_scale->SetBinContent(ibin+NBIN,    eg_elefakepho_syserr_scale->GetBinContent(ibin));
        h_elefakepho_syserr_eleshape->SetBinContent(ibin+NBIN, eg_elefakepho_syserr_eleshape->GetBinContent(ibin));
        h_elefakepho_syserr_jetshape->SetBinContent(ibin+NBIN, eg_elefakepho_syserr_jetshape->GetBinContent(ibin));
        h_elefakepho_syserr_qcdshape->SetBinContent(ibin+NBIN, eg_elefakepho_syserr_qcdshape->GetBinContent(ibin));
        h_elefakepho_syserr_xs->SetBinContent(ibin+NBIN,       eg_elefakepho_syserr_xs->GetBinContent(ibin));
        h_elefakepho_syserr_lumi->SetBinContent(ibin+NBIN,     eg_elefakepho_syserr_lumi->GetBinContent(ibin));
        h_elefakepho_syserr_isr->SetBinContent(ibin+NBIN,      eg_elefakepho_syserr_isr->GetBinContent(ibin));

        h_jetfakepho_norm->SetBinContent(ibin+NBIN,            eg_jetfakepho_norm->GetBinContent(ibin));
        h_jetfakepho_controlsample->SetBinContent(ibin+NBIN,   eg_jetfakepho_controlsample->GetBinContent(ibin));
        h_jetfakepho_transferfactor->SetBinContent(ibin+NBIN,  eg_jetfakepho_transferfactor->GetBinContent(ibin));
        h_jetfakepho_syserr_jes->SetBinContent(ibin+NBIN,      eg_jetfakepho_syserr_jes->GetBinContent(ibin));
        h_jetfakepho_syserr_jer->SetBinContent(ibin+NBIN,      eg_jetfakepho_syserr_jer->GetBinContent(ibin));
        h_jetfakepho_syserr_esf->SetBinContent(ibin+NBIN,      eg_jetfakepho_syserr_esf->GetBinContent(ibin));
        h_jetfakepho_syserr_scale->SetBinContent(ibin+NBIN,    eg_jetfakepho_syserr_scale->GetBinContent(ibin));
        h_jetfakepho_syserr_eleshape->SetBinContent(ibin+NBIN, eg_jetfakepho_syserr_eleshape->GetBinContent(ibin));
        h_jetfakepho_syserr_jetshape->SetBinContent(ibin+NBIN, eg_jetfakepho_syserr_jetshape->GetBinContent(ibin));
        h_jetfakepho_syserr_qcdshape->SetBinContent(ibin+NBIN, eg_jetfakepho_syserr_qcdshape->GetBinContent(ibin));
        h_jetfakepho_syserr_xs->SetBinContent(ibin+NBIN,       eg_jetfakepho_syserr_xs->GetBinContent(ibin));
        h_jetfakepho_syserr_lumi->SetBinContent(ibin+NBIN,     eg_jetfakepho_syserr_lumi->GetBinContent(ibin));
        h_jetfakepho_syserr_isr->SetBinContent(ibin+NBIN,      eg_jetfakepho_syserr_isr->GetBinContent(ibin));

        h_qcdfakelep_norm->SetBinContent(ibin+NBIN,            eg_qcdfakelep_norm->GetBinContent(ibin));
        h_qcdfakelep_controlsample->SetBinContent(ibin+NBIN,   eg_qcdfakelep_controlsample->GetBinContent(ibin));
        h_qcdfakelep_transferfactor->SetBinContent(ibin+NBIN,  eg_qcdfakelep_transferfactor->GetBinContent(ibin));
        h_qcdfakelep_syserr_jes->SetBinContent(ibin+NBIN,      eg_qcdfakelep_syserr_jes->GetBinContent(ibin));
        h_qcdfakelep_syserr_jer->SetBinContent(ibin+NBIN,      eg_qcdfakelep_syserr_jer->GetBinContent(ibin));
        h_qcdfakelep_syserr_esf->SetBinContent(ibin+NBIN,      eg_qcdfakelep_syserr_esf->GetBinContent(ibin));
        h_qcdfakelep_syserr_scale->SetBinContent(ibin+NBIN,    eg_qcdfakelep_syserr_scale->GetBinContent(ibin));
        h_qcdfakelep_syserr_eleshape->SetBinContent(ibin+NBIN, eg_qcdfakelep_syserr_eleshape->GetBinContent(ibin));
        h_qcdfakelep_syserr_jetshape->SetBinContent(ibin+NBIN, eg_qcdfakelep_syserr_jetshape->GetBinContent(ibin));
        h_qcdfakelep_syserr_qcdshape->SetBinContent(ibin+NBIN, eg_qcdfakelep_syserr_qcdshape->GetBinContent(ibin));
        h_qcdfakelep_syserr_xs->SetBinContent(ibin+NBIN,       eg_qcdfakelep_syserr_xs->GetBinContent(ibin));
        h_qcdfakelep_syserr_lumi->SetBinContent(ibin+NBIN,     eg_qcdfakelep_syserr_lumi->GetBinContent(ibin));
        h_qcdfakelep_syserr_isr->SetBinContent(ibin+NBIN,      eg_qcdfakelep_syserr_isr->GetBinContent(ibin));

        h_VGamma_norm->SetBinContent(ibin+NBIN,            eg_VGamma_norm->GetBinContent(ibin));
        h_VGamma_syserr_jes->SetBinContent(ibin+NBIN,      eg_VGamma_syserr_jes->GetBinContent(ibin));
        h_VGamma_syserr_jer->SetBinContent(ibin+NBIN,      eg_VGamma_syserr_jer->GetBinContent(ibin));
        h_VGamma_syserr_esf->SetBinContent(ibin+NBIN,      eg_VGamma_syserr_esf->GetBinContent(ibin));
        h_VGamma_syserr_scale->SetBinContent(ibin+NBIN,    eg_VGamma_syserr_scale->GetBinContent(ibin));
        h_VGamma_syserr_eleshape->SetBinContent(ibin+NBIN, eg_VGamma_syserr_eleshape->GetBinContent(ibin));
        h_VGamma_syserr_jetshape->SetBinContent(ibin+NBIN, eg_VGamma_syserr_jetshape->GetBinContent(ibin));
        h_VGamma_syserr_qcdshape->SetBinContent(ibin+NBIN, eg_VGamma_syserr_qcdshape->GetBinContent(ibin));
        h_VGamma_syserr_xs->SetBinContent(ibin+NBIN,       eg_VGamma_syserr_xs->GetBinContent(ibin));
        h_VGamma_syserr_lumi->SetBinContent(ibin+NBIN,     eg_VGamma_syserr_lumi->GetBinContent(ibin));
        h_VGamma_syserr_isr->SetBinContent(ibin+NBIN,      eg_VGamma_syserr_isr->GetBinContent(ibin));
        h_VGamma_syserr_pu->SetBinContent(ibin+NBIN,       eg_VGamma_syserr_pu->GetBinContent(ibin));
        h_VGamma_syserr_qcdmuF->SetBinContent(ibin+NBIN,   eg_VGamma_syserr_qcdmuF->GetBinContent(ibin));
        h_VGamma_syserr_pdf->SetBinContent(ibin+NBIN,      eg_VGamma_syserr_pdf->GetBinContent(ibin));

        h_rare_norm->SetBinContent(ibin+NBIN,            eg_rare_norm->GetBinContent(ibin));
        h_rare_syserr_jes->SetBinContent(ibin+NBIN,      eg_rare_syserr_jes->GetBinContent(ibin));
        h_rare_syserr_jer->SetBinContent(ibin+NBIN,      eg_rare_syserr_jer->GetBinContent(ibin));
        h_rare_syserr_esf->SetBinContent(ibin+NBIN,      eg_rare_syserr_esf->GetBinContent(ibin));
        h_rare_syserr_scale->SetBinContent(ibin+NBIN,    eg_rare_syserr_scale->GetBinContent(ibin));
        h_rare_syserr_eleshape->SetBinContent(ibin+NBIN, eg_rare_syserr_eleshape->GetBinContent(ibin));
        h_rare_syserr_jetshape->SetBinContent(ibin+NBIN, eg_rare_syserr_jetshape->GetBinContent(ibin));
        h_rare_syserr_qcdshape->SetBinContent(ibin+NBIN, eg_rare_syserr_qcdshape->GetBinContent(ibin));
        h_rare_syserr_xs->SetBinContent(ibin+NBIN,       eg_rare_syserr_xs->GetBinContent(ibin));
        h_rare_syserr_lumi->SetBinContent(ibin+NBIN,     eg_rare_syserr_lumi->GetBinContent(ibin));
        h_rare_syserr_isr->SetBinContent(ibin+NBIN,      eg_rare_syserr_isr->GetBinContent(ibin));
        h_rare_syserr_pu->SetBinContent(ibin+NBIN,       eg_rare_syserr_pu->GetBinContent(ibin));
        h_rare_syserr_qcdmuR->SetBinContent(ibin+NBIN,   eg_rare_syserr_qcdmuR->GetBinContent(ibin));
        h_rare_syserr_qcdmuF->SetBinContent(ibin+NBIN,   eg_rare_syserr_qcdmuF->GetBinContent(ibin));
        h_rare_syserr_pdf->SetBinContent(ibin+NBIN,      eg_rare_syserr_pdf->GetBinContent(ibin));

        // ---- Statistical errors on normalisation and transfer factor histograms ----
        h_elefakepho_norm->SetBinError(ibin,           mg_elefakepho_norm->GetBinError(ibin));
        h_jetfakepho_norm->SetBinError(ibin,           mg_jetfakepho_norm->GetBinError(ibin));
        h_qcdfakelep_norm->SetBinError(ibin,           mg_qcdfakelep_norm->GetBinError(ibin));
        h_VGamma_norm->SetBinError(ibin,               mg_VGamma_norm->GetBinError(ibin));
        h_rare_norm->SetBinError(ibin,                 mg_rare_norm->GetBinError(ibin));
        h_elefakepho_transferfactor->SetBinError(ibin, mg_elefakepho_transferfactor->GetBinError(ibin));
        h_jetfakepho_transferfactor->SetBinError(ibin, mg_jetfakepho_transferfactor->GetBinError(ibin));
        h_qcdfakelep_transferfactor->SetBinError(ibin, mg_qcdfakelep_transferfactor->GetBinError(ibin));

        h_elefakepho_norm->SetBinError(ibin+NBIN,           eg_elefakepho_norm->GetBinError(ibin));
        h_jetfakepho_norm->SetBinError(ibin+NBIN,           eg_jetfakepho_norm->GetBinError(ibin));
        h_qcdfakelep_norm->SetBinError(ibin+NBIN,           eg_qcdfakelep_norm->GetBinError(ibin));
        h_VGamma_norm->SetBinError(ibin+NBIN,               eg_VGamma_norm->GetBinError(ibin));
        h_rare_norm->SetBinError(ibin+NBIN,                 eg_rare_norm->GetBinError(ibin));
        h_elefakepho_transferfactor->SetBinError(ibin+NBIN, eg_elefakepho_transferfactor->GetBinError(ibin));
        h_jetfakepho_transferfactor->SetBinError(ibin+NBIN, eg_jetfakepho_transferfactor->GetBinError(ibin));
        h_qcdfakelep_transferfactor->SetBinError(ibin+NBIN, eg_qcdfakelep_transferfactor->GetBinError(ibin));
    }

    // =========================================================================
    // Pileup reweighting range (for printing/reference, from v1)
    // Not used in any calculation here but kept for documentation purposes.
    // =========================================================================
    float PUmin = 0.0;
    float PUmax;
    if      (RunYear == 2018) PUmax = 26.4;
    else if (RunYear == 2017) PUmax = 19.3;
    else                      PUmax = 6.6;   // 2016 pre/post VFP

    // =========================================================================
    // Write LaTeX signal region table (from v1)
    // Output: logs/SignalRegionTable_<year><VFP>.tex
    // Columns: Bin | channel | MET | HT | pT(gamma) |
    //          e->gamma fake | jet->gamma fake | fake lepton |
    //          W/Z+gamma | EWK | Total bkg | Obs
    // =========================================================================
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

    // Running totals for per-channel yield summary (printed to stdout)
    double bkgContent_elefakepho(0), bkgContent_jetfakepho(0), bkgContent_jetfakelep(0);
    double bkgContent_VGamma(0), bkgContent_rare(0), bkgContent_Total(0), dataEvents(0);
    double bkgContent_elefakepho_mg(0), bkgContent_jetfakepho_mg(0), bkgContent_jetfakelep_mg(0);
    double bkgContent_VGamma_mg(0), bkgContent_rare_mg(0), bkgContent_Total_mg(0), dataEvents_mg(0);
    double VGamma_syserr(0), rare_syserr(0);

    // Map to accumulate systematic uncertainty ranges (min, max) across all bins
    // Used for the summary table printed to stdout
    std::map<std::string, std::pair<double,double>> sysRanges;

    // =========================================================================
    // Main bin loop: compute total background + uncertainties per bin,
    // fill drawing histograms and graphs, write LaTeX table rows.
    // =========================================================================
    for (int ibin = 1; ibin <= 2*NBIN; ibin++) {
        bkgContent = 0; bkgError = 0; statError = 0; totalerror = 0;
        VGamma_syserr = 0; rare_syserr = 0;

        // Total background: sum of all components in this bin
        bkgContent = h_elefakepho_norm->GetBinContent(ibin)
                   + h_jetfakepho_norm->GetBinContent(ibin)
                   + h_qcdfakelep_norm->GetBinContent(ibin)
                   + h_VGamma_norm->GetBinContent(ibin)
                   + h_rare_norm->GetBinContent(ibin);

        // Stacked cumulative sums for drawing (each h_bkg_* is drawn on top of previous)
        double bkgqcdfakelep = h_qcdfakelep_norm->GetBinContent(ibin);
        double bkgelefakepho = bkgqcdfakelep + h_elefakepho_norm->GetBinContent(ibin);
        double bkgjetfakepho = bkgelefakepho  + h_jetfakepho_norm->GetBinContent(ibin);
        double bkgVGamma     = bkgjetfakepho  + h_VGamma_norm->GetBinContent(ibin);

        // ----- Systematic uncertainty (added in quadrature) -----
        // Only positive (upward) systematic shifts are included (convention: 0 = no shift)
        if (h_elefakepho_syserr_eleshape->GetBinContent(ibin) > 0)
            bkgError += pow(h_elefakepho_syserr_eleshape->GetBinContent(ibin), 2);
        if (h_jetfakepho_syserr_jetshape->GetBinContent(ibin) > 0)
            bkgError += pow(h_jetfakepho_syserr_jetshape->GetBinContent(ibin), 2);
        if (h_qcdfakelep_syserr_scale->GetBinContent(ibin) > 0)
            bkgError += pow(h_qcdfakelep_syserr_scale->GetBinContent(ibin), 2);

        if (h_VGamma_syserr_jes->GetBinContent(ibin)   > 0) bkgError += pow(h_VGamma_syserr_jes->GetBinContent(ibin), 2);
        if (h_VGamma_syserr_jer->GetBinContent(ibin)   > 0) bkgError += pow(h_VGamma_syserr_jer->GetBinContent(ibin), 2);
        if (h_VGamma_syserr_esf->GetBinContent(ibin)   > 0) bkgError += pow(h_VGamma_syserr_esf->GetBinContent(ibin), 2);
        if (h_VGamma_syserr_scale->GetBinContent(ibin) > 0) bkgError += pow(h_VGamma_syserr_scale->GetBinContent(ibin), 2);
        if (h_VGamma_syserr_qcdmuF->GetBinContent(ibin)> 0) bkgError += pow(h_VGamma_syserr_qcdmuF->GetBinContent(ibin), 2);
        if (h_VGamma_syserr_pdf->GetBinContent(ibin)   > 0) bkgError += pow(h_VGamma_syserr_pdf->GetBinContent(ibin), 2);

        if (h_rare_syserr_jes->GetBinContent(ibin)    > 0) bkgError += pow(h_rare_syserr_jes->GetBinContent(ibin), 2);
        if (h_rare_syserr_jer->GetBinContent(ibin)    > 0) bkgError += pow(h_rare_syserr_jer->GetBinContent(ibin), 2);
        if (h_rare_syserr_esf->GetBinContent(ibin)    > 0) bkgError += pow(h_rare_syserr_esf->GetBinContent(ibin), 2);
        if (h_rare_syserr_xs->GetBinContent(ibin)     > 0) bkgError += pow(h_rare_syserr_xs->GetBinContent(ibin), 2);
        if (h_rare_syserr_lumi->GetBinContent(ibin)   > 0) bkgError += pow(h_rare_syserr_lumi->GetBinContent(ibin), 2);
        if (h_rare_syserr_qcdmuR->GetBinContent(ibin) > 0) bkgError += pow(h_rare_syserr_qcdmuR->GetBinContent(ibin), 2);
        if (h_rare_syserr_qcdmuF->GetBinContent(ibin) > 0) bkgError += pow(h_rare_syserr_qcdmuF->GetBinContent(ibin), 2);
        if (h_rare_syserr_pdf->GetBinContent(ibin)    > 0) bkgError += pow(h_rare_syserr_pdf->GetBinContent(ibin), 2);

        // Individual component systematic totals (for LaTeX table)
        VGamma_syserr = sqrt(pow(h_VGamma_syserr_jes->GetBinContent(ibin), 2)
                           + pow(h_VGamma_syserr_jer->GetBinContent(ibin), 2)
                           + pow(h_VGamma_syserr_esf->GetBinContent(ibin), 2)
                           + pow(h_VGamma_syserr_scale->GetBinContent(ibin), 2)
                           + pow(h_VGamma_syserr_qcdmuF->GetBinContent(ibin), 2)
                           + pow(h_VGamma_syserr_pdf->GetBinContent(ibin), 2));
        rare_syserr   = sqrt(pow(h_rare_syserr_jes->GetBinContent(ibin), 2)
                           + pow(h_rare_syserr_jer->GetBinContent(ibin), 2)
                           + pow(h_rare_syserr_esf->GetBinContent(ibin), 2)
                           + pow(h_rare_syserr_xs->GetBinContent(ibin), 2)
                           + pow(h_rare_syserr_lumi->GetBinContent(ibin), 2)
                           + pow(h_rare_syserr_qcdmuR->GetBinContent(ibin), 2)
                           + pow(h_rare_syserr_qcdmuF->GetBinContent(ibin), 2)
                           + pow(h_rare_syserr_pdf->GetBinContent(ibin), 2));

        // ----- Statistical uncertainty -----
        if (h_elefakepho_norm->GetBinError(ibin) > 0) statError += pow(h_elefakepho_norm->GetBinError(ibin), 2);
        if (h_jetfakepho_norm->GetBinError(ibin) > 0) statError += pow(h_jetfakepho_norm->GetBinError(ibin), 2);
        if (h_qcdfakelep_norm->GetBinError(ibin) > 0) statError += pow(h_qcdfakelep_norm->GetBinError(ibin), 2);
        if (h_VGamma_norm->GetBinError(ibin)     > 0) statError += pow(h_VGamma_norm->GetBinError(ibin), 2);
        if (h_rare_norm->GetBinError(ibin)       > 0) statError += pow(h_rare_norm->GetBinError(ibin), 2);

        totalerror = sqrt(bkgError + statError);  // combined stat+syst in quadrature
        bkgError   = sqrt(bkgError);              // systematic only (for table)

        // ----- Sanity check: flag any systematic > 100% of the background norm -----
        std::vector<double> sysChecks = {
            (h_elefakepho_norm->GetBinContent(ibin) != 0) ? 100*h_elefakepho_syserr_eleshape->GetBinContent(ibin)/h_elefakepho_norm->GetBinContent(ibin) : 0,
            (h_jetfakepho_norm->GetBinContent(ibin) != 0) ? 100*h_jetfakepho_syserr_jetshape->GetBinContent(ibin)/h_jetfakepho_norm->GetBinContent(ibin) : 0,
            (h_qcdfakelep_norm->GetBinContent(ibin) != 0) ? 100*h_qcdfakelep_syserr_scale->GetBinContent(ibin)/h_qcdfakelep_norm->GetBinContent(ibin) : 0,
            (h_VGamma_norm->GetBinContent(ibin)     != 0) ? 100*h_VGamma_syserr_esf->GetBinContent(ibin)/h_VGamma_norm->GetBinContent(ibin) : 0,
            (h_VGamma_norm->GetBinContent(ibin)     != 0) ? 100*h_VGamma_syserr_jes->GetBinContent(ibin)/h_VGamma_norm->GetBinContent(ibin) : 0,
            (h_rare_norm->GetBinContent(ibin)       != 0) ? 100*h_rare_syserr_jes->GetBinContent(ibin)/h_rare_norm->GetBinContent(ibin) : 0,
        };
        for (float val : sysChecks) {
            if (val > 100.1) {
                std::cerr << "WARNING: >100% uncertainty in bin " << ibin << " (value=" << val << "%)\n";
                break;
            }
        }

        // ----- Accumulate systematic ranges for summary table -----
        // Track min/max percentage uncertainty per source across all bins (v2 approach)
        auto updateRange = [&](const std::string &key, double norm, double sysval) {
            if (norm <= 0) return;
            double pct = round(100.0 * sysval / norm * 10) / 10.0;
            if (pct <= 0) return;
            if (sysRanges.find(key) == sysRanges.end())
                sysRanges[key] = {pct, pct};
            else {
                sysRanges[key].first  = std::min(sysRanges[key].first,  pct);
                sysRanges[key].second = std::max(sysRanges[key].second, pct);
            }
        };
        updateRange("Electron Fake Photon Shape", h_elefakepho_norm->GetBinContent(ibin), h_elefakepho_syserr_eleshape->GetBinContent(ibin));
        updateRange("Jet Fake Photon Shape",      h_jetfakepho_norm->GetBinContent(ibin), h_jetfakepho_syserr_jetshape->GetBinContent(ibin));
        updateRange("QCD Fake Lepton Scale",      h_qcdfakelep_norm->GetBinContent(ibin), h_qcdfakelep_syserr_scale->GetBinContent(ibin));
        updateRange("JEC",      h_VGamma_norm->GetBinContent(ibin), h_VGamma_syserr_jes->GetBinContent(ibin));
        updateRange("JER",      h_VGamma_norm->GetBinContent(ibin), h_VGamma_syserr_jer->GetBinContent(ibin));
        updateRange("ESF",      h_VGamma_norm->GetBinContent(ibin), h_VGamma_syserr_esf->GetBinContent(ibin));
        updateRange("Scale",    h_VGamma_norm->GetBinContent(ibin), h_VGamma_syserr_scale->GetBinContent(ibin));
        updateRange("JEC_rare", h_rare_norm->GetBinContent(ibin),   h_rare_syserr_jes->GetBinContent(ibin));
        updateRange("JER_rare", h_rare_norm->GetBinContent(ibin),   h_rare_syserr_jer->GetBinContent(ibin));
        updateRange("ESF_rare", h_rare_norm->GetBinContent(ibin),   h_rare_syserr_esf->GetBinContent(ibin));
        updateRange("XS",       h_rare_norm->GetBinContent(ibin),   h_rare_syserr_xs->GetBinContent(ibin));
        updateRange("Lumi",     h_rare_norm->GetBinContent(ibin),   h_rare_syserr_lumi->GetBinContent(ibin));
        updateRange("QCDmuF_VGamma", h_VGamma_norm->GetBinContent(ibin), h_VGamma_syserr_qcdmuF->GetBinContent(ibin));
        updateRange("QCDmuR_rare",   h_rare_norm->GetBinContent(ibin),   h_rare_syserr_qcdmuR->GetBinContent(ibin));
        updateRange("QCDmuF_rare",   h_rare_norm->GetBinContent(ibin),   h_rare_syserr_qcdmuF->GetBinContent(ibin));
        updateRange("PDF_VGamma",    h_VGamma_norm->GetBinContent(ibin), h_VGamma_syserr_pdf->GetBinContent(ibin));
        updateRange("PDF_rare",      h_rare_norm->GetBinContent(ibin),   h_rare_syserr_pdf->GetBinContent(ibin));

        // ----- Write LaTeX table row (v1) -----
        TString channel, ptmiss, ht, ptgamma;
        getSignalBinLabels(ibin, channel, ptmiss, ht, ptgamma);
        tex << ibin << " & " << channel << " & " << ptmiss << " & " << ht << " & " << ptgamma << " & "
            << Form("%.2f $\\pm$ %.2f", h_elefakepho_norm->GetBinContent(ibin), max(0.0, h_elefakepho_syserr_eleshape->GetBinContent(ibin))) << " & "
            << Form("%.2f $\\pm$ %.2f", h_jetfakepho_norm->GetBinContent(ibin), max(0.0, h_jetfakepho_syserr_jetshape->GetBinContent(ibin))) << " & "
            << Form("%.2f $\\pm$ %.2f", h_qcdfakelep_norm->GetBinContent(ibin), max(0.0, h_qcdfakelep_syserr_scale->GetBinContent(ibin))) << " & "
            << Form("%.2f $\\pm$ %.2f", h_VGamma_norm->GetBinContent(ibin),     max(0.0, VGamma_syserr)) << " & "
            << Form("%.2f $\\pm$ %.2f", h_rare_norm->GetBinContent(ibin),       max(0.0, rare_syserr)) << " & "
            << Form("%.2f $\\pm$ %.2f", bkgContent, max(0.0, bkgError)) << " & "
            << " \\\\\n";
        if (ibin == NBIN) tex << "\\hline\n";  // divider between mg and eg halves

        // ----- Fill drawing histograms -----
        h_bkg->SetBinContent(ibin, bkgContent);
        h_bkg->SetBinError(ibin, totalerror);
        h_bkg_qcdfakepho->SetBinContent(ibin, bkgqcdfakelep);
        h_bkg_elefakepho->SetBinContent(ibin, bkgelefakepho);
        h_bkg_jetfakepho->SetBinContent(ibin, bkgjetfakepho);
        h_bkg_VGamma->SetBinContent(ibin, bkgVGamma);
        // BUG FIX: x-position of graph points must match bin centres of
        // histograms declared as (NBIN*2, 1, NBIN*2+1).
        // Bin ibin now spans [ibin, ibin+1) so its centre is ibin + 0.5.
        error_bkg->SetPoint(ibin-1, ibin + 0.5, bkgContent);
        error_bkg->SetPointError(ibin-1, 0.5, totalerror);
        error_ratio->SetPoint(ibin-1, ibin + 0.5, 1.0);
        error_ratio->SetPointError(ibin-1, 0.5, totalerror / bkgContent);

        // ----- Accumulate per-channel yield totals -----
        if (ibin <= NBIN) {
            dataEvents_mg        += h_sig->GetBinContent(ibin);
            bkgContent_elefakepho_mg += h_elefakepho_norm->GetBinContent(ibin);
            bkgContent_jetfakepho_mg += h_jetfakepho_norm->GetBinContent(ibin);
            bkgContent_jetfakelep_mg += h_qcdfakelep_norm->GetBinContent(ibin);
            bkgContent_VGamma_mg     += h_VGamma_norm->GetBinContent(ibin);
            bkgContent_rare_mg       += h_rare_norm->GetBinContent(ibin);
            bkgContent_Total_mg      += bkgContent;
        } else {
            dataEvents           += h_sig->GetBinContent(ibin);
            bkgContent_elefakepho += h_elefakepho_norm->GetBinContent(ibin);
            bkgContent_jetfakepho += h_jetfakepho_norm->GetBinContent(ibin);
            bkgContent_jetfakelep += h_qcdfakelep_norm->GetBinContent(ibin);
            bkgContent_VGamma     += h_VGamma_norm->GetBinContent(ibin);
            bkgContent_rare       += h_rare_norm->GetBinContent(ibin);
            bkgContent_Total      += bkgContent;
        }
    }

    tex << "\\hline\n\\end{tabular}}\n\\end{table}\n";
    tex.close();

    // =========================================================================
    // Write SR counts file for simplified likelihood (from v1)
    // Output: SR_counts_SUS-24-011_<year><VFP>.txt
    // =========================================================================
    {
        std::string yearTag = std::to_string(RunYear) + whichVFP;
        std::ofstream srout(Form("SR_counts_SUS-24-011_%s.txt", yearTag.c_str()));
        srout << "\"\"\"\n"
              << "This file contains the observed data, background and uncertainties "
              << "for each signal region. Its purpose is to be imported when computing "
              << "the simplified likelihood\n"
              << "\"\"\"\n\n"
              << "counts = {}\n"
              << "counts[\"cms_sus_24_011\"] = {}\n\n"
              << "#NEED TO ADD THE BIN AS A VALUE\n";

        for (int ibin = 1; ibin <= 2*NBIN; ibin++) {
            double nb_bin = h_elefakepho_norm->GetBinContent(ibin)
                          + h_jetfakepho_norm->GetBinContent(ibin)
                          + h_qcdfakelep_norm->GetBinContent(ibin)
                          + h_VGamma_norm->GetBinContent(ibin)
                          + h_rare_norm->GetBinContent(ibin);

            double bkgErr2 = 0;
            auto addSys = [&](TH1D *h) { if (h->GetBinContent(ibin) > 0) bkgErr2 += pow(h->GetBinContent(ibin), 2); };
            addSys(h_elefakepho_syserr_eleshape);
            addSys(h_jetfakepho_syserr_jetshape);
            addSys(h_qcdfakelep_syserr_scale);
            addSys(h_VGamma_syserr_jes);  addSys(h_VGamma_syserr_jer);  addSys(h_VGamma_syserr_esf);
            addSys(h_VGamma_syserr_scale); addSys(h_VGamma_syserr_qcdmuF); addSys(h_VGamma_syserr_pdf);
            addSys(h_rare_syserr_jes);    addSys(h_rare_syserr_jer);    addSys(h_rare_syserr_esf);
            addSys(h_rare_syserr_xs);     addSys(h_rare_syserr_lumi);
            addSys(h_rare_syserr_qcdmuR); addSys(h_rare_syserr_qcdmuF); addSys(h_rare_syserr_pdf);

            double statErr2 = 0;
            auto addStat = [&](TH1D *h) { if (h->GetBinError(ibin) > 0) statErr2 += pow(h->GetBinError(ibin), 2); };
            addStat(h_elefakepho_norm); addStat(h_jetfakepho_norm);
            addStat(h_qcdfakelep_norm); addStat(h_VGamma_norm); addStat(h_rare_norm);

            double totalErr = sqrt(bkgErr2 + statErr2);
            int    nobs     = (int)h_sig->GetBinContent(ibin);

            srout << "counts[\"cms_sus_24_011\"][" << ibin << "] = {"
                  << "\"n_obs\": " << nobs << ", "
                  << "\"n_bkg\": " << Form("%.4f", nb_bin) << ", "
                  << "\"uncertainty\": " << Form("%.4f", totalErr) << "}\n";
        }
        srout.close();
    }

    // =========================================================================
    // Print per-channel yield summary to stdout (from v1)
    // =========================================================================
    cout << "\nPer-bin background breakdown:\n";
    for (int ibin = 1; ibin <= 2*NBIN; ibin++) {
        cout << std::fixed << std::setprecision(2)
             << ibin << "\t"
             << h_elefakepho_norm->GetBinContent(ibin) << "\t"
             << h_jetfakepho_norm->GetBinContent(ibin) << "\t"
             << h_qcdfakelep_norm->GetBinContent(ibin) << "\t"
             << h_VGamma_norm->GetBinContent(ibin) << "\t"
             << h_rare_norm->GetBinContent(ibin) << "\n";
    }
    cout << "\n--- mu+gamma channel ---\n";
    cout << "e->gamma fake  : " << bkgContent_elefakepho_mg << "\n";
    cout << "jet->gamma fake: " << bkgContent_jetfakepho_mg << "\n";
    cout << "fake lepton    : " << bkgContent_jetfakelep_mg << "\n";
    cout << "V+gamma        : " << bkgContent_VGamma_mg << "\n";
    cout << "Electroweak    : " << bkgContent_rare_mg << "\n";
    cout << "Total SM bkg   : " << bkgContent_Total_mg << "  |  Data: " << dataEvents_mg << "\n";
    cout << "Discrepancy    : " << 100*(bkgContent_Total_mg-dataEvents_mg)/dataEvents_mg << " %\n";
    cout << "\n--- e+gamma channel ---\n";
    cout << "e->gamma fake  : " << bkgContent_elefakepho << "\n";
    cout << "jet->gamma fake: " << bkgContent_jetfakepho << "\n";
    cout << "fake lepton    : " << bkgContent_jetfakelep << "\n";
    cout << "V+gamma        : " << bkgContent_VGamma << "\n";
    cout << "Electroweak    : " << bkgContent_rare << "\n";
    cout << "Total SM bkg   : " << bkgContent_Total << "  |  Data: " << dataEvents << "\n";
    cout << "Discrepancy    : " << 100*(bkgContent_Total-dataEvents)/dataEvents << " %\n";

    // =========================================================================
    // Print LaTeX systematic uncertainty summary table to stdout (v2 style)
    // Rows: uncertainty source; Columns: background component
    // =========================================================================
    std::cout << "\n\\begin{table}[htbp]\n  \\centering\n  \\resizebox{\\textwidth}{!}{%\n"
              << "  \\begin{tabular}{lccccc}\n  \\hline\n"
              << "  Uncertainty source & $e\\to\\gamma$ fake & jet$\\to\\gamma$ fake "
              << "& fake lep & $W/Z\\gamma$ & EWK \\\\\n  \\hline\n";

    auto rangeStr = [&](const std::string &key) -> std::string {
        if (sysRanges.find(key) == sysRanges.end()) return "-";
        return to_string_with_precision(sysRanges[key].first) + "-" + to_string_with_precision(sysRanges[key].second);
    };
    std::map<std::string, std::vector<std::string>> tableData = {
        {"Jet energy scale",        {"-","-","-", rangeStr("JEC"),  rangeStr("JEC_rare")}},
        {"Jet energy resolution",   {"-","-","-", rangeStr("JER"),  rangeStr("JER_rare")}},
        {"ID/trigger ESF",          {"-","-","-", rangeStr("ESF"),  rangeStr("ESF_rare")}},
        {"e-fake-photon shape",     {rangeStr("Electron Fake Photon Shape"),"-","-","-","-"}},
        {"jet-fake-photon shape",   {"-", rangeStr("Jet Fake Photon Shape"),"-","-","-"}},
        {"Normalization scale",     {"-","-", rangeStr("QCD Fake Lepton Scale"), rangeStr("Scale"),"-"}},
        {"Integrated luminosity",   {"-","-","-","-", to_string_with_precision(sysRanges.count("Lumi") ? sysRanges["Lumi"].first : 0)}},
        {"Cross section",           {"-","-","-","-", to_string_with_precision(sysRanges.count("XS")   ? sysRanges["XS"].first   : 0)}},
        {"QCD scale",               {"-","-","-", rangeStr("QCDmuF_VGamma"),
                                     to_string_with_precision(std::min(sysRanges.count("QCDmuR_rare") ? sysRanges["QCDmuR_rare"].first : 0,
                                                                       sysRanges.count("QCDmuF_rare") ? sysRanges["QCDmuF_rare"].first : 0)) + "-" +
                                     to_string_with_precision(std::max(sysRanges.count("QCDmuR_rare") ? sysRanges["QCDmuR_rare"].second : 0,
                                                                       sysRanges.count("QCDmuF_rare") ? sysRanges["QCDmuF_rare"].second : 0))}},
        {"PDF",                     {"-","-","-", rangeStr("PDF_VGamma"), rangeStr("PDF_rare")}},
    };
    for (auto &row : tableData) {
        std::cout << "  " << row.first;
        for (const auto &col : row.second) std::cout << " & " << col;
        std::cout << " \\\\\n";
    }
    std::cout << "  \\hline\n  \\end{tabular}}\n"
              << "  \\caption{Systematic uncertainties of the SM backgrounds.}\n"
              << "  \\label{table:ch4-systematic}\n\\end{table}\n";

    // =========================================================================
    // Save all combined systematic histograms to output ROOT file (from v1)
    // Output: /uscms_data/d3/tmishra/Output/SignalSystematic_<year><VFP>.root
    // =========================================================================
    TFile *outputfile = TFile::Open(
        Form("/uscms_data/d3/tmishra/Output/SignalSystematic_%d%s.root", RunYear, whichVFP.c_str()),
        "RECREATE");
    outputfile->cd();
    // -- e->gamma fake --
    h_elefakepho_norm->Write();           h_elefakepho_controlsample->Write();
    h_elefakepho_transferfactor->Write(); h_elefakepho_syserr_jes->Write();
    h_elefakepho_syserr_jer->Write();     h_elefakepho_syserr_esf->Write();
    h_elefakepho_syserr_scale->Write();   h_elefakepho_syserr_eleshape->Write();
    h_elefakepho_syserr_jetshape->Write();h_elefakepho_syserr_qcdshape->Write();
    h_elefakepho_syserr_xs->Write();      h_elefakepho_syserr_lumi->Write();
    h_elefakepho_syserr_isr->Write();
    // -- jet->gamma fake --
    h_jetfakepho_norm->Write();           h_jetfakepho_controlsample->Write();
    h_jetfakepho_transferfactor->Write(); h_jetfakepho_syserr_jes->Write();
    h_jetfakepho_syserr_jer->Write();     h_jetfakepho_syserr_esf->Write();
    h_jetfakepho_syserr_scale->Write();   h_jetfakepho_syserr_eleshape->Write();
    h_jetfakepho_syserr_jetshape->Write();h_jetfakepho_syserr_qcdshape->Write();
    h_jetfakepho_syserr_xs->Write();      h_jetfakepho_syserr_lumi->Write();
    h_jetfakepho_syserr_isr->Write();
    // -- QCD fake lepton --
    h_qcdfakelep_norm->Write();           h_qcdfakelep_controlsample->Write();
    h_qcdfakelep_transferfactor->Write(); h_qcdfakelep_syserr_jes->Write();
    h_qcdfakelep_syserr_jer->Write();     h_qcdfakelep_syserr_esf->Write();
    h_qcdfakelep_syserr_scale->Write();   h_qcdfakelep_syserr_eleshape->Write();
    h_qcdfakelep_syserr_jetshape->Write();h_qcdfakelep_syserr_qcdshape->Write();
    h_qcdfakelep_syserr_xs->Write();      h_qcdfakelep_syserr_lumi->Write();
    h_qcdfakelep_syserr_isr->Write();
    // -- V+gamma --
    h_VGamma_norm->Write();               h_VGamma_syserr_jes->Write();
    h_VGamma_syserr_jer->Write();         h_VGamma_syserr_esf->Write();
    h_VGamma_syserr_scale->Write();       h_VGamma_syserr_eleshape->Write();
    h_VGamma_syserr_jetshape->Write();    h_VGamma_syserr_qcdshape->Write();
    h_VGamma_syserr_xs->Write();          h_VGamma_syserr_lumi->Write();
    h_VGamma_syserr_isr->Write();         h_VGamma_syserr_pu->Write();
    h_VGamma_syserr_qcdmuF->Write();      h_VGamma_syserr_pdf->Write();
    // -- Rare EWK --
    h_rare_norm->Write();                 h_rare_syserr_jes->Write();
    h_rare_syserr_jer->Write();           h_rare_syserr_esf->Write();
    h_rare_syserr_scale->Write();         h_rare_syserr_eleshape->Write();
    h_rare_syserr_jetshape->Write();      h_rare_syserr_qcdshape->Write();
    h_rare_syserr_xs->Write();            h_rare_syserr_lumi->Write();
    h_rare_syserr_isr->Write();           h_rare_syserr_pu->Write();
    h_rare_syserr_qcdmuR->Write();        h_rare_syserr_qcdmuF->Write();
    h_rare_syserr_pdf->Write();
    // -- Data --
    h_sig->Write();
    outputfile->Close();

    // =========================================================================
    // Canvas and pad setup
    // Two-pad layout:
    //   pad1 (top, 65%): signal region event count plot (log scale)
    //   pad2 (bot, 35%): data/background ratio
    // =========================================================================
    TCanvas *can = new TCanvas("can", "", 1200, 800);
    gStyle->SetPadLeftMargin(0.10);
    can->SetLeftMargin(0.10);
    can->SetRightMargin(0.04);
    can->SetTopMargin(0.08);
    can->SetBottomMargin(0.1);
    can->cd();

    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
    setTopPad(pad1);
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

    const int col_ewk = TColor::GetColor("#9c9ca1");
    const int col_elefake = TColor::GetColor("#f89c20");
    const int col_jetfake = TColor::GetColor("#e42536");
    const int col_qcdfake = TColor::GetColor("#5790fc");
    const int col_vgamma  = TColor::GetColor("#964a8b");

    // Background colors
    h_bkg->SetFillColor(col_ewk); h_bkg->SetLineColor(col_ewk);
    h_bkg_elefakepho->SetFillColor(col_elefake); h_bkg_elefakepho->SetLineColor(col_elefake);
    h_bkg_jetfakepho->SetFillColor(col_jetfake); h_bkg_jetfakepho->SetLineColor(col_jetfake);
    h_bkg_qcdfakepho->SetFillColor(col_qcdfake); h_bkg_qcdfakepho->SetLineColor(col_qcdfake);
    h_bkg_VGamma->SetFillColor(col_vgamma); h_bkg_VGamma->SetLineColor(col_vgamma);

    // ----- Draw background stack first, then SUSY signals, then data (v2 draw order) -----
    h_bkg->Draw("hist same");         // total background (rare EWK base)
    h_bkg_VGamma->Draw("hist same");  // V+gamma on top
    h_bkg_jetfakepho->Draw("hist same");
    h_bkg_elefakepho->Draw("hist same");
    h_bkg_qcdfakepho->Draw("hist same");

    const int col_T5Wg    = TColor::GetColor("#832db6");
    const int col_TChiWG  = TColor::GetColor("#92dadd");

    p_t5wg_signal->SetLineColor(col_T5Wg);         // blue
    p_t5wg_signal->SetLineWidth(3);
    p_t5wg_signal->SetFillStyle(0);
    //p_t5wg_signal->Draw("HIST same");

    p_tchiwg_signal->SetLineColor(col_TChiWG);      // brown/olive
    p_tchiwg_signal->SetLineStyle(2);       // dashed
    p_tchiwg_signal->SetLineWidth(3);
    p_tchiwg_signal->SetFillStyle(0);
    //p_tchiwg_signal->Draw("HIST same");

    // =========================================================================
    // Build the observed data points with explicit Garwood/Poisson asymmetric
    // error bars (CMS recommendation, see getGarwoodInterval() above), rather
    // than relying on TH1::SetBinErrorOption(kPoisson) + Draw("E1P"). This
    // ensures (a) bin contents are correctly rounded to integer counts before
    // computing the interval, and (b) the exact same interval definition is
    // reused for the ratio-panel data errors below, so the two panels are
    // guaranteed consistent with each other (and with the corresponding
    // postfit plot, once that script is updated the same way).
    // =========================================================================
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
    g_data->Draw("P same");

    // ----- Legend -----
    TLegend *leg = new TLegend(0.15, 0.82, 0.92, 0.9);
    leg->SetNColumns(7);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    h_bkg->SetMarkerSize(0);
    h_bkg_elefakepho->SetMarkerSize(0);
    h_bkg_jetfakepho->SetMarkerSize(0);
    h_bkg_qcdfakepho->SetMarkerSize(0);
    h_bkg_VGamma->SetMarkerSize(0);
    error_bkg->SetMarkerSize(0);
    error_bkg->SetLineWidth(0);
    leg->AddEntry(g_data,           "Data",                      "ep");
    leg->AddEntry(h_bkg,           "t#bar{t}#gamma / WW#gamma / WZ#gamma");
    leg->AddEntry(h_bkg_VGamma,    "W#gamma / Z#gamma");
    //leg->AddEntry(p_t5wg_signal,   "T5Wg (1700,1000)"); 
    leg->AddEntry(h_bkg_qcdfakepho,"fake leptons");
    //leg->AddEntry(p_tchiwg_signal, "TChiWg (800)"); 
    leg->AddEntry(h_bkg_jetfakepho,"j->#gamma fake");
    leg->AddEntry(h_bkg_elefakepho,"e->#gamma fake");
    leg->AddEntry(error_bkg,       "Unc");
    leg->Draw("same");

    // ----- Total uncertainty band -----
    error_bkg->SetFillColor(12);
    error_bkg->SetFillStyle(3345);
    error_bkg->Draw("E2 same");

    // ----- Dividing lines between kinematic regions -----
    // BUG FIX: line x-positions shift by +1 to match the new histogram
    // binning (NBIN*2, 1, NBIN*2+1) where the boundary between bin 9 and
    // bin 10 is now at x=10, between bin 18 and 19 at x=19, etc.
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

    // ----- Kinematic region labels -----
    TLatex *latex = new TLatex();
    latex->SetTextSize(0.04);
    // BUG FIX: label x-positions shift by +1 to match new bin layout.
    // Bin 1 centre is now 1.5, bin 9 right-edge is 10, etc.
    latex->DrawLatex(2.5,  50000, "#mu#gamma, p_{T}^{#gamma} < 200 GeV");
    latex->DrawLatex(11.5, 50000, "#mu#gamma, p_{T}^{#gamma} > 200 GeV");
    latex->DrawLatex(20.5, 50000, "e#gamma, p_{T}^{#gamma} < 200 GeV");
    latex->DrawLatex(29.5, 50000, "e#gamma, p_{T}^{#gamma} > 200 GeV");
    latex->DrawLatex(4,  15000, "p^{miss}_{T} (GeV)");
    latex->DrawLatex(13,  15000, "p^{miss}_{T} (GeV)");
    latex->DrawLatex(22,  15000, "p^{miss}_{T} (GeV)");
    latex->DrawLatex(31,  15000, "p^{miss}_{T} (GeV)");
    // MET bin labels (repeated for each ptgamma block)
    latex->DrawLatex(2,  5000, "< 200");   latex->DrawLatex(4.2,    5000, "200,400"); latex->DrawLatex(7.4,  5000, "> 400");
    latex->DrawLatex(11, 5000, "< 200");   latex->DrawLatex(13.2,   5000, "200,400"); latex->DrawLatex(16.4, 5000, "> 400");
    latex->DrawLatex(20, 5000, "< 200");   latex->DrawLatex(22.2,   5000, "200,400"); latex->DrawLatex(25.4, 5000, "> 400");
    latex->DrawLatex(29, 5000, "< 200");   latex->DrawLatex(31.2,   5000, "200,400"); latex->DrawLatex(34.4, 5000, "> 400");

    // ----- CMS luminosity label (all years + Run2 combined) -----
    if      (RunYear == 2016 && preVFP == 1)  CMS_lumi(pad1, 1, 1, 11);  // 2016 preVFP
    else if (RunYear == 2016 && preVFP == 0)  CMS_lumi(pad1, 2, 1, 11);  // 2016 postVFP
    else if (RunYear == 2017)                  CMS_lumi(pad1, 3, 1, 11);  // 2017
    else if (RunYear == 2018)                  CMS_lumi(pad1, 4, 1, 11);  // 2018
    else if (RunYear == 678)                   CMS_lumi(pad1, 6, 1, 11);  // full Run 2 (v2 addition)

    pad1->RedrawAxis();  // redraw axes on top of histogram fills (v2 addition)

    // =========================================================================
    // Ratio pad (pad2): Data / Background
    // Uses a TH1F frame with numbered bin labels (v2 approach, cleaner than v1).
    // Uncertainty on ratio is propagated from totalerror / bkg.
    // Confidence interval for data uses the CMS-recommended Garwood
    // (Poisson, gamma-quantile-based) method -- see getGarwoodInterval()
    // and https://twiki.cern.ch/twiki/bin/view/CMS/PoissonErrorBars
    // =========================================================================
    can->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.35);
    pad2->SetBottomMargin(0.3);
    pad2->Draw();
    pad2->cd();

    // Create frame histogram with numerical bin labels 1..36
    // BUG FIX: was (36, 0, 36) so bin i had centre i-0.5, misaligning with the
    // upper pad where histograms are centred at integer values after the
    // declarations were changed to (NBIN*2, 1, NBIN*2+1).
    // Now (36, 1, 37): bin i spans [i, i+1) with centre i+0.5 — but since
    // the histogram bins are labelled 1..36 via SetBinLabel, the label "1"
    // appears at bin 1 whose centre is now 1.5, matching the upper-pad bin 1.
    TH1F *ratioframe = new TH1F("frame", ";Search region; Data / Pred.", 36, 1, 37);
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

    // Compute Garwood/Poisson confidence intervals for the data/bkg ratio,
    // using the SAME helper (and hence the SAME interval) as the upper-panel
    // data points above, so the two panels are guaranteed consistent.
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

    flatratio->Draw("same");

    // =========================================================================
    // Save canvas
    // =========================================================================
    can->Update();
    can->SaveAs(Form("/uscms_data/d3/tmishra/Output/signalCount_%d%s.png", RunYear, whichVFP.c_str()));
    can->SaveAs(Form("/uscms_data/d3/tmishra/Output/signalCount_%d%s.pdf", RunYear, whichVFP.c_str()));

    // =========================================================================
    // Validation test: cross-check gamma_quantile against TH1::kPoisson errors
    // (kept from both versions; L and U lines are commented for brevity)
    // =========================================================================
    TH1D *test_h = new TH1D("test_h", "", 18, 0, 18);
    test_h->SetBinErrorOption(TH1::kPoisson);
    for (unsigned ibin = 1; ibin <= 18; ibin++) test_h->SetBinContent(ibin, ibin-1);
    // Uncomment to print cross-check values:
    // for (unsigned ibin = 1; ibin <= 18; ibin++) {
    //     int N = test_h->GetBinContent(ibin);
    //     double L = (N==0)? 0: ROOT::Math::gamma_quantile(alpha/2, N, 1.);
    //     double U = ROOT::Math::gamma_quantile(1-alpha/2, N+1, 1.);
    //     cout << ibin << " N=" << N << " L=" << L << " U=" << U
    //          << " ROOT_low=" << test_h->GetBinErrorLow(ibin)
    //          << " ROOT_up=" << test_h->GetBinErrorUp(ibin) << endl;
    // }
}
