// =============================================================================
// plot_invmass.C
//
// Stacked data / prediction plot in bins of the invariant mass m(l,gamma) of the
// lepton + photon pair of the events selected in the signal region (SR).
// Analogue of plot_eventct.C, which does the same in the signal-region bins.
//
//   eg channel -> m(e,gamma)       mg channel -> m(mu,gamma)
//
// Inputs  (written by the updated pred_*.C scripts, same files plot_eventct.C reads):
//   signalTree_<ch>_eleBkg_<year>.root   e->gamma fakes       (toy-based shape error)
//   signalTree_<ch>_jetbkg_<year>.root   jet->gamma fakes     (toy-based shape error)
//   signalTree_<ch>_qcd_<year>.root      fake leptons         (norm-scale error)
//   signalTree_<ch>_VGBkg_<year>.root    W/Z + gamma          (JES, JER, ESF, scale, QCD muF, PDF)
//   signalTree_<ch>_rareBkg_<year>.root  ttgamma/WWgamma/...  (JES, JER, ESF, xs, lumi, QCD muR/muF, PDF)
//   signalTree_<ch>_signal_<year>.root   observed data
//   with <ch> = egamma (eg) or mg.
//
// Usage (from the same directory / with the same environment as plot_eventct.C):
//   root -l -b -q 'plot_invmass.C(1)'     // e+gamma channel   (default)
//   root -l -b -q 'plot_invmass.C(2)'     // mu+gamma channel
//   root -l -b -q 'plot_invmass.C(0)'     // both, as two separate plots
//
// The total uncertainty band is built exactly like in plot_eventct.C (quadrature
// sum of the same sources + the MC/proxy statistical errors), with ONE difference:
// the QCD-muF / QCD-muR / PDF histograms are stored by the pred_*.C scripts as
// FRACTIONAL variations (|kappa-1|), so here they are multiplied by the yield of
// the process before being added in quadrature (plot_eventct.C adds the raw
// fraction to the absolute yield uncertainties, which under-/over-states them).
//
// Output: invMassCount_<ch>_<year>.{png,pdf}, InvMassSystematic_<ch>_<year>.root
// =============================================================================

// ----- Standard library includes -----
#include <sstream>
#include <iomanip>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <initializer_list>

// ----- ROOT includes -----
#include <TROOT.h>
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"
#include "TColor.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "Math/QuantFuncMathCore.h"

// ----- Analysis framework includes -----
// Provides: setTDRStyle(), SetSignalConfig(), CMS_lumi(), setTopPad(),
//           and global variables RunYear, preVFP
#include "../../include/analysis_commoncode.h"
// Mass binning shared with the pred_*.C scripts
#include "invmass_binning.h"


// Same switch as in plot_eventct.C: treat 2016 pre+post VFP as one dataset.
// (names carry an "_invmass" tag so this file can be loaded in the same ROOT
//  session as plot_eventct.C without redefinition clashes)
bool total16_invmass = false;


// CMS-recommended Garwood (Poisson) interval for an observed count N, identical
// to getGarwoodInterval() in plot_eventct.C:
// https://twiki.cern.ch/twiki/bin/view/CMS/PoissonErrorBars
void getGarwoodIntervalInvMass(int N, double &L, double &U, double CL = 0.6827)
{
    double alpha = 1.0 - CL;
    L = (N == 0) ? 0.0 : ROOT::Math::gamma_quantile(alpha / 2.0,       N,   1.0);
    U =                  ROOT::Math::gamma_quantile(1.0 - alpha / 2.0, N+1, 1.0);
}


// Fetch a histogram from a file and check that it exists and has the mass
// binning of invmass_binning.h (i.e. that the pred_*.C scripts and this plot
// script were built from the same header).
TH1D *getInvMassHist(TFile *f, const std::string &name)
{
    TH1D *h = (TH1D*)f->Get(name.c_str());
    if (!h)
        throw std::runtime_error("histogram '" + name + "' not found in " + f->GetName() +
                                 " -- re-run the updated pred_*.C script that produces it");
    if (h->GetNbinsX() != nInvMassBins ||
        TMath::Abs(h->GetXaxis()->GetXmin() - invMassMin) > 1e-6 ||
        TMath::Abs(h->GetXaxis()->GetXmax() - invMassMax) > 1e-6)
        throw std::runtime_error("histogram '" + name + "' in " + f->GetName() +
                                 " does not have the binning of invmass_binning.h -- re-run the pred_*.C scripts");
    // pred_sig.C fills the data histogram without foldInvMass(): fold under/overflow here
    // (no-op for the background histograms, which are already folded when filled)
    const int n = h->GetNbinsX();
    h->SetBinContent(1, h->GetBinContent(0) + h->GetBinContent(1));
    h->SetBinContent(n, h->GetBinContent(n) + h->GetBinContent(n + 1));
    h->SetBinContent(0, 0.); h->SetBinContent(n + 1, 0.);
    return h;
}


// =============================================================================
// Main function
//   channel : 1 = e+gamma (default), 2 = mu+gamma, 0 = both (two separate plots)
// =============================================================================
void plot_invmass(int channel = 1)
{
    if (channel == 0) { plot_invmass(1); plot_invmass(2); return; }
    if (channel != 1 && channel != 2) {
        std::cerr << "plot_invmass: channel must be 0 (both), 1 (e+gamma) or 2 (mu+gamma)" << std::endl;
        return;
    }

    // ----- ROOT style setup (same as plot_eventct.C) -----
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    setTDRStyle();
    SetSignalConfig();                              // sets RunYear, preVFP, ... from config
    gStyle->SetErrorX(0);                           // no horizontal error bars on data points
    gSystem->Load("../../lib/libAnaClasses.so");

    // ----- VFP suffix for the 2016 file names -----
    std::string whichVFP;
    if      (RunYear == 2016 && preVFP == 1) whichVFP = "preVFP";
    else if (RunYear == 2016 && preVFP == 0) whichVFP = "postVFP";
    else                                     whichVFP = "";
    if (RunYear == 2016 && total16_invmass == true) whichVFP = "";

    // ----- Channel-dependent names -----
    const std::string chFile  = (channel == 1) ? "egamma" : "mg";   // tag in the file names
    const std::string chPre   = (channel == 1) ? "eg"     : "mg";   // prefix of the histogram names
    const TString     chLabel = (channel == 1) ? "e#gamma" : "#mu#gamma";
    const TString     xTitle  = (channel == 1) ? "m_{e#gamma} (GeV)" : "m_{#mu#gamma} (GeV)";
    const std::string yearTag = std::to_string(RunYear) + whichVFP;
    const std::string uniq    = chPre + "_" + yearTag;              // keeps ROOT object names unique
    const std::string outDir  = "/uscms_data/d3/tmishra/Output/";
    cout << "Channel: " << chFile << "\tYear: " << RunYear << "\t" << whichVFP << endl;

    // =========================================================================
    // Open the inputs
    // =========================================================================
    auto openFile = [&](const std::string &tag) -> TFile* {
        std::string name = outDir + "signalTree_" + chFile + "_" + tag + "_" + yearTag + ".root";
        TFile *f = TFile::Open(name.c_str());
        if (!f || f->IsZombie()) throw std::runtime_error("cannot open " + name);
        return f;
    };
    TFile *f_ele  = openFile("eleBkg");
    TFile *f_jet  = openFile("jetbkg");
    TFile *f_qcd  = openFile("qcd");
    TFile *f_VG   = openFile("VGBkg");
    TFile *f_rare = openFile("rareBkg");
    TFile *f_sig  = openFile("signal");

    const std::string P = chPre + "_";

    // ----- observed data -----
    TH1D *d_mass = getInvMassHist(f_sig, "p_invmass");

    // ----- nominal predictions (bin content = yield, bin error = statistical) -----
    TH1D *ele_norm  = getInvMassHist(f_ele,  P + "elefakepho_invmass_norm");
    TH1D *jet_norm  = getInvMassHist(f_jet,  P + "jetfakepho_invmass_norm");
    TH1D *qcd_norm  = getInvMassHist(f_qcd,  P + "qcdfakelep_invmass_norm");
    TH1D *VG_norm   = getInvMassHist(f_VG,   P + "VGamma_invmass_norm");
    TH1D *rare_norm = getInvMassHist(f_rare, P + "rare_invmass_norm");

    // ----- systematic uncertainties, same set as plot_eventct.C uses for the total -----
    // data-driven fakes: shape / scale error (absolute yield, -1 = not available)
    TH1D *ele_sys   = getInvMassHist(f_ele, P + "elefakepho_invmass_syserr_e_to_pho");
    TH1D *jet_sys   = getInvMassHist(f_jet, P + "jetfakepho_invmass_syserr_j_to_pho");
    TH1D *qcd_sys   = getInvMassHist(f_qcd, P + "qcdfakelep_invmass_syserr_scale");
    // V+gamma: absolute: jes, jer, esf, scale; fractional: qcdmuF, pdf
    TH1D *VG_jes    = getInvMassHist(f_VG,  P + "VGamma_invmass_syserr_jes");
    TH1D *VG_jer    = getInvMassHist(f_VG,  P + "VGamma_invmass_syserr_jer");
    TH1D *VG_esf    = getInvMassHist(f_VG,  P + "VGamma_invmass_syserr_esf");
    TH1D *VG_scale  = getInvMassHist(f_VG,  P + "VGamma_invmass_syserr_scale");
    TH1D *VG_muF    = getInvMassHist(f_VG,  P + "VGamma_invmass_syserr_qcdmuF");
    TH1D *VG_pdf    = getInvMassHist(f_VG,  P + "VGamma_invmass_syserr_pdf");
    // rare: absolute: jes, jer, esf, xs, lumi; fractional: qcdmuR, qcdmuF, pdf
    TH1D *rare_jes  = getInvMassHist(f_rare, P + "rare_invmass_syserr_jes");
    TH1D *rare_jer  = getInvMassHist(f_rare, P + "rare_invmass_syserr_jer");
    TH1D *rare_esf  = getInvMassHist(f_rare, P + "rare_invmass_syserr_esf");
    TH1D *rare_xs   = getInvMassHist(f_rare, P + "rare_invmass_syserr_xs");
    TH1D *rare_lumi = getInvMassHist(f_rare, P + "rare_invmass_syserr_lumi");
    TH1D *rare_muR  = getInvMassHist(f_rare, P + "rare_invmass_syserr_qcdmuR");
    TH1D *rare_muF  = getInvMassHist(f_rare, P + "rare_invmass_syserr_qcdmuF");
    TH1D *rare_pdf  = getInvMassHist(f_rare, P + "rare_invmass_syserr_pdf");

    // =========================================================================
    // Combine: yields, stacked histograms and uncertainties per mass bin
    //   h_<proc>        : yield of each process, bin error = stat (+) syst of that process
    //   h_bkg           : total background, bin error = total stat (+) syst
    //   h_bkg_qcdfakepho .. h_bkg_VGamma : cumulative sums used for the stacked drawing
    //                                      (same order as plot_eventct.C:
    //                                       fake lep, e->g, jet->g, V+gamma, rare on top)
    // =========================================================================
    auto newH = [&](const std::string &n) -> TH1D* {
        TH1D *h = new TH1D((n + "_" + uniq).c_str(), "", nInvMassBins, invMassMin, invMassMax);
        h->SetDirectory(0);     // keep it alive independent of the input files
        return h;
    };
    TH1D *h_data = newH("h_data");
    TH1D *h_ele  = newH("h_elefakepho");
    TH1D *h_jet  = newH("h_jetfakepho");
    TH1D *h_qcd  = newH("h_qcdfakelep");
    TH1D *h_VG   = newH("h_VGamma");
    TH1D *h_rare = newH("h_rare");
    TH1D *h_bkg            = newH("h_bkg");
    TH1D *h_bkg_qcdfakepho = newH("h_bkg_qcdfakepho");
    TH1D *h_bkg_elefakepho = newH("h_bkg_elefakepho");
    TH1D *h_bkg_jetfakepho = newH("h_bkg_jetfakepho");
    TH1D *h_bkg_VGamma     = newH("h_bkg_VGamma");

    TGraphErrors      *error_bkg   = new TGraphErrors(nInvMassBins);
    TGraphErrors      *error_ratio = new TGraphErrors(nInvMassBins);
    TGraphAsymmErrors *ratio       = new TGraphAsymmErrors();
    TGraphAsymmErrors *g_data      = new TGraphAsymmErrors(nInvMassBins);

    // systematic value of a histogram bin; -1 ("not available") counts as 0
    auto pos = [](TH1D *h, int i) -> double { double v = h->GetBinContent(i); return v > 0 ? v : 0.0; };
    auto sq  = [](double x) -> double { return x * x; };

    double sum_ele(0), sum_jet(0), sum_qcd(0), sum_VG(0), sum_rare(0), sum_bkg(0), sum_data(0);
    std::cout << "\nm(l,gamma) bins:  e->g   jet->g   fake-lep   V+g   rare   |  total +- stat +- syst  |  data\n";

    for (int ibin = 1; ibin <= nInvMassBins; ibin++) {
        const double nEle  = ele_norm ->GetBinContent(ibin);
        const double nJet  = jet_norm ->GetBinContent(ibin);
        const double nQcd  = qcd_norm ->GetBinContent(ibin);
        const double nVG   = VG_norm  ->GetBinContent(ibin);
        const double nRare = rare_norm->GetBinContent(ibin);
        const double bkgContent = nEle + nJet + nQcd + nVG + nRare;

        // ----- systematic uncertainty per process (quadrature sum of the sources) -----
        const double sysEle2  = sq(pos(ele_sys, ibin));
        const double sysJet2  = sq(pos(jet_sys, ibin));
        const double sysQcd2  = sq(pos(qcd_sys, ibin));
        double sysVG2 = 0;
        for (TH1D *h : {VG_jes, VG_jer, VG_esf, VG_scale}) sysVG2 += sq(pos(h, ibin));
        for (TH1D *h : {VG_muF, VG_pdf})                   sysVG2 += sq(pos(h, ibin) * nVG);     // fractional
        double sysRare2 = 0;
        for (TH1D *h : {rare_jes, rare_jer, rare_esf, rare_xs, rare_lumi}) sysRare2 += sq(pos(h, ibin));
        for (TH1D *h : {rare_muR, rare_muF, rare_pdf})                     sysRare2 += sq(pos(h, ibin) * nRare); // fractional

        // ----- statistical uncertainty per process -----
        const double statEle  = ele_norm ->GetBinError(ibin);
        const double statJet  = jet_norm ->GetBinError(ibin);
        const double statQcd  = qcd_norm ->GetBinError(ibin);
        const double statVG   = VG_norm  ->GetBinError(ibin);
        const double statRare = rare_norm->GetBinError(ibin);

        const double sysTot2  = sysEle2 + sysJet2 + sysQcd2 + sysVG2 + sysRare2;
        const double statTot2 = sq(statEle) + sq(statJet) + sq(statQcd) + sq(statVG) + sq(statRare);
        const double totalerror = std::sqrt(sysTot2 + statTot2);

        // ----- per-process histograms (also written to the output file) -----
        h_ele ->SetBinContent(ibin, nEle);  h_ele ->SetBinError(ibin, std::sqrt(sysEle2  + sq(statEle)));
        h_jet ->SetBinContent(ibin, nJet);  h_jet ->SetBinError(ibin, std::sqrt(sysJet2  + sq(statJet)));
        h_qcd ->SetBinContent(ibin, nQcd);  h_qcd ->SetBinError(ibin, std::sqrt(sysQcd2  + sq(statQcd)));
        h_VG  ->SetBinContent(ibin, nVG);   h_VG  ->SetBinError(ibin, std::sqrt(sysVG2   + sq(statVG)));
        h_rare->SetBinContent(ibin, nRare); h_rare->SetBinError(ibin, std::sqrt(sysRare2 + sq(statRare)));

        // ----- cumulative stack (each h_bkg_* is drawn on top of the previous ones) -----
        const double cumQcd = nQcd;
        const double cumEle = cumQcd + nEle;
        const double cumJet = cumEle + nJet;
        const double cumVG  = cumJet + nVG;
        h_bkg->SetBinContent(ibin, bkgContent);
        h_bkg->SetBinError(ibin, totalerror);
        h_bkg_qcdfakepho->SetBinContent(ibin, cumQcd);
        h_bkg_elefakepho->SetBinContent(ibin, cumEle);
        h_bkg_jetfakepho->SetBinContent(ibin, cumJet);
        h_bkg_VGamma    ->SetBinContent(ibin, cumVG);

        // ----- data -----
        const int N = TMath::Nint(d_mass->GetBinContent(ibin));   // observed integer count
        h_data->SetBinContent(ibin, N);
        double L, U;
        getGarwoodIntervalInvMass(N, L, U);
        const double x = h_bkg->GetBinCenter(ibin);
        const double halfw = 0.5 * h_bkg->GetBinWidth(ibin);
        g_data->SetPoint(ibin-1, x, N);
        g_data->SetPointEXlow(ibin-1,  0);
        g_data->SetPointEXhigh(ibin-1, 0);
        g_data->SetPointEYlow(ibin-1,  N - L);
        g_data->SetPointEYhigh(ibin-1, U - N);

        // ----- uncertainty band (upper pad: absolute, ratio pad: relative) -----
        error_bkg->SetPoint(ibin-1, x, bkgContent);
        error_bkg->SetPointError(ibin-1, halfw, totalerror);
        error_ratio->SetPoint(ibin-1, x, 1.0);
        error_ratio->SetPointError(ibin-1, halfw, bkgContent > 0 ? totalerror / bkgContent : 0.0);

        sum_ele += nEle; sum_jet += nJet; sum_qcd += nQcd; sum_VG += nVG; sum_rare += nRare;
        sum_bkg += bkgContent; sum_data += N;

        std::cout << std::fixed << std::setprecision(2)
                  << Form("%4.0f-%-4.0f", h_bkg->GetBinLowEdge(ibin), h_bkg->GetBinLowEdge(ibin) + h_bkg->GetBinWidth(ibin))
                  << "  " << nEle << "  " << nJet << "  " << nQcd << "  " << nVG << "  " << nRare
                  << "  |  " << bkgContent << " +- " << std::sqrt(statTot2) << " +- " << std::sqrt(sysTot2)
                  << "  |  " << N << "\n";
    }
    std::cout << std::fixed << std::setprecision(2)
              << "Sum over mass bins:  e->g " << sum_ele << "  jet->g " << sum_jet << "  fake-lep " << sum_qcd
              << "  V+g " << sum_VG << "  rare " << sum_rare << "  |  total " << sum_bkg << "  |  data " << sum_data << "\n";

    // ----- closure check: the sum over mass bins must reproduce the SR-bin yields -----
    // (the signal-bin histograms live in the same input files). Small differences are
    // expected for V+gamma / rare in bins where the low-stat JES/JER fix of the pred_*.C
    // scripts replaced the nominal yield.
    {
        auto srSum = [](TFile *f, const std::string &name) -> double {
            TH1D *h = (TH1D*)f->Get(name.c_str());
            return h ? h->Integral() : -1.0;
        };
        std::cout << "Closure, sum over SR bins (for comparison):  e->g " << srSum(f_ele,  P + "elefakepho_norm")
                  << "  jet->g "   << srSum(f_jet,  P + "jetfakepho_norm")
                  << "  fake-lep " << srSum(f_qcd,  P + "qcdfakelep_norm")
                  << "  V+g "      << srSum(f_VG,   P + "VGamma_norm")
                  << "  rare "     << srSum(f_rare, P + "rare_norm")
                  << "  |  data "  << srSum(f_sig,  "p_eventcount") << "\n";
    }

    // =========================================================================
    // Save the combined histograms (yields + total errors) for later use
    // =========================================================================
    {
        TFile *fout = TFile::Open((outDir + "InvMassSystematic_" + chFile + "_" + yearTag + ".root").c_str(), "RECREATE");
        fout->cd();
        h_data->Write("h_data");
        h_bkg ->Write("h_bkg_total");
        h_ele ->Write("h_elefakepho");
        h_jet ->Write("h_jetfakepho");
        h_qcd ->Write("h_qcdfakelep");
        h_VG  ->Write("h_VGamma");
        h_rare->Write("h_rare");
        fout->Close();
    }

    // =========================================================================
    // Canvas and pads: pad1 (top 65%) stacked yields, pad2 (bottom 35%) data / pred
    // =========================================================================
    TCanvas *can = new TCanvas(("can_" + uniq).c_str(), "", 1200, 800);
    gStyle->SetPadLeftMargin(0.10);
    can->SetLeftMargin(0.10);
    can->SetRightMargin(0.04);
    can->SetTopMargin(0.08);
    can->SetBottomMargin(0.1);
    can->cd();

    TPad *pad1 = new TPad(("pad1_" + uniq).c_str(), "pad1", 0, 0.35, 1, 1.0);
    setTopPad(pad1);
    pad1->Draw();
    pad1->cd();
    gPad->SetLogy();

    // ----- Axis ranges -----
    double ymaxData = 0, ymaxBkg = 0;
    for (int ibin = 1; ibin <= nInvMassBins; ibin++) {
        ymaxData = std::max(ymaxData, h_data->GetBinContent(ibin) + g_data->GetErrorYhigh(ibin-1));
        ymaxBkg  = std::max(ymaxBkg,  h_bkg ->GetBinContent(ibin) + h_bkg->GetBinError(ibin));
    }
    const double ymax = 200.0 * std::max(1.0, std::max(ymaxData, ymaxBkg));  // head-room for legend + labels
    h_bkg->GetYaxis()->SetRangeUser(0.05, ymax);
    h_bkg->GetXaxis()->SetRangeUser(invMassMin, invMassMax);
    h_bkg->GetXaxis()->SetLabelSize(0);       // x labels are shown on the ratio pad
    h_bkg->GetYaxis()->SetTitle("Events / bin");
    h_bkg->GetYaxis()->SetTitleOffset(1.4);

    // ----- Background stack colours (same as plot_eventct.C) -----
    const int col_ewk     = TColor::GetColor("#9c9ca1");
    const int col_elefake = TColor::GetColor("#f89c20");
    const int col_jetfake = TColor::GetColor("#e42536");
    const int col_qcdfake = TColor::GetColor("#5790fc");
    const int col_vgamma  = TColor::GetColor("#964a8b");

    h_bkg->SetFillColor(col_ewk);                  h_bkg->SetLineColor(col_ewk);
    h_bkg_elefakepho->SetFillColor(col_elefake);   h_bkg_elefakepho->SetLineColor(col_elefake);
    h_bkg_jetfakepho->SetFillColor(col_jetfake);   h_bkg_jetfakepho->SetLineColor(col_jetfake);
    h_bkg_qcdfakepho->SetFillColor(col_qcdfake);   h_bkg_qcdfakepho->SetLineColor(col_qcdfake);
    h_bkg_VGamma->SetFillColor(col_vgamma);        h_bkg_VGamma->SetLineColor(col_vgamma);

    // ----- Draw: total first (rare EWK is what remains visible on top), then the cumulative sums -----
    h_bkg->Draw("hist");
    h_bkg_VGamma->Draw("hist same");
    h_bkg_jetfakepho->Draw("hist same");
    h_bkg_elefakepho->Draw("hist same");
    h_bkg_qcdfakepho->Draw("hist same");

    // ----- Data points with Garwood asymmetric errors, drawn on top -----
    g_data->SetMarkerStyle(20);
    g_data->SetMarkerSize(1.0);
    g_data->SetLineWidth(1);
    g_data->SetLineColor(kBlack);
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
    leg->AddEntry(g_data,           "Data",                        "ep");
    leg->AddEntry(h_bkg,            "t#bar{t}#gamma / WW#gamma / WZ#gamma");
    leg->AddEntry(h_bkg_VGamma,     "W#gamma / Z#gamma");
    leg->AddEntry(h_bkg_qcdfakepho, "misid. leptons");
    leg->AddEntry(h_bkg_jetfakepho, "j#rightarrow#gamma misid.");
    leg->AddEntry(h_bkg_elefakepho, "e#rightarrow#gamma misid.");
    leg->AddEntry(error_bkg,        "Unc");
    leg->Draw("same");

    // ----- Total uncertainty band -----
    error_bkg->SetFillColor(12);
    error_bkg->SetFillStyle(3345);
    error_bkg->Draw("E2 same");

    // ----- Channel / selection label -----
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.04);
    // region label from the config cuts (the inputs are whatever pred_*.C wrote with the same config)
    TString region;
    if (highMET > 0) region = Form("control region, p_{T}^{miss} < %d GeV", highMET);
    else             region = Form("p_{T}^{miss} > %d GeV, m_{T} > %d GeV", lowMET, lowMt);
    latex->DrawLatex(0.14, 0.76, (chLabel + ", " + region + " (last bin includes overflow)").Data());

    // ----- CMS luminosity label (same convention as plot_eventct.C) -----
    if      (RunYear == 2016 && preVFP == 1)  CMS_lumi(pad1, 1, 1, 11);  // 2016 preVFP
    else if (RunYear == 2016 && preVFP == 0)  CMS_lumi(pad1, 2, 1, 11);  // 2016 postVFP
    else if (RunYear == 2017)                  CMS_lumi(pad1, 3, 1, 11);  // 2017
    else if (RunYear == 2018)                  CMS_lumi(pad1, 4, 1, 11);  // 2018
    else if (RunYear == 678)                   CMS_lumi(pad1, 6, 1, 11);  // full Run 2

    pad1->RedrawAxis();

    // =========================================================================
    // Ratio pad: Data / Prediction, Garwood errors on data (same interval as above),
    // relative total uncertainty band on the prediction
    // =========================================================================
    can->cd();
    TPad *pad2 = new TPad(("pad2_" + uniq).c_str(), "pad2", 0, 0, 1, 0.35);
    pad2->SetBottomMargin(0.3);
    pad2->Draw();
    pad2->cd();

    TH1F *ratioframe = new TH1F(("frame_" + uniq).c_str(), (";" + xTitle + "; Data / Pred.").Data(), nInvMassBins, invMassMin, invMassMax);
    ratioframe->GetYaxis()->SetTitleOffset(1.4);
    ratioframe->GetYaxis()->SetRangeUser(0, 2.0);
    ratioframe->Draw();

    TLine *flatratio = new TLine(invMassMin, 1, invMassMax, 1);

    int ipt = 0;
    for (int ibin = 1; ibin <= nInvMassBins; ibin++) {
        const int    N      = TMath::Nint(h_data->GetBinContent(ibin));
        const double bkgval = h_bkg->GetBinContent(ibin);
        if (bkgval <= 0) continue;                       // no prediction in this bin
        double L, U;
        getGarwoodIntervalInvMass(N, L, U);
        ratio->SetPoint(ipt, h_bkg->GetBinCenter(ibin), (double)N / bkgval);
        ratio->SetPointEXlow(ipt, 0);
        ratio->SetPointEXhigh(ipt, 0);
        ratio->SetPointEYlow( ipt, (N - L) / bkgval);
        ratio->SetPointEYhigh(ipt, (U - N) / bkgval);
        ipt++;
    }
    ratio->SetMarkerStyle(20);
    ratio->SetLineColor(kBlack);
    ratio->Draw("P same");

    error_ratio->SetFillColor(12);
    error_ratio->SetFillStyle(3345);
    error_ratio->Draw("E2 same");

    flatratio->Draw("same");

    // =========================================================================
    // Save canvas
    // =========================================================================
    can->Update();
    can->SaveAs((outDir + "invMassCount_" + chFile + "_" + yearTag + ".png").c_str());
    can->SaveAs((outDir + "invMassCount_" + chFile + "_" + yearTag + ".pdf").c_str());
}
