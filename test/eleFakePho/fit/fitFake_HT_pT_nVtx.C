#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "../../../include/tdrstyle.C"
#include "TROOT.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooAbsReal.h"
#include "RooKeysPdf.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooNumIntConfig.h"
#include "RooFFTConvPdf.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooHistPdf.h"
#include "RooMCStudy.h"
#include "RooChi2MCSModule.h"
#include "RooFitResult.h"
#include "TFitResult.h"
#include "TVirtualFitter.h"
#include "RooMultiVarGaussian.h"
#include "../../../include/RooCMSShape.h"
#include "../../../include/RooDCBShape.h"
#include "../../../include/RooUserPoly.h"
// #include "../../../include/analysis_fakes.h"
#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_mcData.h"
#define NTOY 100

#define MAXVTX 48

int MINVTX = 18;
bool doEB = true;
const char*processName ="Data"; bool doDrellYan = false;

bool RUN_VS_NVTX_BY_HT   = false;   
bool RUN_VS_NVTX_BY_PT   = true;  
bool RUN_VS_PT_BY_NVTX   = true;  
bool RUN_NOMINAL_VS_PT   = true;   // fake rate vs pT,   no nVtx sub-range (fr_bothcount_pt)
bool RUN_NOMINAL_VS_NVTX = true;   // fake rate vs nVtx, no pT sub-range   (fr_bothcount_vtx)

// Consistency-check overlays: sub-ranged curves (points only, no error bars)
// drawn on top of the nominal curve (light-grey band in the background).
bool RUN_PT_CONSISTENCY   = true;  // pT dependence:  nVtx8-18/18-28/28-38/>38 + nominal   (5 sets)
bool RUN_NVTX_CONSISTENCY = true;  // nVtx dependence: pT30-60/60-120/>120 + nominal        (4 sets)

// -----------------------------------------------------------------------------
// Which functional form to fit with. nVtx-axis plots keep the original straight
// line (pol1). pT-axis plots use a dedicated turn-on-style parameterization:
//   fakerate(pT) = (par[0]*pT + par[1])^par[2]   (slope, constant, index)
// -----------------------------------------------------------------------------
enum class FitFunctionType { kPol1, kPtDependence };

double fakerate_ptDependence(double *x, double *par)
{
    double slope    = par[0];
    double constant = par[1];
    double index    = par[2];
    double pt = TMath::Max(x[0], 0.000001);
    double arg = slope*pt + constant;
    if(arg <= 0) return 0.0;
    return pow(arg, index);
}

// -----------------------------------------------------------------------------
// Build a TGraphErrors of fake rate vs the binning variable from already-parsed
// den/num arrays (with their errors) and fit it. Shared by the plain
// (statistical-only) reader and the systematic-band reader below.
// -----------------------------------------------------------------------------
void BuildGraphAndFit(const std::vector<double> &x,
                       const std::vector<double> &den, const std::vector<double> &denErr,
                       const std::vector<double> &num, const std::vector<double> &numErr,
                       double fitXMin, double fitXMax,
                       TGraphErrors *&gr, TF1 *&fit, int color, const char* fitname,
                       FitFunctionType fitFuncType, int RunYear,
                       TFitResultPtr *outResult = nullptr)
{
    int nBins = (int)x.size();

    gr = new TGraphErrors(nBins);
    for(int i=0;i<nBins;i++){
        double fr = num[i]/den[i];
        double er = sqrt( pow(numErr[i]/den[i],2)
                        + pow(num[i]*denErr[i]/(den[i]*den[i]),2) );

        // Bin center (not left edge): midpoint to the next bin's low edge, or
        // to fitXMax for the last (overflow) bin -- matches fitFakeFunc.C.
        double xvalue = (i+1 < nBins) ? (x[i]+x[i+1])/2.0 : (fitXMax + x[i])/2.0;
        double xerror = (i+1 < nBins) ? (x[i+1]-x[i])/2.0 : (fitXMax - x[i])/2.0;

        gr->SetPoint(i,xvalue,fr);
        gr->SetPointError(i,xerror,er);
    }

    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(color);
    gr->SetLineColor(color);
    gr->SetLineWidth(2);

    // ---------- Fit ----------
    if(fitFuncType == FitFunctionType::kPol1){
        fit = new TF1(fitname,"pol1",fitXMin,fitXMax);
        fit->SetLineColor(color);
        fit->SetParameters(0.02,0.0);

        TFitResultPtr r = gr->Fit(fit,"S");
        if(outResult) *outResult = r;
    }
    else { // FitFunctionType::kPtDependence
        fit = new TF1(fitname, fakerate_ptDependence, fitXMin, fitXMax, 3);
        fit->SetLineColor(color);
        fit->SetParNames("slope","constant","index");

        TVirtualFitter::SetMaxIterations(1000000);

        fit->SetParameter(0, 1.0e+04);
        fit->SetParLimits(0, 1.001e+02, 1.25601e+04);
        fit->SetParameter(2, -0.4);
        fit->SetParLimits(2, -0.4, -0.1);

        // "constant" (par[1]) starting point/limits are year-dependent,
        // matching the working reference fit in fitFakeFunc.C (Data only).
        if(RunYear==2016){
            fit->SetParameter(1, -10000);
            fit->SetParLimits(1, -15000.0, -5000);
        } else if(RunYear==2017){
            fit->SetParameter(1, -20000);
            fit->SetParLimits(1, -39070.0, -10000);
        } else if(RunYear==2018){
            fit->SetParameter(1, -30000);
            fit->SetParLimits(1, -38685.0, 0);
        } else {
            fit->SetParameter(1, -20000);
            fit->SetParLimits(1, -39070.0, 0);
        }

        TFitResultPtr r = gr->Fit(fit,"R S");
        if(outResult) *outResult = r;
    }

    double chi2ndf= (fit->GetNDF()>0) ? fit->GetChisquare()/fit->GetNDF() : -1;

    if(fitFuncType == FitFunctionType::kPol1){
        std::cout<<"Intercept p0 = "<<fit->GetParameter(0)<<" +/- "<<fit->GetParError(0)<<std::endl;
        std::cout<<"Slope     p1 = "<<fit->GetParameter(1)<<" +/- "<<fit->GetParError(1)<<std::endl;
    }
    else {
        std::cout<<"slope    = "<<fit->GetParameter(0)<<" +/- "<<fit->GetParError(0)<<std::endl;
        std::cout<<"constant = "<<fit->GetParameter(1)<<" +/- "<<fit->GetParError(1)<<std::endl;
        std::cout<<"index    = "<<fit->GetParameter(2)<<" +/- "<<fit->GetParError(2)<<std::endl;
    }
    std::cout<<"chi2/ndf     = "<<chi2ndf<<std::endl;
}

// -----------------------------------------------------------------------------
// Plain reader: statistical-only errors (fitrms from the one input file).
// Used for the HT/pT/nVtx cross-binning comparisons, which don't have DY/expo
// systematic-variant files produced per cross-bin sub-range.
// -----------------------------------------------------------------------------
bool ReadFakeRateFile(const std::string &path, double fitXMin, double fitXMax,
                       TGraphErrors *&gr, TF1 *&fit, int color, const char* fitname,
                       FitFunctionType fitFuncType, int RunYear)
{
    std::ifstream fin(path);
    if(!fin.is_open()){
        std::cout<<"Cannot open "<<path<<std::endl;
        return false;
    }

    std::vector<double> x, den, num, denErr, numErr;
    std::string bintype,numtype;
    float lowcut, signal1, error1, fitmean, fitrms;

    std::string line;
    int nLines=0;
    while(std::getline(fin,line)) nLines++;
    fin.clear(); fin.seekg(0);

    int nBins = nLines/2;

    for(int i=0;i<nBins;i++){
        fin>>bintype>>numtype>>lowcut>>signal1>>error1>>fitmean>>fitrms;
        den.push_back(fitmean);
        denErr.push_back(fitrms);
        x.push_back(lowcut);
    }
    for(int i=0;i<nBins;i++){
        fin>>bintype>>numtype>>lowcut>>signal1>>error1>>fitmean>>fitrms;
        num.push_back(fitmean);
        numErr.push_back(fitrms);
    }

    std::cout<<"===================================="<<std::endl;
    BuildGraphAndFit(x, den, denErr, num, numErr, fitXMin, fitXMax,
                      gr, fit, color, fitname, fitFuncType, RunYear);
    return true;
}

// -----------------------------------------------------------------------------
// Systematic-band reader: combines the nominal (Bw-ker) result with a DY-ker
// variant and a Bw-expo variant, taking the larger |variant - nominal| per bin
// as the systematic, added in quadrature with the statistical fitrms --
// matches fitFakeFunc.C's pt_denerror/pt_numerror (and vtx_* ) construction.
// Falls back to statistical-only errors if the variant files aren't found.
// -----------------------------------------------------------------------------
bool ReadFakeRateFileWithSyst(const std::string &nominalPath,
                               const std::string &dyPath,
                               const std::string &polPath,
                               double fitXMin, double fitXMax,
                               TGraphErrors *&gr, TF1 *&fit, int color, const char* fitname,
                               FitFunctionType fitFuncType, int RunYear,
                               TFitResultPtr *outResult = nullptr)
{
    std::ifstream finNom(nominalPath);
    if(!finNom.is_open()){
        std::cout<<"Cannot open "<<nominalPath<<std::endl;
        return false;
    }
    std::ifstream finDY(dyPath);
    std::ifstream finPol(polPath);
    bool haveSyst = finDY.is_open() && finPol.is_open();
    if(!haveSyst){
        std::cout<<"Warning: systematic variant file(s) not found ("<<dyPath<<", "<<polPath
                  <<"), falling back to statistical-only errors for "<<nominalPath<<std::endl;
    }

    std::string bintype,numtype;
    float lowcut, signal1, error1, fitmean, fitrms;
    float DYmean, DYrms, Polmean, Polrms;

    std::string line;
    int nLines=0;
    while(std::getline(finNom,line)) nLines++;
    finNom.clear(); finNom.seekg(0);

    int nBins = nLines/2;

    std::vector<double> x, den, denErr, num, numErr;

    for(int i=0;i<nBins;i++){
        finNom>>bintype>>numtype>>lowcut>>signal1>>error1>>fitmean>>fitrms;
        x.push_back(lowcut);
        den.push_back(fitmean);
        if(haveSyst){
            finDY>>bintype>>numtype>>lowcut>>signal1>>error1>>DYmean>>DYrms;
            finPol>>bintype>>numtype>>lowcut>>signal1>>error1>>Polmean>>Polrms;
            double sysdiff = std::max(fabs(DYmean-fitmean), fabs(Polmean-fitmean));
            denErr.push_back(sqrt(fitrms*fitrms + sysdiff*sysdiff));
        } else {
            denErr.push_back(fitrms);
        }
    }
    for(int i=0;i<nBins;i++){
        finNom>>bintype>>numtype>>lowcut>>signal1>>error1>>fitmean>>fitrms;
        num.push_back(fitmean);
        if(haveSyst){
            finDY>>bintype>>numtype>>lowcut>>signal1>>error1>>DYmean>>DYrms;
            finPol>>bintype>>numtype>>lowcut>>signal1>>error1>>Polmean>>Polrms;
            double sysdiff = std::max(fabs(DYmean-fitmean), fabs(Polmean-fitmean));
            numErr.push_back(sqrt(fitrms*fitrms + sysdiff*sysdiff));
        } else {
            numErr.push_back(fitrms);
        }
    }

    std::cout<<"===================================="<<std::endl;
    std::cout<<(haveSyst ? "Using DY/expo systematic band for " : "Statistical-only for ")<<nominalPath<<std::endl;
    BuildGraphAndFit(x, den, denErr, num, numErr, fitXMin, fitXMax,
                      gr, fit, color, fitname, fitFuncType, RunYear, outResult);
    return true;
}

// -----------------------------------------------------------------------------
// Build a full comparison plot (several graphs+fits overlaid) from a list of
// result files, and save it to a .png. `outSuffix` is used both to keep ROOT
// object names unique and to tag the output filename.
// -----------------------------------------------------------------------------
void PlotFakeRateComparison(
    const std::vector<std::string> &files,
    const std::vector<std::string> &labels,
    const std::vector<int> &colors,
    const char* xTitle, double xAxisMin, double xAxisMax, int nBinsAxis,
    double fitXMin, double fitXMax, double yMax,
    const char* outSuffix, int RunYear, const std::string &whichVFP,
    FitFunctionType fitFuncType,
    const std::vector<std::string> &dyFiles = {},
    const std::vector<std::string> &polFiles = {})
{
    // Use the DY/expo systematic band only if a matching variant file was
    // supplied for every entry in `files`; otherwise fall back to the plain,
    // statistical-only reader (this is the case for Modes 1-3, which don't
    // have DY/expo variants produced per cross-binning sub-range).
    bool useSyst = (dyFiles.size() == files.size()) && (polFiles.size() == files.size());

    std::vector<TGraphErrors*> graphs;
    std::vector<TF1*> fits;

    for(size_t ifile=0; ifile<files.size(); ifile++){
        std::cout<<"===================================="<<std::endl;
        std::cout<<"Fit results : "<<labels[ifile]<<std::endl;

        TGraphErrors *gr=nullptr;
        TF1 *fit=nullptr;
        bool ok;
        if(useSyst){
            ok = ReadFakeRateFileWithSyst(files[ifile], dyFiles[ifile], polFiles[ifile],
                                           fitXMin, fitXMax, gr, fit,
                                           colors[ifile], Form("fit_%s_%zu",outSuffix,ifile),
                                           fitFuncType, RunYear);
        } else {
            ok = ReadFakeRateFile(files[ifile], fitXMin, fitXMax, gr, fit,
                                   colors[ifile], Form("fit_%s_%zu",outSuffix,ifile),
                                   fitFuncType, RunYear);
        }
        if(!ok) continue;

        graphs.push_back(gr);
        fits.push_back(fit);
    }

    if(graphs.empty()){
        std::cout<<"No usable files found for "<<outSuffix<<", skipping this plot."<<std::endl;
        return;
    }

    // ---------- Plot ----------
    TCanvas *c = new TCanvas(Form("c_%s",outSuffix),"",700,600);

    TH1D *dummy = new TH1D(Form("dummy_%s",outSuffix), Form(";%s;fake rate",xTitle),
                            nBinsAxis, xAxisMin, xAxisMax);
    dummy->SetMaximum(yMax);
    dummy->Draw();

    TLegend *leg = new TLegend(0.45,0.7,0.78,0.88);

    for(size_t i=0; i<graphs.size(); i++){
        graphs[i]->Draw("EP same");
        fits[i]->Draw("same");
        leg->AddEntry(graphs[i],labels[i].c_str(),"pe");
    }

    leg->Draw();

    c->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_%s_%d%s.png",
                   outSuffix, RunYear, whichVFP.c_str()));
}

// -----------------------------------------------------------------------------
// Consistency-check overlay: draws the nominal fake rate (statistical + DY/expo
// systematic band) as a light-grey filled band in the background, then overlays
// the sub-range (nVtx- or pT-binned) curves as bare points on top -- no
// horizontal bin-width error bars and no vertical (statistical) error bars on
// those points, so the eye compares shapes/normalization against the nominal
// band rather than getting distracted by per-point error bars.
// -----------------------------------------------------------------------------
void PlotFakeRateConsistency(
    const std::vector<std::string> &files,
    const std::vector<std::string> &labels,
    const std::vector<int> &colors,
    const std::string &nominalFile,
    const std::string &nominalDyFile,
    const std::string &nominalPolFile,
    const char* xTitle, double xAxisMin, double xAxisMax, int nBinsAxis,
    double fitXMin, double fitXMax, double yMax,
    const char* outSuffix, int RunYear, const std::string &whichVFP,
    FitFunctionType fitFuncType)
{
    // ---------- Nominal (background band) ----------
    TGraphErrors *grNom=nullptr;
    TF1 *fitNom=nullptr;
    TFitResultPtr fitResultNom;
    bool okNom = ReadFakeRateFileWithSyst(nominalFile, nominalDyFile, nominalPolFile,
                                           fitXMin, fitXMax, grNom, fitNom,
                                           kRed, Form("fit_%s_nom",outSuffix),
                                           fitFuncType, RunYear, &fitResultNom);

    // Nominal is represented by only the fit curve and its systematic
    // uncertainty band (no data points), both in red -- matching
    // fitFakeFunc.C's "Fit result" / "Syst. Unc." style. The band itself is
    // built the same way fitFakeFunc.C builds fr_pt_sigmaband/fr_vtx_sigmaband:
    // toy pseudo-experiments drawn from the fit's parameter covariance matrix
    // (RooMultiVarGaussian), evaluated on a fine x grid, band half-width =
    // (max-min)/2 of the toys at each point.
    TGraphErrors *grNomBand = nullptr;
    if(okNom && fitResultNom.Get() != nullptr){
        int nPar = (fitFuncType==FitFunctionType::kPtDependence) ? 3 : 2;

        TMatrixDSym cov = fitResultNom->GetCovarianceMatrix();
        TVectorD mean(nPar);
        for(int p=0;p<nPar;p++) mean(p) = fitResultNom->Parameter(p);

        RooArgList parList;
        std::vector<RooRealVar*> vars;
        for(int p=0;p<nPar;p++){
            double err = fitResultNom->ParError(p);
            vars.push_back(new RooRealVar(Form("toypar_%s_%d",outSuffix,p), "",
                                           mean(p)-err, mean(p)+err));
            parList.add(*vars.back());
        }

        RooMultiVarGaussian mvg(Form("mvg_%s",outSuffix), "", parList, mean, cov);
        RooDataSet *toyData = mvg.generate(RooArgSet(parList), NTOY);

        const int nSmooth = 200;
        double xLo = fitXMin, xHi = fitXMax;
        std::vector<double> loVal(nSmooth, 1e30), hiVal(nSmooth, -1e30);

        for(int itoy=0; itoy<NTOY; itoy++){
            double p[3] = {0,0,0};
            for(int pi=0; pi<nPar; pi++)
                p[pi] = toyData->get(itoy)->getRealValue(vars[pi]->GetName());

            TF1 *toyF;
            if(fitFuncType==FitFunctionType::kPtDependence){
                toyF = new TF1(Form("toyF_%s_%d",outSuffix,itoy), fakerate_ptDependence, xLo, xHi, 3);
                toyF->SetParameters(p[0],p[1],p[2]);
            } else {
                toyF = new TF1(Form("toyF_%s_%d",outSuffix,itoy), "pol1", xLo, xHi);
                toyF->SetParameters(p[0],p[1]);
            }
            for(int j=0;j<nSmooth;j++){
                double xx = xLo + (xHi-xLo)*j/(nSmooth-1);
                double val = toyF->Eval(xx);
                if(val<loVal[j]) loVal[j]=val;
                if(val>hiVal[j]) hiVal[j]=val;
            }
            delete toyF;
        }

        grNomBand = new TGraphErrors(nSmooth);
        for(int j=0;j<nSmooth;j++){
            double xx = xLo + (xHi-xLo)*j/(nSmooth-1);
            grNomBand->SetPoint(j, xx, fitNom->Eval(xx));
            grNomBand->SetPointError(j, 0, 0.5*(hiVal[j]-loVal[j]));
        }
        grNomBand->SetFillColorAlpha(kRed, 0.30);
        grNomBand->SetFillStyle(1001);
        grNomBand->SetLineColor(kRed);
        grNomBand->SetLineWidth(0);
        grNomBand->SetMarkerSize(0);

        fitNom->SetLineColor(kRed);
        fitNom->SetLineStyle(1);
        fitNom->SetLineWidth(2);

        delete toyData;
        for(auto v : vars) delete v;
    }

    // ---------- Sub-range curves (points only, no error bars) ----------
    // Built as plain TGraphs (x,y only -- no error arrays at all) so no
    // vertical or horizontal error bar can be drawn regardless of option.
    std::vector<TGraph*> drawPts;
    std::vector<TF1*> fits;

    for(size_t ifile=0; ifile<files.size(); ifile++){
        std::cout<<"===================================="<<std::endl;
        std::cout<<"Fit results : "<<labels[ifile]<<std::endl;

        TGraphErrors *gr=nullptr;
        TF1 *fit=nullptr;
        bool ok = ReadFakeRateFile(files[ifile], fitXMin, fitXMax, gr, fit,
                                    colors[ifile], Form("fit_%s_%zu",outSuffix,ifile),
                                    fitFuncType, RunYear);
        if(!ok) continue;

        TGraph *pts = new TGraph(gr->GetN(), gr->GetX(), gr->GetY());
        pts->SetMarkerStyle(20);
        pts->SetMarkerColor(colors[ifile]);
        pts->SetMarkerSize(1.1);

        drawPts.push_back(pts);
        fits.push_back(fit);
    }

    if(drawPts.empty() && !okNom){
        std::cout<<"No usable files found for "<<outSuffix<<", skipping this plot."<<std::endl;
        return;
    }

    // ---------- Plot ----------
    TCanvas *c = new TCanvas(Form("c_%s",outSuffix),"",700,600);

    TH1D *dummy = new TH1D(Form("dummy_%s",outSuffix), Form(";%s;fake rate",xTitle),
                            nBinsAxis, xAxisMin, xAxisMax);
    dummy->SetMaximum(yMax);
    dummy->Draw();

    TLegend *leg = new TLegend(0.45,0.7,0.78,0.88);

    // Nominal drawn first (red band + red fit line) so it sits in the
    // background -- only the fit and its systematic band are shown for
    // nominal, no data points.
    if(okNom && grNomBand){
        grNomBand->Draw("e3 same");
        fitNom->Draw("same");
        leg->AddEntry(fitNom,"Fit result","l");
        leg->AddEntry(grNomBand,"Syst. Unc.","f");
    }

    // Sub-range curves drawn on top, points only -- plain TGraph with "P":
    // no horizontal bin-width bars, no vertical error bars.
    for(size_t i=0; i<drawPts.size(); i++){
        drawPts[i]->Draw("P same");
        fits[i]->Draw("same");
        leg->AddEntry(drawPts[i],labels[i].c_str(),"p");
    }

    leg->Draw();

    c->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_%s_%d%s.png",
                   outSuffix, RunYear, whichVFP.c_str()));
}

void fitFake_HT_pT_nVtx(int RunYear, bool preVFP)
{
    gROOT->SetBatch(1);
    setTDRStyle();
    gStyle->SetOptStat(0);

    std::string whichVFP="";
    if(RunYear==2016 && preVFP)  whichVFP="preVFP";
    if(RunYear==2016 && !preVFP) whichVFP="postVFP";

    // =========================================================================
    // Mode 1: fake rate vs nVtx, in bins of HT  (original/default behavior)
    // Reads DATAResult<year><VFP>_HT<range>/EleFakeRate-DATA-Bw-ker-vtx-60-120.txt
    // =========================================================================
    if(RUN_VS_NVTX_BY_HT){
        std::vector<std::string> dirs   = {"HT0-100","HT100-400","HT-gt400"};
        std::vector<std::string> labels = {"HT 0-100","HT 100-400","HT >400"};
        std::vector<int> colors         = {kBlue,kRed,kGreen+2};

        std::vector<std::string> files;
        for(auto &d : dirs)
            files.push_back(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s_%s/EleFakeRate-DATA-Bw-ker-vtx-60-120.txt",
                                  RunYear, whichVFP.c_str(), d.c_str()));

        PlotFakeRateComparison(files, labels, colors,
                                "nVtx", 0, MAXVTX, MAXVTX,
                                MINVTX, 46, 0.07,
                                "vtx_HTcompare", RunYear, whichVFP,
                                FitFunctionType::kPol1);
    }

    // =========================================================================
    // Mode 2: fake rate vs nVtx, in bins of pT
    // Reads DATAResult<year><VFP>_pT<range>/EleFakeRate-DATA-Bw-ker-vtx-60-120.txt
    // =========================================================================
    if(RUN_VS_NVTX_BY_PT){
        std::vector<std::string> dirs   = {"pT30-60","pT60-120","pT-gt120"};
        std::vector<std::string> labels = {"p_{T} 30-60","p_{T} 60-120","p_{T} >120"};
        std::vector<int> colors         = {kBlue,kRed,kGreen+2};

        std::vector<std::string> files;
        for(auto &d : dirs)
            files.push_back(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s_%s/EleFakeRate-DATA-Bw-ker-vtx-60-120.txt",
                                  RunYear, whichVFP.c_str(), d.c_str()));

        PlotFakeRateComparison(files, labels, colors,
                                "nVtx", 0, MAXVTX, MAXVTX,
                                MINVTX, 46, 0.07,
                                "vtx_pTcompare", RunYear, whichVFP,
                                FitFunctionType::kPol1);
    }

    // =========================================================================
    // Mode 3: fake rate vs pT, in bins of nVtx
    // Reads DATAResult<year><VFP>_nVtx<range>/EleFakeRate-DATA-Bw-ker-pt-60-120.txt
    // Fit range matches the working reference (fitFakeFunc.C): 30-200 GeV.
    // =========================================================================
    if(RUN_VS_PT_BY_NVTX){
        std::vector<std::string> dirs   = {"nVtx8-18","nVtx18-28","nVtx28-38","nVtx-gt38"};
        std::vector<std::string> labels = {"nVtx 8-18","nVtx 18-28","nVtx 28-38","nVtx >38"};
        std::vector<int> colors         = {kBlue,kRed,kGreen+2,kMagenta+1};

        std::vector<std::string> files;
        for(auto &d : dirs)
            files.push_back(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s_%s/EleFakeRate-DATA-Bw-ker-pt-60-120.txt",
                                  RunYear, whichVFP.c_str(), d.c_str()));

        const double PTMIN = 30, PTMAX = 200;
        PlotFakeRateComparison(files, labels, colors,
                                "p_{T} [GeV]", PTMIN, PTMAX, 34,
                                PTMIN, PTMAX, 0.07,
                                "pt_nVtxcompare", RunYear, whichVFP,
                                FitFunctionType::kPtDependence);
    }

    // =========================================================================
    // Mode 4: nominal fake rate vs pT, no nVtx sub-range restriction
    // (fr_bothcount_pt in fitFakeFunc.C). Reads the plain, non-cross-binned
    // DATAResult<year><VFP>/EleFakeRate-DATA-Bw-ker-pt-60-120.txt, plus DY-ker
    // and Bw-expo variants to build the systematic band on den/num, exactly
    // as fitFakeFunc.C does for pt_denerror/pt_numerror.
    //
    // NOTE: fitFakeFunc.C's own Pt_Polfile actually re-opens the DY-ker file a
    // second time (likely a copy/paste slip), unlike the eta/vtx blocks which
    // correctly use a Bw-expo variant for the third file. Here we use Bw-expo
    // for pT too, mirroring the eta/vtx pattern. If you specifically want the
    // literal (buggy) reference behavior reproduced, change polFilesPT below
    // to point at the same "-DY-ker-" file as dyFilesPT.
    // =========================================================================
    if(RUN_NOMINAL_VS_PT){
        std::vector<std::string> files  = {
            Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-ker-pt-60-120.txt",
                 RunYear, whichVFP.c_str())
        };
        std::vector<std::string> dyFilesPT  = {
            Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-DY-ker-pt-60-120.txt",
                 RunYear, whichVFP.c_str())
        };
        std::vector<std::string> polFilesPT = {
            Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-expo-pt-60-120.txt",
                 RunYear, whichVFP.c_str())
        };
        std::vector<std::string> labels = {"Data"};
        std::vector<int> colors         = {kBlue};

        const double PTMIN = 30, PTMAX = 200;
        PlotFakeRateComparison(files, labels, colors,
                                "p_{T} [GeV]", PTMIN, PTMAX, 34,
                                PTMIN, PTMAX, 0.07,
                                "pt_nominal", RunYear, whichVFP,
                                FitFunctionType::kPtDependence,
                                dyFilesPT, polFilesPT);
    }

    // =========================================================================
    // Mode 5: nominal fake rate vs nVtx, no pT sub-range restriction
    // (fr_bothcount_vtx in fitFakeFunc.C). Reads the plain, non-cross-binned
    // DATAResult<year><VFP>/EleFakeRate-DATA-Bw-ker-vtx-60-120.txt, plus
    // DY-ker and Bw-expo variants for the systematic band (matches reference).
    // =========================================================================
    if(RUN_NOMINAL_VS_NVTX){
        std::vector<std::string> files  = {
            Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-ker-vtx-60-120.txt",
                 RunYear, whichVFP.c_str())
        };
        std::vector<std::string> dyFilesVtx  = {
            Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-DY-ker-vtx-60-120.txt",
                 RunYear, whichVFP.c_str())
        };
        std::vector<std::string> polFilesVtx = {
            Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-expo-vtx-60-120.txt",
                 RunYear, whichVFP.c_str())
        };
        std::vector<std::string> labels = {"Data"};
        std::vector<int> colors         = {kBlue};

        PlotFakeRateComparison(files, labels, colors,
                                "nVtx", 0, MAXVTX, MAXVTX,
                                MINVTX, 46, 0.07,
                                "vtx_nominal", RunYear, whichVFP,
                                FitFunctionType::kPol1,
                                dyFilesVtx, polFilesVtx);
    }

    // =========================================================================
    // Mode 6: consistency check -- fake rate vs pT, nVtx sub-ranges overlaid on
    // the nominal pT dependence (5 sets: nVtx8-18, nVtx18-28, nVtx28-38,
    // nVtx>38, + nominal). Nominal drawn as a light-grey band in the
    // background; the nVtx-ranged points are drawn with no horizontal
    // bin-width or vertical error bars.
    // =========================================================================
    if(RUN_PT_CONSISTENCY){
        std::vector<std::string> dirs   = {"nVtx8-18","nVtx18-28","nVtx28-38","nVtx-gt38"};
        std::vector<std::string> labels = {"nVtx 8-18","nVtx 18-28","nVtx 28-38","nVtx >38"};
        std::vector<int> colors         = {kBlue,kRed,kGreen+2,kMagenta+1};

        std::vector<std::string> files;
        for(auto &d : dirs)
            files.push_back(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s_%s/EleFakeRate-DATA-Bw-ker-pt-60-120.txt",
                                  RunYear, whichVFP.c_str(), d.c_str()));

        std::string nominalFile = Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-ker-pt-60-120.txt",
                                        RunYear, whichVFP.c_str());
        std::string nominalDyFile = Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-DY-ker-pt-60-120.txt",
                                          RunYear, whichVFP.c_str());
        std::string nominalPolFile = Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-expo-pt-60-120.txt",
                                           RunYear, whichVFP.c_str());

        const double PTMIN = 30, PTMAX = 200;
        PlotFakeRateConsistency(files, labels, colors,
                                 nominalFile, nominalDyFile, nominalPolFile,
                                 "p_{T} [GeV]", PTMIN, PTMAX, 34,
                                 PTMIN, PTMAX, 0.07,
                                 "pt_consistency", RunYear, whichVFP,
                                 FitFunctionType::kPtDependence);
    }

    // =========================================================================
    // Mode 7: consistency check -- fake rate vs nVtx, pT sub-ranges overlaid on
    // the nominal nVtx dependence (4 sets: pT30-60, pT60-120, pT>120, +
    // nominal). Nominal drawn as a light-grey band in the background; the
    // pT-ranged points are drawn with no horizontal bin-width or vertical
    // error bars.
    // =========================================================================
    if(RUN_NVTX_CONSISTENCY){
        std::vector<std::string> dirs   = {"pT30-60","pT60-120","pT-gt120"};
        std::vector<std::string> labels = {"p_{T} 30-60","p_{T} 60-120","p_{T} >120"};
        std::vector<int> colors         = {kBlue,kRed,kGreen+2};

        std::vector<std::string> files;
        for(auto &d : dirs)
            files.push_back(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s_%s/EleFakeRate-DATA-Bw-ker-vtx-60-120.txt",
                                  RunYear, whichVFP.c_str(), d.c_str()));

        std::string nominalFile = Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-ker-vtx-60-120.txt",
                                        RunYear, whichVFP.c_str());
        std::string nominalDyFile = Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-DY-ker-vtx-60-120.txt",
                                          RunYear, whichVFP.c_str());
        std::string nominalPolFile = Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-expo-vtx-60-120.txt",
                                           RunYear, whichVFP.c_str());

        PlotFakeRateConsistency(files, labels, colors,
                                 nominalFile, nominalDyFile, nominalPolFile,
                                 "nVtx", 0, MAXVTX, MAXVTX,
                                 MINVTX, 46, 0.07,
                                 "vtx_consistency", RunYear, whichVFP,
                                 FitFunctionType::kPol1);
    }
}
