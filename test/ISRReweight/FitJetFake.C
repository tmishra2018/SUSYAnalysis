#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <ctime>
#include <cmath>

#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TH2F.h"
#include "TMath.h"
#include "TSystem.h"
#include "TFractionFitter.h"
#include "TError.h"
#include "TStyle.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooAbsReal.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"

#include "../../include/tdrstyle.C"
#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_jet.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_fakes.h"

bool useMC     = true;
bool doIterate = true;

// Data eras (same as run_Reweight_step1.sh) merged into one fake-rate era
// (2017: MuonEG lumi covers D-F only, DoubleEG lumi covers B-F)
static std::vector<std::string> getEras(int RunYear, bool preVFP, const std::string& channel){
    if(RunYear == 2016 && preVFP)  return {"B-ver1","B-ver2","C","D","E","F_preVFP"};
    if(RunYear == 2016 && !preVFP) return {"F_postVFP","G","H"};
    if(RunYear == 2017)            return (channel == "mg") ? std::vector<std::string>{"D","E","F"}
                                                            : std::vector<std::string>{"B","C","D","E","F"};
    return {"A","B","C","D"};
}

// ---------------------------------------------------------------------------
// Returns 0 on success, 1 on empty input: root -q exits with this value,
// and runFit_step2.sh (set -e) stops on a non-zero exit.
int FitJetFake(float lowercut, float uppercut,
               int detType, int RunYear, bool preVFP,
               const char* channelArg = "mg")
{
    setTDRStyle();
    time_t now = time(0);
    // TFractionFitter repeats "Did not find solution for t_i" per bin and
    // iteration; at Warning level step2.log reached ~0.5 GB and filled the
    // home quota. Keep errors only.
    gErrorIgnoreLevel = kError;

    const std::string channel = channelArg;  // "mg" or "egamma"
    const std::string det     = (detType == 1) ? "EB" : "EE";
    const std::string vfpTag  = (RunYear == 2016) ? (preVFP ? "preVFP" : "postVFP") : "";
    const std::string eraTag  = Form("%d%s", RunYear, vfpTag.c_str());

    std::cout << "[FitJetFake] pT=[" << lowercut << "," << uppercut << ")"
              << "  det=" << det << "  channel=" << channel
              << "  era=" << eraTag << "\n";

    gSystem->Load("../../lib/libAnaClasses.so");

    // ---- Output fake rate text file (channel + detector + era tagged) ------
    //   JetFakeRate-ISR-{EB|EE}_{mg|egamma}_{year}{VFP}.txt, read by plot_ISRweight.C
    // One line per pT bin, appended; runFit_step2.sh removes the file first.
    // Columns: ptLow ptHigh fakerate total_error syst_error mc_fakerate
    std::string outname = Form("/eos/uscms/store/user/tmishra/jetfakepho/files/JetFakeRate-ISR-%s_%s_%s.txt",
                               det.c_str(), channel.c_str(), eraTag.c_str());

    std::ofstream myfile(outname.c_str(), std::ios_base::app);

    // ---- pT bin label strings ---------------------------------------------
    char lowername[16], uppername[16];
    sprintf(lowername, "%d", (int)lowercut);
    if(uppercut < 1000) sprintf(uppername, "%d", (int)uppercut);
    else                sprintf(uppername, "Inf");

    // ---- sigma_ietaieta sideband cuts (EB / EE) ---------------------------
    double StandardCut(0), StandardIso(0);
    // Fall17v2 cut-based loose photon ID (AN Sec. 3.1, src/analysis_photon.cc),
    // same as ../jetFakePho/fit/FitJetFake.C:
    //   EB sigmaIetaIeta < 0.0106, Iso_h+- < 1.694;  EE 0.0272, 2.089
    Double_t SigmaCutLower_EB[] = {0.0106,0.0107,0.0108,0.0109,0.0110,
                                   0.0111,0.0112,0.0113,0.0114,0.0115};
    Double_t SigmaCutUpper_EB[] = {0.0140,0.0145,0.0150,0.0155,0.0160,
                                   0.0165,0.0170,0.0175,0.0180,0.0185};
    Double_t SigmaCutLower_EE[] = {0.0272,0.0273,0.0274,0.0275,0.0276,
                                   0.0277,0.0278,0.0279,0.0280,0.0281};
    Double_t SigmaCutUpper_EE[] = {0.035,0.036,0.037,0.038,0.039,0.04};

    std::vector<double> SigmaCutLower, SigmaCutUpper;
    if(detType == 1){
        StandardCut = 0.0106; StandardIso = 1.694;
        for(auto v : SigmaCutLower_EB) SigmaCutLower.push_back(v);
        for(auto v : SigmaCutUpper_EB) SigmaCutUpper.push_back(v);
    } else {
        StandardCut = 0.0272; StandardIso = 2.089;
        for(auto v : SigmaCutLower_EE) SigmaCutLower.push_back(v);
        for(auto v : SigmaCutUpper_EE) SigmaCutUpper.push_back(v);
    }
    unsigned nUpper = SigmaCutUpper.size();
    unsigned nLower = SigmaCutLower.size();

    // ---- Histograms -------------------------------------------------------
    std::ostringstream hname;
    hname.str(""); hname << "ChIsoTar-pt-" << lowername << "-" << uppername;
    TH1D* h_target = new TH1D(hname.str().c_str(), ";Iso_{h^{#pm}} (GeV);", 20, 0.0, 20);

    hname.str(""); hname << "mcChIso-pt-" << lowername << "-" << uppername;
    TH1D* mc_sig   = new TH1D(hname.str().c_str(), ";Iso_{h^{#pm}} (GeV);", 20, 0.0, 20);

    hname.str(""); hname << "mcChIso-fake-pt-" << lowername << "-" << uppername;
    TH1D* mc_fake  = new TH1D(hname.str().c_str(), ";Iso_{h^{#pm}} (GeV);", 20, 0.0, 20);

    TH1D* h_bg[nLower][nUpper];
    TH1D* mc_sbcontamination[nLower][nUpper];
    float templateCorrFactor[nLower][nUpper];
    TObjArray* templateContainer[nLower][nUpper];
    TFractionFitter* fitter[nLower][nUpper];
    TH1D* result[nLower][nUpper];
    TH1* mc_predict[nLower][nUpper];

    for(unsigned iU(0); iU < nUpper; ++iU){
        for(unsigned iL(0); iL < nLower; ++iL){
            templateContainer[iL][iU] = new TObjArray(2);
            hname.str("");
            hname << "ChIsoHad-pt-" << lowername << "-" << uppername
                  << "_" << SigmaCutLower[iL] << "_" << SigmaCutUpper[iU];
            h_bg[iL][iU] = new TH1D(hname.str().c_str(), hname.str().c_str(), 20, 0.0, 20);
            hname.str("");
            hname << "mcChIso-sideband-pt-" << lowername << "-" << uppername
                  << "_" << SigmaCutLower[iL] << "_" << SigmaCutUpper[iU];
            mc_sbcontamination[iL][iU] = new TH1D(hname.str().c_str(), hname.str().c_str(), 20, 0.0, 20);
        }
    }

    // ---- MC signal template (GJet MC) -------------------------------------
    // True-photon ChIso template, same for both channels;
    // made with ../jetFakePho/src/analysis_mgGJet.C
    std::string mcGJetPath = Form(
        "/eos/uscms/store/user/tmishra/jetfakepho/files/"
        "plot_hadron_mgGJet_%s.root",
        eraTag.c_str());

    TChain* mctree = new TChain("egTree");
    mctree->Add(mcGJetPath.c_str());

    float mc_phoEt(0), mc_phoEta(0), mc_phoPhi(0),
          mc_phoSigma(0), mc_phoChIso(0);
    std::vector<int>   *mc_mcPID=0, *mc_mcMomPID=0;
    std::vector<float> *mc_mcEta=0, *mc_mcPhi=0, *mc_mcPt=0;

    mctree->SetBranchAddress("phoEt",    &mc_phoEt);
    mctree->SetBranchAddress("phoEta",   &mc_phoEta);
    mctree->SetBranchAddress("phoPhi",   &mc_phoPhi);
    mctree->SetBranchAddress("phoSigma", &mc_phoSigma);
    mctree->SetBranchAddress("phoChIso", &mc_phoChIso);
    if(useMC){
        // GJet MC may not have mcPID branches in the hadron proxy tree;
        // guard with TBranch existence check
        if(mctree->GetBranch("mcPID")){
            mctree->SetBranchAddress("mcPID",    &mc_mcPID);
            mctree->SetBranchAddress("mcEta",    &mc_mcEta);
            mctree->SetBranchAddress("mcPhi",    &mc_mcPhi);
            mctree->SetBranchAddress("mcPt",     &mc_mcPt);
            mctree->SetBranchAddress("mcMomPID", &mc_mcMomPID);
        } else {
            useMC = false;  // fallback: no MC truth available
        }
    }

    // ---- Data target region (signalTree from analysis_Hadron) -------------
    // signalTree contains phoEt, phoSigma, phoChIso, dilepMass in the
    // Z(ll)+gamma control region (80 < m(ll) < 100 GeV).
    // Per-era files written by analysis_Hadron.C, chained here.
    TChain* egtree = new TChain("signalTree");
    std::string dataPath;
    for(const auto& era : getEras(RunYear, preVFP, channel)){
        std::string f = Form(
            "/eos/uscms/store/user/tmishra/jetfakepho/files/"
            "plot_hadron_ISR_%s_%d%s.root",
            channel.c_str(), RunYear, era.c_str());
        egtree->Add(f.c_str());
        dataPath += f + " ";
    }
    if(mctree->GetEntries() == 0 || egtree->GetEntries() == 0){
        std::cerr << "[FitJetFake] ERROR: empty input (MC " << mctree->GetEntries()
                  << ", data " << egtree->GetEntries() << ") for " << channel
                  << " " << eraTag << "\n";
        return 1;
    }

    float phoEt(0), phoEta(0), phoSigma(0), phoChIso(0), llmass(0);
    egtree->SetBranchAddress("phoEt",     &phoEt);
    egtree->SetBranchAddress("phoEta",    &phoEta);
    egtree->SetBranchAddress("phoSigma",  &phoSigma);
    egtree->SetBranchAddress("phoChIso",  &phoChIso);
    egtree->SetBranchAddress("dilepMass", &llmass);

    // ---- Fill MC histograms -----------------------------------------------
    double nMCtarget(0), nMCfaketarget(0), nsignalregion(0);
    double nsideband[nLower][nUpper], nsbPassIso[nLower][nUpper];
    for(unsigned ii(0); ii < nLower; ++ii)
        for(unsigned jj(0); jj < nUpper; ++jj)
            nsideband[ii][jj] = nsbPassIso[ii][jj] = 0;

    for(unsigned ievt(0); ievt < (unsigned)mctree->GetEntries(); ++ievt){
        mctree->GetEntry(ievt);
        if(mc_phoEt < lowercut || mc_phoEt > uppercut) continue;
        if(detType == 1 && fabs(mc_phoEta) > 1.4442) continue;
        else if(detType == 2 && (fabs(mc_phoEta) < 1.56 || fabs(mc_phoEta) > 2.4)) continue;

        bool isFake(false);
        if(useMC && mc_mcPID && !mc_mcPID->empty()){
            unsigned mcIdx(0); float mindR(0.7);
            for(unsigned ii(0); ii < mc_mcPID->size(); ++ii){
                float dR = DeltaR(mc_phoEta, mc_phoPhi, (*mc_mcEta)[ii], (*mc_mcPhi)[ii]);
                float dE = fabs(mc_phoEt - (*mc_mcPt)[ii]) / mc_phoEt;
                if(dR < mindR && dE < 0.7){ mindR = dR; mcIdx = ii; }
            }
            if(mindR < 0.7)
                isFake = isHad(fabs((*mc_mcPID)[mcIdx]), fabs((*mc_mcMomPID)[mcIdx]));
        }

        if(mc_phoSigma <= StandardCut && isFake)  mc_fake->Fill(mc_phoChIso);
        if(mc_phoSigma <= StandardCut && !isFake){
            mc_sig->Fill(mc_phoChIso);
            nsignalregion += 1;
            if(mc_phoChIso < StandardIso) nMCtarget += 1;
        }
        if(mc_phoSigma <= StandardCut && isFake && mc_phoChIso < StandardIso)
            nMCfaketarget += 1;

        for(unsigned iU(0); iU < nUpper; ++iU)
            for(unsigned iL(0); iL < nLower; ++iL)
                if(mc_phoSigma > SigmaCutLower[iL] && mc_phoSigma < SigmaCutUpper[iU] && !isFake){
                    mc_sbcontamination[iL][iU]->Fill(mc_phoChIso);
                    nsideband[iL][iU] += 1;
                    if(mc_phoChIso < StandardIso) nsbPassIso[iL][iU] += 1;
                }
    }

    double mcfakerate = (nMCtarget + nMCfaketarget) > 0
                      ? nMCfaketarget / (nMCtarget + nMCfaketarget)
                      : 0.0;

    // ---- Fill data histograms --------------------------------------------
    double nden(0), totalden(0);
    double nnum[nLower][nUpper], totalnum[nLower][nUpper];
    for(unsigned ii(0); ii < nLower; ++ii)
        for(unsigned jj(0); jj < nUpper; ++jj)
            nnum[ii][jj] = totalnum[ii][jj] = 0;

    for(unsigned ievt(0); ievt < (unsigned)egtree->GetEntries(); ++ievt){
        egtree->GetEntry(ievt);
        if(llmass < 80 || llmass > 100) continue;
        if(phoEt < lowercut || phoEt > uppercut) continue;
        if(detType == 1 && fabs(phoEta) > 1.4442) continue;
        else if(detType == 2 && (fabs(phoEta) < 1.56 || fabs(phoEta) > 2.4)) continue;

        if(phoSigma <= StandardCut){
            h_target->Fill(phoChIso);
            totalden += 1;
            if(phoChIso < StandardIso) nden += 1;
        } else {
            for(unsigned iU(0); iU < nUpper; ++iU)
                for(unsigned iL(0); iL < nLower; ++iL)
                    if(phoSigma > SigmaCutLower[iL] && phoSigma < SigmaCutUpper[iU]){
                        h_bg[iL][iU]->Fill(phoChIso);
                        totalnum[iL][iU] += 1;
                        if(phoChIso < StandardIso) nnum[iL][iU] += 1;
                    }
        }
    }

    // ---- Template fit (same iterative TFractionFitter as original) --------
    double hadfrac(0), fittingError(0);
    double lowestfrac(1), highestfrac(0);

    TH1D* fracHad1D = new TH1D(Form("fracHad1D_%s_%s",lowername,uppername),
                                "", 500, 0, 1.0);

    for(unsigned iU(0); iU < nUpper; ++iU){
        for(unsigned iL(0); iL < nLower; ++iL){
            double iterUnc(1), lastPurity(1), phoPurity(0);
            double fracPassSigma = (mc_sbcontamination[iL][iU]->GetEntries() > 0)
                                 ? nsbPassIso[iL][iU] / mc_sbcontamination[iL][iU]->GetEntries()
                                 : 0;
            templateCorrFactor[iL][iU] = 0.0;
            fitter[iL][iU] = nullptr;

            // the purity can oscillate in low-stat bins (EE 200-250 never
            // converged and grew to 11 GB), so cap the iterations
            const int maxIter = 50;
            int nIter = 0;
            while(iterUnc > 0.01 && nIter++ < maxIter){
                TH1D* hist_new = (TH1D*)h_bg[iL][iU]->Clone();
                if(mc_sbcontamination[iL][iU]->GetEntries() > 0)
                    hist_new->Add(mc_sbcontamination[iL][iU],
                                  -1.0*templateCorrFactor[iL][iU]/mc_sbcontamination[iL][iU]->GetEntries());
                hist_new->Sumw2();
                // subtracting the signal contamination can leave negative bins,
                // which TFractionFitter rejects (CheckConsistency error)
                for(int ib = 1; ib <= hist_new->GetNbinsX(); ++ib)
                    if(hist_new->GetBinContent(ib) < 0) hist_new->SetBinContent(ib, 0);
                templateContainer[iL][iU]->Add(mc_sig);
                templateContainer[iL][iU]->Add(hist_new);

                delete fitter[iL][iU];
                fitter[iL][iU] = new TFractionFitter(h_target, templateContainer[iL][iU], "Q");
                fitter[iL][iU]->Constrain(1, 0.0, 1.0);
                Int_t status = fitter[iL][iU]->Fit();
                double tmpSig(0), tmpErr(0);
                if(status == 0){
                    fitter[iL][iU]->GetResult(0, tmpSig, tmpErr);
                    phoPurity = tmpSig;
                }
                iterUnc  = fabs(lastPurity - phoPurity) / std::max(phoPurity, 1e-9);
                lastPurity = phoPurity;
                templateCorrFactor[iL][iU] = totalden * phoPurity
                                           * nsideband[iL][iU]
                                           / std::max(nsignalregion, 1.0);
                templateContainer[iL][iU]->Clear();
                hist_new->Delete();
            }
            // a non-converged sideband has no meaningful purity: leave it out of
            // the result and of the sideband spread (systematic)
            if(iterUnc > 0.01){
                std::cerr << "[FitJetFake] WARNING: purity not converged after " << maxIter
                          << " iterations (sideband " << iL << "," << iU << ", last change "
                          << iterUnc << "), skipped\n";
                continue;
            }

            double fracBkg(0), errBkg(0), fracSig(0), errSig(0);
            fitter[iL][iU]->GetResult(0, fracSig, errSig);
            fitter[iL][iU]->GetResult(1, fracBkg, errBkg);

            double den    = (totalden > 0) ? nden / totalden : 0;
            double denerr = (nden > 0 && totalden > 0)
                          ? (1./sqrt(nden) + 1./sqrt(totalden)) * den : 0;
            double denom  = totalnum[iL][iU] - templateCorrFactor[iL][iU];
            double num    = (denom > 0)
                          ? (nnum[iL][iU] - templateCorrFactor[iL][iU]*fracPassSigma) / denom
                          : 0;
            double numerr = (nnum[iL][iU] > 0 && denom > 0)
                          ? (1./sqrt(nnum[iL][iU]) + 1./sqrt(denom)) * num : 0;

            double thisHadFrac = (den > 0 && fracBkg > 0) ? num * fracBkg / den : 0;

            if(SigmaCutLower[iL] == StandardCut && iU == 0){
                hadfrac      = thisHadFrac;
                fittingError = (num > 0 && den > 0 && fracBkg > 0)
                             ? (numerr/num + denerr/den + errBkg/fracBkg) * thisHadFrac
                             : 0;
            }
            if(fracBkg > 0 && fracBkg < 1) fracHad1D->Fill(thisHadFrac);
            if(thisHadFrac < lowestfrac)   lowestfrac  = thisHadFrac;
            if(thisHadFrac > highestfrac)  highestfrac = thisHadFrac;
        }
    }

    if(hadfrac < 1e-4) hadfrac = fracHad1D->GetMean();
    double systErr = highestfrac - lowestfrac;
    double totalErr = sqrt(fittingError*fittingError + systErr*systErr);

    // ---- Write to output file --------------------------------------------
    // numeric columns only: plot_ISRweight.C parses six numbers per line
    myfile << lowercut << " " << uppercut << " "
           << hadfrac   << " " << totalErr  << " " << systErr << " "
           << mcfakerate << "\n";
    myfile.close();

    std::cout << "[FitJetFake] input data: " << dataPath << " (" << ctime(&now) << ")";
    std::cout << "[FitJetFake] " << channel << " " << det << " " << eraTag
              << " pT=[" << lowername << "," << uppername << ")"
              << "  fakerate=" << hadfrac
              << " +/- " << totalErr << "\n";

    delete fracHad1D;
    delete mctree;
    delete egtree;
    return 0;
}
