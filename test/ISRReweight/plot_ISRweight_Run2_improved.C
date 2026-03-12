// plot_ISRweight_Run2_improved.C
// ================================
// Derives ISR reweighting factors with the following improvements:
//
//   Rec.1: Separate derivation per channel (mg=MuonEG, eg=DoubleEG)
//   Rec.2: Use gen-level bosonPt (MC branch "bosonPt") as reweighting
//          variable instead of reco jet pT when useGenPt=true
//   Rec.3: Separate weight files per process (ZGToLLG, WGToLNuG)
//   Rec.4: Rare backgrounds (TTG/WWG/WZG) rejected at entry — no
//          weight derived for them
//   Rec.5: Normalization derived SEPARATELY from shape correction:
//          Step A — shape-only ratio from normalized data/MC
//          Step B — normalization = data_integral / shape-corrected_MC
//
// Branch names confirmed from ZTree.Print() on actual files:
//   Data tree:  ISRJetPt  (reco vector-sum jet pT)
//               NO bosonPt branch
//   MC tree:    JetPt     (reco vector-sum jet pT — different name!)
//               bosonPt   (gen boson pT — NOT bosonPt_gen)
//               nEvts, crosssection, MCweight, PUweight, mcType
//
// Output:
//   ISRweights_{channel}_{procTag}_{yearTag}.txt   — weight table
//   ISRweights_{channel}_{procTag}_{yearTag}.root  — histograms
//   ISRweights_{channel}_{procTag}_{yearTag}_derivation.pdf
//
// Usage:
//   root -b -q "plot_ISRweight_Run2_improved.C+(2018,false,\"mg\",\"ZGToLLG\",true)"
//   root -b -q "plot_ISRweight_Run2_improved.C+(2018,false,\"eg\",\"ZGToLLG\",true)"
//   root -b -q "plot_ISRweight_Run2_improved.C+(2018,false,\"mg\",\"WGToLNuG\",true)"
//   root -b -q "plot_ISRweight_Run2_improved.C+(2018,false,\"eg\",\"WGToLNuG\",true)"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLorentzVector.h"


#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_jet.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_fakes.h"
#include "../../include/analysis_scalefactor.h"
#include "../../include/tdrstyle.C"

void plot_ISRweight_Run2_improved(
    int         RunYear  = 2018,
    bool        preVFP   = false,
    std::string channel  = "mg",       // "mg" (MuonEG) or "eg" (DoubleEG)
    std::string procTag  = "ZGToLLG",  // "ZGToLLG" or "WGToLNuG"
    bool        useGenPt = true        // Rec.2: true=gen bosonPt, false=reco JetPt
){
    setTDRStyle();
    gROOT->SetBatch(kTRUE);
    gSystem->Load("../../lib/libAnaClasses.so");

    // ---- Rec.4: block rare backgrounds ------------------------------------
    if(procTag=="TTGJets" || procTag=="WWG" || procTag=="WZG"){
        std::cout << "[ERROR] Rec.4: ISR weights must NOT be derived for rare"
                  << " backgrounds (" << procTag << ").\n"
                  << "  Assign 100% normalization uncertainty instead." << std::endl;
        return;
    }

    // ---- year tag ----------------------------------------------------------
    std::string vfpTag = "";
    if(RunYear==2016 &&  preVFP) vfpTag = "preVFP";
    if(RunYear==2016 && !preVFP) vfpTag = "postVFP";
    std::string yearTag = std::to_string(RunYear) + vfpTag;

    std::cout << "\n=================================================" << std::endl;
    std::cout << " ISR weight derivation" << std::endl;
    std::cout << " channel=" << channel << "  proc=" << procTag
              << "  year=" << yearTag
              << "  useGenPt=" << (useGenPt?"true":"false") << std::endl;
    std::cout << "=================================================" << std::endl;

    // ---- Binning (same as original plot_ISRweight_Run2.C) ------------------
    const int    nBins   = 7;
    Double_t     ptBins[nBins+1] = {0,50,100,150,200,250,300,800};

    // ---- Output files ------------------------------------------------------
    std::string outBase = Form("ISRweights_%s_%s_%s",
                                channel.c_str(), procTag.c_str(), yearTag.c_str());
    TFile *outROOT = TFile::Open((outBase+".root").c_str(), "RECREATE");
    std::ofstream outTXT(outBase+".txt");
    outTXT << "# ISR reweight factors" << std::endl;
    outTXT << "# channel=" << channel << "  procTag=" << procTag
           << "  year="    << yearTag << "  useGenPt=" << useGenPt << std::endl;
    outTXT << "# ptLow  ptHigh  reweightF_shapeOnly  reweightF_shapeAndNorm  stat_unc" << std::endl;

    esfScaleFactor objectESF;

    // ========================================================================
    // PART 1 — Load DATA ZTree
    // Rec.1: channel-specific data file
    // Data branch for reco jet pT: "ISRJetPt"  (confirmed from ZTree.Print())
    // ========================================================================
    TH1D *h_data = new TH1D("h_data", "data ISRJetPt", nBins, ptBins);
    h_data->Sumw2();

    // Rec.1: channel-specific file name
    std::string dataFile;
    if(channel == "mg")
        dataFile = Form("/eos/uscms/store/user/tmishra/ISRweighting/"
                        "resTree_ISR_data_%s.root", yearTag.c_str());
    else
        dataFile = Form("/eos/uscms/store/user/tmishra/ISRweighting/"
                        "resTree_ISR_data_eg_%s.root", yearTag.c_str());

    TChain *dataTree = new TChain("ZTree");
    dataTree->Add(dataFile.c_str());
    std::cout << "\nData file : " << dataFile
              << "\n  entries = " << dataTree->GetEntries() << std::endl;

    // Confirmed data branch names
    float d_phoEt(0), d_phoEta(0), d_dRPhoLep(0), d_dilepMass(0), d_ISRJetPt(0);
    dataTree->SetBranchAddress("phoEt",     &d_phoEt);
    dataTree->SetBranchAddress("phoEta",    &d_phoEta);
    dataTree->SetBranchAddress("dRPhoLep",  &d_dRPhoLep);
    dataTree->SetBranchAddress("dilepMass", &d_dilepMass);
    dataTree->SetBranchAddress("ISRJetPt",  &d_ISRJetPt);  // confirmed name

    std::string frBase = "/eos/uscms/store/user/tmishra/jetfakepho/files/"
                         "JetFakeRate-ISR-";

    for(unsigned ievt=0; ievt<(unsigned)dataTree->GetEntries(); ievt++){
        dataTree->GetEntry(ievt);
        if(d_dRPhoLep  < 0.8)           continue;
        if(d_dilepMass < 80 || d_dilepMass > 100) continue;

        // j->gamma fake rate subtraction weight (same as original)
        double weight = 1.0;
        std::string frFile = frBase;
        frFile += (fabs(d_phoEta)<1.4442) ? "EB_" : "EE_";
        frFile += yearTag + ".txt";

        std::ifstream jf(frFile.c_str());
        if(jf.is_open()){
            float pl, ph, fr, er, sy, tf;
            for(int i=0; i<9; i++){
                jf >> pl >> ph >> fr >> er >> sy >> tf;
                if(std::isnan(fr) || fr<=1e-2 || fr>=1.0) fr=0.03;
                if(d_phoEt > pl && d_phoEt < ph) weight = 1.0-fr;
            }
            jf.close();
        }

        double ptVar = (d_ISRJetPt > 799) ? 799 : d_ISRJetPt;
        h_data->Fill(ptVar, weight);
    }
    std::cout << "  Data entries passing cuts: "
              << (int)h_data->GetEntries() << std::endl;

    // ========================================================================
    // PART 2 — Load MC ZTree
    // Rec.3: process-specific MC file
    // MC branch for reco jet pT: "JetPt"   (confirmed — different from data!)
    // MC branch for gen boson pT: "bosonPt" (confirmed — not "bosonPt_gen")
    // ========================================================================
    TH1D *h_mc = new TH1D("h_mc", "MC ISR pt", nBins, ptBins);
    h_mc->Sumw2();

    std::string mcFile = Form(
        "/eos/uscms/store/user/tmishra/ISRweighting/"
        "resTree_ISR_%s_%s.root", procTag.c_str(), yearTag.c_str());

    TChain *mcTree = new TChain("ZTree");
    mcTree->Add(mcFile.c_str());
    std::cout << "\nMC file   : " << mcFile
              << "\n  entries = " << mcTree->GetEntries() << std::endl;

    // Rec.2: check which branches actually exist in this MC file
    bool hasBosonPt = (mcTree->GetBranch("bosonPt")     != nullptr);
    bool hasJetPt   = (mcTree->GetBranch("JetPt")       != nullptr);
    bool hasISRJetPt= (mcTree->GetBranch("ISRJetPt")    != nullptr);

    std::cout << "  MC branches found:"
              << "  bosonPt="  << hasBosonPt
              << "  JetPt="    << hasJetPt
              << "  ISRJetPt=" << hasISRJetPt << std::endl;

    // Rec.2: even if bosonPt branch exists, verify it is not all-zeros.
    // A compression factor >>10 for a float branch means a constant (zero) fill —
    // this happens when the MC was produced with the old analysis_ISRMC.C which
    // searched MCData (status=1 only) and never found the intermediate Z/W boson.
    bool bosonPtUsable = false;
    if(hasBosonPt && useGenPt){
        float probe_bosonPt = 0;
        mcTree->SetBranchAddress("bosonPt", &probe_bosonPt);
        int nProbe = (int)std::min((Long64_t)500, mcTree->GetEntries());
        int nNonZero = 0;
        for(int ip=0; ip<nProbe; ip++){
            mcTree->GetEntry(ip);
            if(probe_bosonPt > 0) nNonZero++;
        }
        mcTree->ResetBranchAddresses();
        bosonPtUsable = (nNonZero > nProbe/5); // >20% non-zero = properly filled
        if(!bosonPtUsable){
            std::cout << "  [WARN] bosonPt branch exists but appears all-zeros "
                      << "(" << nNonZero << "/" << nProbe << " non-zero).\n"
                      << "  [WARN] MC was produced with old analysis_ISRMC.C.\n"
                      << "  [WARN] Re-run analysis_ISRMC_perprocess.C (fixed bosonPt loop)\n"
                      << "  [WARN] to regenerate this MC file with correct bosonPt.\n"
                      << "  [WARN] Falling back to reco JetPt for this run." << std::endl;
        } else {
            std::cout << "  bosonPt branch OK (" << nNonZero << "/" << nProbe
                      << " non-zero)" << std::endl;
        }
    }

    // Decide reweighting variable for MC
    // Priority: gen bosonPt (if usable) > reco JetPt > reco ISRJetPt
    bool doGenPt = (useGenPt && hasBosonPt && bosonPtUsable);
    if(useGenPt && !doGenPt)
        std::cout << "  [INFO] useGenPt requested but bosonPt not usable "
                  << "— falling back to reco pT" << std::endl;

    std::string mcPtBranch = doGenPt ? "bosonPt" : (hasJetPt ? "JetPt" : "ISRJetPt");
    std::cout << "  MC reweight variable: " << mcPtBranch << std::endl;

    // Confirmed MC branch names
    float mc_MCweight(1), mc_PUweight(1);
    float mc_phoEt(0),  mc_phoEta(0);
    float mc_lepPt(0),  mc_lepEta(0);
    float mc_dRPhoLep(0), mc_dilepMass(0);
    float mc_JetPt(0),    mc_ISRJetPt_mc(0), mc_bosonPt(0);

    mcTree->SetBranchAddress("MCweight",   &mc_MCweight);
    mcTree->SetBranchAddress("PUweight",   &mc_PUweight);
    mcTree->SetBranchAddress("phoEt",      &mc_phoEt);
    mcTree->SetBranchAddress("phoEta",     &mc_phoEta);
    mcTree->SetBranchAddress("lepPt",      &mc_lepPt);
    mcTree->SetBranchAddress("lepEta",     &mc_lepEta);
    mcTree->SetBranchAddress("dRPhoLep",   &mc_dRPhoLep);
    mcTree->SetBranchAddress("dilepMass",  &mc_dilepMass);

    if(hasJetPt)    mcTree->SetBranchAddress("JetPt",    &mc_JetPt);
    if(hasISRJetPt) mcTree->SetBranchAddress("ISRJetPt", &mc_ISRJetPt_mc);
    if(hasBosonPt)  mcTree->SetBranchAddress("bosonPt",  &mc_bosonPt);

    for(unsigned ievt=0; ievt<(unsigned)mcTree->GetEntries(); ievt++){
        mcTree->GetEntry(ievt);
        if(mc_dRPhoLep  < 0.8)                         continue;
        if(mc_dilepMass < 80 || mc_dilepMass > 100)    continue;

        double sf = objectESF.getMuonESF(mc_lepPt, mc_lepEta)
                  * objectESF.getPhotonESF(mc_phoEt, mc_phoEta)
                  * objectESF.getMuonEGTRGESF(mc_phoEt, mc_lepPt);
        double w  = mc_MCweight * mc_PUweight * sf;

        // Rec.2: use gen bosonPt if available, else reco JetPt
        double ptVar = 0;
        if(doGenPt)         ptVar = mc_bosonPt;
        else if(hasJetPt)   ptVar = mc_JetPt;
        else                ptVar = mc_ISRJetPt_mc;
        if(ptVar > 799)     ptVar = 799;

        h_mc->Fill(ptVar, w);
    }
    std::cout << "  MC entries passing cuts (weighted): "
              << h_mc->GetEntries() << std::endl;

    // ========================================================================
    // PART 3 — Rec.5: Shape correction THEN normalization separately
    // ========================================================================

    // Step A: normalize both to unit area, take bin-by-bin ratio = shape weight
    TH1D *h_data_norm = (TH1D*)h_data->Clone("h_data_norm");
    TH1D *h_mc_norm   = (TH1D*)h_mc->Clone("h_mc_norm");
    if(h_data_norm->Integral()>0) h_data_norm->Scale(1.0/h_data_norm->Integral());
    if(h_mc_norm->Integral()>0)   h_mc_norm->Scale(1.0/h_mc_norm->Integral());

    TH1D *h_shape = (TH1D*)h_data_norm->Clone("h_shapeRatio");
    h_shape->SetTitle("Shape-only weight (norm. data / norm. MC)");
    h_shape->Divide(h_mc_norm);

    // protect against zero/negative MC bins
    for(int ib=1; ib<=nBins; ib++){
        if(h_shape->GetBinContent(ib) <= 0){
            std::cout << "  [WARN] bin " << ib << " shape weight <= 0, setting to 1" << std::endl;
            h_shape->SetBinContent(ib, 1.0);
            h_shape->SetBinError(ib, 0.5);
        }
    }

    // Step B: apply shape weights to MC, then derive normalization
    // Rec.5: normFactor is derived AFTER shape correction, kept separate
    TH1D *h_mc_afterShape = (TH1D*)h_mc->Clone("h_mc_afterShape");
    for(int ib=1; ib<=nBins; ib++)
        h_mc_afterShape->SetBinContent(ib,
            h_mc->GetBinContent(ib) * h_shape->GetBinContent(ib));

    double normFactor = (h_mc_afterShape->Integral() > 0)
                      ? h_data->Integral() / h_mc_afterShape->Integral()
                      : 1.0;
    std::cout << "\n  Rec.5 — normalization factor (separate from shape): "
              << normFactor << std::endl;

    // Combined (shape * norm) histogram
    TH1D *h_total = (TH1D*)h_shape->Clone("h_totalRatio");
    h_total->SetTitle("Shape+Norm weight");
    h_total->Scale(normFactor);

    // ========================================================================
    // PART 4 — Write output text file
    // ========================================================================
    outTXT << "# normFactor (Rec.5, separate from shape): " << normFactor << std::endl;

    std::cout << "\n  Bin weights:" << std::endl;
    std::cout << Form("  %-8s %-8s %-20s %-20s %-10s",
                      "ptLow","ptHigh","wShape","wTotal","stat_unc") << std::endl;
    std::cout << "  " << std::string(68,'-') << std::endl;

    for(int ib=1; ib<=nBins; ib++){
        double ptLow  = h_shape->GetBinLowEdge(ib);
        double ptHigh = h_shape->GetBinLowEdge(ib+1);
        double wShape = h_shape->GetBinContent(ib);
        double wTotal = h_total->GetBinContent(ib);
        double unc    = h_shape->GetBinError(ib);

        std::cout << Form("  %-8.0f %-8.0f %-20.6f %-20.6f %-10.6f",
                          ptLow, ptHigh, wShape, wTotal, unc) << std::endl;
        outTXT    << ptLow  << "  " << ptHigh << "  "
                  << wShape << "  " << wTotal << "  " << unc << std::endl;
    }
    outTXT << "normFactor  " << normFactor << std::endl;
    outTXT.close();
    std::cout << "\n  Written: " << outBase << ".txt" << std::endl;

    // ========================================================================
    // PART 5 — Diagnostic plots
    // ========================================================================
    TCanvas *c1 = new TCanvas("c1","ISR derivation",1400,500);
    c1->Divide(3);

    // Panel 1: data vs MC (MC scaled to data integral)
    c1->cd(1);
    gPad->SetLogy();
    TH1D *h_mc_scaled = (TH1D*)h_mc->Clone("h_mc_scaled");
    if(h_mc_scaled->Integral()>0)
        h_mc_scaled->Scale(h_data->Integral()/h_mc_scaled->Integral());
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.9);
    h_data->SetTitle(Form("%s %s %s;pT (GeV);Events",
                          channel.c_str(), procTag.c_str(), yearTag.c_str()));
    h_data->Draw("EP");
    h_mc_scaled->SetLineColor(kRed);
    h_mc_scaled->SetLineWidth(2);
    h_mc_scaled->Draw("hist same");
    TLegend *leg1 = new TLegend(0.45,0.72,0.88,0.88);
    leg1->AddEntry(h_data,     "data (j#rightarrow#gamma subtracted)","ep");
    leg1->AddEntry(h_mc_scaled,Form("MC %s (norm to data)",procTag.c_str()),"l");
    leg1->SetTextSize(0.038);
    leg1->Draw();
    TLatex lat;
    lat.SetNDC(); lat.SetTextSize(0.035);
    lat.DrawLatex(0.15,0.92,Form("Rec.2: var=%s", mcPtBranch.c_str()));

    // Panel 2: shape-only weight
    c1->cd(2);
    h_shape->SetTitle("Shape-only weight (Rec.5 Step A);pT (GeV);w_{shape}");
    h_shape->SetLineColor(kBlue);
    h_shape->SetMarkerStyle(20);
    h_shape->SetMarkerColor(kBlue);
    h_shape->GetYaxis()->SetRangeUser(0.3, 2.0);
    h_shape->Draw("EP");
    TLine *l1 = new TLine(ptBins[0],1,ptBins[nBins],1);
    l1->SetLineColor(kRed); l1->SetLineStyle(2); l1->SetLineWidth(2);
    l1->Draw();
    lat.DrawLatex(0.15,0.92,"Rec.5: shape only (norm. removed)");

    // Panel 3: shape+norm weight
    c1->cd(3);
    h_total->SetTitle(Form("Shape+Norm weight (normF=%.4f);pT (GeV);w_{total}",
                           normFactor));
    h_total->SetLineColor(kGreen+2);
    h_total->SetMarkerStyle(20);
    h_total->SetMarkerColor(kGreen+2);
    h_total->GetYaxis()->SetRangeUser(0.3, 2.0);
    h_total->Draw("EP");
    TLine *l2 = new TLine(ptBins[0],1,ptBins[nBins],1);
    l2->SetLineColor(kRed); l2->SetLineStyle(2); l2->SetLineWidth(2);
    l2->Draw();
    lat.DrawLatex(0.15,0.92,Form("Rec.5: shape #times norm (%.4f)",normFactor));

    c1->SaveAs((outBase+"_derivation.pdf").c_str());
    std::cout << "  Written: " << outBase << "_derivation.pdf" << std::endl;

    // Save histograms
    outROOT->cd();
    h_data->Write("h_data");
    h_mc->Write("h_mc");
    h_mc_scaled->Write("h_mc_scaledToData");
    h_shape->Write("h_shapeWeight");
    h_total->Write("h_totalWeight");
    outROOT->Close();
    std::cout << "  Written: " << outBase << ".root\n" << std::endl;

    std::cout << "Rec.4 reminder: TTGJets/WWG/WZG get ISRWeight=1 "
              << "with 100% uncertainty in apply_ISRWeight.C" << std::endl;
}
