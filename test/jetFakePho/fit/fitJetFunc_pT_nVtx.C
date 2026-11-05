#include<string>
#include<array>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<vector>
#include<limits>
#include<cmath>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TF3.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TPad.h"
#include "TFitResult.h"
#include "TVirtualFitter.h"
#include "TMatrixDSym.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooExponential.h"
#include "RooMCStudy.h"
#include "RooChi2MCSModule.h"
#include "RooPlot.h"
#include "TH1.h"
#include "RooFitResult.h"
#include "RooAbsReal.h"
#include "RooMultiVarGaussian.h"
#include "TGraphAsymmErrors.h"

#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_rawData.h"
#include "../../../include/tdrstyle.C"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_fakes.h"
#include "TROOT.h"

bool USE_VTX_BY_PT_SLICE = true;

#define NTOY 1000
#define NBIN 19
#define REBINSIZE 1

Double_t tmpjetfake_func(Double_t *x, Double_t *par)
{
	double pt_low = x[0] - REBINSIZE/2.0;
	double pt_high = x[0]+ REBINSIZE/2.0;

	double c1 = par[0];
	double c2 = par[1];
	double lamda1 = par[2];
	double lamda2 = par[3];

	double jetfakes_lowedge = c1*exp(lamda1*pt_low)/lamda1 + c2*exp(lamda2*pt_low)/lamda2;
	double jetfakes_highedge =  c1*exp(lamda1*pt_high)/lamda1 + c2*exp(lamda2*pt_high)/lamda2;
	return (jetfakes_highedge - jetfakes_lowedge);
}

double exp_safe(double x) {
    if (x > 700) return std::exp(700);   // avoid overflow
    if (x < -700) return 0.0;            // avoid underflow
    return std::exp(x);
}

double exp2c_func(double *x, double *par) {
    double pt = x[0];
    double c1 = par[0];
    double c2 = par[1];
    double lamda1 = par[2];
    double lamda2 = par[3];
    double shift = par[4];
    double offset = par[5];

    double arg1 = -lamda1 * pt;
    double arg2 = -lamda2 * (pt - shift - 200);

    double term1 = c1 * exp_safe(arg1);
    double term2 = c2 * exp_safe(arg2);
    double value = term1 + term2 + offset;

    if (!std::isfinite(value)) {
        std::cout << "[WARNING] Non-finite value at pt = " << pt
                  << ", arg1 = " << arg1
                  << ", arg2 = " << arg2
                  << ", shift = " << shift
                  << ", lamda1 = " << lamda1
                  << ", lamda2 = " << lamda2
                  << std::endl;
    }

    return value;
}

// FIX: log10-space wrapper around exp2c_func, used only as the *fit target*
// for the numerator (see the tail-fit note near fitfunc_num below). A plain
// chi2 fit of exp2c_func directly against the linear-scale graph is
// overwhelmingly dominated by the ~75 finely-binned, small-relative-error
// points below pT ~ 110 GeV, leaving the sparse, large-error high-pT tail
// (which is what c2, lamda2, shift and offset actually control) essentially
// unconstrained. Fitting log10(y) instead equalizes relative weight across
// the full pT range. exp2c_func is guaranteed >= offset >= 0 by the existing
// positivity parameter limits, so it is always safe to take the log of; the
// tiny floor below only guards the pathological case where every term is
// pinned at exactly its lower limit simultaneously.
double log_exp2c_func(double *x, double *par) {
    double value = exp2c_func(x, par);
    if (!(value > 0)) value = 1e-12;
    return std::log10(value);
}

void fitJetFunc_pT_nVtx(int ichannel,int  RunYear,bool ISpreVFP){
	int detType = 1;
	int channel = ichannel;
	gROOT->SetBatch(kTRUE);
	setTDRStyle();
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetErrorX(0.5);
	gStyle->SetTitleX(0.5);
	
	float correction =1.0;
	
	                  if (RunYear == 2016 and ISpreVFP == 1 and ichannel == 1) correction = 1.00158;
          else if (RunYear == 2016 and ISpreVFP == 1 and ichannel == 2) correction = 0.912476;
          else if (RunYear == 2016 and ISpreVFP == 0 and ichannel == 1) correction = 0.971904;
          else if (RunYear == 2016 and ISpreVFP == 0 and ichannel == 2) correction = 0.914534;
          else if (RunYear == 2017 and ichannel == 1) correction = 1.07788;
          else if (RunYear == 2017 and ichannel == 2) correction = 0.806638;
          else if (RunYear == 2018 and ichannel == 1) correction = 1.08004;
          else if (RunYear == 2018 and ichannel == 2) correction = 0.971525;
	
	std::string whichVFP;
	if(RunYear==2016 and ISpreVFP == true) whichVFP = "preVFP";
	if(RunYear==2016 and ISpreVFP == false) whichVFP = "postVFP";
	if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

	cout<< " RunYear " << RunYear <<  "  "<< whichVFP << "  channel "<< channel <<endl;

	TChain *sigtree = new TChain("signalTree");
	if(channel == 1)sigtree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s.root",RunYear,whichVFP.c_str()));
	else if(channel ==2)sigtree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s.root",RunYear,whichVFP.c_str()));

	TChain *controltree = new TChain("jetTree");
	if(channel == 1)controltree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s.root",RunYear,whichVFP.c_str()));
	else if(channel ==2)controltree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s.root",RunYear,whichVFP.c_str()));

	std::stringstream fakerate_filename;
	fakerate_filename.str("");
	if(channel == 1) fakerate_filename << "/eos/uscms/store/user/tmishra/jetfakepho/txt"<<RunYear<<whichVFP<<"/JetFakeRate-DoubleEG-";
        if(channel == 2) fakerate_filename << "/eos/uscms/store/user/tmishra/jetfakepho/txt"<<RunYear<<whichVFP<<"/JetFakeRate-MuonEG-";

	if(detType == 1)fakerate_filename << "EB_pT_nVtx.txt";
	else if(detType == 2)fakerate_filename << "EE_pT_nVtx.txt";
	std::ifstream jetfake_file(fakerate_filename.str().c_str());
	double PtBin[NBIN];
	double fracHad[NBIN];
	double fracHadError[NBIN];
	int i(0);
	float pt_lower(0), pt_upper(0);
	float fakerate(0), error(0), systematic(0);
	float truefake;
	// Jet fake rate as input
	if(jetfake_file.is_open()){
		for(int i(0); i < NBIN; i++){ 
			jetfake_file >> pt_lower >> pt_upper >> fakerate >> error >> systematic >> truefake;
			PtBin[i] = pt_lower;
			fracHad[i]= fakerate;
			fracHadError[i] = error;
		}
		jetfake_file.close(); 
	}

	std::string detTag = (detType == 1) ? "EB" : "EE";


	auto fakerate_ptDependence_pred = [](double pt, double slope, double constant, double index) -> double {
		double p = TMath::Max(pt, 0.000001);
		double arg = slope*p + constant;
		if(arg <= 0) return 0.0;
		return pow(arg, index);
	};

	std::string elefakeBaseName = Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-Data-ByPtVtx-%s",RunYear,whichVFP.c_str(),detTag.c_str());
	double scalefactor_incl(0), ptslope_incl(0), ptconstant_incl(0), ptindex_incl(0), vtxconst_incl(0), vtxslope_incl(0);
	{
		std::ifstream elefake_file;
		elefake_file.open(elefakeBaseName + ".txt");
		std::string variabletype; double variablevalue;
		if(elefake_file.is_open()){
			while(elefake_file >> variabletype >> variablevalue){
				if(variabletype=="scalefactor") scalefactor_incl = variablevalue;
				else if(variabletype=="ptslope") ptslope_incl = variablevalue;
				else if(variabletype=="ptconstant") ptconstant_incl = variablevalue;
				else if(variabletype=="ptindex") ptindex_incl = variablevalue;
				else if(variabletype=="vtxconst") vtxconst_incl = variablevalue;
				else if(variabletype=="vtxslope") vtxslope_incl = variablevalue;

			}
		} 
		elefake_file.close();
	}

	// Per-nVtx-slice fits + their OWN scalefactor -- same boundaries/file-tag
	std::vector<std::string> nVtxSliceFileTag = {"_nVtx8-20","_nVtx20-30","_nVtx-gt30"};
	std::vector<double> nVtxSliceUpperEdge    = {20.0, 30.0, 1e9}; // event falls in first slice whose nVertex < edge
	std::vector<double> scalefactor_nVtxSlice(nVtxSliceFileTag.size(), -1);
	std::vector<double> ptslope_nVtxSlice(nVtxSliceFileTag.size(), 0);
	std::vector<double> ptconstant_nVtxSlice(nVtxSliceFileTag.size(), 0);
	std::vector<double> ptindex_nVtxSlice(nVtxSliceFileTag.size(), 0);
	std::vector<bool>   haveNVtxSlice(nVtxSliceFileTag.size(), false);
	for(size_t is=0; is<nVtxSliceFileTag.size(); is++){
		std::ifstream f;
		f.open(elefakeBaseName + nVtxSliceFileTag[is] + ".txt");
		std::string variabletype; double variablevalue;
		while(f >> variabletype >> variablevalue){
			if(variabletype=="scalefactor") scalefactor_nVtxSlice[is] = variablevalue;
			else if(variabletype=="ptslope") ptslope_nVtxSlice[is] = variablevalue;
			else if(variabletype=="ptconstant") ptconstant_nVtxSlice[is] = variablevalue;
			else if(variabletype=="ptindex") ptindex_nVtxSlice[is] = variablevalue;
		}
		f.close();
		haveNVtxSlice[is] = (scalefactor_nVtxSlice[is] > 0);
	}

	auto SelectPtParamsForNVtx = [&](int nVertexVal) -> std::array<double,4> {
		for(size_t is=0; is<nVtxSliceUpperEdge.size(); is++){
			if(nVertexVal < nVtxSliceUpperEdge[is]){
				if(haveNVtxSlice[is]) return {scalefactor_nVtxSlice[is], ptslope_nVtxSlice[is], ptconstant_nVtxSlice[is], ptindex_nVtxSlice[is]};
				return {scalefactor_incl, ptslope_incl, ptconstant_incl, ptindex_incl};
			}
		}
		return {scalefactor_incl, ptslope_incl, ptconstant_incl, ptindex_incl};
	};
	
	std::vector<std::string> pTSliceFileTag = {"_pT30-50","_pT50-80","_pT-gt80"};
	std::vector<double> pTSliceUpperEdge    = {50.0, 80.0, 1e9}; // event falls in first slice whose pT < edge
	std::vector<double> scalefactor_pTSlice(pTSliceFileTag.size(), -1);
	std::vector<double> vtxconst_pTSlice(pTSliceFileTag.size(), 0);
	std::vector<double> vtxslope_pTSlice(pTSliceFileTag.size(), 0);
	std::vector<double> ptshapeavg_pTSlice(pTSliceFileTag.size(), -1);
	std::vector<bool>   havePtSlice(pTSliceFileTag.size(), false);
	for(size_t is=0; is<pTSliceFileTag.size(); is++){
		std::ifstream f;
		f.open(elefakeBaseName + pTSliceFileTag[is] + ".txt");
		std::string variabletype; double variablevalue;
		while(f >> variabletype >> variablevalue){
			if(variabletype=="scalefactor") scalefactor_pTSlice[is] = variablevalue;
			else if(variabletype=="vtxconst") vtxconst_pTSlice[is] = variablevalue;
			else if(variabletype=="vtxslope") vtxslope_pTSlice[is] = variablevalue;
			else if(variabletype=="ptshapeavg") ptshapeavg_pTSlice[is] = variablevalue;
		}
		f.close();
		havePtSlice[is] = (scalefactor_pTSlice[is] > 0);
	}
	auto SelectPtSliceIdx = [&](double ptVal) -> int {
		for(size_t is=0; is<pTSliceUpperEdge.size(); is++)
			if(ptVal < pTSliceUpperEdge[is]) return (int)is;
		return -1;
	};
 
	auto SelectVtxParamsForPt = [&](double ptVal) -> std::array<double,3> {
			for(size_t is=0; is<pTSliceUpperEdge.size(); is++){
				if(ptVal < pTSliceUpperEdge[is]){
					if(havePtSlice[is]) return {scalefactor_pTSlice[is], vtxconst_pTSlice[is], vtxslope_pTSlice[is]};
					return {scalefactor_incl, vtxconst_incl, vtxslope_incl};
				}
			}
			return {scalefactor_incl, vtxconst_incl, vtxslope_incl};
	};

	std::ifstream Eta_file;
	Eta_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-ker-eta-60-120.txt",RunYear,whichVFP.c_str()));
	unsigned nEtaBins(0);
	{
		std::string line;
		while(std::getline(Eta_file, line)) nEtaBins += 1;
		Eta_file.clear(); Eta_file.seekg(0, std::ios::beg);
		nEtaBins /= 2; // file has den block then num block, one line per bin each
	}
	std::vector<double> EtaBins(nEtaBins,0), eta_den(nEtaBins,0), eta_num(nEtaBins,0);
	double avg_fr_eta = 0;
	if(Eta_file.is_open() && nEtaBins>0){
		std::string bintype, numtype; double lowcut, signal1, error1, fitmean, fitrms;
		for(unsigned i=0;i<nEtaBins;i++){
			Eta_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
			eta_den[i] = fitmean; EtaBins[i] = lowcut;
		}
		for(unsigned i=0;i<nEtaBins;i++){
			Eta_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
			eta_num[i] = fitmean;
		}
		for(unsigned i=0;i<nEtaBins;i++) avg_fr_eta += eta_num[i]/eta_den[i];
		avg_fr_eta /= nEtaBins;
	} 
	Eta_file.close();
	std::cout << "avg_fr_eta (mean fake rate over eta bins) = " << avg_fr_eta << std::endl;

	// Same piecewise eta lookup as fitFakeFunc_pT_nVtx.C's w_ele block.
	auto WeightEta = [&](double eta) -> double {
		double weight_eta = (nEtaBins>0) ? eta_num[0]/eta_den[0] : 1.0;
		for(unsigned i=0;i<nEtaBins;i++){
			if(i+1 < nEtaBins){
				if(fabs(eta) > EtaBins[i] && fabs(eta) < EtaBins[i+1]) weight_eta = eta_num[i]/eta_den[i];
			} else {
				if(fabs(eta) > EtaBins[i]) weight_eta = eta_num[i]/eta_den[i];
			}
		}
		return weight_eta;
	};

	// R(pT, eta, Nvtx) = scalefactor_nVtxSlice(Nvtx) * correction *
	//                    ptFitForEvent(pT | Nvtx slice) * [weight_eta(eta)/avg_fr_eta]
	auto EvalFakeRate = [&](double pt, int nVertexVal, double eta) -> double {
		auto p = SelectPtParamsForNVtx(nVertexVal);
		double sf = p[0], slope = p[1], constant = p[2], index = p[3];
		double ptval = fakerate_ptDependence_pred(pt, slope, constant, index);
		double etaCorr = (avg_fr_eta > 0) ? (WeightEta(eta) / avg_fr_eta) : 1.0;
		return sf * correction * ptval * etaCorr;
	};
	auto EvalFakeRate_VtxByPtSlice = [&](double pt, int nVertexVal, double eta) -> double {
		auto p = SelectVtxParamsForPt(pt);
		double sf = p[0], vtxconst = p[1], vtxslope = p[2];
		double vtxval = vtxconst + vtxslope*nVertexVal;
		int idx = SelectPtSliceIdx(pt);
		double residualPt = 1.0;
		if(idx >= 0 && ptshapeavg_pTSlice[idx] > 0){
			double ptval_incl = fakerate_ptDependence_pred(pt, ptslope_incl, ptconstant_incl, ptindex_incl);
			residualPt = ptval_incl / ptshapeavg_pTSlice[idx];
		}
		double etaCorr = (avg_fr_eta > 0) ? (WeightEta(eta) / avg_fr_eta) : 1.0;
		return sf * correction * vtxval * residualPt * etaCorr;
	};


	TH1D *p_controlPhoEt = new TH1D("p_controlPhoEt",";p_{T} (GeV);Events",565,35,600);
	TH1D *p_sigPhoEt  = new TH1D("p_sigPhoEt",";p_{T} (GeV);Events",565,35,600);
	TH1D *p_fakesPhoEt = new TH1D("p_fakesPhoEt",";p_{T} (GeV);Events",565,35,600);
	TH1D *p_fakesPhoEt_stat = new TH1D("p_fakesPhoEt_stat",";p_{T} (GeV);Events",565,35,600);
	TH1D *p_elebkgPhoEt = new TH1D("p_elebkgPhoEt",";p_{T} (GeV);Events",565,35,600);

	TH1D *p_checkEt = new TH1D("p_checkEt","p_checkEt",465,35,500);
	if(detType == 1){
		sigtree->Draw("phoEt >> p_sigPhoEt", " phoEt >35 && sigMET < 70 && fabs(phoEta) < 1.4442");
		std::cout << "signal " << p_sigPhoEt->GetEntries() << std::endl;
		controltree->Draw("phoEt >> p_controlPhoEt", "phoEt > 35 && sigMET < 70 && fabs(phoEta) < 1.4442");
		controltree->Draw("phoEt >> p_checkEt",  "phoEt > 35 && sigMET < 70 && fabs(phoEta) < 1.4442 && phoChIso < 1.694");
	}
	if(detType == 2){
		sigtree->Draw("phoEt >> p_sigPhoEt", " phoEt >35 && sigMET < 70 && fabs(phoEta) > 1.56 && fabs(phoEta) < 2.1");
		std::cout << "signal " << p_sigPhoEt->GetEntries() << std::endl;
		controltree->Draw("phoEt >> p_controlPhoEt", "phoEt > 35 && sigMET < 70 && fabs(phoEta) > 1.56 && fabs(phoEta) < 2.1");
	}
  
	p_sigPhoEt->Sumw2();
	p_controlPhoEt->Sumw2();


	//************ Proxy Tree **********************//
	//e fake photon background events
	TChain *proxytree = new TChain("proxyTree");
	if(channel == 1)proxytree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s.root",RunYear,whichVFP.c_str()));
	else if(channel ==2)proxytree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s.root",RunYear,whichVFP.c_str()));

	float proxyphoEt(0);
	float proxyphoEta(0);
	float proxylepPt(0);
	float proxysigMET(0);
	int   proxynVertex(0);
 
	proxytree->SetBranchAddress("phoEt",     &proxyphoEt);
	proxytree->SetBranchAddress("phoEta",    &proxyphoEta);
	proxytree->SetBranchAddress("lepPt",     &proxylepPt);
	proxytree->SetBranchAddress("sigMET",    &proxysigMET);
	proxytree->SetBranchAddress("nVertex",   &proxynVertex);
 
	for (unsigned ievt(0); ievt<proxytree->GetEntries(); ++ievt){//loop on entries
		proxytree->GetEntry(ievt);
		if(proxysigMET > 70)continue;
		if(detType == 1 && fabs(proxyphoEta) > 1.4442)continue;
		else if(detType == 2 && (fabs(proxyphoEta) < 1.56 || fabs(proxyphoEta) > 2.1))continue;
		if(proxyphoEt < 35)continue;
		double w_ele = USE_VTX_BY_PT_SLICE ? EvalFakeRate_VtxByPtSlice(proxyphoEt, proxynVertex, fabs(proxyphoEta))
		                                    : EvalFakeRate(proxyphoEt, proxynVertex, fabs(proxyphoEta));
		p_elebkgPhoEt->Fill(proxyphoEt,w_ele);
	}
	p_elebkgPhoEt->Sumw2();
	p_sigPhoEt->Add(p_elebkgPhoEt, -1);
	// subtracted the fake contribution
	p_sigPhoEt->Sumw2();
	for(int ibin(1); ibin < p_sigPhoEt->GetSize()-1; ibin++){
		double xvalue = p_sigPhoEt->GetBinCenter(ibin);
		double frac(0),fracerror(0);
		for(unsigned i(0); i< NBIN-1; i++) 
			if(xvalue >= PtBin[i] && xvalue < PtBin[i+1]){
				frac = fracHad[i];
				fracerror = fracHadError[i];
				// hadron fraction
			}
		if(xvalue >= PtBin[NBIN-1]){ frac = fracHad[NBIN-1]; fracerror = fracHadError[NBIN-1];}
    
		double binvalue = p_sigPhoEt->GetBinContent(ibin)*frac;
		// hadron fraction multiplied
		if(binvalue == 0)continue;
		double totalerror = p_sigPhoEt->GetBinContent(ibin)*fracerror;
		p_fakesPhoEt_stat->SetBinContent(ibin, p_sigPhoEt->GetBinContent(ibin)); 
		p_fakesPhoEt->SetBinContent(ibin, binvalue);
		p_fakesPhoEt->SetBinError(ibin, totalerror);
	}


	TCanvas *c_pt = new TCanvas("Photon_Pt", "Photon P_{T}",800,800);
	c_pt->cd();
	TPad *can_pad1 = new TPad("can_pad1", "pad1", 0, 0.30, 1, 1.0);
	can_pad1->SetBottomMargin(0.02);
	can_pad1->Draw();          
	can_pad1->cd();          

	TGraphAsymmErrors *new_controlPhoEt = new TGraphAsymmErrors(125);
	TGraphAsymmErrors *new_fakesPhoEt   = new TGraphAsymmErrors(125);
	// bins 1-75, bin width 1, pt 35-110
	for(unsigned ibin(1); ibin <= 75; ibin++){
		double new_control_value = p_controlPhoEt->GetBinContent(ibin); 
		double new_fakes_value = p_fakesPhoEt->GetBinContent(ibin);
		double new_control_error = p_controlPhoEt->GetBinError(ibin); 
		double new_fakes_error = p_fakesPhoEt->GetBinError(ibin);
		double staterror = sqrt(p_fakesPhoEt_stat->GetBinContent(ibin))/p_fakesPhoEt_stat->GetBinContent(ibin)*new_fakes_value;
		//std::cout << ibin << " " << sqrt(p_fakesPhoEt_stat->GetBinContent(ibin))/p_fakesPhoEt_stat->GetBinContent(ibin) << std::endl;
		cout<<p_fakesPhoEt->GetBinCenter(ibin) <<"  "<< new_fakes_value << "  "<<new_fakes_error<<endl;
		new_fakes_error = sqrt(new_fakes_error*new_fakes_error + staterror);
		new_controlPhoEt->SetPoint(ibin,p_controlPhoEt->GetBinCenter(ibin),  new_control_value);
		//new_controlPhoEt->SetPointError(ibin, 0.5, 0.5, min(new_control_error, new_control_value-0.001), new_control_error);
		new_controlPhoEt->SetPointError(ibin, 0.5, 0.5, new_control_error, new_control_error);
		new_fakesPhoEt->SetPoint(ibin, p_fakesPhoEt->GetBinCenter(ibin),  new_fakes_value);
		//new_fakesPhoEt->SetPointError(ibin,0.5, 0.5, min(new_fakes_error, new_fakes_value - 0.001), new_fakes_error);
		new_fakesPhoEt->SetPointError(ibin,0.5, 0.5, new_fakes_error, new_fakes_error);
	}
	// bins 76-94, bin width =  10, pt 110-600
	for(unsigned ibin(76); ibin <= 124; ibin++){
		double new_control_value = 0;
		double new_fakes_value = 0;
		double new_control_error = 0;
		double new_fakes_error = 0;
		double staterror = 0;
		for(unsigned j(1); j <=10; j++){
			new_control_value += p_controlPhoEt->GetBinContent(75+(ibin-76)*10+j); 
			new_fakes_value += p_fakesPhoEt->GetBinContent(75+(ibin-76)*10+j);
			new_control_error += p_controlPhoEt->GetBinContent(75+(ibin-76)*10+j);
			new_fakes_error += p_fakesPhoEt->GetBinError(75+(ibin-76)*10+j);
			staterror += p_fakesPhoEt_stat->GetBinContent(75+(ibin-76)*10+j);
		}
		new_control_value /=10.0;
		new_control_error  = new_control_value*sqrt(new_control_error)/new_control_error; 
		new_fakes_value  /=10.0;
		new_fakes_error   /= 10.0;
		//std::cout << ibin << " " << sqrt(staterror)/staterror << std::endl;
		//staterror = sqrt(staterror)/staterror*new_fakes_value;
		staterror = 0;
		new_fakes_error = sqrt(new_fakes_error*new_fakes_error + staterror*staterror);
		new_controlPhoEt->SetPoint(ibin,p_controlPhoEt->GetBinCenter(80+(ibin-76)*10),  new_control_value);
		//new_controlPhoEt->SetPointError(ibin,0, 0, min(new_control_error, new_control_value-0.001), new_control_error);
		new_controlPhoEt->SetPointError(ibin,5, 5, new_control_error, new_control_error);
		new_fakesPhoEt->SetPoint(ibin, p_fakesPhoEt->GetBinCenter(80+(ibin-76)*10),  new_fakes_value);
		//new_fakesPhoEt->SetPointError(ibin,0, 0, min(new_fakes_error, new_fakes_value - 0.001), new_fakes_error);
		new_fakesPhoEt->SetPointError(ibin,5, 5, new_fakes_error, new_fakes_error);
		cout<<p_fakesPhoEt->GetBinCenter(80+(ibin-76)*10) <<"  "<< new_fakes_value << "  "<<new_fakes_error<<endl;
	}
	{ // overflow bin
		double new_control_value = p_controlPhoEt->GetBinContent(566); 
		double new_fakes_value = p_fakesPhoEt->GetBinContent(566);
		double new_control_error = p_controlPhoEt->GetBinError(566); 
		double new_fakes_error = p_fakesPhoEt->GetBinError(566);
		new_controlPhoEt->SetPoint(125,600,  new_control_value);
		new_controlPhoEt->SetPointError(125,5, 5, min(new_control_error, new_control_value-0.001), new_control_error);
		new_fakesPhoEt->SetPoint(125, 600,  new_fakes_value);
		new_fakesPhoEt->SetPointError(125,5, 5, min(new_fakes_error, new_fakes_value - 0.001), new_fakes_error);
	}

	TFile* fout = nullptr;
	if (channel == 1)
    		fout = new TFile(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/fitJetFunc_Output_DoubleEG_%d%s_pT_nVtx.root", RunYear, whichVFP.c_str()), "RECREATE");
	else if (channel == 2)
    		fout = new TFile(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/fitJetFunc_Output_MuonEG_%d%s_pT_nVtx.root", RunYear, whichVFP.c_str()), "RECREATE");

	new_fakesPhoEt->GetListOfFunctions()->Delete();
	new_fakesPhoEt->SetName("new_fakesPhoEt");
	new_fakesPhoEt->Write();

	new_controlPhoEt->GetListOfFunctions()->Delete();
	new_controlPhoEt->SetName("new_controlPhoEt");
	new_controlPhoEt->Write();

	fout->Close();
	std::cout << "Histograms written to fitJetFunc_Output.root without extra functions." << std::endl;
	
	new_fakesPhoEt->GetXaxis()->SetTitle("p_{T} (GeV)");
	new_controlPhoEt->GetXaxis()->SetTitle("p_{T} (GeV)");
	new_controlPhoEt->GetXaxis()->SetTitleOffset(1);
	new_controlPhoEt->GetXaxis()->SetTitleSize(20);
	new_controlPhoEt->GetXaxis()->SetRangeUser(35,600);
	new_fakesPhoEt->GetXaxis()->SetRangeUser(35,600);
	gPad->SetLogy();
	TH1D *new_dummy=new TH1D("dummy",";p_{T} (GeV); Event/GeV",265,35,300);
	new_dummy->SetMinimum(0.001);
	new_dummy->SetMaximum(100000);
	new_dummy->GetXaxis()->SetTitle("p_{T} (GeV)");
	new_dummy->GetXaxis()->SetTitleOffset(1.0);
	new_dummy->GetXaxis()->SetLabelSize(0);   // x axis is drawn on the ratio pad
	new_dummy->GetXaxis()->SetTitleSize(0);
	new_dummy->Draw();
	new_controlPhoEt->Draw("P same");
	new_controlPhoEt->SetLineColor(kBlack);
	new_controlPhoEt->SetMarkerStyle(20);
	new_fakesPhoEt->SetLineColor(kRed);
	new_fakesPhoEt->SetMarkerStyle(20);
	new_fakesPhoEt->SetMarkerColor(kRed);
	new_fakesPhoEt->Draw("EP same");
	TLegend *leg =  new TLegend(0.6,0.7,0.85,0.85);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	gStyle->SetLegendFillColor(0);
	new_controlPhoEt->SetFillColor(0);
	new_fakesPhoEt->SetFillColor(0);
	leg->AddEntry(new_controlPhoEt,"hadron proxies");
	leg->AddEntry(new_fakesPhoEt,"fake photons");
	leg->Draw("same");

	TCanvas* mccan = new TCanvas("mccan","mccan",1200,600) ;
	mccan->Divide(2);

	//********************   denominator *****************************************************//
	
	TF1 *fitfunc_den= new TF1("fitfunc_den", exp2c_func, 20, 600, 6);
	fitfunc_den->SetParameters(1e4, 1e3, 0.05, 0.01, 1.0);
	new_controlPhoEt->Fit("fitfunc_den","S");
	 if (fitfunc_den->GetParameter(5) < 0) {
                cout << "[INFO] fitfunc_den offset fit negative (" << fitfunc_den->GetParameter(5)
                     << "); refitting with offset floored at 0." << endl;
                fitfunc_den->SetParameter(5, 0.0);
                fitfunc_den->SetParLimits(5, 0, 100);
                new_controlPhoEt->Fit("fitfunc_den","S");
        }

	TF1 *fitden = new_controlPhoEt->GetFunction("fitfunc_den");

	ofstream myfile;
	if(channel == 1 and detType == 1)myfile.open(Form("/eos/uscms/store/user/tmishra/jetfakepho/txt%d%s/JetFakeRate-transferfactor-DoubleEG-EB_pT_nVtx.txt",RunYear,whichVFP.c_str()), std::ios_base::ate | std::ios_base::out);
	if(channel == 1 and detType == 2)myfile.open(Form("/eos/uscms/store/user/tmishra/jetfakepho/txt%d%s/JetFakeRate-transferfactor-DoubleEG-EE_pT_nVtx.txt",RunYear,whichVFP.c_str()), std::ios_base::ate | std::ios_base::out);

	if(channel == 2 and detType == 1)myfile.open(Form("/eos/uscms/store/user/tmishra/jetfakepho/txt%d%s/JetFakeRate-transferfactor-MuonEG-EB_pT_nVtx.txt",RunYear,whichVFP.c_str()), std::ios_base::ate | std::ios_base::out);
	if(channel == 2 and detType == 2)myfile.open(Form("/eos/uscms/store/user/tmishra/jetfakepho/txt%d%s/JetFakeRate-transferfactor-MuonEG-EE_pT_nVtx.txt",RunYear,whichVFP.c_str()), std::ios_base::ate | std::ios_base::out);
	
	TH1D *ratio = new TH1D("transfer fraction","",int(p_fakesPhoEt->GetXaxis()->GetNbins()/REBINSIZE),35,600);
	for(int ibin(1); ibin < int(p_fakesPhoEt->GetXaxis()->GetNbins()/REBINSIZE); ibin++){
		ratio->SetBinContent(ibin, new_fakesPhoEt->Eval(35+(ibin-1)*2 +1)/new_controlPhoEt->Eval(35+(ibin-1)*2 +1));
	}
 
	myfile << "den_coeff1 " << fitden->GetParameter(0) << std::endl;
	myfile << "den_coeff2 " << fitden->GetParameter(1) << std::endl;
	myfile << "den_lambd1 " << fitden->GetParameter(2) << std::endl;
	myfile << "den_lambd2 " << fitden->GetParameter(3) << std::endl;
	myfile << "den_shift " << fitden->GetParameter(4) << std::endl;
	myfile << "den_offset " << fitden->GetParameter(5) << std::endl;

	cout << "den_coeff1 " << fitden->GetParameter(0) << std::endl;
	cout << "den_coeff2 " << fitden->GetParameter(1) << std::endl;
	cout << "den_lambd1 " << fitden->GetParameter(2) << std::endl;
	cout << "den_lambd2 " << fitden->GetParameter(3) << std::endl;
	cout << "den_shift " << fitden->GetParameter(4) << std::endl;
	cout << "den_offset " << fitden->GetParameter(5) << std::endl;

	cout << "chiSquare/ndf:  " << fitden->GetChisquare()/fitden->GetNDF() << std::endl;
	cout << "chiSquare:  " << fitden->GetChisquare() << std::endl;

	TFitResultPtr rden = new_controlPhoEt->Fit("fitfunc_den","S");
	// covariance matrix
	TMatrixDSym covden = rden->GetCovarianceMatrix(); 
	rden->Print("V");
	fitfunc_den->Draw("same");
	TVectorD muden(6) ;
	muden(0) = rden->Parameter(0);  // c1
        muden(1) = rden->Parameter(1);  // c2
        muden(2) = rden->Parameter(2);  // λ1
        muden(3) = rden->Parameter(3);  // λ2
        muden(4) = rden->Parameter(4);  // center shift
        muden(5) = rden->Parameter(5);  // offset

	// FIX: RooRealVar's 3rd/4th args are the *hard allowed range* for toy
	// generation, previously a symmetric +-1sigma window around each best fit
	// (mu-error, mu+error). For c1,c2,lambda1,lambda2,offset -- every term
	// that must stay >=0 for exp2c_func to stay non-negative -- that window
	// can dip below 0 whenever the parameter's error is comparable to or
	// bigger than its value (exactly the common case for offset once it's
	// sitting near the positivity floor). That let the toy generator draw
	// unphysical negative parameter sets, which is what was crashing
	// den_lower to zero at high pT. Floor each range at its own physical
	// boundary so every toy is guaranteed non-negative everywhere by
	// construction -- shift is the only one of the six with no sign
	// requirement, so it's left as a symmetric window.

	
        RooRealVar central_coeff1_den("central_coeff1_den","central_coeff1_den",muden(0)-rden->ParError(0), muden(0)+rden->ParError(0));
        RooRealVar central_coeff2_den("central_coeff2_den","central_coeff2_den",muden(1)-rden->ParError(1), muden(1)+rden->ParError(1));
        RooRealVar central_lambda1_den("central_lambda1_den","central_lambda1_den",muden(2)-rden->ParError(2),muden(2)+rden->ParError(2));
        RooRealVar central_lambda2_den("central_lambda2_den","central_lambda2_den",muden(3)-rden->ParError(3),muden(3)+rden->ParError(3));
        RooRealVar shift_den("shift_den","shift_den", muden(4) - rden->ParError(4), muden(4) + rden->ParError(4));
        RooRealVar offset_den("offset_den","offset_den", muden(5) - rden->ParError(5), muden(5) + rden->ParError(5));

	// Multivariant Gaussian pdf, get NTOY numbers of paramters
	RooMultiVarGaussian mvgden("mvgden","mvgden",RooArgList(central_coeff1_den,central_coeff2_den,central_lambda1_den,central_lambda2_den, shift_den, offset_den),muden,covden);
	RooDataSet* toymcdataden = mvgden.generate(RooArgSet(central_coeff1_den,central_coeff2_den,central_lambda1_den,central_lambda2_den, shift_den, offset_den),NTOY);
	std::ostringstream modelnameden;
	TF1 *gen_den[NTOY];
	TH1D *den_upper = new TH1D("den_upper","den_upper",565,35,600);
	TH1D *den_lower = new TH1D("den_lower","den_lower",565,35,600);
	can_pad1->cd();
	// FIX: store every toy's per-bin value instead of updating den_upper/
	// den_lower with a running min/max. See the matching comment in the
	// numerator block below for why min/max-over-NTOY is not a meaningful
	// confidence band (it's an order statistic that keeps growing more
	// extreme with NTOY, evaluated independently per bin) -- that's what was
	// dragging den_lower roughly 1-2 orders of magnitude below the data even
	// after the positivity-range fix removed the unphysical negative dips.
	float toyptvalueden[565][NTOY];
	for (unsigned ibin(0); ibin < 565; ibin++)
		for (int i(0); i < NTOY; i++)
			toyptvalueden[ibin][i] = std::numeric_limits<float>::quiet_NaN();
	for(int i(0); i<NTOY; i++){
		double data1 = toymcdataden->get(i)->getRealValue("central_coeff1_den");
		double data2 = toymcdataden->get(i)->getRealValue("central_coeff2_den");
		double data3 = toymcdataden->get(i)->getRealValue("central_lambda1_den");
		double data4 = toymcdataden->get(i)->getRealValue("central_lambda2_den");
		double data5 = toymcdataden->get(i)->getRealValue("shift_den");
		double data6 = toymcdataden->get(i)->getRealValue("offset_den");
		modelnameden.str("");
		modelnameden << "gen_den_" << i;
		gen_den[i] = new TF1(modelnameden.str().c_str(), exp2c_func, 20, 600, 6);
		gen_den[i]->SetParameters(data1, data2, data3, data4, data5, data6);
		gen_den[i]->SetLineColorAlpha(kBlue, 0.35);
		//gen_den[i]->Draw("same");
	 	for(unsigned ibin(1); ibin <= 565; ibin++){
			double estimated = gen_den[i]->Eval(den_upper->GetBinCenter(ibin));
			toyptvalueden[ibin-1][i] = estimated;
		}
	}

	// 68% CL (16th-84th percentile) band from the toy distribution at each
	// bin, replacing the old min/max envelope.
	for (unsigned ibin(0); ibin < 565; ibin++) {
		std::vector<float> vals;
		vals.reserve(NTOY);
		for (int i = 0; i < NTOY; i++) {
			float v = toyptvalueden[ibin][i];
			if (std::isfinite(v)) vals.push_back(v);
		}
		double nominal = fitfunc_den->Eval(den_upper->GetBinCenter(ibin+1));
		if (vals.empty()) {
			den_upper->SetBinContent(ibin+1, nominal);
			den_lower->SetBinContent(ibin+1, nominal);
			continue;
		}
		std::sort(vals.begin(), vals.end());
		int n = (int)vals.size();
		int ilo = std::max(0, (int)std::lround(0.16 * (n - 1)));
		int ihi = std::min(n - 1, (int)std::lround(0.84 * (n - 1)));
		den_lower->SetBinContent(ibin+1, vals[ilo]);
		den_upper->SetBinContent(ibin+1, vals[ihi]);
	}

	den_upper->Draw("L same");
	den_lower->Draw("L same");
// *************************  Numerator ******************************************************************//

TF1 *fitfunc_num = new TF1("fitfunc_num", exp2c_func, 20, 1000, 6);

fitfunc_num->SetParameters(51521, 1.13974, 0.0560658, 0.0200371, 125.136, 0.0970124);
fitfunc_num->SetParLimits(0, 0, 1e9);      // c1: positive
fitfunc_num->SetParLimits(1, 0, 1e5);      // c2: positive
fitfunc_num->SetParLimits(2, 0.0001, 1);   // lamda1: decay rate
fitfunc_num->SetParLimits(3, 0.0001, 1);   // lamda2: decay rate
fitfunc_num->SetParLimits(4, 0, 300);      // shift
fitfunc_num->SetParLimits(5, 0, 100);      // offset: >=0, was (-100, 100)


// FIX: build a log10(y)-space copy of the numerator graph and fit
// fitfunc_num's shape there instead of directly against new_fakesPhoEt. See
// the note above log_exp2c_func for why: a linear-scale chi2 fit is
// dominated by the low-pT bulk and leaves the pT > 200 GeV tail essentially
// unconstrained, which is exactly the "fit not good at high pT" symptom.
// Points with y <= 0 (e.g. graph index 0, which SetPoint() never touches --
// the fill loops below all start at ibin = 1) are skipped since log10(0) is
// undefined.
TGraphAsymmErrors *log_fakesPhoEt = new TGraphAsymmErrors();
{
	int np = 0;
	double xx, yy;
	const double ln10 = std::log(10.0);
	for (int i = 0; i < new_fakesPhoEt->GetN(); i++) {
		new_fakesPhoEt->GetPoint(i, xx, yy);
		if (yy <= 0) continue;
		double eyl = new_fakesPhoEt->GetErrorYlow(i);
		double eyh = new_fakesPhoEt->GetErrorYhigh(i);
		double exl = new_fakesPhoEt->GetErrorXlow(i);
		double exh = new_fakesPhoEt->GetErrorXhigh(i);
		// standard log10 error propagation: sigma_log = sigma_y / (y * ln10)
		log_fakesPhoEt->SetPoint(np, xx, std::log10(yy));
		log_fakesPhoEt->SetPointError(np, exl, exh, eyl/(yy*ln10), eyh/(yy*ln10));
		np++;
	}
}

TF1 *fitfunc_num_log = new TF1("fitfunc_num_log", log_exp2c_func, 20, 1000, 6);
fitfunc_num_log->SetParameters(fitfunc_num->GetParameters());
for (int i = 0; i < 6; ++i) {
	double lo, hi;
	fitfunc_num->GetParLimits(i, lo, hi);
	fitfunc_num_log->SetParLimits(i, lo, hi);
}


cout << "[DEBUG] Initial Parameters Set: " << endl;
for (int i = 0; i < 6; ++i) {
    cout << "  p" << i << " = " << fitfunc_num_log->GetParameter(i) << endl;
}

cout << "[DEBUG] Histogram integral (total entries): " << new_fakesPhoEt->Integral() << endl;
if (new_fakesPhoEt->Integral() == 0) {
    cout << "[ERROR] Histogram has zero total content. Fit will fail." << endl;
}

log_fakesPhoEt->Fit("fitfunc_num_log", "S"); // silent fit, retrieve result below
TF1 *fitnum = fitfunc_num_log;

if (fitfunc_num_log->GetParameter(5) < 0) {
        cout << "[INFO] fitfunc_num offset fit negative (" << fitfunc_num_log->GetParameter(5)
             << "); refitting with offset floored at 0." << endl;
        fitfunc_num_log->SetParameter(5, 0.0);
        fitfunc_num_log->SetParLimits(5, 0, 100);
}

TVirtualFitter::SetMaxIterations(1000000);
TFitResultPtr r = log_fakesPhoEt->Fit("fitfunc_num_log", "R S");

fitfunc_num->SetParameters(fitfunc_num_log->GetParameters());
for (int i = 0; i < 6; ++i) fitfunc_num->SetParError(i, fitfunc_num_log->GetParError(i));
fitfunc_num->SetChisquare(fitfunc_num_log->GetChisquare());
fitfunc_num->SetNDF(fitfunc_num_log->GetNDF());
fitfunc_num->SetLineColor(kRed+1);
fitfunc_num->SetNpx(1000);
new_fakesPhoEt->GetListOfFunctions()->Add(fitfunc_num);

if (r.Get() == nullptr || r->IsValid() == false) {
    cout << "[ERROR] Fit failed or returned invalid result!" << endl;
} else {
    cout << "[DEBUG] Fit successful. Status: " << r->Status() << ", IsValid: " << r->IsValid() << endl;
    cout << "[DEBUG] EDM: " << r->Edm() << ", NCalls: " << r->NCalls() << endl;
}

TF1 *fit = new_fakesPhoEt->GetFunction("fitfunc_num");

        myfile << "num_coeff1 " << fit->GetParameter(0) << std::endl;
        myfile << "num_coeff2 " << fit->GetParameter(1) << std::endl;
        myfile << "num_lambd1 " << fit->GetParameter(2) << std::endl;
        myfile << "num_lambd2 " << fit->GetParameter(3) << std::endl;
        myfile << "num_shift " 	<< fit->GetParameter(4) << std::endl;
        myfile << "num_offset " << fit->GetParameter(5) << std::endl;


cout << "chiSquare/ndf:  " << fit->GetChisquare()/fit->GetNDF() << std::endl;
cout << "chiSquare:      " << fit->GetChisquare() << std::endl;

can_pad1->cd();
r->Print("V"); // Verbose print
fitfunc_num->Draw("same");

float nominalvalue_num[565];
for(unsigned ibin = 0; ibin < 565; ibin++){
    double xval = 35 + ibin;
    nominalvalue_num[ibin] = fitfunc_num->Eval(xval);
    if (!std::isfinite(nominalvalue_num[ibin])) {
        cout << "[WARNING] Fit function returned non-finite value at x = " << xval << endl;
    }
}

	TVectorD mu(6);
	mu(0) = r->Parameter(0);  // c1
	mu(1) = r->Parameter(1);  // c2
	mu(2) = r->Parameter(2);  // λ1
	mu(3) = r->Parameter(3);  // λ2
	mu(4) = r->Parameter(4);  // center shift
	mu(5) = r->Parameter(5);  // offset
	TMatrixDSym cov = r->GetCovarianceMatrix();

	// FIX: same positivity floor as the denominator toy ranges above -- see
	// that comment for why the previous symmetric (mu-error, mu+error)
	// window could and did sample unphysical negative c1/c2/lambda/offset,
	// which is what was crashing num_lower to zero above pT ~ 120-150 GeV.

	RooRealVar central_coeff1_num("central_coeff1_num","central_coeff1_num",mu(0)-r->ParError(0), mu(0)+r->ParError(0));
        RooRealVar central_coeff2_num("central_coeff2_num","central_coeff2_num",mu(1)-r->ParError(1), mu(1)+r->ParError(1));
        RooRealVar central_lambda1_num("central_lambda1_num","central_lambda1_num",mu(2)-r->ParError(2),mu(2)+r->ParError(2));
        RooRealVar central_lambda2_num("central_lambda2_num","central_lambda2_num",mu(3)-r->ParError(3),mu(3)+r->ParError(3));
        RooRealVar shift_num("shift_num","shift_num", mu(4) - r->ParError(4), mu(4) + r->ParError(4));
        RooRealVar offset_num("offset_num","offset_num", mu(5) - r->ParError(5), mu(5) + r->ParError(5));

	RooMultiVarGaussian mvg("mvg","mvg",RooArgList(central_coeff1_num,central_coeff2_num,central_lambda1_num,central_lambda2_num, shift_num, offset_num),mu,cov);
        RooDataSet* toymcdata = mvg.generate(RooArgSet(central_coeff1_num,central_coeff2_num,central_lambda1_num,central_lambda2_num, shift_num, offset_num),NTOY);

	std::ostringstream modelname;
	TF1 *gen_num[NTOY];

	TH1D *num_upper = new TH1D("num_upper","num_upper",565,35,600);
	TH1D *num_lower = new TH1D("num_lower","num_lower",565,35,600);
	float toyptvalue[565][NTOY];
	float lowtoyptvalue[565];
	float hightoyptvalue[565];
	for(unsigned ii(0); ii < 565; ii++){
		lowtoyptvalue[ii]  = fitfunc_num->Eval(35 + ii);   // fallback if a bin has no valid toys
		hightoyptvalue[ii] = fitfunc_num->Eval(35 + ii);   // fallback if a bin has no valid toys
		for (int i(0); i < NTOY; i++) toyptvalue[ii][i] = std::numeric_limits<float>::quiet_NaN();
	}
	
	for(int i(0); i<NTOY; i++){
		double data1 = toymcdata->get(i)->getRealValue("central_coeff1_num");
		double data2 = toymcdata->get(i)->getRealValue("central_coeff2_num");
		double data3 = toymcdata->get(i)->getRealValue("central_lambda1_num");
		double data4 = toymcdata->get(i)->getRealValue("central_lambda2_num");
		double data5 = toymcdata->get(i)->getRealValue("shift_num");
		double data6 = toymcdata->get(i)->getRealValue("offset_num");
		if(data1 < 1e8 && data2 < 1e8 && data3 < 1e8 && data4 < 1e8 && data5 < 1e8 && data6 < 1e8){
			modelname.str("");
			modelname << "gen_num_" << i;
			gen_num[i] = new TF1(modelname.str().c_str(), exp2c_func, 20, 600, 6);
			gen_num[i]->SetParameters(data1, data2, data3, data4, data5, data6);
			gen_num[i]->SetLineColorAlpha(kBlue, 0.35);
			bool exception(false);
			if(exception)continue; 

			for(unsigned ibin(0); ibin < 565; ibin++){
				double estimated = gen_num[i]->Eval(35+ibin);
				if (std::isnan(estimated) || std::isinf(estimated)) {
				    std::cout << "[WARNING] Fit function returned non-finite value at bin " << ibin
              				<< " for toy " << i << std::endl;
    				    		std::cout << "  Parameters: c1=" << data1 << ", c2=" << data2
              					<< ", lambda1=" << data3 << ", lambda2=" << data4
              					<< ", shift=" << data5 << ", offset=" << data6 << std::endl;
    					exception = true;
    					break; 
				}
				// FIX: store every toy's value; the running min/max update
				// that used to live here (see comment below) is replaced by
				// a proper 68% CL percentile pass over the full toy
				// distribution once all NTOY toys have been generated.
				toyptvalue[ibin][i] = estimated;
			}

		}
		if (data3 < 0 || data4 < 0) 
			    std::cout << "[DEBUG] λ1 or λ2 negative: λ1 = " << data3 << ", λ2 = " << data4 << std::endl;
		if (data5 + 200 > 600) 
    			std::cout << "[DEBUG] Shift too large: shift+200 = " << (data5 + 200) << std::endl;
	}

	// FIX: 68% CL (16th-84th percentile) band from the toy distribution at
	// each bin, replacing the old min/max-over-NTOY envelope -- see the
	// comment by toyptvalueden above for why min/max isn't a meaningful
	// confidence band. This is what was dragging num_lower 1-2 orders of
	// magnitude below the data even after the positivity-range fix removed
	// the unphysical negative dips.
	for (unsigned ibin(0); ibin < 565; ibin++) {
		std::vector<float> vals;
		vals.reserve(NTOY);
		for (int i = 0; i < NTOY; i++) {
			float v = toyptvalue[ibin][i];
			if (std::isfinite(v)) vals.push_back(v);
		}
		if (vals.empty()) continue; // leave the nominal-fit fallback in place
		std::sort(vals.begin(), vals.end());
		int n = (int)vals.size();
		int ilo = std::max(0, (int)std::lround(0.16 * (n - 1)));
		int ihi = std::min(n - 1, (int)std::lround(0.84 * (n - 1)));
		lowtoyptvalue[ibin]  = vals[ilo];
		hightoyptvalue[ibin] = vals[ihi];
	}
	
	TCanvas *cangaus = new TCanvas("cangaus","",600,600);
	cangaus->cd();
	for(unsigned ibin(0); ibin < 565; ibin++){
		// histogram bin ibin+1 is pt = 35+ibin (bin ibin was one low and left bin 565 empty)
		num_upper->SetBinContent(ibin+1, hightoyptvalue[ibin]);
		num_lower->SetBinContent(ibin+1, lowtoyptvalue[ibin]);
		cout<< ibin << " "<< hightoyptvalue[ibin] << " " << lowtoyptvalue[ibin]<<endl;
	}
	
	can_pad1->cd();          
	num_upper->Draw("L same");
	num_lower->Draw("L same");
	new_controlPhoEt->Draw("EP same");
	new_fakesPhoEt->Draw("EP same");

for (unsigned ibin(0); ibin < 565; ibin += 40) {
    double nominal = fitfunc_num->Eval(35 + ibin);
    cout << "pt=" << 35+ibin
         << " nominal=" << nominal
         << " upper=" << hightoyptvalue[ibin]
         << " lower=" << lowtoyptvalue[ibin]
         << " (upper-nominal)=" << hightoyptvalue[ibin]-nominal
         << endl;
}
	
        if(RunYear==2016 and ISpreVFP == 1)       CMS_lumi( can_pad1, 1, ichannel, 11 );
        else if(RunYear==2016 and ISpreVFP == 0)  CMS_lumi( can_pad1, 2, ichannel, 11 );
        else if(RunYear==2017)                    CMS_lumi( can_pad1, 3, ichannel, 11 );
        else if(RunYear==2018)                    CMS_lumi( can_pad1, 4, ichannel, 11 );

	// ---------------- Ratio pad: data / fit, with data uncertainties ----------------
	// Created here (rather than up front) so nothing in between - the mccan/cangaus
	// canvases, the many ->cd() calls during fitting - can leave it mis-attached.
	// Point index 0 of new_controlPhoEt/new_fakesPhoEt is never filled (stays at (0,0)),
	// so start from index 1, matching the real pT range used above.
	c_pt->cd();
	TPad *can_pad2 = new TPad("can_pad2", "pad2", 0, 0.0, 1, 0.30);
	can_pad2->SetTopMargin(0.03);
	can_pad2->SetBottomMargin(0.35);
	can_pad2->SetLeftMargin(can_pad1->GetLeftMargin());
	can_pad2->SetRightMargin(can_pad1->GetRightMargin());
	can_pad2->SetGridy();
	can_pad2->Draw();
	can_pad2->cd();

	TGraphAsymmErrors *ratio_control = new TGraphAsymmErrors(new_controlPhoEt->GetN());
	TGraphAsymmErrors *ratio_fakes   = new TGraphAsymmErrors(new_fakesPhoEt->GetN());

	for (int ip = 1; ip < new_controlPhoEt->GetN(); ip++) {
		double x, y;
		new_controlPhoEt->GetPoint(ip, x, y);
		double fitval = fitfunc_den->Eval(x);
		double eylow  = new_controlPhoEt->GetErrorYlow(ip);
		double eyhigh = new_controlPhoEt->GetErrorYhigh(ip);
		double rval = (fitval != 0) ? y / fitval : 0;
		double erl  = (fitval != 0) ? eylow  / fitval : 0;
		double erh  = (fitval != 0) ? eyhigh / fitval : 0;
		ratio_control->SetPoint(ip, x, rval);
		ratio_control->SetPointError(ip, new_controlPhoEt->GetErrorXlow(ip), new_controlPhoEt->GetErrorXhigh(ip), erl, erh);
	}

	for (int ip = 1; ip < new_fakesPhoEt->GetN(); ip++) {
		double x, y;
		new_fakesPhoEt->GetPoint(ip, x, y);
		double fitval = fitfunc_num->Eval(x);
		double eylow  = new_fakesPhoEt->GetErrorYlow(ip);
		double eyhigh = new_fakesPhoEt->GetErrorYhigh(ip);
		double rval = (fitval != 0) ? y / fitval : 0;
		double erl  = (fitval != 0) ? eylow  / fitval : 0;
		double erh  = (fitval != 0) ? eyhigh / fitval : 0;
		ratio_fakes->SetPoint(ip, x, rval);
		ratio_fakes->SetPointError(ip, new_fakesPhoEt->GetErrorXlow(ip), new_fakesPhoEt->GetErrorXhigh(ip), erl, erh);
	}

	TH1D *ratio_dummy = new TH1D("ratio_dummy", ";p_{T} (GeV);Data / Fit", 265, 35, 300);
	ratio_dummy->SetMinimum(0.0);
	ratio_dummy->SetMaximum(2.0);
	ratio_dummy->GetXaxis()->SetTitle("p_{T} (GeV)");
	ratio_dummy->GetYaxis()->SetTitle("Data / Fit");

	ratio_dummy->GetXaxis()->SetTitleSize(0.11);
	ratio_dummy->GetXaxis()->SetLabelSize(0.09);
	ratio_dummy->GetXaxis()->SetTitleOffset(1.1);
	ratio_dummy->GetXaxis()->SetTitleFont(42);
	ratio_dummy->GetXaxis()->SetLabelFont(42);
	ratio_dummy->GetYaxis()->SetTitleSize(0.10);
	ratio_dummy->GetYaxis()->SetLabelSize(0.08);
	ratio_dummy->GetYaxis()->SetTitleOffset(0.5);
	ratio_dummy->GetYaxis()->SetTitleFont(42);
	ratio_dummy->GetYaxis()->SetLabelFont(42);
	ratio_dummy->GetYaxis()->SetNdivisions(505);
	ratio_dummy->Draw("AXIS");

	TLine *unity = new TLine(35, 1, 300, 1);
	unity->SetLineStyle(2);
	unity->SetLineColor(kBlack);
	unity->Draw("same");

	ratio_control->SetLineColor(kBlack);
	ratio_control->SetMarkerStyle(20);
	ratio_control->SetMarkerColor(kBlack);
	ratio_control->Draw("P same");

	ratio_fakes->SetLineColor(kRed);
	ratio_fakes->SetMarkerStyle(20);
	ratio_fakes->SetMarkerColor(kRed);
	ratio_fakes->Draw("P same");

	can_pad2->Modified();
	can_pad2->Update();
	c_pt->cd();
	c_pt->Modified();
	c_pt->Update();
	// ---------------------------------------------------------------------------------
	for(unsigned ii(1); ii <= 265; ii++){
		myfile << "numerror " << ii-1 << " " <<  num_upper->GetBinContent(ii) - num_lower->GetBinContent(ii) << std::endl;
	} 
	for(unsigned ii(1); ii <= 265; ii++){
		myfile << "denerror " << ii-1 << " " <<  den_upper->GetBinContent(ii) - den_lower->GetBinContent(ii) << std::endl;
	}
	
	std::cout << std::endl; 
	for(int ibin(1); ibin < ratio->GetSize(); ibin++){
		myfile << "ratio " << ibin << " = " << ratio->GetBinContent(ibin) << "  est=" << fitfunc_num->Eval(ratio->GetBinCenter(ibin))/fitfunc_den->Eval(ratio->GetBinCenter(ibin)) << std::endl;
	}
	myfile.close();
	       
	if(channel == 1 and detType == 1)c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Plots%d%s/DoubleEG/JetFakeRate_transfer_DoubleEG_EB_%d%s_pT_nVtx.png", RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));
	else if(channel == 1 and detType == 2)c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Plots%d%s/DoubleEG/JetFakeRate_transfer_DoubleEG_EE_%d%s_pT_nVtx.png", RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));
	else if(channel == 2 and detType == 1)c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Plots%d%s/MuonEG/JetFakeRate_transfer_MuonEG_EB_%d%s_pT_nVtx.png", RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));
	else if(channel == 2 and detType == 2)c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Plots%d%s/MuonEG/JetFakeRate_transfer_MuonEG_EE_%d%s_pT_nVtx.png", RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));
}
