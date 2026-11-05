#include<string>
#include<array>
#include "../../include/analysis_commoncode.h"

#define NTOY 1000
bool useGaussFit=false;

bool USE_VTX_BY_PT_SLICE = true;

void analysis_eleBkg_pT_nVtx(){

	SetRunConfig();
	setTDRStyle();

  	gSystem->Load("../../lib/libAnaClasses.so");
  	int channelType = ichannel; // eg = 1; mg =2;

  	std::string whichVFP;
  	if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
  	if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
  	if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

	float correction =1.0;
	
	  if (RunYear == 2016 and preVFP == 1 and ichannel == 1) correction = 1.00158;
          else if (RunYear == 2016 and preVFP == 1 and ichannel == 2) correction = 0.912476;
          else if (RunYear == 2016 and preVFP == 0 and ichannel == 1) correction = 0.971904;
          else if (RunYear == 2016 and preVFP == 0 and ichannel == 2) correction = 0.914534;
          else if (RunYear == 2017 and ichannel == 1) correction = 1.07788;
          else if (RunYear == 2017 and ichannel == 2) correction = 0.806638;
          else if (RunYear == 2018 and ichannel == 1) correction = 1.08004;
          else if (RunYear == 2018 and ichannel == 2) correction = 0.971525;
	
  	/************************************/
	/*	double normfactor = par[0]; */  
  	/*	double slope = par[1];				*/
  	/*	double constant = par[2];			*/
  	/*	double index = par[3];				*/
  	/*	double coeff = par[4]; 				*/
	/*	double vtx_constant = par[5];	*/
	/*	double vtx_slope = par[6];		*/
  	/**********************************/
// =========================================================================
// rate lookup: R(pT,eta,Nvtx) = scalefactor_nVtxSlice * correction *
// ptFitForEvent(pT) * [weight_eta / avg_fr_eta], where ptFitForEvent(pT) is
// the pT-dependence curve fit independently to the matching nVtx slice
// (8-20 / 20-30 / >30), NOT the old separable
// R = scalefactor*correction*f(pT)*f(Nvtx)*f(eta) that the shared
// fakerate_func()/TF3 h_nominal_fakerate machinery in analysis_fakes.h
// implements. Same reading/weight-construction logic as

	auto fakerate_ptDependence_pred = [](double pt, double slope, double constant, double index) -> double {
		double p = TMath::Max(pt, 0.000001);
		double arg = slope*p + constant;
		if(arg <= 0) return 0.0;
		return pow(arg, index);
	};

	std::string elefakeBaseName = Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-Data-ByPtVtx-EB",RunYear,whichVFP.c_str());

	// Inclusive fit: used as the fallback pT curve + scalefactor for any
	// event whose nVtx slice fit didn't converge (same fallback fit script
	// uses via SelectPtFitForNVtx -> fit_fakerate_pt).
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

	// Returns {scalefactor, slope, constant, index} for a given nVertex,
	// falling back to the inclusive fit if that slice's fit didn't converge
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

	// Returns {scalefactor, vtxconst, vtxslope} for a given pT, falling back
	// to the inclusive fit if that pT slice's file/fit wasn't available.
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
	} else {
		std::cout << "WARNING: could not open/parse eta file -- eta shape correction will be a no-op (=1)" << std::endl;
	}
	Eta_file.close();
	std::cout << "avg_fr_eta (mean fake rate over eta bins) = " << avg_fr_eta << std::endl;

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

	// R(pT, eta, Nvtx) = scalefactor_pTSlice(pT) * correction *
	//                    [vtxconst_pTSlice(pT) + vtxslope_pTSlice(pT)*Nvtx] * [weight_eta(eta)/avg_fr_eta]
	// This is the requested method: fake rate dependence on nVtx, fit
	// independently per pT slice, rather than the reverse above.
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

	// Factorized-model toys (ToyFakeRate_Data_EB.txt): used only when USE_VTX_BY_PT_SLICE is false.
	// The VtxByPtSlice toys are read below and evaluated by ToyFakeRate().
	TF3 *h_toymc_fakerate[NTOY];
	std::ostringstream funcname;
	std::ifstream elefake_toyfile;
	elefake_toyfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/ToyFakeRate_Data_EB.txt",RunYear,whichVFP.c_str()));
	if(elefake_toyfile.is_open()){
		double toy_scalefactor, toy_ptslope, toy_ptconstant, toy_ptindex, toy_vtxconst, toy_vtxslope;
  	for(int i(0); i<NTOY; i++){ 
			elefake_toyfile >> toy_scalefactor >> toy_ptslope >> toy_ptconstant >> toy_ptindex >> toy_vtxconst >> toy_vtxslope;
			funcname.str("");
			funcname << "h_toymc_fakerate_" << i;
			h_toymc_fakerate[i] = new TF3(funcname.str().c_str(), fakerate_func,10,1000,0,100,0,1.5,7);
			h_toymc_fakerate[i]->SetParameters(toy_scalefactor*correction, toy_ptslope, toy_ptconstant, toy_ptindex, 1, toy_vtxconst, toy_vtxslope); 
  	}
	}
	elefake_toyfile.close();

	// Toys for the VtxByPtSlice model (written by eleFakePho/fit/fitFakeFunc_pT_nVtx.C):
	// per line [scalefactor vtxconst vtxslope ptshapeavg] x 3 pT slices, then ptslope ptconstant ptindex.
	// Used instead of the factorized-model toys above when USE_VTX_BY_PT_SLICE is set.
	std::vector<std::array<double,15> > toySlicePar;
	if(USE_VTX_BY_PT_SLICE){
		std::ifstream toySliceFile(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/ToyFakeRate_Data_EB_VtxByPtSlice.txt",RunYear,whichVFP.c_str()));
		std::array<double,15> p;
		while(toySlicePar.size() < NTOY){
			bool okRead = true;
			for(int k=0; k<15 && okRead; k++) okRead = static_cast<bool>(toySliceFile >> p[k]);
			if(!okRead) break;
			toySlicePar.push_back(p);
		}
		if(toySlicePar.size() < NTOY){
			std::cout << "ERROR: only " << toySlicePar.size() << " VtxByPtSlice toys found (need " << NTOY
			          << "); rerun eleFakePho/fit/runFakeFunc.sh" << std::endl;
			return;
		}
	}
	// toy fake rate: same formula as EvalFakeRate_VtxByPtSlice with the toy parameters
	auto ToyFakeRate = [&](int it, double pt, int nVertexVal, double eta) -> double {
		if(!USE_VTX_BY_PT_SLICE) return h_toymc_fakerate[it]->Eval(pt, nVertexVal, eta);
		const std::array<double,15> &p = toySlicePar[it];
		int s = SelectPtSliceIdx(pt);
		if(s < 0) s = 2;
		double sf = p[4*s], vtxconst = p[4*s+1], vtxslope = p[4*s+2], ptavg = p[4*s+3];
		double residualPt = (ptavg > 0) ? fakerate_ptDependence_pred(pt, p[12], p[13], p[14]) / ptavg : 1.0;
		double etaCorr = (avg_fr_eta > 0) ? (WeightEta(eta) / avg_fr_eta) : 1.0;
		return sf * correction * (vtxconst + vtxslope*nVertexVal) * residualPt * etaCorr;
	};

	//*********** histo list ********************** 
	std::ostringstream histname;
	TH1D *p_PhoEt = new TH1D("p_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *p_LepPt = new TH1D("p_LepPt","p_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *p_MET = new TH1D("p_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *p_Mt = new TH1D("p_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *p_HT = new TH1D("p_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *p_PhoEta = new TH1D("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *p_LepEta = new TH1D("p_LepEta","p_LepEta",60,-3,3);
	TH1D *p_dPhiEleMET = new TH1D("p_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *p_PU = new TH1D("p_PU","",100,0,100);
	TH1D *p_nJet = new TH1D("p_nJet","p_nJet",10,0,10);
	TH1D *p_nBJet = new TH1D("p_nBJet","p_nBJet",5,0,5);

	TH1D *p_LepPt_TT = new TH1D("p_LepPt_TT","p_LepPt",nBkgPtBins,bkgPtBins);
        TH1D *p_nJet_TT = new TH1D("p_nJet_TT","p_nJet",10,0,10);
        TH1D *p_nBJet_TT = new TH1D("p_nBJet_TT","p_nBJet",5,0,5);
	TH1D *p_PhoEt_TT = new TH1D("p_PhoEt_TT","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *p_MET_TT = new TH1D("p_MET_TT","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *p_Mt_TT = new TH1D("p_Mt_TT","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *p_HT_TT = new TH1D("p_HT_TT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *p_dPhiEleMET_TT = new TH1D("p_dPhiEleMET_TT","dPhiEleMET",32,0,3.2);

	const int nHTSlices = 3;
        double HTsliceMin[nHTSlices] = {0.,   100., 400.};
        double HTsliceMax[nHTSlices] = {100., 400., 1e6};
        TH2D *p_MET_vs_HTslice = new TH2D("p_MET_vs_HTslice", "MET in HT slices; MET (GeV); HT slice", nBkgMETBins, bkgMETBins, nHTSlices, 0, nHTSlices);
        p_MET_vs_HTslice->GetYaxis()->SetBinLabel(1,"HT: 0-100");
        p_MET_vs_HTslice->GetYaxis()->SetBinLabel(2,"HT: 100-400");
        p_MET_vs_HTslice->GetYaxis()->SetBinLabel(3,"HT: >400");

	const int nMETSlices = 3;
        double METsliceMin[nMETSlices] = {120., 200., 400.};
        double METsliceMax[nMETSlices] = {200., 400., 1e6};
        TH2D *p_HT_vs_METslice = new TH2D("p_HT_vs_METslice", "HT in MET slices; H_{T} (GeV); MET slice", nBkgHTBins, bkgHTBins, nMETSlices, 0, nMETSlices);
        p_HT_vs_METslice->GetYaxis()->SetBinLabel(1,"MET: 120-200");
        p_HT_vs_METslice->GetYaxis()->SetBinLabel(2,"MET: 200-400");
        p_HT_vs_METslice->GetYaxis()->SetBinLabel(3,"MET: >400");

	TH1D *toy_PhoEt[NTOY];
	TH1D *toy_LepPt[NTOY];
	TH1D *toy_HT[NTOY];
	TH1D *toy_MET[NTOY];
	TH1D *toy_MET_HT[NTOY];
	TH1D *toy_HT_MET[NTOY];
	TH1D *toy_Mt[NTOY];
	TH1D *toy_dPhiEleMET[NTOY];

	TH1D *toy_PhoEt_TT[NTOY];
	TH1D *toy_HT_TT[NTOY];
	TH1D *toy_MET_TT[NTOY];
	TH1D *toy_Mt_TT[NTOY];
	TH1D *toy_dPhiEleMET_TT[NTOY];
	for(unsigned ih(0); ih < NTOY; ih++){
		histname.str("");
		histname << "toy_PhoEt_ " << ih;
		toy_PhoEt[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgEtBins,bkgEtBins);
		histname.str("");
		histname << "toy_LepPt_" << ih;
		toy_LepPt[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgPtBins,bkgPtBins);
		histname.str("");
		histname << "toy_MET_ " << ih;
		toy_MET[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgMETBins, bkgMETBins);
		histname.str("");
		histname << "toy_MET_HT_ " << ih;
		toy_MET_HT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(), 3*nBkgMETBins, 0, 3*nBkgMETBins);
		histname.str("");
		histname << "toy_HT_MET_ " << ih;
		toy_HT_MET[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(), 3*nBkgHTBins, 0, 3*nBkgHTBins);
		histname.str("");
		histname << "toy_Mt_ " << ih;
		toy_Mt[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgMtBins,bkgMtBins);
		histname.str("");
		histname << "toy_HT_ " << ih;
		toy_HT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgHTBins, bkgHTBins);
		histname.str("");
		histname << "toy_eledPhiEleMET_" << ih;
		toy_dPhiEleMET[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),32,0,3.2);
		histname.str("");
		histname << "toy_PhoEt_TT_ " << ih;
		toy_PhoEt_TT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgEtBins,bkgEtBins);
		histname.str("");
		histname << "toy_MET_TT_" << ih;
		toy_MET_TT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgMETBins, bkgMETBins);
		histname.str("");
		histname << "toy_HT_TT_" << ih;
		toy_HT_TT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgHTBins, bkgHTBins);
		histname.str("");
		histname << "toy_Mt_TT_" << ih;
		toy_Mt_TT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgMtBins,bkgMtBins);
		histname.str("");
		histname << "toy_eledPhiEleMET_TT_" << ih;
		toy_dPhiEleMET_TT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),32,0,3.2);
	}
	//************ Proxy Tree **********************//
	// background estimated from data, with proxyTree
	TChain *proxytree = new TChain("proxyTree");
                if(channelType==1)proxytree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s.root",RunYear,whichVFP.c_str()));
                if(channelType==2)proxytree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s.root",RunYear,whichVFP.c_str()));
	float phoEt(0);
	float phoEta(0);
	float phoPhi(0);
	float lepPt(0);
	float lepEta(0);
	float lepPhi(0);
	float sigMT(0);
	float sigMET(0);
	float sigMETPhi(0);
	float dPhiLepMET(0);
	int   nVertex(0);
	float dRPhoLep(0);
	float HT(0);
	float nJetFloat(0);
	int nJetInt(0);

	int   nBJet(0);	


	proxytree->SetBranchAddress("phoEt",     &phoEt);
	proxytree->SetBranchAddress("phoEta",    &phoEta);
	proxytree->SetBranchAddress("phoPhi",    &phoPhi);
	proxytree->SetBranchAddress("lepPt",     &lepPt);
	proxytree->SetBranchAddress("lepEta",    &lepEta);
	proxytree->SetBranchAddress("lepPhi",    &lepPhi);
	proxytree->SetBranchAddress("sigMT",     &sigMT);
	proxytree->SetBranchAddress("sigMET",    &sigMET);
	proxytree->SetBranchAddress("sigMETPhi", &sigMETPhi);
	proxytree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
	proxytree->SetBranchAddress("nVertex",   &nVertex);
	proxytree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
	proxytree->SetBranchAddress("HT",        &HT);
	proxytree->SetBranchAddress("nBJet",     &nBJet);
	if (channelType == 1) proxytree->SetBranchAddress("nJet", &nJetFloat);
	else proxytree->SetBranchAddress("nJet", &nJetInt);

	for (unsigned ievt(0); ievt<proxytree->GetEntries(); ++ievt){//loop on entries
		proxytree->GetEntry(ievt);
		p_PU->Fill(nVertex);
		/** cut flow *****/
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
		if(sigMET < lowMET)continue;
		if(highMET > 0 && sigMET > highMET)continue;
		if(sigMT < lowMt)continue;
		if(highMt > 0 && sigMT > highMt)continue;
		if(lepPt < lowPt)continue;
		if(highPt > 0 && lepPt > highPt)continue;
		if(HT < lowHT)continue;
                if(highHT > 0 && HT > highHT)continue;
                if(phoEt < lowPhoEt)continue;
                if(highPhoEt > 0 && phoEt > highPhoEt)continue;
		// different MET, MT and lep pT cuts
		double w_ele = USE_VTX_BY_PT_SLICE ? EvalFakeRate_VtxByPtSlice(phoEt, nVertex, fabs(phoEta))
		                                    : EvalFakeRate(phoEt, nVertex, fabs(phoEta));
		// fake rate weight
		p_PhoEt->Fill(phoEt,w_ele);
		p_PhoEta->Fill(phoEta, w_ele);
		p_MET->Fill(sigMET, w_ele);
		p_Mt->Fill(sigMT, w_ele);
		p_HT->Fill(HT, w_ele);
		p_LepPt->Fill(lepPt, w_ele);
		p_LepEta->Fill(lepEta, w_ele);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), w_ele);
		if (channelType == 1) p_nJet->Fill(nJetFloat, w_ele);
                if (channelType == 2) p_nJet->Fill(nJetInt, w_ele);

		//p_nJet->Fill(nJet, w_ele);
		p_nBJet->Fill(nBJet, w_ele);
		for(unsigned it(0); it < NTOY; it++){
			double toy_ele = ToyFakeRate(it, phoEt, nVertex, fabs(phoEta));
			toy_PhoEt[it]->Fill(phoEt,toy_ele);
			toy_MET[it]->Fill(sigMET, toy_ele);
			toy_Mt[it]->Fill(sigMT, toy_ele);
			toy_HT[it]->Fill(HT, toy_ele);
			toy_LepPt[it]->Fill(lepPt, toy_ele);
			toy_dPhiEleMET[it]->Fill(fabs(dPhiLepMET), toy_ele);
		}

		// ttbar events when nBJets >= 1
		if(nBJet >= 1){
			p_PhoEt_TT->Fill(phoEt,  w_ele);
			p_MET_TT->Fill(sigMET,  w_ele);
			p_Mt_TT->Fill(sigMT,  w_ele);
			p_HT_TT->Fill(HT,  w_ele);
			p_LepPt_TT->Fill(lepPt,  w_ele);
                        if (channelType == 1) p_nJet_TT->Fill(nJetFloat,  w_ele);
                        if (channelType == 2) p_nJet_TT->Fill(nJetInt,  w_ele);
                        p_nBJet_TT->Fill(nBJet, w_ele);
			p_dPhiEleMET_TT->Fill(fabs(dPhiLepMET), w_ele);
			for(unsigned it(0); it < NTOY; it++){
				double toy_ele = ToyFakeRate(it, phoEt, nVertex, fabs(phoEta));
				toy_PhoEt_TT[it]->Fill(phoEt,toy_ele);
				toy_MET_TT[it]->Fill(sigMET, toy_ele);
				toy_HT_TT[it]->Fill(HT, toy_ele);
				toy_Mt_TT[it]->Fill(sigMT, toy_ele);
				toy_dPhiEleMET_TT[it]->Fill(fabs(dPhiLepMET), toy_ele);
			}
		}
		if(anatype != 0){
                        int htslice = -1;
                        if(HT < 100)        htslice = 0;
                        else if(HT < 400)   htslice = 1;
                        else                htslice = 2;
                        p_MET_vs_HTslice->Fill(sigMET, htslice + 0.5, w_ele);}

		if (anatype != 0){
                        int metslice = -1;
                        if      (sigMET >= 120 && sigMET < 200) metslice = 0;
                        else if (sigMET >= 200 && sigMET < 400) metslice = 1;
                        else if (sigMET >= 400)                 metslice = 2;
                        p_HT_vs_METslice->Fill(HT, metslice + 0.5, w_ele);}
	}
        
	TH1D *h_MET_HT0 = p_MET_vs_HTslice->ProjectionX("h_MET_HT0",1,1);
        TH1D *h_MET_HT1 = p_MET_vs_HTslice->ProjectionX("h_MET_HT1",2,2);
        TH1D *h_MET_HT2 = p_MET_vs_HTslice->ProjectionX("h_MET_HT2",3,3);
        int nTotBins = nHTSlices * nBkgMETBins;
        TH1D *h_MET_HT_concat = new TH1D("h_MET_HT_concat", "MET in HT slices; MET / HT slice index;Events", nTotBins, 0, nTotBins);
        TH1D *hists[3] = {h_MET_HT0, h_MET_HT1, h_MET_HT2};
        for(int i=0;i<3;i++){
                for(int b=1;b<=nBkgMETBins;b++){
                        int newBin = i*nBkgMETBins + b;
                        h_MET_HT_concat->SetBinContent(newBin, hists[i]->GetBinContent(b));
                        h_MET_HT_concat->SetBinError  (newBin, hists[i]->GetBinError(b));
                }
        }
        delete h_MET_HT0;
        delete h_MET_HT1;
        delete h_MET_HT2;
        delete p_MET_vs_HTslice;

	TH1D *h_HT_MET0 = p_HT_vs_METslice->ProjectionX("h_HT_MET0",1,1);
        TH1D *h_HT_MET1 = p_HT_vs_METslice->ProjectionX("h_HT_MET1",2,2);
        TH1D *h_HT_MET2 = p_HT_vs_METslice->ProjectionX("h_HT_MET2",3,3);
        nTotBins = nMETSlices * nBkgHTBins;
        TH1D *h_HT_MET_concat = new TH1D("h_HT_MET_concat","HT in MET slices; HT / MET slice index; Events",nTotBins, 0, nTotBins);
        TH1D *p_HT_MET[3] = {h_HT_MET0, h_HT_MET1, h_HT_MET2};
        for (int i = 0; i < 3; i++) {
                for (int b = 1; b <= nBkgHTBins; b++) {
                        int newBin = i*nBkgHTBins + b;
                        h_HT_MET_concat->SetBinContent(newBin, p_HT_MET[i]->GetBinContent(b));
                        h_HT_MET_concat->SetBinError  (newBin, p_HT_MET[i]->GetBinError(b));
                }
        }
        delete h_HT_MET0;
        delete h_HT_MET1;
        delete h_HT_MET2;
        delete p_HT_vs_METslice;


	// Error is assigned with Gauss fit
	std::vector<double> toyvec; 

	 for(int i=0;i<3;i++){
                for(int b=1;b<=nBkgMETBins;b++){
                        int newBin = i*nBkgMETBins + b;
                	toyvec.clear();
                	toyvec.push_back(h_MET_HT_concat->GetBinContent(newBin));
                	for(unsigned it(0); it < NTOY; it++)
				toyvec.push_back(toy_MET_HT[it]->GetBinContent(newBin));
                	double syserr = calcToyError( toyvec, useGaussFit, channelType);
                	double totalerror = sqrt(syserr*syserr + h_MET_HT_concat->GetBinError(newBin)*h_MET_HT_concat->GetBinError(newBin));
                	h_MET_HT_concat->SetBinError(newBin, totalerror);
		}
        }

 	for(int i=0;i<3;i++){
    		for(int b=1;b<=nBkgHTBins;b++){
        		int newBin = i*nBkgHTBins + b;
        		toyvec.clear();
        		toyvec.push_back(h_HT_MET_concat->GetBinContent(newBin));
        		for(unsigned it(0); it < NTOY; it++)
            			toyvec.push_back(toy_HT_MET[it]->GetBinContent(newBin));
        		double syserr = calcToyError( toyvec, useGaussFit, channelType);
        		double totalerror = sqrt(syserr*syserr +
            		h_HT_MET_concat->GetBinError(newBin)*h_HT_MET_concat->GetBinError(newBin));
        		h_HT_MET_concat->SetBinError(newBin, totalerror);
    		}
	}

	for(int ibin(1); ibin < p_PhoEt->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_PhoEt->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_PhoEt[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); // systematic error from Toy study
		double totalerror = sqrt(syserr*syserr + p_PhoEt->GetBinError(ibin)*p_PhoEt->GetBinError(ibin)); // stat + syst error
		p_PhoEt->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_LepPt->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_LepPt->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_LepPt[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_LepPt->GetBinError(ibin)*p_LepPt->GetBinError(ibin));
		p_LepPt->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_MET->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_MET->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_MET[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_MET->GetBinError(ibin)*p_MET->GetBinError(ibin));
		p_MET->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_Mt->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_Mt->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_Mt[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_Mt->GetBinError(ibin)*p_Mt->GetBinError(ibin));
		p_Mt->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_HT->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_HT->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_HT[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_HT->GetBinError(ibin)*p_HT->GetBinError(ibin));
		p_HT->SetBinError(ibin, totalerror);
	}
		
	for(int ibin(1); ibin < p_PhoEt_TT->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_PhoEt_TT->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_PhoEt_TT[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_PhoEt_TT->GetBinError(ibin)*p_PhoEt_TT->GetBinError(ibin));
		p_PhoEt_TT->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_MET_TT->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_MET_TT->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_MET_TT[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_MET_TT->GetBinError(ibin)*p_MET_TT->GetBinError(ibin));
		p_MET_TT->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_HT_TT->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_HT_TT->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_HT_TT[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_HT_TT->GetBinError(ibin)*p_HT_TT->GetBinError(ibin));
		p_HT_TT->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_Mt_TT->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_Mt_TT->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_Mt_TT[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_Mt_TT->GetBinError(ibin)*p_Mt_TT->GetBinError(ibin));
		p_Mt_TT->SetBinError(ibin, totalerror);
	}
	for(int ibin(1); ibin < p_dPhiEleMET_TT->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_dPhiEleMET_TT->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_dPhiEleMET_TT[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_dPhiEleMET_TT->GetBinError(ibin)*p_dPhiEleMET_TT->GetBinError(ibin));
		p_dPhiEleMET_TT->SetBinError(ibin, totalerror);
	}
	std::ostringstream outputname;
	outputname << "/uscms_data/d3/tmishra/Background/";
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_eleBkg";
	else if(channelType==2)outputname << "mg_eleBkg";
	if(anatype ==0 or anatype ==1)	outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt;
	bool htRestricted = !(lowHT <= 0 && highHT < 0);
        if((anatype ==0 or anatype ==1) and htRestricted) outputname << "_HT" << lowHT << "_" << highHT;
        if((anatype ==0 or anatype ==1) and !(lowPhoEt <= 0 && highPhoEt < 0)) outputname << "_pho" << lowPhoEt << "_" << highPhoEt;  // photon-pT split (SF binning study)
	outputname <<"_" << RunYear<<whichVFP <<".root";

	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
	p_PhoEt->Write();
	p_PhoEta->Write();
	p_LepPt->Write();
	p_LepEta->Write();
	p_MET->Write();
	h_MET_HT_concat->Write();
	h_HT_MET_concat->Write();
	p_Mt->Write();
	p_HT->Write();
	p_dPhiEleMET->Write();
	p_PU->Write();
	p_nJet->Write();
	p_nBJet->Write();
	p_PhoEt_TT->Write();
	p_dPhiEleMET_TT->Write();
	p_MET_TT->Write();
	p_Mt_TT->Write();
	p_HT_TT->Write();
	p_LepPt_TT->Write();
	p_nJet_TT->Write();
	p_nBJet_TT->Write();
	for(unsigned it(0); it < NTOY; it++){
		toy_dPhiEleMET[it]->Write();
	}
	outputfile->Write();
	outputfile->Close();
}
