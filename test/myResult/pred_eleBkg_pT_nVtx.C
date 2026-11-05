#include "../../include/analysis_commoncode.h"
#include <array>
#include "invmass_binning.h"
#define NTOY 1000
bool useGaussFit=false;
//bool dopostfit=true;

// true  -> use the fake rate's nVtx dependence, fit independently per pT
//          slice (30-50 / 50-80 / >80) -- reads
//          EleFakeRate-Data-ByPtVtx-EB_pT<range>.txt (scalefactor/vtxconst/
// false -> use the fake rate's pT dependence, fit independently per nVtx
//          slice (8-20 / 20-30 / >30) -- the method this script used before,
//          kept available for cross-checks.
bool USE_VTX_BY_PT_SLICE = true;

void pred_eleBkg_pT_nVtx(){
	SetSignalConfig();
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);
	setTDRStyle();
  	gSystem->Load("../../lib/libAnaClasses.so");
  	int channelType = ichannel; // eg = 1; mg =2;
        std::string whichVFP;	
        if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
        else if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
        else whichVFP = "";

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
		} else {
			std::cout << "WARNING: could not open inclusive fake rate file " << (elefakeBaseName+".txt") << std::endl;
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
		if(!f.is_open()){
			std::cout << "WARNING: could not open " << (elefakeBaseName+nVtxSliceFileTag[is]+".txt")
			          << " -- events in this nVtx slice will fall back to the inclusive fit" << std::endl;
			continue;
		}
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

// -----------------------------------------------------------------------------
// Per-pT-slice nVtx fits + their OWN scalefactor -- reads the
// EleFakeRate-Data-ByPtVtx-EB_pT<range>.txt files written by
// params are written there, since pT is the fixed slice selector for these
// files rather than a fit variable -- nVtx and eta are kept at their full/
// inclusive range, only pT is restricted). This is the REVERSE slicing
// direction from nVtxSliceFileTag above: here we pick a pT slice from the
// event's own pT, then evaluate that slice's independent nVtx dependence.
// -----------------------------------------------------------------------------
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
		if(!f.is_open()){
			std::cout << "WARNING: could not open " << (elefakeBaseName+pTSliceFileTag[is]+".txt")
			          << " -- events in this pT slice will fall back to the inclusive fit" << std::endl;
			continue;
		}
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
	TH1D *p_PhoEt = new TH1D("p_PhoEt","; p_{T}^{#gamma} (GeV);",nSigEtBins,sigEtBins);
	TH1D *p_LepPt = new TH1D("p_LepPt","p_LepPt",nSigPtBins,sigPtBins);
	TH1D *p_MET = new TH1D("p_MET","; p_{T}^{miss} (GeV);",nSigMETBins, sigMETBins);
	TH1D *p_Mt = new TH1D("p_Mt","; M_{T} (GeV);",nSigMtBins,sigMtBins);
	TH1D *p_HT = new TH1D("p_HT","HT; H_{T} (GeV);",nSigHTBins, sigHTBins); 
	TH1D *p_PhoEta = new TH1D("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *p_LepEta = new TH1D("p_LepEta","p_LepEta",60,-3,3);
	TH1D *p_dPhiEleMET = new TH1D("p_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *p_PU = new TH1D("p_PU","",100,0,100);
	TH1D *p_nJet = new TH1D("p_nJet","p_nJet",10,0,10);

	double predict_lowEt(0), predict_highEt(0), proxy_lowEt(0), proxy_highEt(0);
	double toy_predict_lowEt[NTOY];
	double toy_predict_highEt[NTOY];
	for(unsigned i(0); i<NTOY; i++){
		toy_predict_lowEt[i] = 0;
	  toy_predict_highEt[i] = 0;	
	}

	TH1D *h_elefakepho_norm = 0;
	TH1D *h_elefakepho_controlsample = 0;
	TH1D *h_elefakepho_transferfactor = 0;
	TH1D *h_elefakepho_syserr_jes = 0;
	TH1D *h_elefakepho_syserr_jer = 0;
	TH1D *h_elefakepho_syserr_esf = 0;
	TH1D *h_elefakepho_syserr_scale = 0;
	TH1D *h_elefakepho_syserr_e_to_pho = 0;
	TH1D *h_elefakepho_syserr_j_to_pho = 0;
	TH1D *h_elefakepho_syserr_j_to_lep = 0;
	TH1D *h_elefakepho_syserr_xs = 0;
	TH1D *h_elefakepho_syserr_lumi = 0;
	TH1D *h_elefakepho_syserr_isr = 0;
	if(channelType==1){
		h_elefakepho_norm            = new TH1D("eg_elefakepho_norm","eventcount",NBIN,0,NBIN);
		h_elefakepho_controlsample   = new TH1D("eg_elefakepho_controlsample","",NBIN,0,NBIN);
		h_elefakepho_transferfactor  = new TH1D("eg_elefakepho_transferfactor","",NBIN,0,NBIN);
		h_elefakepho_syserr_jes      = new TH1D("eg_elefakepho_syserr_jes","",NBIN,0,NBIN);	
		h_elefakepho_syserr_jer      = new TH1D("eg_elefakepho_syserr_jer","",NBIN,0,NBIN);	
		h_elefakepho_syserr_esf      = new TH1D("eg_elefakepho_syserr_esf","",NBIN,0,NBIN);	
		h_elefakepho_syserr_scale    = new TH1D("eg_elefakepho_syserr_scale","",NBIN,0,NBIN);	
		h_elefakepho_syserr_e_to_pho = new TH1D("eg_elefakepho_syserr_e_to_pho","",NBIN,0,NBIN);	
		h_elefakepho_syserr_j_to_pho = new TH1D("eg_elefakepho_syserr_j_to_pho","",NBIN,0,NBIN);	
		h_elefakepho_syserr_j_to_lep = new TH1D("eg_elefakepho_syserr_j_to_lep","",NBIN,0,NBIN);	
		h_elefakepho_syserr_xs       = new TH1D("eg_elefakepho_syserr_xs","",NBIN,0,NBIN);	
		h_elefakepho_syserr_lumi     = new TH1D("eg_elefakepho_syserr_lumi","",NBIN,0,NBIN);	
		h_elefakepho_syserr_isr      = new TH1D("eg_elefakepho_syserr_isr","",NBIN,0,NBIN); 
	} 
	else if(channelType==2){
		h_elefakepho_norm            = new TH1D("mg_elefakepho_norm","eventcount",NBIN,0,NBIN);
		h_elefakepho_controlsample   = new TH1D("mg_elefakepho_controlsample","",NBIN,0,NBIN);
		h_elefakepho_transferfactor  = new TH1D("mg_elefakepho_transferfactor","",NBIN,0,NBIN);
		h_elefakepho_syserr_jes      = new TH1D("mg_elefakepho_syserr_jes","",NBIN,0,NBIN);	
		h_elefakepho_syserr_jer      = new TH1D("mg_elefakepho_syserr_jer","",NBIN,0,NBIN);	
		h_elefakepho_syserr_esf      = new TH1D("mg_elefakepho_syserr_esf","",NBIN,0,NBIN);	
		h_elefakepho_syserr_scale    = new TH1D("mg_elefakepho_syserr_scale","",NBIN,0,NBIN);	
		h_elefakepho_syserr_e_to_pho = new TH1D("mg_elefakepho_syserr_e_to_pho","",NBIN,0,NBIN);	
		h_elefakepho_syserr_j_to_pho = new TH1D("mg_elefakepho_syserr_j_to_pho","",NBIN,0,NBIN);	
		h_elefakepho_syserr_j_to_lep = new TH1D("mg_elefakepho_syserr_j_to_lep","",NBIN,0,NBIN);	
		h_elefakepho_syserr_xs       = new TH1D("mg_elefakepho_syserr_xs","",NBIN,0,NBIN);	
		h_elefakepho_syserr_lumi     = new TH1D("mg_elefakepho_syserr_lumi","",NBIN,0,NBIN);	
		h_elefakepho_syserr_isr      = new TH1D("mg_elefakepho_syserr_isr","",NBIN,0,NBIN); 
	} 
	h_elefakepho_controlsample->SetBinErrorOption(TH1::kPoisson);

	TH1D *toy_PhoEt[NTOY];
	TH1D *toy_LepPt[NTOY];
	TH1D *toy_HT[NTOY];
	TH1D *toy_MET[NTOY];
	TH1D *toy_Mt[NTOY];
	TH1D *toy_dPhiEleMET[NTOY];
	TH1D *toy_eventcount[NTOY];
	for(unsigned ih(0); ih < NTOY; ih++){
		histname.str("");
		histname << "toy_PhoEt_ " << ih;
		toy_PhoEt[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nSigEtBins,sigEtBins);
		histname.str("");
		histname << "toy_LepPt_" << ih;
		toy_LepPt[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nSigPtBins,sigPtBins);
		histname.str("");
		histname << "toy_MET_ " << ih;
		toy_MET[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nSigMETBins, sigMETBins);
		histname.str("");
		histname << "toy_Mt_ " << ih;
		toy_Mt[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nSigMtBins,sigMtBins);
		histname.str("");
		histname << "toy_HT_ " << ih;
		toy_HT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nSigHTBins, sigHTBins);
		histname.str("");
		histname << "toy_eledPhiEleMET_" << ih;
		toy_dPhiEleMET[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),32,0,3.2);
		histname.str("");
		histname << "toy_eventcount_" << ih;
		toy_eventcount[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),NBIN,0,NBIN);
	}

	// nominal yield + 1000 toys of the fake-rate parameters, filled for the same
	// events as h_elefakepho_norm / toy_eventcount (SR bin found)
	const std::string chPrefix = (channelType==1) ? "eg" : "mg";
	TH1D *h_elefakepho_invmass_norm          = new TH1D((chPrefix+"_elefakepho_invmass_norm").c_str(),"m_{l#gamma};m_{l#gamma} (GeV);Events / bin",nInvMassBins,invMassMin,invMassMax);
	TH1D *h_elefakepho_invmass_controlsample = new TH1D((chPrefix+"_elefakepho_invmass_controlsample").c_str(),"",nInvMassBins,invMassMin,invMassMax);
	TH1D *h_elefakepho_invmass_syserr_e_to_pho = new TH1D((chPrefix+"_elefakepho_invmass_syserr_e_to_pho").c_str(),"",nInvMassBins,invMassMin,invMassMax);
	h_elefakepho_invmass_norm->Sumw2();
	h_elefakepho_invmass_controlsample->SetBinErrorOption(TH1::kPoisson);
	TH1D *toy_invmass[NTOY];
	for(unsigned ih(0); ih < NTOY; ih++){
		histname.str("");
		histname << "toy_invmass_" << ih;
		toy_invmass[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nInvMassBins,invMassMin,invMassMax);
	}
	//************ Proxy Tree *********************
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
	float invmass(0);
	
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
	if (channelType == 1) proxytree->SetBranchAddress("nJet", &nJetFloat);
        else proxytree->SetBranchAddress("nJet", &nJetInt);
	const char *invmassBranch = proxytree->GetBranch("eg_invmass") ? "eg_invmass" : (proxytree->GetBranch("mg_invmass") ? "mg_invmass" : "invmass");
        proxytree->SetBranchAddress(invmassBranch, &invmass);
	
	for (unsigned ievt(0); ievt<proxytree->GetEntries(); ++ievt){//loop on entries
		proxytree->GetEntry(ievt);
		
		p_PU->Fill(nVertex);
		/** cut flow *****/
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
		if(sigMT < lowMt)continue;
		if(highMt > 0 && sigMT > highMt)continue;
		if(lepPt < lowPt)continue;
		if(highPt > 0 && lepPt > highPt)continue;

		double w_ele = USE_VTX_BY_PT_SLICE ? EvalFakeRate_VtxByPtSlice(phoEt, nVertex, fabs(phoEta))
		                                    : EvalFakeRate(phoEt, nVertex, fabs(phoEta));
		p_MET->Fill(sigMET, w_ele);

		for(unsigned it(0); it < NTOY; it++){
			double toy_ele = ToyFakeRate(it, phoEt, nVertex, fabs(phoEta));
			toy_MET[it]->Fill(sigMET, toy_ele);
		}

		/** cut flow *****/
		if(sigMET < lowMET)continue;
		if(highMET > 0 && sigMET > highMET)continue;

		w_ele = USE_VTX_BY_PT_SLICE ? EvalFakeRate_VtxByPtSlice(phoEt, nVertex, fabs(phoEta))
		                            : EvalFakeRate(phoEt, nVertex, fabs(phoEta));
		p_PhoEt->Fill(phoEt,w_ele);
		p_PhoEta->Fill(phoEta, w_ele);
		p_Mt->Fill(sigMT, w_ele);
		p_HT->Fill(HT, w_ele);
		p_LepPt->Fill(lepPt, w_ele);
		p_LepEta->Fill(lepEta, w_ele);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), w_ele);
		if (channelType == 1) p_nJet->Fill(nJetFloat, w_ele);
                if (channelType == 2) p_nJet->Fill(nJetInt, w_ele);

		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET, HT, phoEt);
		float m_lg = invmass;
		//if(SigBinIndex >=0){
			h_elefakepho_norm->Fill( SigBinIndex, w_ele);
			h_elefakepho_controlsample->Fill( SigBinIndex );
			h_elefakepho_invmass_norm->Fill( m_lg, w_ele);
			h_elefakepho_invmass_controlsample->Fill( m_lg );
			if(phoEt < 200){
				predict_lowEt += w_ele;
				proxy_lowEt   += 1;
			}
			else{
				predict_highEt += w_ele;
				proxy_highEt   += 1;
			}
		//}	

		for(unsigned it(0); it < NTOY; it++){
			double toy_ele = ToyFakeRate(it, phoEt, nVertex, fabs(phoEta));
			toy_PhoEt[it]->Fill(phoEt,toy_ele);
			toy_Mt[it]->Fill(sigMT, toy_ele);
			toy_HT[it]->Fill(HT, toy_ele);
			toy_LepPt[it]->Fill(lepPt, toy_ele);
			toy_dPhiEleMET[it]->Fill(fabs(dPhiLepMET), toy_ele);
			if( SigBinIndex >= 0){
				toy_eventcount[it]->Fill(SigBinIndex, toy_ele);
				toy_invmass[it]->Fill(m_lg, toy_ele);
				if(phoEt < 200){
					toy_predict_lowEt[it] += toy_ele;
				}
				else{
					toy_predict_highEt[it] += toy_ele;
				}
			}
		}
	}

	std::vector<double> toyvec; 
	for(int ibin(1); ibin < p_PhoEt->GetSize(); ibin++){
		toyvec.clear();
		toyvec.push_back(p_PhoEt->GetBinContent(ibin));
		for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_PhoEt[it]->GetBinContent(ibin));
		double syserr = calcToyError( toyvec, useGaussFit, channelType); 
		double totalerror = sqrt(syserr*syserr + p_PhoEt->GetBinError(ibin)*p_PhoEt->GetBinError(ibin));
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
	for(int ibin(1); ibin < h_elefakepho_controlsample->GetSize(); ibin++){
		if(h_elefakepho_controlsample->GetBinContent(ibin) > 0){
			toyvec.clear();
			toyvec.push_back(h_elefakepho_norm->GetBinContent(ibin));
			for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_eventcount[it]->GetBinContent(ibin));
			double syserr = calcToyError( toyvec, useGaussFit, channelType);

			h_elefakepho_transferfactor->SetBinContent(ibin,h_elefakepho_norm->GetBinContent(ibin)/h_elefakepho_controlsample->GetBinContent(ibin)); 
			h_elefakepho_transferfactor->SetBinError(ibin,syserr/h_elefakepho_norm->GetBinContent(ibin));
			h_elefakepho_syserr_jes->SetBinContent(ibin, -1);
			h_elefakepho_syserr_jer->SetBinContent(ibin, -1);
			h_elefakepho_syserr_esf->SetBinContent(ibin, -1);
			h_elefakepho_syserr_scale->SetBinContent(ibin, -1);
			h_elefakepho_syserr_e_to_pho->SetBinContent(ibin, syserr);
			h_elefakepho_syserr_j_to_pho->SetBinContent(ibin, -1);
			h_elefakepho_syserr_j_to_lep->SetBinContent(ibin, -1);
			h_elefakepho_syserr_xs->SetBinContent(ibin, -1);     
			h_elefakepho_syserr_lumi->SetBinContent(ibin, -1);      
			h_elefakepho_syserr_isr->SetBinContent(ibin, -1);      
		}
		else{	
			double sf(1);
			if(ibin <= 9)sf = predict_lowEt/proxy_lowEt;
			else sf = predict_highEt/proxy_highEt;
			toyvec.clear();
			toyvec.push_back(sf);
			for(unsigned it(0); it < NTOY; it++)toyvec.push_back( ibin <= 9? toy_predict_lowEt[it]/proxy_lowEt: toy_predict_highEt[it]/proxy_highEt);
			double syserr = calcToyError( toyvec, useGaussFit, channelType);
			std::cout << "syserr " << ibin << " " << syserr << " proxy " << proxy_highEt << std::endl;

			h_elefakepho_transferfactor->SetBinContent(ibin, sf); 
			h_elefakepho_transferfactor->SetBinError(ibin,   syserr);
			h_elefakepho_syserr_jes->SetBinContent(ibin, -1);
			h_elefakepho_syserr_jer->SetBinContent(ibin, -1);
			h_elefakepho_syserr_esf->SetBinContent(ibin, -1);
			h_elefakepho_syserr_scale->SetBinContent(ibin, -1);
			h_elefakepho_syserr_e_to_pho->SetBinContent(ibin, -1); 
			h_elefakepho_syserr_j_to_pho->SetBinContent(ibin, -1);
			h_elefakepho_syserr_j_to_lep->SetBinContent(ibin, -1);
			h_elefakepho_syserr_xs->SetBinContent(ibin, -1);     
			h_elefakepho_syserr_lumi->SetBinContent(ibin, -1);      
			h_elefakepho_syserr_isr->SetBinContent(ibin, -1);      
		}
	}
		
	// per-mass-bin e->gamma systematic: spread of the 1000 fake-rate toys, same
	// calcToyError() call as for the signal bins. Bins without any proxy event
	// carry no prediction (and no error), flagged with -1 like the SR bins.
	for(int ibin(1); ibin <= nInvMassBins; ibin++){
		if(h_elefakepho_invmass_controlsample->GetBinContent(ibin) > 0){
			toyvec.clear();
			toyvec.push_back(h_elefakepho_invmass_norm->GetBinContent(ibin));
			for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_invmass[it]->GetBinContent(ibin));
			double syserr = calcToyError( toyvec, useGaussFit, channelType);
			h_elefakepho_invmass_syserr_e_to_pho->SetBinContent(ibin, syserr);
		}
		else h_elefakepho_invmass_syserr_e_to_pho->SetBinContent(ibin, -1);
	}

	std::ostringstream outputname;
	outputname << "/uscms_data/d3/tmishra/Output/";
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_eleBkg";
	else if(channelType==2)outputname << "mg_eleBkg";
	if(anatype ==0)outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt;
	outputname << "_" <<RunYear<< whichVFP<<".root";

	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
	p_PhoEt->Write();
	p_PhoEta->Write();
	p_LepPt->Write();
	p_LepEta->Write();
	p_MET->Write();
	p_Mt->Write();
	p_HT->Write();
	p_dPhiEleMET->Write();
	p_PU->Write();
	p_nJet->Write();
	h_elefakepho_norm->Write();
	h_elefakepho_controlsample->Write();
	h_elefakepho_transferfactor->Write();
	h_elefakepho_syserr_jes->Write();       
	h_elefakepho_syserr_jer->Write();       
	h_elefakepho_syserr_esf->Write();       
	h_elefakepho_syserr_scale->Write();     
	h_elefakepho_syserr_e_to_pho->Write();  
	h_elefakepho_syserr_j_to_pho->Write(); 
	h_elefakepho_syserr_j_to_lep->Write(); 
	h_elefakepho_syserr_xs->Write();        
	h_elefakepho_syserr_lumi->Write();      
	h_elefakepho_syserr_isr->Write();
	h_elefakepho_invmass_norm->Write();
	h_elefakepho_invmass_controlsample->Write();
	h_elefakepho_invmass_syserr_e_to_pho->Write();
	for(unsigned it(0); it < NTOY; it++){
	//	toy_PhoEt[it]->Write();
	//	toy_MET[it]->Write();
	//	toy_Mt[it]->Write();
	//	toy_LepPt[it]->Write();
	//	toy_HT[it]->Write();
	//	toy_dPhiEleMET[it]->Write();
	}
	outputfile->Write();
	outputfile->Close();
}
