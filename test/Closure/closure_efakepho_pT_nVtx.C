	// Also check mcfakerate_func() -- only used for the toy-MC placeholder band now, see NOTE below

	#include "../../include/analysis_commoncode.h"
	#include "TArrow.h"
	#include <array>

	//#define NTOY 1000
	#define NTOY 1000   // same as the production prediction (pred_eleBkg_pT_nVtx.C); was 50
	bool useGaussFit=false;

	// true  -> use the fake rate's nVtx dependence, fit independently per pT
	//          slice (30-50 / 50-80 / >80) -- reads
	//          EleFakeRate-DrellYan-ByPtVtx-EB_pT<range>.txt (scalefactor/
	// false -> use the fake rate's pT dependence, fit independently per nVtx
	//          slice (8-20 / 20-30 / >30) -- the method this script used
	//          before, kept available for cross-checks.
	bool USE_VTX_BY_PT_SLICE = true;

	#define MAXET 199
	#define MAXMT 399
	#define MAXMET 399
	#define MAXHT 399

	bool isHardLepton(int momID){
                if((fabs(momID) >= 0 && fabs(momID) <= 6) || momID==21 || fabs(momID) == 15 || fabs(momID) == 999 || fabs(momID) == 23 || fabs(momID) == 24)return true;
                else return false;
	}

	void closure_efakepho_pT_nVtx(int ichannel, int Year, bool ISpreVFP){
	
	  float correction = 1.0;

	  if (Year == 2016 and ISpreVFP == 1 and ichannel == 1) correction = 1.00158;
	  else if (Year == 2016 and ISpreVFP == 1 and ichannel == 2) correction = 0.912476;
	  else if (Year == 2016 and ISpreVFP == 0 and ichannel == 1) correction = 0.971904;
	  else if (Year == 2016 and ISpreVFP == 0 and ichannel == 2) correction = 0.914534;
	  else if (Year == 2017 and ichannel == 1) correction = 1.07788;
	  else if (Year == 2017 and ichannel == 2) correction = 0.806638;
	  else if (Year == 2018 and ichannel == 1) correction = 1.08004;
	  else if (Year == 2018 and ichannel == 2) correction = 0.971525;
	  gROOT->SetBatch(kTRUE);
	  setTDRStyle();
	  gStyle->SetTitleXOffset(2.5);
	  gSystem->Load("../../lib/libAnaClasses.so");
	  int channelType = ichannel; // eg = 1; mg =2;

		std::string whichVFP;
		if(Year==2016 and ISpreVFP == true) whichVFP = "preVFP";
		if(Year==2016 and ISpreVFP == false) whichVFP = "postVFP";
		if(Year==2017 or  Year == 2018) whichVFP = "";
		// =========================================================================
		// ptFitForEvent(pT) * [weight_eta / avg_fr_eta], where ptFitForEvent(pT) is
		// the pT-dependence curve fit independently to the matching nVtx slice
		// (8-20 / 20-30 / >30), NOT the old separable
		// mcfakerate_func()/TF3 h_nominal_fakerate machinery in analysis_fakes.h
		// implements. Reads the DrellYan-based nVtx-slice output files (same
		// naming convention as the inclusive EleFakeRate-DrellYan-ByPtVtx-EB.txt
		// above, just with a "_nVtx<range>" tag), and rebuilds
		// weight_eta/avg_fr_eta from the same raw eta-bin file
		// piece.
		// =========================================================================

		auto fakerate_ptDependence_pred = [](double pt, double slope, double constant, double index) -> double {
			double p = TMath::Max(pt, 0.000001);
			double arg = slope*p + constant;
			if(arg <= 0) return 0.0;
			return pow(arg, index);
		};

		std::string elefakeBaseName = Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/EleFakeRate-DrellYan-ByPtVtx-EB",Year,whichVFP.c_str());

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

		// Raw per-eta-bin numerator/denominator counts -- SAME file
		// independently so weight_eta and avg_fr_eta are rebuilt exactly the same
		// way rather than depending on whatever eta table is (separately, and
		// possibly stale) hardcoded inside analysis_fakes.h's mcfakerate_func.
		std::ifstream Eta_file;
		Eta_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DYResult%d%s/EleFakeRate-DY-Bw-ker-eta-60-120.txt",Year,whichVFP.c_str()));
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
			std::cout << "WARNING: could not open/parse eta file -- eta shape will be a no-op (=1)" << std::endl;
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

		TF3 *h_toymc_fakerate[NTOY];
		std::ostringstream funcname;
		std::ifstream elefake_toyfile(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/ToyFakeRate_DrellYan_EB.txt",Year,whichVFP.c_str()));

		// toy fake rate byPtVtx estimated from DrellYan (old inclusive scheme, placeholder -- see NOTE above)

		if(!USE_VTX_BY_PT_SLICE && elefake_toyfile.is_open()){   // old inclusive-model toys: only for the cross-check mode
  	for(int i(0); i<NTOY; i++){ 
				double toy_scalefactor, toy_ptslope, toy_ptconstant, toy_ptindex, toy_vtxconst, toy_vtxslope;
				elefake_toyfile >> toy_scalefactor >> toy_ptslope >> toy_ptconstant >> toy_ptindex >>  toy_vtxconst >> toy_vtxslope;
				funcname.str("");
				funcname << "h_toymc_fakerate_" << i;
				h_toymc_fakerate[i] = new TF3(funcname.str().c_str(), mcfakerate_func,10,1000,0,100,0,1.5,7);
				h_toymc_fakerate[i]->SetParameters(toy_scalefactor*correction, toy_ptslope, toy_ptconstant, toy_ptindex, 1, toy_vtxconst, toy_vtxslope); 
  	  }
		}
		elefake_toyfile.close();

		// Toys of the SAME model as the nominal prediction (VtxByPtSlice), as in myResult/pred_eleBkg_pT_nVtx.C.
		// Per line: [scalefactor vtxconst vtxslope ptshapeavg] x 3 pT slices, then ptslope ptconstant ptindex
		// (written by eleFakePho/fit/fitFakeFunc_pT_nVtx.C with useDY=true). Before 2026-11-03 the band used the
		// old inclusive-model toys (ToyFakeRate_DrellYan_EB.txt, degenerate pT slope/constant), so it covered the canvas.
		std::vector<std::array<double,15> > toySlicePar;
		if(USE_VTX_BY_PT_SLICE){
			std::ifstream toySliceFile(Form("/eos/uscms/store/user/tmishra/elefakepho/DrellYanResult%d%s/ToyFakeRate_DrellYan_EB_VtxByPtSlice.txt",Year,whichVFP.c_str()));
			std::array<double,15> tp;
			while(toySlicePar.size() < NTOY){
				bool okRead = true;
				for(int k=0; k<15 && okRead; k++) okRead = static_cast<bool>(toySliceFile >> tp[k]);
				if(!okRead) break;
				toySlicePar.push_back(tp);
			}
			if(toySlicePar.size() < NTOY){
				std::cout << "ERROR: only " << toySlicePar.size() << " DY VtxByPtSlice toys (need " << NTOY << ")" << std::endl;
				return;
			}
		}
		auto ToyFakeRate = [&](int it, double pt, int nVertexVal, double eta) -> double {
			if(!USE_VTX_BY_PT_SLICE) return h_toymc_fakerate[it]->Eval(pt, nVertexVal, eta);
			const std::array<double,15> &tp = toySlicePar[it];
			int sl = SelectPtSliceIdx(pt);
			if(sl < 0) sl = 2;
			double sf = tp[4*sl], vtxconst = tp[4*sl+1], vtxslope = tp[4*sl+2], ptavg = tp[4*sl+3];
			double residualPt = (ptavg > 0) ? fakerate_ptDependence_pred(pt, tp[12], tp[13], tp[14]) / ptavg : 1.0;
			double etaCorr = (avg_fr_eta > 0) ? (WeightEta(eta) / avg_fr_eta) : 1.0;
			return sf * correction * (vtxconst + vtxslope*nVertexVal) * residualPt * etaCorr;
		};


	// Signal Tree //
	//*********** histo list **********************//
	TH1D *p_PhoEt = new TH1D("p_PhoEt",";p_{T}^{#gamma} (GeV);Events / bin",nBkgEtBins,bkgEtBins);
	TH1D *p_LepPt = new TH1D("p_LepPt","",nBkgPtBins,bkgPtBins);
	TH1D *p_MET = new TH1D("p_MET","; p_{T}^{miss} (GeV);Events / bin",nBkgMETBins, bkgMETBins);
	TH1D *p_Mt = new TH1D("p_Mt","; M_{T} (GeV);Events / bin",nBkgMtBins,bkgMtBins);
	TH1D *p_HT = new TH1D("p_HT","; HT (GeV);Events / bin",nBkgHTBins, bkgHTBins); 
	TH1D *p_dPhiEleMET = new TH1D("p_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *p_PhoEta = new TH1D("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *p_LepEta = new TH1D("p_LepEta","p_LepEta",60,-3,3);
	TH1D *p_eventcount = new TH1D("p_eventcount","eventcount",9,0,9);
//	TH1D *p_nJet = new TH1D("p_nJet","p_nJet",10,0,10);
	//************ Signal Tree **********************//
	TChain *sigtree = new TChain("signalTree");
	// processes contribute to electron fake photon background <= direct signal event
	if(channelType==1)sigtree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_DYJetsToLL_%d%s.root",Year,whichVFP.c_str()));
	if(channelType==1)sigtree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_TTJets_%d%s.root",Year,whichVFP.c_str()));
	if(channelType==1)sigtree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_WW_%d%s.root",Year,whichVFP.c_str()));
	if(channelType==1)sigtree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_WZ_%d%s.root",Year,whichVFP.c_str()));

	if(channelType==2)sigtree->Add(Form("/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_DYJetsToLL_%d%s.root",Year,whichVFP.c_str()));
	if(channelType==2)sigtree->Add(Form("/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_TTJets_%d%s.root",Year,whichVFP.c_str()));
	if(channelType==2)sigtree->Add(Form("/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_WW_%d%s.root",Year,whichVFP.c_str()));
	if(channelType==2)sigtree->Add(Form("/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_WZ_%d%s.root",Year,whichVFP.c_str()));

	float crosssection(0);
	float ntotalevent(0);
	float phoEt(0);
	float phoEta(0);
	float phoPhi(0);
	float lepPt(0);
	float lepEta(0);
	float lepPhi(0);
	float sigMT(0);
	float sigMET(0);
	float dPhiLepMET(0);
	float sigMETPhi(0);
	float HT(0);
	int   nVertex(0);
	float dRPhoLep(0);
//	float nJet(0);
  	std::vector<int>   *mcPID=0;
  	std::vector<float> *mcEta=0;
  	std::vector<float> *mcPhi=0;
  	std::vector<float> *mcPt=0;
  	std::vector<int>   *mcMomPID=0;
  	std::vector<int>   *mcGMomPID=0;

	sigtree->SetBranchAddress("crosssection",&crosssection);
	sigtree->SetBranchAddress("ntotalevent", &ntotalevent);
	sigtree->SetBranchAddress("phoEt",     &phoEt);
	sigtree->SetBranchAddress("phoEta",    &phoEta);
	sigtree->SetBranchAddress("phoPhi",    &phoPhi);
	sigtree->SetBranchAddress("lepPt",     &lepPt);
	sigtree->SetBranchAddress("lepEta",    &lepEta);
	sigtree->SetBranchAddress("lepPhi",    &lepPhi);
	sigtree->SetBranchAddress("sigMT",     &sigMT);
	sigtree->SetBranchAddress("sigMET",    &sigMET);
  	sigtree->SetBranchAddress("HT",        &HT);
	sigtree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
	sigtree->SetBranchAddress("sigMETPhi", &sigMETPhi);
	sigtree->SetBranchAddress("nVertex",   &nVertex);
	sigtree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
//	sigtree->SetBranchAddress("nJet",      &nJet);
  	sigtree->SetBranchAddress("mcPID",     &mcPID);
  	sigtree->SetBranchAddress("mcEta",     &mcEta);
  	sigtree->SetBranchAddress("mcPhi",     &mcPhi);
  	sigtree->SetBranchAddress("mcPt",      &mcPt);
  	sigtree->SetBranchAddress("mcMomPID",  &mcMomPID);
  	sigtree->SetBranchAddress("mcGMomPID", &mcGMomPID);

	for (unsigned ievt(0); ievt<sigtree->GetEntries(); ++ievt){//loop on entries
		sigtree->GetEntry(ievt);
		// cross sectional weight
		double weight = 1;
		if(channelType == 1){
                        if(Year == 2016 and ISpreVFP == 1)             weight = lumi_2016preVFP_DoubleEG*1000*crosssection/ntotalevent;
                        else if(Year == 2016 and ISpreVFP == 0)        weight = lumi_2016postVFP_DoubleEG*1000*crosssection/ntotalevent;
                        else if(Year == 2017)                        weight = lumi_2017_DoubleEG*1000*crosssection/ntotalevent;
                        else if(Year == 2018)                        weight = lumi_2018_DoubleEG*1000*crosssection/ntotalevent;}

                else if(channelType == 2){
                        if(Year == 2016 and ISpreVFP == 1)             weight = lumi_2016preVFP_MuonEG*1000*crosssection/ntotalevent;
                        else if(Year == 2016 and ISpreVFP == 0)        weight = lumi_2016postVFP_MuonEG*1000*crosssection/ntotalevent;
                        else if(Year == 2017)                        weight = lumi_2017_MuonEG*1000*crosssection/ntotalevent;
                        else if(Year == 2018)                        weight = lumi_2018_MuonEG*1000*crosssection/ntotalevent;}
		//if(ievt==0) cout<<Year<< " "<<crosssection<<" "<<ntotalevent<<" "<<weight<<endl;
		//cout<<"weight is "<<weight<<endl;
		/** cut flow *****/
		if ((ichannel == 1 && (phoEt < 35 || lepPt < 25)) || (ichannel == 2 && (phoEt < 35 || lepPt < 20))) continue;
		if(fabs(phoEta) > 1.4442 || fabs(lepEta) > 2.5)continue;
//		if(sigMT < 100)continue;
//		if(sigMET < 120)continue;
//		if(sigMET > 70)continue;

		if (Year == 2017 && (phoEta > 1.0 && phoEta < 1.5 && phoPhi > 2.7)) continue; 
		if (Year == 2017 && (lepEta > 1.0 && lepEta < 1.5 && lepPhi > 2.7)) continue; 

		if (Year == 2018 && (phoEta > 0.3 && phoEta < 1.2 && phoPhi > 0.4 && phoPhi < 0.8)) continue; 
		if (Year == 2018 && (lepEta > 0.3 && lepEta < 1.2 && lepPhi > 0.4 && lepPhi < 0.8)) continue; 
		

		bool isFakePho(false);
		double mindR(0.3),deltaE(1);
		double minEledR(3.0);
		unsigned matchIndex(0);
		for(unsigned iMC(0); iMC < mcPID->size(); iMC++){
			double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], phoEta, phoPhi);
			double dE = fabs((*mcPt)[iMC] - phoEt)/phoEt;
			//if(dR < mindR && dE < 0.5){mindR=dR; matchIndex=iMC;deltaE = dE;}
			if(dR < mindR){mindR=dR; matchIndex=iMC;deltaE = dE;}
			if(fabs((*mcPID)[iMC]) == 11 && dR < minEledR)minEledR = dR;
			// find the mc particle (and electron) closest to the signal photon
		}
		if(mindR < 0.1){
			if(((*mcPID)[matchIndex] == 11 || (*mcPID)[matchIndex] == -11))isFakePho = true;
			// If dR(ele,photon) < 0.1, It is a ele-fake photon
		}
		if(minEledR < 0.02)isFakePho = true;

		if(!isFakePho){
			//std::cout << "event " << ievt << std::endl; 
			for(unsigned iMC(0); iMC < mcPID->size(); iMC++){
				double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], phoEta, phoPhi);
				double dE = fabs((*mcPt)[iMC] - phoEt)/phoEt;
				//std::cout << fabs((*mcPID)[iMC]) << " dR " << dR << " dE " << dE << std::endl;
			}
		}
	
		if(!isFakePho)continue;

		// only XSec wt
		p_PhoEt->Fill(phoEt, weight);
		p_PhoEta->Fill(phoEta, weight);
		p_LepPt->Fill(lepPt, weight);
		p_LepEta->Fill(lepEta, weight);
		p_MET->Fill(sigMET, weight);
		p_Mt->Fill(sigMT, weight);
		p_HT->Fill(HT, weight);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), weight);
//		p_nJet->Fill(nJet, weight);
		
	}        
	
	// checking the ele-fake-pho events from direct simulation for each years. 
	TFile *outputFile;
	if (channelType == 1) 
    		outputFile = new TFile(Form("/eos/uscms/store/user/tmishra/elefakepho/Closure/Ele_fake_pho_Direct-Simulation_eg_%d%s.root", Year, whichVFP.c_str()), "RECREATE");

	if (channelType == 2) 
    		outputFile = new TFile(Form("/eos/uscms/store/user/tmishra/elefakepho/Closure/Ele_fake_pho_Direct-Simulation_mg_%d%s.root", Year, whichVFP.c_str()), "RECREATE");


	p_PhoEt->Write();
	p_PhoEta->Write();
	p_LepPt->Write();
	p_LepEta->Write();
	p_MET->Write();
	p_Mt->Write();
	p_HT->Write();
	p_dPhiEleMET->Write();
//	p_nJet->Write();

	outputFile->Close();
	delete outputFile;

	// Fake Tree //
	//*********** histo list **********************//
	std::ostringstream histname;
	TH1D *pred_PhoEt = new TH1D("pred_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *pred_LepPt = new TH1D("pred_LepPt","pred_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *pred_MET = new TH1D("pred_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *pred_Mt = new TH1D("pred_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *pred_HT = new TH1D("pred_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *pred_PhoEta = new TH1D("pred_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *pred_LepEta = new TH1D("pred_LepEta","pred_LepEta",60,-3,3);
	TH1D *pred_dPhiEleMET = new TH1D("pred_dPhiEleMET","dPhiEleMET",32,0,3.2); 
//	TH1D *pred_nJet = new TH1D("pred_nJet","pred_nJet",10,0,10);

	TH1D *DY_PhoEt = new TH1D("DY_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *DY_LepPt = new TH1D("DY_LepPt","DY_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *DY_MET = new TH1D("DY_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *DY_Mt = new TH1D("DY_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *DY_HT = new TH1D("DY_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *DY_PhoEta = new TH1D("DY_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *DY_LepEta = new TH1D("DY_LepEta","DY_LepEta",60,-3,3);
	TH1D *DY_dPhiEleMET = new TH1D("DY_dPhiEleMET","dPhiEleMET",32,0,3.2); 
//	TH1D *DY_nJet = new TH1D("DY_nJet","DY_nJet",10,0,10);

	TH1D *toy_PhoEt[NTOY];
	TH1D *toy_LepPt[NTOY];
	TH1D *toy_HT[NTOY];
	TH1D *toy_MET[NTOY];
	TH1D *toy_Mt[NTOY];
	TH1D *toy_dPhiEleMET[NTOY];
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
		histname << "toy_Mt_ " << ih;
		toy_Mt[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgMtBins,bkgMtBins);
		histname.str("");
		histname << "toy_HT_ " << ih;
		toy_HT[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),nBkgHTBins, bkgHTBins);
		histname.str("");
		histname << "toy_eledPhiEleMET_" << ih;
		toy_dPhiEleMET[ih] = new TH1D(histname.str().c_str(), histname.str().c_str(),32,0,3.2);
	}
	//************ Proxy Tree ********************** 
	TChain *proxytree = new TChain("proxyTree");
	// DY majorly contribute to electron fake photon background, proxy event
	if(channelType==1)proxytree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_DYJetsToLL_%d%s.root",Year,whichVFP.c_str()));
	if(channelType==2)proxytree->Add(Form("/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_DYJetsToLL_%d%s.root",Year,whichVFP.c_str()));

	float proxycrosssection(0);
	float proxyntotalevent(0);
	float proxyphoEt(0);
	float proxyphoEta(0);
	float proxyphoPhi(0);
	float proxylepPt(0);
	float proxylepEta(0);
	float proxylepPhi(0);
	float proxysigMT(0);
	float proxysigMET(0);
	float proxysigMETPhi(0);
	float proxydPhiLepMET(0);
	int   proxynVertex(0);
	float proxydRPhoLep(0);
	float proxyHT(0);
//	float proxynJet(0);

	std::vector<int>   *proxy_mcPID=0;
	std::vector<float> *proxy_mcEta=0;
	std::vector<float> *proxy_mcPhi=0;
	std::vector<float> *proxy_mcPt=0;
	std::vector<int>   *proxy_mcMomPID=0;
	std::vector<int>   *proxy_mcGMomPID=0;

	proxytree->SetBranchAddress("crosssection",&proxycrosssection);
	proxytree->SetBranchAddress("ntotalevent", &proxyntotalevent);
	proxytree->SetBranchAddress("phoEt",     	 &proxyphoEt);
	proxytree->SetBranchAddress("phoEta",    	 &proxyphoEta);
	proxytree->SetBranchAddress("phoPhi",    	 &proxyphoPhi);
	proxytree->SetBranchAddress("lepPt",     	 &proxylepPt);
	proxytree->SetBranchAddress("lepEta",    	 &proxylepEta);
	proxytree->SetBranchAddress("lepPhi",    	 &proxylepPhi);
	proxytree->SetBranchAddress("sigMT",     	 &proxysigMT);
	proxytree->SetBranchAddress("sigMET",    	 &proxysigMET);
	proxytree->SetBranchAddress("sigMETPhi", 	 &proxysigMETPhi);
	proxytree->SetBranchAddress("dPhiLepMET",	 &proxydPhiLepMET);
	proxytree->SetBranchAddress("nVertex",   	 &proxynVertex);
	proxytree->SetBranchAddress("dRPhoLep",  	 &proxydRPhoLep);
	proxytree->SetBranchAddress("HT",        	 &proxyHT);
//	proxytree->SetBranchAddress("nJet",      	 &proxynJet);
	proxytree->SetBranchAddress("mcPID",       &proxy_mcPID);
	proxytree->SetBranchAddress("mcEta",       &proxy_mcEta);
	proxytree->SetBranchAddress("mcPhi",       &proxy_mcPhi);
	proxytree->SetBranchAddress("mcPt",        &proxy_mcPt);
	proxytree->SetBranchAddress("mcMomPID",    &proxy_mcMomPID);
	proxytree->SetBranchAddress("mcGMomPID",   &proxy_mcGMomPID);
	
	// ---- Fake-lepton fraction counters for the DY proxy (ele->photon) sample ----
	//double proxy_wele_total   = 0.0;  // sum of ele-fake weights for all proxy events
	//double proxy_wele_fakelep = 0.0;  // sum of ele-fake weights for events with a fake lepton
	//double proxy_n_total      = 0.0;  // unweighted event count (for cross-check)
	//double proxy_n_fakelep    = 0.0;  // unweighted fake-lepton event count

	cout<<"Year "<<Year<<endl<<endl;
	for (unsigned ievt(0); ievt<proxytree->GetEntries(); ++ievt){//loop on entries
		proxytree->GetEntry(ievt);
		
		double weight = 1;
		if(channelType == 1){
                        if(Year == 2016 and ISpreVFP == 1)             weight = lumi_2016preVFP_DoubleEG*1000*proxycrosssection/proxyntotalevent;
                        else if(Year == 2016 and ISpreVFP == 0)        weight = lumi_2016postVFP_DoubleEG*1000*proxycrosssection/proxyntotalevent;
                        else if(Year == 2017)                        weight = lumi_2017_DoubleEG*1000*proxycrosssection/proxyntotalevent;
                        else if(Year == 2018)                        weight = lumi_2018_DoubleEG*1000*proxycrosssection/proxyntotalevent;}

                else if(channelType == 2){
                        if(Year == 2016 and ISpreVFP == 1)             weight = lumi_2016preVFP_MuonEG*1000*proxycrosssection/proxyntotalevent;
                        else if(Year == 2016 and ISpreVFP == 0)        weight = lumi_2016postVFP_MuonEG*1000*proxycrosssection/proxyntotalevent;
                        else if(Year == 2017)                        weight = lumi_2017_MuonEG*1000*proxycrosssection/proxyntotalevent;
                        else if(Year == 2018)                        weight = lumi_2018_MuonEG*1000*proxycrosssection/proxyntotalevent;}
		
		if(ievt==0) cout<<"DY cross-section "<<proxycrosssection<<endl;

		/** cut flow *****/
		if ((ichannel == 1 && (proxyphoEt < 35 || proxylepPt < 25)) || (ichannel == 2 && (proxyphoEt < 35 || proxylepPt < 20))) continue;
		if(fabs(proxyphoEta) > 1.4442 || fabs(proxylepEta) > 2.5)continue;

		double w_ele = USE_VTX_BY_PT_SLICE ? EvalFakeRate_VtxByPtSlice(proxyphoEt, proxynVertex, fabs(proxyphoEta))
		                                    : EvalFakeRate(proxyphoEt, proxynVertex, fabs(proxyphoEta));
		w_ele = w_ele*weight;
		// weight is XSec weight times eleFakePho weight

//		if(proxysigMT < 100)continue;
//              if(proxysigMET < 120)continue;
//		if(proxysigMET > 70)continue;

		if (Year == 2017 && (proxyphoEta > 1.0 && proxyphoEta < 1.5 && proxyphoPhi > 2.7)) continue;
                if (Year == 2017 && (proxylepEta > 1.0 && proxylepEta < 1.5 && proxylepPhi > 2.7)) continue;

                if (Year == 2018 && (proxyphoEta > 0.3 && proxyphoEta < 1.2 && proxyphoPhi > 0.4 && proxyphoPhi < 0.8)) continue;
                if (Year == 2018 && (proxylepEta > 0.3 && proxylepEta < 1.2 && proxylepPhi > 0.4 && proxylepPhi < 0.8)) continue;

		// ---- Classify whether the lepton in this proxy event is a fake ----
		// A lepton is considered prompt (real) if it MC-matches within dR<0.1
		// to a generator-level lepton (|PID|=11 or 13) whose mother is a
		// W (24), Z (23), or tau (15). Otherwise it is a fake lepton.
		/*bool isRealLep = false;
		double lepMindR = 999.0;
		unsigned lepMatchIdx = 0;
		for(unsigned iMC(0); iMC < proxy_mcPID->size(); iMC++){
			int absPID = fabs((*proxy_mcPID)[iMC]);
			if(absPID != 11 && absPID != 13)continue; // only match to gen leptons
			double dR = DeltaR((*proxy_mcEta)[iMC], (*proxy_mcPhi)[iMC], proxylepEta, proxylepPhi);
			if(dR < lepMindR){ lepMindR = dR; lepMatchIdx = iMC; }
		}
		if(lepMindR < 0.1){
                        if(isHardLepton((*proxy_mcMomPID)[lepMatchIdx])) isRealLep = true;
                }

		bool isFakeLep = !isRealLep;

		// Accumulate fake-lepton fraction counters
		proxy_wele_total   += w_ele;
		proxy_n_total      += 1.0;
		if(isFakeLep){
			proxy_wele_fakelep += w_ele;
			proxy_n_fakelep    += 1.0;
		}*/

		pred_PhoEt->Fill(proxyphoEt,w_ele);
		pred_PhoEta->Fill(proxyphoEta, w_ele);
		pred_MET->Fill(proxysigMET, w_ele);
		pred_Mt->Fill(proxysigMT, w_ele);
		pred_HT->Fill(proxyHT, w_ele);
		pred_LepPt->Fill(proxylepPt, w_ele);
		pred_LepEta->Fill(proxylepEta, w_ele);
		pred_dPhiEleMET->Fill(fabs(proxydPhiLepMET), w_ele);
//		pred_nJet->Fill(proxynJet, w_ele);

		DY_PhoEt->Fill(proxyphoEt,w_ele);
		DY_PhoEta->Fill(proxyphoEta, w_ele);
		DY_MET->Fill(proxysigMET, w_ele);
		DY_Mt->Fill(proxysigMT, w_ele);
		DY_HT->Fill(proxyHT, w_ele);
		DY_LepPt->Fill(proxylepPt, w_ele);
		DY_LepEta->Fill(proxylepEta, w_ele);
		DY_dPhiEleMET->Fill(fabs(proxydPhiLepMET), w_ele);
//		DY_nJet->Fill(proxynJet, w_ele);

		for(unsigned it(0); it < NTOY; it++){
			// here is NTOY fake rates
			double toy_ele = ToyFakeRate(it, proxyphoEt, proxynVertex, fabs(proxyphoEta));
			toy_ele = toy_ele*weight;
			toy_PhoEt[it]->Fill(proxyphoEt,toy_ele);
			toy_MET[it]->Fill(proxysigMET, toy_ele);
			toy_Mt[it]->Fill(proxysigMT, toy_ele);
			toy_HT[it]->Fill(proxyHT, toy_ele);
			toy_LepPt[it]->Fill(proxylepPt, toy_ele);
			toy_dPhiEleMET[it]->Fill(fabs(proxydPhiLepMET), toy_ele);
		}
	}

	// ---- Print fake-lepton fraction for the DY proxy (ele->photon) sample ----
	/*double fakeLepFrac_wele  = (proxy_wele_total  > 0) ? proxy_wele_fakelep / proxy_wele_total  : 0.0;
	double fakeLepFrac_unwtd = (proxy_n_total      > 0) ? proxy_n_fakelep    / proxy_n_total     : 0.0;
	std::cout << "============================================================" << std::endl;
	std::cout << " Ele->photon proxy sample (DY): fake-lepton fraction" << std::endl;
	std::cout << "  Weighted  : " << proxy_wele_fakelep << " / " << proxy_wele_total
	          << "  = " << fakeLepFrac_wele*100.0 << " %" << std::endl;
	std::cout << "  Unweighted: " << proxy_n_fakelep    << " / " << proxy_n_total
	          << "  = " << fakeLepFrac_unwtd*100.0 << " %" << std::endl;
	std::cout << "============================================================" << std::endl;*/
		
	//************ Proxy Tree **********************
	TChain *raretree = new TChain("proxyTree");
	// TTJets, WW, WZ rarely contribute to electron fake photon background, proxy event
	if(channelType==1)raretree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_TTJets_%d%s.root",Year,whichVFP.c_str()));
	if(channelType==1)raretree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_WW_%d%s.root",Year,whichVFP.c_str()));
	if(channelType==1)raretree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_WZ_%d%s.root",Year,whichVFP.c_str()));

	if(channelType==2)raretree->Add(Form("/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_TTJets_%d%s.root",Year,whichVFP.c_str()));
	if(channelType==2)raretree->Add(Form("/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_WW_%d%s.root",Year,whichVFP.c_str()));
	if(channelType==2)raretree->Add(Form("/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_WZ_%d%s.root",Year,whichVFP.c_str()));

	float rarecrosssection(0);
	float rarentotalevent(0);
	float rarephoEt(0);
	float rarephoEta(0);
	float rarephoPhi(0);
	float rarelepPt(0);
	float rarelepEta(0);
	float rarelepPhi(0);
	float raresigMT(0);
	float raresigMET(0);
	float raresigMETPhi(0);
	float raredPhiLepMET(0);
	int   rarenVertex(0);
	float raredRPhoLep(0);
	float rareHT(0);
//	float rarenJet(0);
	// MC truth branches for lepton fake classification in rare proxy sample
	std::vector<int>   *rare_mcPID=0;
	std::vector<float> *rare_mcEta=0;
	std::vector<float> *rare_mcPhi=0;
	std::vector<float> *rare_mcPt=0;
	std::vector<int>   *rare_mcMomPID=0;
	std::vector<int>   *rare_mcGMomPID=0;
		
	raretree->SetBranchAddress("crosssection",&rarecrosssection);
	raretree->SetBranchAddress("ntotalevent", &rarentotalevent);
	raretree->SetBranchAddress("phoEt",     	 &rarephoEt);
	raretree->SetBranchAddress("phoEta",    	 &rarephoEta);
	raretree->SetBranchAddress("phoPhi",    	 &rarephoPhi);
	raretree->SetBranchAddress("lepPt",     	 &rarelepPt);
	raretree->SetBranchAddress("lepEta",    	 &rarelepEta);
	raretree->SetBranchAddress("lepPhi",    	 &rarelepPhi);
	raretree->SetBranchAddress("sigMT",     	 &raresigMT);
	raretree->SetBranchAddress("sigMET",    	 &raresigMET);
	raretree->SetBranchAddress("sigMETPhi", 	 &raresigMETPhi);
	raretree->SetBranchAddress("dPhiLepMET",	 &raredPhiLepMET);
	raretree->SetBranchAddress("nVertex",   	 &rarenVertex);
	raretree->SetBranchAddress("dRPhoLep",  	 &raredRPhoLep);
	raretree->SetBranchAddress("HT",        	 &rareHT);
//	raretree->SetBranchAddress("nJet",      	 &rarenJet);
	raretree->SetBranchAddress("mcPID",        &rare_mcPID);
	raretree->SetBranchAddress("mcEta",        &rare_mcEta);
	raretree->SetBranchAddress("mcPhi",        &rare_mcPhi);
	raretree->SetBranchAddress("mcPt",         &rare_mcPt);
	raretree->SetBranchAddress("mcMomPID",     &rare_mcMomPID);
	raretree->SetBranchAddress("mcGMomPID",    &rare_mcGMomPID);

	// ---- Fake-lepton fraction counters for the rare proxy (TT/WW/WZ) sample ----
	double rare_wele_total   = 0.0;
	double rare_wele_fakelep = 0.0;
	double rare_n_total      = 0.0;
	double rare_n_fakelep    = 0.0;

	for (unsigned ievt(0); ievt<raretree->GetEntries(); ++ievt){//loop on entries
		raretree->GetEntry(ievt);
		
		double weight = 1;
		if(channelType == 1){
                        if(Year == 2016 and ISpreVFP == 1)             weight = lumi_2016preVFP_DoubleEG*1000*rarecrosssection/rarentotalevent;
                        else if(Year == 2016 and ISpreVFP == 0)        weight = lumi_2016postVFP_DoubleEG*1000*rarecrosssection/rarentotalevent;
                        else if(Year == 2017)                        weight = lumi_2017_DoubleEG*1000*rarecrosssection/rarentotalevent;
                        else if(Year == 2018)                        weight = lumi_2018_DoubleEG*1000*rarecrosssection/rarentotalevent;}

                else if(channelType == 2){
                        if(Year == 2016 and ISpreVFP == 1)             weight = lumi_2016preVFP_MuonEG*1000*rarecrosssection/rarentotalevent;
                        else if(Year == 2016 and ISpreVFP == 0)        weight = lumi_2016postVFP_MuonEG*1000*rarecrosssection/rarentotalevent;
                        else if(Year == 2017)                        weight = lumi_2017_MuonEG*1000*rarecrosssection/rarentotalevent;
                        else if(Year == 2018)                        weight = lumi_2018_MuonEG*1000*rarecrosssection/rarentotalevent;}
		
		if(ievt==0) cout<<"tt cross-section "<<rarecrosssection<<endl<<endl;
//		if(raresigMT < 100)continue;
//              if(raresigMET < 120)continue;
//		if(raresigMET > 70)continue;
		/** cut flow *****/
		if ((ichannel == 1 && (rarephoEt < 35 || rarelepPt < 25)) || (ichannel == 2 && (rarephoEt < 35 || rarelepPt < 20))) continue;
		if(fabs(rarephoEta) > 1.4442 || fabs(rarelepEta) > 2.5)continue;
		
		if (Year == 2017 && (rarephoEta > 1.0 && rarephoEta < 1.5 && rarephoPhi > 2.7)) continue;
                if (Year == 2017 && (rarelepEta > 1.0 && rarelepEta < 1.5 && rarelepPhi > 2.7)) continue;

                if (Year == 2018 && (rarephoEta > 0.3 && rarephoEta < 1.2 && rarephoPhi > 0.4 && rarephoPhi < 0.8)) continue;
                if (Year == 2018 && (rarelepEta > 0.3 && rarelepEta < 1.2 && rarelepPhi > 0.4 && rarelepPhi < 0.8)) continue;

		double w_ele = USE_VTX_BY_PT_SLICE ? EvalFakeRate_VtxByPtSlice(rarephoEt, rarenVertex, fabs(rarephoEta))
		                                    : EvalFakeRate(rarephoEt, rarenVertex, fabs(rarephoEta));
		// ele fake photon weight * XSec weight
		w_ele = w_ele*weight;

		// ---- Classify lepton in rare proxy event ----
		/*bool isRealLep_rare = false;
		double lepMindR_rare = 999.0;
		unsigned lepMatchIdx_rare = 0;
		for(unsigned iMC(0); iMC < rare_mcPID->size(); iMC++){
			int absPID = fabs((*rare_mcPID)[iMC]);
			if(absPID != 11 && absPID != 13)continue;
			double dR = DeltaR((*rare_mcEta)[iMC], (*rare_mcPhi)[iMC], rarelepEta, rarelepPhi);
			if(dR < lepMindR_rare){ lepMindR_rare = dR; lepMatchIdx_rare = iMC; }
		}
		if(lepMindR_rare < 0.1){
                        if(isHardLepton((*rare_mcMomPID)[lepMatchIdx_rare])) isRealLep_rare = true;
                }

		bool isFakeLep_rare = !isRealLep_rare;

		rare_wele_total   += w_ele;
		rare_n_total      += 1.0;
		if(isFakeLep_rare){
			rare_wele_fakelep += w_ele;
			rare_n_fakelep    += 1.0;
		}*/

		// predicted electron fake background using proxy events
		pred_PhoEt->Fill(rarephoEt,w_ele);
		pred_PhoEta->Fill(rarephoEta, w_ele);
		pred_MET->Fill(raresigMET, w_ele);
		pred_Mt->Fill(raresigMT, w_ele);
		pred_HT->Fill(rareHT, w_ele);
		pred_LepPt->Fill(rarelepPt, w_ele);
		pred_LepEta->Fill(rarelepEta, w_ele);
		pred_dPhiEleMET->Fill(fabs(raredPhiLepMET), w_ele);
	//	pred_nJet->Fill(rarenJet, w_ele);
		for(unsigned it(0); it < NTOY; it++){
			double toy_ele = ToyFakeRate(it, rarephoEt, rarenVertex, fabs(rarephoEta));
			toy_ele = toy_ele*weight;
			toy_PhoEt[it]->Fill(rarephoEt,toy_ele);
			toy_MET[it]->Fill(raresigMET, toy_ele);
			toy_Mt[it]->Fill(raresigMT, toy_ele);
			toy_HT[it]->Fill(rareHT, toy_ele);
			toy_LepPt[it]->Fill(rarelepPt, toy_ele);
			toy_dPhiEleMET[it]->Fill(fabs(raredPhiLepMET), toy_ele);
		}
	}

	// ---- Print fake-lepton fraction for the rare proxy (TT/WW/WZ) sample ----
	/*double fakeLepFrac_rare_wele  = (rare_wele_total  > 0) ? rare_wele_fakelep / rare_wele_total  : 0.0;
	double fakeLepFrac_rare_unwtd = (rare_n_total      > 0) ? rare_n_fakelep    / rare_n_total     : 0.0;
	std::cout << "============================================================" << std::endl;
	std::cout << " Rare proxy sample TT/WW/WZ (ele->photon): fake-lepton fraction" << std::endl;
	std::cout << "  Weighted  : " << rare_wele_fakelep << " / " << rare_wele_total
	          << "  = " << fakeLepFrac_rare_wele*100.0 << " %" << std::endl;
	std::cout << "  Unweighted: " << rare_n_fakelep    << " / " << rare_n_total
	          << "  = " << fakeLepFrac_rare_unwtd*100.0 << " %" << std::endl;
	std::cout << "============================================================" << std::endl;

	// ---- Combined fake-lepton fraction across all proxy samples ----
	double combined_wele_total   = proxy_wele_total   + rare_wele_total;
	double combined_wele_fakelep = proxy_wele_fakelep + rare_wele_fakelep;
	double combined_n_total      = proxy_n_total      + rare_n_total;
	double combined_n_fakelep    = proxy_n_fakelep    + rare_n_fakelep;
	double fakeLepFrac_combined_wele  = (combined_wele_total > 0) ? combined_wele_fakelep / combined_wele_total  : 0.0;
	double fakeLepFrac_combined_unwtd = (combined_n_total    > 0) ? combined_n_fakelep    / combined_n_total     : 0.0;
	std::cout << "============================================================" << std::endl;
	std::cout << " Combined proxy sample (ele->photon): fake-lepton fraction" << std::endl;
	std::cout << "  Weighted  : " << combined_wele_fakelep << " / " << combined_wele_total
	          << "  = " << fakeLepFrac_combined_wele*100.0 << " %" << std::endl;
	std::cout << "  Unweighted: " << combined_n_fakelep    << " / " << combined_n_total
	          << "  = " << fakeLepFrac_combined_unwtd*100.0 << " %" << std::endl;
	std::cout << "============================================================" << std::endl;*/

// total errors in bins of PhoEt, LepPt, MET, Mt, dPhiEleMET, HT
	std::vector<double> toyvec; 
	for(int ibin(1); ibin < pred_PhoEt->GetSize(); ibin++){
    toyvec.clear();
    toyvec.push_back(pred_PhoEt->GetBinContent(ibin));
    for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_PhoEt[it]->GetBinContent(ibin));
	// systematic error from predicted bkg and toys with Gauss fit
    		double syserr = calcToyError( toyvec, useGaussFit, channelType);
		// Total stat and syst error
		double totalerror = sqrt(syserr*syserr + pred_PhoEt->GetBinError(ibin)*pred_PhoEt->GetBinError(ibin));
		pred_PhoEt->SetBinError(ibin, totalerror);
		//pred_PhoEt->SetBinError(ibin, pred_PhoEt->GetBinError(ibin));
	}
	for(int ibin(1); ibin < pred_LepPt->GetSize(); ibin++){
    toyvec.clear();
    toyvec.push_back(pred_LepPt->GetBinContent(ibin));
    for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_LepPt[it]->GetBinContent(ibin));
    double syserr = calcToyError( toyvec, useGaussFit, channelType);

		double totalerror = sqrt(syserr*syserr + pred_LepPt->GetBinError(ibin)*pred_LepPt->GetBinError(ibin));
		pred_LepPt->SetBinError(ibin, totalerror);
		//pred_LepPt->SetBinError(ibin, pred_LepPt->GetBinError(ibin));

	}
	for(int ibin(1); ibin < pred_MET->GetSize(); ibin++){
    toyvec.clear();
    toyvec.push_back(pred_MET->GetBinContent(ibin));
    for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_MET[it]->GetBinContent(ibin));
    double syserr = calcToyError( toyvec, useGaussFit, channelType);

		double totalerror = sqrt(syserr*syserr + pred_MET->GetBinError(ibin)*pred_MET->GetBinError(ibin));
		pred_MET->SetBinError(ibin, totalerror);
		//pred_MET->SetBinError(ibin, pred_MET->GetBinError(ibin));
	}
	for(int ibin(1); ibin < pred_Mt->GetSize(); ibin++){
    toyvec.clear();
    toyvec.push_back(pred_Mt->GetBinContent(ibin));
    for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_Mt[it]->GetBinContent(ibin));
    double syserr = calcToyError( toyvec, useGaussFit, channelType);

		double totalerror = sqrt(syserr*syserr + pred_Mt->GetBinError(ibin)*pred_Mt->GetBinError(ibin));
		pred_Mt->SetBinError(ibin, totalerror);
		//pred_Mt->SetBinError(ibin, pred_Mt->GetBinError(ibin));
	}
	for(int ibin(1); ibin < pred_dPhiEleMET->GetSize(); ibin++){
    toyvec.clear();
    toyvec.push_back(pred_dPhiEleMET->GetBinContent(ibin));
    for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_dPhiEleMET[it]->GetBinContent(ibin));
    double syserr = calcToyError( toyvec, useGaussFit, channelType);

		double totalerror = sqrt(syserr*syserr + pred_dPhiEleMET->GetBinError(ibin)*pred_dPhiEleMET->GetBinError(ibin));
		pred_dPhiEleMET->SetBinError(ibin, totalerror);
		//pred_dPhiEleMET->SetBinError(ibin, pred_dPhiEleMET->GetBinError(ibin));
	}
	for(int ibin(1); ibin < pred_HT->GetSize(); ibin++){
    toyvec.clear();
    toyvec.push_back(pred_HT->GetBinContent(ibin));
    for(unsigned it(0); it < NTOY; it++)toyvec.push_back(toy_HT[it]->GetBinContent(ibin));
    double syserr = calcToyError( toyvec, useGaussFit, channelType);

		double totalerror = sqrt(syserr*syserr + pred_HT->GetBinError(ibin)*pred_HT->GetBinError(ibin));
		pred_HT->SetBinError(ibin, totalerror);
		//pred_HT->SetBinError(ibin, pred_HT->GetBinError(ibin));
	}
		
	TGraphErrors *error_PhoEt = new TGraphErrors(nBkgEtBins);
	TGraphErrors *error_MET = new TGraphErrors(nBkgMETBins); 
	TGraphErrors *error_Mt = new TGraphErrors(nBkgMtBins); 
	TGraphErrors *error_HT = new TGraphErrors(nBkgHTBins);
	TGraphErrors *ratioerror_PhoEt = new TGraphErrors(nBkgEtBins);
	TGraphErrors *ratioerror_MET = new TGraphErrors(nBkgMETBins); 
	TGraphErrors *ratioerror_Mt = new TGraphErrors(nBkgMtBins); 
	TGraphErrors *ratioerror_HT = new TGraphErrors(nBkgHTBins);
// from signal events
	TH1F *ratio=(TH1F*)p_PhoEt->Clone("transfer factor");
	TH1F *ratio_met=(TH1F*)p_MET->Clone("transfer factor");
	TH1F *ratio_mt=(TH1F*)p_Mt->Clone("transfer factor");
	TH1F *ratio_HT=(TH1F*)p_HT->Clone("transfer factor");

	TCanvas *c_pt = new TCanvas("Photon_Pt", "Photon P_{T}",600,600);  
	setCanvas(c_pt); 
	c_pt->cd();                                                             
	TPad *pt_pad1 = new TPad("pt_pad1", "pt_pad1", 0, 0.35, 1, 1.0);        
	setTopPad(pt_pad1); 
	pt_pad1->SetBottomMargin(0);
	pt_pad1->Draw();                                                        
	pt_pad1->cd();                                                          
	gPad->SetLogy();
	p_PhoEt->SetMaximum(200*p_PhoEt->GetBinContent(p_PhoEt->GetMaximumBin())); 
	p_PhoEt->GetXaxis()->SetRangeUser(35,200);                              
	p_PhoEt->GetXaxis()->SetLabelSize(0);
	p_PhoEt->Draw();                                                        
	p_PhoEt->SetLineColor(kBlack);                                          
	p_PhoEt->SetMarkerStyle(20);                                            
	pred_PhoEt->SetFillStyle(1001);                                         
	pred_PhoEt->SetLineColor(kRed);                                         
	pred_PhoEt->SetFillColor(kRed);                                         
	DY_PhoEt->SetFillStyle(1001);                                           
	DY_PhoEt->SetLineColor(kYellow);
	DY_PhoEt->SetFillColor(kYellow);
// from proxy events, predicted backgrounds
	pred_PhoEt->Draw("hist same");
	DY_PhoEt->Draw("hist same");
	TLegend *leg =  new TLegend(0.5,0.55,0.9,0.8);
	leg->SetFillStyle(0);
	gStyle->SetLegendBorderSize(1);
	gStyle->SetLegendFillColor(0);
	pred_PhoEt->SetMarkerColor(kRed);
	DY_PhoEt->SetMarkerColor(kYellow);
	error_PhoEt->SetMarkerSize(0);
	error_PhoEt->SetLineWidth(0);
	leg->AddEntry(p_PhoEt,"Direct simulation", "pel");
//	leg->AddEntry(pred_PhoEt,"t#bar{t}");
	leg->AddEntry(pred_PhoEt,"t#bar{t} / WW / WZ");
	leg->AddEntry(DY_PhoEt,"Drell#scale[0.5]{ }#minusYan");

	leg->AddEntry(error_PhoEt, "Total uncertainty");
	leg->Draw("same");
  	TLatex chantex;
  	chantex.SetNDC();
  	chantex.SetTextFont(42);
  	chantex.SetTextSize(0.07);    
  	if(channelType==1) chantex.DrawLatex(0.58,0.82," e + #gamma");
  	if(channelType==2) chantex.DrawLatex(0.58,0.82," #mu + #gamma");
 	gPad->RedrawAxis();
	if(Year==2016 and ISpreVFP == 1)       CMS_lumi( pt_pad1, 1,ichannel, 11 );
        else if(Year==2016 and ISpreVFP == 0)  CMS_lumi( pt_pad1, 2,ichannel, 11 );
        else if(Year==2017)                  CMS_lumi( pt_pad1, 3,ichannel, 11 );
        else if(Year==2018)                  CMS_lumi( pt_pad1, 4,ichannel, 11 );
		
	p_PhoEt->Draw("same");
	for(int ibin(1); ibin < pred_PhoEt->GetSize(); ibin++){
		error_PhoEt->SetPoint(ibin-1,pred_PhoEt->GetBinCenter(ibin), pred_PhoEt->GetBinContent(ibin));
		// syst + stat error
		float prederror = pred_PhoEt->GetBinError(ibin);
		error_PhoEt->SetPointError(ibin-1,(pred_PhoEt->GetBinLowEdge(ibin+1)-pred_PhoEt->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_PhoEt->SetPoint(ibin-1,pred_PhoEt->GetBinCenter(ibin), 1); 
		ratioerror_PhoEt->SetPointError(ibin-1,(pred_PhoEt->GetBinLowEdge(ibin+1)-pred_PhoEt->GetBinLowEdge(ibin))/2, prederror/pred_PhoEt->GetBinContent(ibin)); 
	}
	error_PhoEt->SetFillColor(12);
  	error_PhoEt->SetFillStyle(3345);
	error_PhoEt->Draw("E2 same");

	c_pt->cd();
	TPad *pt_pad2 = new TPad("pt_pad2", "pt_pad2", 0, 0, 1, 0.35);
	pt_pad2->SetTopMargin(0);
	pt_pad2->SetBottomMargin(0.3);
	pt_pad2->Draw();
	pt_pad2->cd();
  TLine *flatratio = new TLine(35,1,200,1);
	ratio->GetXaxis()->SetRangeUser(35,200);
	ratio->GetYaxis()->SetRangeUser(0.4,1.7);
	ratio->GetYaxis()->SetNdivisions(504);
	ratio->SetMarkerStyle(20);
	ratio->SetLineColor(kBlack);
	// from signal events - predicted events
	ratio->Divide(pred_PhoEt);
	
	ratio->SetTitle("");
	ratio->GetXaxis()->SetTitleOffset(0.9);
	ratio->GetXaxis()->SetTitle("p_{T}^{#gamma} (GeV)");
	ratio->GetYaxis()->SetTitle("#frac{Simulation}{Prediction} ");
	ratio->GetYaxis()->SetNdivisions(504);
	ratio->Draw();
	ratioerror_PhoEt->SetFillColor(12);
	ratioerror_PhoEt->SetFillStyle(3345);
	ratioerror_PhoEt->Draw("E2 same");
	ratio->Draw("same");
	flatratio->Draw("same");
	TLine *ratioValue_PhoEt = new TLine(35,p_PhoEt->Integral()/pred_PhoEt->Integral(),200,p_PhoEt->Integral()/pred_PhoEt->Integral());
        ratioValue_PhoEt->SetLineColor(kRed);
        ratioValue_PhoEt->Draw("same");
	
	if(channelType==1) c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Closure/closure_elefakepho_PhotonEt_eg_%d%s_pT_nVtx.png",Year,whichVFP.c_str()));
	if(channelType==2) c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Closure/closure_elefakepho_PhotonEt_mg_%d%s_pT_nVtx.png",Year,whichVFP.c_str()));


// similar for other variables
// ******** MET ************************//
	TCanvas *c_met = new TCanvas("MET", "MET",600,600);
	c_met->cd();
	TPad *met_pad1 = new TPad("met_pad1", "met_pad1", 0, 0.35, 1, 1.0);
	setTopPad(met_pad1); 
	met_pad1->SetBottomMargin(0);
	met_pad1->Draw();  
	met_pad1->cd();  
	gPad->SetLogy();
	cout<<"Direct simulation "<< p_MET->Integral()<<endl;
	cout<<"DY " <<DY_MET->Integral()<<endl;
	cout<<"DY + t#bar{t} + WW + WZ "<<pred_MET->Integral()<<endl;
	cout<<"Ratio "<<p_MET->Integral()/pred_MET->Integral()<<endl;
	p_MET->SetMaximum(1000*p_MET->GetBinContent(p_MET->GetMaximumBin())); 
	p_MET->SetMinimum(0.05);
	p_MET->GetXaxis()->SetRangeUser(0,400);
	p_MET->GetXaxis()->SetLabelSize(0);
	p_MET->SetLineColor(1);
	p_MET->SetMarkerStyle(20);
	p_MET->Draw("P");
	DY_MET->SetFillStyle(1001);
	DY_MET->SetLineColor(kYellow-4);
	DY_MET->SetFillColor(kYellow-4);
	pred_MET->SetFillStyle(1001);
	pred_MET->SetLineColor(kRed);
	pred_MET->SetFillColor(kRed);
	for(int ibin(1); ibin < pred_MET->GetSize(); ibin++){
		float prederror = pred_MET->GetBinError(ibin);
		// Total stat and syst error
		error_MET->SetPoint(ibin-1,pred_MET->GetBinCenter(ibin), pred_MET->GetBinContent(ibin));
		error_MET->SetPointError(ibin-1,(pred_MET->GetBinLowEdge(ibin+1)-pred_MET->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_MET->SetPoint(ibin-1,pred_MET->GetBinCenter(ibin), 1); 
		ratioerror_MET->SetPointError(ibin-1,(pred_MET->GetBinLowEdge(ibin+1)-pred_MET->GetBinLowEdge(ibin))/2, prederror/pred_MET->GetBinContent(ibin)); 
	}
	pred_MET->Draw("hist same");
	DY_MET->Draw("hist same");
  	error_MET->SetFillColor(12);
  	error_MET->SetFillStyle(3345);
	error_MET->Draw("E2 same");
	leg->Draw("same");
	p_MET->Draw("E same");
  	if(channelType==1) chantex.DrawLatex(0.58,0.82," e + #gamma");
  	if(channelType==2) chantex.DrawLatex(0.58,0.82," #mu + #gamma");
	//TLine *line_met = new TLine(70,0,70,10000);
	TLine *line_met = new TLine(70,0.05,70,1000*p_MET->GetBinContent(p_MET->GetMaximumBin()));
	line_met->SetLineStyle(2);
	line_met->Draw("same");
	TLatex* latex = new TLatex();
	latex->SetTextSize(0.05);
	latex->DrawLatex(20, 60000,"control");
	latex->DrawLatex(20, 30000,"region");
 	gPad->RedrawAxis();
	if(Year==2016 and ISpreVFP == 1)       CMS_lumi( met_pad1, 1,ichannel, 11 );
        else if(Year==2016 and ISpreVFP == 0)  CMS_lumi( met_pad1, 2,ichannel, 11 );
        else if(Year==2017)                  CMS_lumi( met_pad1, 3,ichannel, 11 );
        else if(Year==2018)                  CMS_lumi( met_pad1, 4,ichannel, 11 );

	c_met->cd();
	TPad *met_pad2 = new TPad("met_pad2", "met_pad2", 0, 0, 1, 0.35);
	met_pad2->SetTopMargin(0);
	met_pad2->SetBottomMargin(0.3);
	met_pad2->Draw();
	met_pad2->cd();
  TLine *flatratio_met = new TLine(0,1,400,1);
	ratio_met->GetXaxis()->SetRangeUser(0,400);
	ratio_met->GetYaxis()->SetNdivisions(504);
	ratio_met->SetLineColor(kBlack);
	ratio_met->SetMarkerStyle(20);
	ratio_met->Divide(pred_MET);
	ratio_met->SetTitle("");
	ratio_met->GetXaxis()->SetTitleOffset(0.9);
	ratio_met->GetXaxis()->SetTitle("p_{T}^{miss} (GeV)");
	ratio_met->GetYaxis()->SetTitle("#frac{Simulation}{Prediction} ");
	ratio_met->GetYaxis()->SetRangeUser(0.4,1.7);
	ratio_met->GetYaxis()->SetNdivisions(504);
	ratio_met->Draw();
	ratioerror_MET->SetFillColor(12);
	ratioerror_MET->SetFillStyle(3345);
	ratioerror_MET->Draw("E2 same");
	ratio_met->Draw("same");
	flatratio_met->Draw("same");
	
	TLine *ratioValue_met = new TLine(0,p_MET->Integral()/pred_MET->Integral(),400,p_MET->Integral()/pred_MET->Integral());
        ratioValue_met->SetLineColor(kRed);
        ratioValue_met->Draw("same");

	TLine *line_met_ratio = new TLine(70,0.4,70,1.7);
	line_met_ratio->SetLineStyle(2);
	line_met_ratio->Draw("same");
 	gPad->RedrawAxis();
	if(channelType==1) c_met->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Closure/closure_elefakepho_MET_eg_%d%s_pT_nVtx.png",Year,whichVFP.c_str()));
	if(channelType==2) c_met->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Closure/closure_elefakepho_MET_mg_%d%s_pT_nVtx.png",Year,whichVFP.c_str()));

// ******** Mt ************************
	TCanvas *c_mt = new TCanvas("Mt", "Mt",600,600);
	c_mt->cd();
	TPad *mt_pad1 = new TPad("mt_pad1", "mt_pad1", 0, 0.35, 1, 1.0);
	setTopPad(mt_pad1); 
	mt_pad1->SetBottomMargin(0);
	mt_pad1->Draw();  
	mt_pad1->cd();  
	gPad->SetLogy();
//	p_Mt->SetMinimum(0.05);
	p_Mt->SetMaximum(100*p_Mt->GetBinContent(p_Mt->GetMaximumBin())); 
	p_Mt->GetXaxis()->SetRangeUser(0,400);
	p_Mt->GetXaxis()->SetLabelSize(0);
	p_Mt->SetLineColor(1);
	p_Mt->SetMarkerStyle(20);
	p_Mt->Draw("P");
	DY_Mt->SetFillStyle(1001);
	DY_Mt->SetLineColor(kYellow-4);
	DY_Mt->SetFillColor(kYellow-4);
	pred_Mt->SetFillStyle(1001);
	pred_Mt->SetLineColor(kRed);
	pred_Mt->SetFillColor(kRed);
	for(int ibin(1); ibin < pred_Mt->GetSize(); ibin++){
		error_Mt->SetPoint(ibin-1,pred_Mt->GetBinCenter(ibin), pred_Mt->GetBinContent(ibin));
		float prederror = pred_Mt->GetBinError(ibin);
		error_Mt->SetPointError(ibin-1,(pred_Mt->GetBinLowEdge(ibin+1)-pred_Mt->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_Mt->SetPoint(ibin-1,pred_Mt->GetBinCenter(ibin), 1); 
		ratioerror_Mt->SetPointError(ibin-1,(pred_Mt->GetBinLowEdge(ibin+1)-pred_Mt->GetBinLowEdge(ibin))/2, prederror/pred_Mt->GetBinContent(ibin)); 
	}
	pred_Mt->Draw("hist same");
	DY_Mt->Draw("hist same");
  error_Mt->SetFillColor(12);
  error_Mt->SetFillStyle(3345);
	error_Mt->Draw("E2 same");
	leg->Draw("same");
	p_Mt->Draw("E same");
  	if(channelType==1) chantex.DrawLatex(0.58,0.82," e + #gamma");
  	if(channelType==2) chantex.DrawLatex(0.58,0.82," #mu + #gamma");
 	gPad->RedrawAxis();
	if(Year==2016 and ISpreVFP == 1)       CMS_lumi( mt_pad1, 1,ichannel, 11 );
        else if(Year==2016 and ISpreVFP == 0)  CMS_lumi( mt_pad1, 2,ichannel, 11 );
        else if(Year==2017)                  CMS_lumi( mt_pad1, 3,ichannel, 11 );
        else if(Year==2018)                  CMS_lumi( mt_pad1, 4,ichannel, 11 );

	c_mt->cd();
	TPad *mt_pad2 = new TPad("mt_pad2", "mt_pad2", 0, 0, 1, 0.35);
	mt_pad2->SetTopMargin(0);
	mt_pad2->SetBottomMargin(0.3);
	mt_pad2->Draw();
	mt_pad2->cd();
  	TLine *flatratio_mt = new TLine(0,1,400,1);
	ratio_mt->SetMarkerStyle(20);
	ratio_mt->SetLineColor(kBlack);
	ratio_mt->GetXaxis()->SetRangeUser(0,400);
	ratio_mt->GetYaxis()->SetRangeUser(0.4,1.7);
	ratio_mt->GetYaxis()->SetNdivisions(504);
	ratio_mt->Divide(pred_Mt);
	ratio_mt->SetTitle("");
	ratio_mt->GetXaxis()->SetTitleOffset(0.9);
	ratio_mt->GetXaxis()->SetTitle("M_{T} (GeV)");
	ratio_mt->GetYaxis()->SetTitle("#frac{Simulation}{Prediction} ");
	ratio_mt->Draw();
	ratioerror_Mt->SetFillColor(12);
	ratioerror_Mt->SetFillStyle(3345);
	ratioerror_Mt->Draw("E2 same");
	ratio_mt->Draw("same");
	flatratio_mt->Draw("same");

	TLine *ratioValue_mt = new TLine(0,p_Mt->Integral()/pred_Mt->Integral(),400,p_Mt->Integral()/pred_Mt->Integral());
        ratioValue_mt->SetLineColor(kRed);
        ratioValue_mt->Draw("same");

	if(channelType==1) c_mt->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Closure/closure_elefakepho_MT_eg_%d%s_pT_nVtx.png",Year,whichVFP.c_str()));
	if(channelType==2) c_mt->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Closure/closure_elefakepho_MT_mg_%d%s_pT_nVtx.png",Year,whichVFP.c_str()));

// ******** HT ************************ 
	TCanvas *c_HT = new TCanvas("HT", "HT",600,600);
	c_HT->cd();
	TPad *HT_pad1 = new TPad("HT_pad1", "HT_pad1", 0, 0.35, 1, 1.0);
	setTopPad(HT_pad1); 
	HT_pad1->SetBottomMargin(0);
	HT_pad1->Draw();  
	HT_pad1->cd();  
	gPad->SetLogy();
	p_HT->SetMaximum(100*p_HT->GetBinContent(p_HT->GetMaximumBin())); 
	p_HT->GetXaxis()->SetRangeUser(0,400);
	p_HT->GetXaxis()->SetLabelSize(0);

//	p_HT->SetMinimum(5);
	p_HT->SetLineColor(1);
	p_HT->SetMarkerStyle(20);
	p_HT->Draw("P");
	DY_HT->SetFillStyle(1001);
	DY_HT->SetLineColor(kYellow-4);
	DY_HT->SetFillColor(kYellow-4);
	pred_HT->SetFillStyle(1001);
	pred_HT->SetLineColor(kRed);
	pred_HT->SetFillColor(kRed);
	for(int ibin(1); ibin < pred_HT->GetSize(); ibin++){
		error_HT->SetPoint(ibin-1,pred_HT->GetBinCenter(ibin), pred_HT->GetBinContent(ibin));
		float prederror = pred_HT->GetBinError(ibin);
		error_HT->SetPointError(ibin-1,(pred_HT->GetBinLowEdge(ibin+1)-pred_HT->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_HT->SetPoint(ibin-1,pred_HT->GetBinCenter(ibin), 1); 
		ratioerror_HT->SetPointError(ibin-1,(pred_HT->GetBinLowEdge(ibin+1)-pred_HT->GetBinLowEdge(ibin))/2, prederror/pred_HT->GetBinContent(ibin)); 
	}
	pred_HT->Draw("hist same");
	DY_HT->Draw("hist same");
	error_HT->SetFillColor(12);
  	error_HT->SetFillStyle(3345);
	error_HT->Draw("E2 same");
	leg->Draw("same");
	p_HT->Draw("E same");
	  if(channelType==1) chantex.DrawLatex(0.58,0.82," e + #gamma");
	  if(channelType==2) chantex.DrawLatex(0.58,0.82," #mu + #gamma");
 	gPad->RedrawAxis();
	if(Year==2016 and ISpreVFP == 1)       CMS_lumi( HT_pad1, 1, ichannel, 11 );
        else if(Year==2016 and ISpreVFP == 0)  CMS_lumi( HT_pad1, 2, ichannel, 11 );
        else if(Year==2017)                  CMS_lumi( HT_pad1, 3, ichannel, 11 );
        else if(Year==2018)                  CMS_lumi( HT_pad1, 4, ichannel, 11 );

	c_HT->cd();
	TPad *HT_pad2 = new TPad("HT_pad2", "HT_pad2", 0, 0, 1, 0.35);
	HT_pad2->SetTopMargin(0);
	HT_pad2->SetBottomMargin(0.3);
	HT_pad2->Draw();
	HT_pad2->cd();
  	TLine *flatratio_HT = new TLine(0,1,420,1);
	ratio_HT->SetMarkerStyle(20);
	ratio_HT->SetLineColor(kBlack);
	ratio_HT->GetXaxis()->SetRangeUser(0,420);
	ratio_HT->GetYaxis()->SetRangeUser(0.4,1.7);
	ratio_HT->GetYaxis()->SetNdivisions(504);
	ratio_HT->Divide(pred_HT);
	ratio_HT->SetTitle("");
	ratio_HT->GetXaxis()->SetTitleOffset(0.9);
	ratio_HT->GetXaxis()->SetTitle("H_{T} (GeV)");
	ratio_HT->GetYaxis()->SetTitle("#frac{Simulation}{Prediction} ");
	ratio_HT->Draw();
	ratioerror_HT->SetFillColor(12);
	ratioerror_HT->SetFillStyle(3345);
	ratioerror_HT->Draw("E2 same");
	flatratio_HT->Draw("same");
	TLine *ratioValue_ht = new TLine(0,p_HT->Integral()/pred_HT->Integral(),420,p_HT->Integral()/pred_HT->Integral());
        ratioValue_ht->SetLineColor(kRed);
        ratioValue_ht->Draw("same");
	if(channelType==1) c_HT->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Closure/closure_elefakepho_HT_eg_%d%s_pT_nVtx.png",Year,whichVFP.c_str()));
	if(channelType==2) c_HT->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Closure/closure_elefakepho_HT_mg_%d%s_pT_nVtx.png",Year,whichVFP.c_str()));

}
