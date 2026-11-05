// g++ `root-config --cflags` ../../lib/libAnaClasses.so plot_ISRweight.C -o plot_ISRweight.exe `root-config --libs`
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<ctime>
#include<vector>
#include<cstdlib>
#include<iomanip>

#include "TROOT.h"
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
#include "TLatex.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TGraphErrors.h"

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

// ../../include/analysis_scalefactor.h
// change the file /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/SUSYAnalysis/test/Background/BkgPredConfig.txt
//bool doEB=false;

namespace {
struct ISRBinDef {
	double low;
	double high;
};

static std::vector<ISRBinDef> GetISRBinDefs(){
	return { {0,50}, {50,100}, {100,150}, {150,200}, {200,250}, {250,300}, {300,1.0e9} };
}

static std::string GetEraTag(int RunYear, bool preVFP){
	std::ostringstream os;
	os << RunYear;
	if(RunYear == 2016) os << (preVFP ? "preVFP" : "postVFP");
	return os.str();
}

static std::string GetISRWeightFileName(int RunYear, bool preVFP, string whichChannel){
	std::ostringstream os;
	os << "/eos/uscms/store/user/tmishra/ISRweighting/ISRweights_" << GetEraTag(RunYear, preVFP) << "_" << whichChannel << ".txt";
	return os.str();
}

static bool WriteISRWeights(const std::string& path, const std::vector<double>& weights, double normalization){
	std::ofstream out(path.c_str());
	if(!out.is_open()) return false;
	auto bins = GetISRBinDefs();
	out << "# normalization " << normalization << "\n";
	out << "# ptLow ptHigh weight\n";
	for(size_t i=0; i<bins.size() && i<weights.size(); ++i){
		out << bins[i].low << " " << bins[i].high << " " << std::setprecision(8) << weights[i] << "\n";
	}
	out.close();
	return true;
}

static std::vector<double> ReadISRWeights(const std::string& path, double &normalization){
	std::ifstream in(path.c_str());
	std::vector<double> weights(GetISRBinDefs().size(),1.0);
	normalization = 1.0;
	if(!in.is_open()) return weights;
	std::string line;
	int ib(0);
	while(std::getline(in,line)){
		if(line.empty()) continue;
		if(line.rfind("# normalization",0)==0){
			std::istringstream hs(line);
			std::string hash, key;
			double nval(1.0);
			if(hs >> hash >> key >> nval) normalization = nval;
			continue;
		}
		if(line[0]=='#') continue;
		std::istringstream ss(line);
		double low(0), high(0), w(1.0);
		if(!(ss>>low>>high>>w)) continue;
		if(ib < (int)weights.size()) weights[ib] = w;
		++ib;
	}
	in.close();
	return weights;
}

struct FakeRateBin {
	double low;
	double high;
	double rate;
};

static std::string GetJetFakeRateFileName(const std::string& det, const std::string& fakeChannel, int RunYear, bool preVFP){
	std::ostringstream os;
	os << "/eos/uscms/store/user/tmishra/jetfakepho/files/JetFakeRate-ISR-" << det << "_" << fakeChannel << "_" << GetEraTag(RunYear, preVFP) << ".txt";
	return os.str();
}

// Fixed choices of the 2016 analysis (mengleis .../test/Reweight/plot_ISRweight.C),
// where the fits run out of statistics:
const double kFakeNoCorrPt   = 250.;  // photon pT >= 250: no fake subtraction (weight 1)
const double kFakeEEFixedLow = 150.;  // EE, 150 <= pT < 250: fixed fake rate
const double kFakeEEFixed    = 0.03;
const double kFakeFallback   = 0.03;  // failed fit (rate NaN, <= 0, >= 1, or non-finite error)

// One bin per line: ptLow ptHigh fakerate [total_error syst_error mc_fakerate]
static std::vector<FakeRateBin> ReadJetFakeRate(const std::string& path){
	std::vector<FakeRateBin> bins;
	std::ifstream in(path.c_str());
	if(!in.is_open()){
		std::cout << "[WARN] missing jet fake rate file " << path << ", no fake subtraction" << std::endl;
		return bins;
	}
	std::string line;
	while(std::getline(in,line)){
		if(line.empty() || line[0]=='#') continue;
		std::istringstream ss(line);
		FakeRateBin fb;
		std::string errStr;
		if(!(ss >> fb.low >> fb.high >> fb.rate)) continue;
		ss >> errStr;  // total error; "inf"/"nan" marks a failed fit
		bool badErr = (!errStr.empty() && (errStr.find("inf") != std::string::npos || errStr.find("nan") != std::string::npos));
		if(std::isnan(fb.rate) || fb.rate <= 0 || fb.rate >= 1.0 || badErr){
			std::cout << "[WARN] failed fake-rate fit in " << path << " bin [" << fb.low << "," << fb.high
			          << "): rate " << fb.rate << " err " << errStr << " -> " << kFakeFallback << std::endl;
			fb.rate = kFakeFallback;
		}
		bins.push_back(fb);
	}
	in.close();
	std::cout << "[INFO] read " << bins.size() << " jet fake rate bins from " << path << std::endl;
	return bins;
}

// Same object/trigger SFs as Background/analysis_VGBkg.C: channel 1 = egamma, 2 = mgamma
static double GetChannelSF(esfScaleFactor& esf, int channel, float lepPt, float lepEta, float phoEt, float phoEta){
	if(channel == 1)
		return esf.getElectronESF(lepPt,lepEta)*esf.getPhotonESF(phoEt,phoEta)*esf.getegPhotonTRGESF(phoEt,phoEta)*esf.getElectronTRGESF(lepPt,lepEta);
	return esf.getMuonESF(lepPt,lepEta)*esf.getPhotonESF(phoEt,phoEta)*esf.getMuonEGTRGESF(phoEt, lepPt);
}

static double ISRWeightFromPt(double pt, const std::vector<double>& weights, double normalization){
	auto bins = GetISRBinDefs();
	for(size_t i=0; i<bins.size() && i<weights.size(); ++i){
		if(pt >= bins[i].low && pt < bins[i].high) return weights[i]*normalization;
	}
	return 1.0;
}
}

static double SafeRatioISR(double num, double den){
	return den>0 ? num/den : 1.0;
}


void plot_ISRweight(int RunYear=2016, bool preVFP=true, bool applyISR=false, int channel=1){//main

	setTDRStyle();
	TH1::AddDirectory(kFALSE);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);
	gStyle->SetErrorX(0.5);

  	gSystem->Load("../../lib/libAnaClasses.so");
	
	std::string whichVFP;
  	if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
  	if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
  	if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

	std::string whichChannel;
        if(channel == 1) whichChannel = "egamma";
        if(channel == 2) whichChannel = "mgamma";

	gROOT->SetBatch(kTRUE);
	esfScaleFactor  objectESF;
	const std::string isrWgtFile = GetISRWeightFileName(RunYear, preVFP, whichChannel);
	std::vector<double> loadedISRWeights;
	double loadedNormalization(1.0);
	if(applyISR) loadedISRWeights = ReadISRWeights(isrWgtFile, loadedNormalization);

	Double_t plotEtBins[]={0,50,100,150,200,250,300,800};
	TH1F *p_phoEt_data     = new TH1F("p_phoEt_data","",7,plotEtBins);
	TH1F *p_JetPt_data     = new TH1F("p_JetPt_data","",7,plotEtBins);
	TH1F *p_JetPt_total    = new TH1F("p_JetPt_total","",7,plotEtBins);
	TH1F *p_llmass_data    = new TH1F("p_llmass_data","",100,30,130);
	TH1F *p_phoEt_data_highEt   = new TH1F("p_phoEt_data_highEt","",7,plotEtBins);
	TH1F *p_JetPt_data_highEt   = new TH1F("p_JetPt_data_highEt","",7,plotEtBins);
	TH1F *p_llmass_data_highEt  = new TH1F("p_llmass_data_highEt","",100,30,130);

	TH1F *p_phoEt_ZG     = new TH1F("p_phoEt_ZG","",7,plotEtBins);
	TH1F *p_JetPt_ZG     = new TH1F("p_JetPt_ZG","",7,plotEtBins);
	TH1F *p_llmass_ZG    = new TH1F("p_llmass_ZG","",100,30,130);
	TH1F *p_phoEt_ZG_highEt   = new TH1F("p_phoEt_ZG_highEt","",7,plotEtBins);
	TH1F *p_JetPt_ZG_highEt   = new TH1F("p_JetPt_ZG_highEt","",7,plotEtBins);
	TH1F *p_llmass_ZG_highEt  = new TH1F("p_llmass_ZG_highEt","",100,30,130);

	TH1F *p_phoEt_rare     = new TH1F("p_phoEt_rare","",7,plotEtBins);
	TH1F *p_JetPt_rare     = new TH1F("p_JetPt_rare","",7,plotEtBins);
	TH1F *p_llmass_rare    = new TH1F("p_llmass_rare","",100,30,130);
	TH1F *p_phoEt_rare_highEt   = new TH1F("p_phoEt_rare_highEt","",7,plotEtBins);
	TH1F *p_JetPt_rare_highEt   = new TH1F("p_JetPt_rare_highEt","",7,plotEtBins);
	TH1F *p_llmass_rare_highEt  = new TH1F("p_llmass_rare_highEt","",100,30,130);

	TProfile *p_scalefactor  = new TProfile("p_scalefactor","p_scalefactor",50,-2.5,2.5);

	//************ Signal Tree **********************//
	TChain *tree = new TChain("ZTree");
	// using analysis_ISR.C
  	tree->Add(Form("/eos/uscms/store/user/tmishra/ISRweighting/%s/resTree_ISR_data_%d%s.root", whichChannel.c_str(), RunYear,whichVFP.c_str()));
  	float phoEt(0);
  	float phoEta(0);
  	float phoPhi(0);
  	float lepPt(0);
  	float lepEta(0);
  	float lepPhi(0);
  	float sigMT(0);
  	float sigMET(0);
  	float dPhiLepMET(0);
  	float dRPhoLep(0);
  	float HT(0);
  	float nJet(0);
	float threeMass(0);
	float dilepMass(0);
	float JetPt=0;

  	tree->SetBranchAddress("phoEt",     &phoEt);
  	tree->SetBranchAddress("phoEta",    &phoEta);
  	tree->SetBranchAddress("phoPhi",    &phoPhi);
  	tree->SetBranchAddress("lepPt",     &lepPt);
  	tree->SetBranchAddress("lepEta",    &lepEta);
  	tree->SetBranchAddress("lepPhi",    &lepPhi);
  	tree->SetBranchAddress("sigMT",     &sigMT);
  	tree->SetBranchAddress("sigMET",    &sigMET);
  	tree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
  	tree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
  	tree->SetBranchAddress("HT",        &HT);
  	tree->SetBranchAddress("nJet",      &nJet);
	tree->SetBranchAddress("threeMass", &threeMass);
	tree->SetBranchAddress("dilepMass", &dilepMass);
	tree->SetBranchAddress("ISRJetPt",  &JetPt);

	// jet->photon fake rates from runFit_step2.sh (FitJetFake.C), per channel and era
	const std::string fakeChannel = (channel == 1) ? "egamma" : "mg";
	const std::vector<FakeRateBin> fakeBinsEB = ReadJetFakeRate(GetJetFakeRateFileName("EB", fakeChannel, RunYear, preVFP));
	const std::vector<FakeRateBin> fakeBinsEE = ReadJetFakeRate(GetJetFakeRateFileName("EE", fakeChannel, RunYear, preVFP));

  	for(unsigned ievt(0); ievt<tree->GetEntries(); ++ievt){//loop on entries
		tree->GetEntry(ievt);
		// always set 0 as doEB=0, so both EB and EE considered
		//if(doEB && fabs(phoEta) > 1.4442)continue;
		
		//else if(!doEB && (fabs(phoEta) < 1.56 || fabs(phoEta) > 2.4))continue;

		if(phoEt > 799)phoEt = 799;
		if(JetPt > 799)JetPt = 799;
		if(dRPhoLep < 0.8)continue;
		if(dilepMass < 80 || dilepMass > 100)continue;
		double weight = 1;

		// subtract jet->photon fakes: weight = 1 - fakerate(phoEt, EB/EE), with the 2016 choices
		const bool isEB = (fabs(phoEta) < 1.4442);
		const std::vector<FakeRateBin>& fakeBins = isEB ? fakeBinsEB : fakeBinsEE;
		if(phoEt >= kFakeNoCorrPt) weight = 1;
		else if(!isEB && phoEt >= kFakeEEFixedLow) weight = 1 - kFakeEEFixed;
		else{
			for(const auto& fb : fakeBins){
				if(phoEt >= fb.low && phoEt < fb.high){ weight = 1 - fb.rate; break; }
			}
		}

	 	p_JetPt_data->Fill(JetPt, weight);
		p_phoEt_data->Fill(phoEt, weight); 
		p_llmass_data->Fill(dilepMass, weight);

		//if(phoEt > 145){
		//	p_llmass_data_highEt->Fill(dilepMass, weight);
		//	p_JetPt_data_highEt->Fill(JetPt, weight);
		//	p_phoEt_data_highEt->Fill(phoEt, weight);
		//}
	}//loop on  events


	double reweightF(1);
        double Normalization=1;
	float ISRWeight(0);

	//************ Signal Tree **********************//
  	TChain *ZGtree = new TChain("ZTree");
	// using analysis_ISRMC.C
	ZGtree->Add(Form("/eos/uscms/store/user/tmishra/ISRweighting/%s/resTree_ISR_ZGToLLG_%d%s.root", whichChannel.c_str(), RunYear,whichVFP.c_str()));
	float ZG_MCweight(0);
	float ZG_PUweight(0);
  	float ZG_phoEt(0);
  	float ZG_phoEta(0);
  	float ZG_phoPhi(0);
  	float ZG_lepPt(0);
  	float ZG_lepEta(0);
  	float ZG_lepPhi(0);
	float ZG_bosonPt(0);
  	float ZG_sigMT(0);
  	float ZG_sigMET(0);
  	float ZG_dPhiLepMET(0);
  	float ZG_dRPhoLep(0);
  	float ZG_HT(0);
  	float ZG_nJet(0);
	float ZG_JetPt(0);
	float ZG_threeMass(0);
	float ZG_dilepMass(0);
  	std::vector<int>   *ZG_mcPID=0;
  	std::vector<float> *ZG_mcEta=0;
  	std::vector<float> *ZG_mcPhi=0;
  	std::vector<float> *ZG_mcPt=0;
  	std::vector<int>   *ZG_mcMomPID=0;
  	std::vector<int>   *ZG_mcGMomPID=0;

	// MC weight is Lumi weight
	ZGtree->SetBranchAddress("MCweight",  &ZG_MCweight); 
	ZGtree->SetBranchAddress("PUweight",  &ZG_PUweight);
  	ZGtree->SetBranchAddress("phoEt",     &ZG_phoEt);
  	ZGtree->SetBranchAddress("phoEta",    &ZG_phoEta);
  	ZGtree->SetBranchAddress("phoPhi",    &ZG_phoPhi);
  	ZGtree->SetBranchAddress("lepPt",     &ZG_lepPt);
  	ZGtree->SetBranchAddress("lepEta",    &ZG_lepEta);
  	ZGtree->SetBranchAddress("lepPhi",    &ZG_lepPhi);
	ZGtree->SetBranchAddress("bosonPt",   &ZG_bosonPt);
  	ZGtree->SetBranchAddress("sigMT",     &ZG_sigMT);
  	ZGtree->SetBranchAddress("sigMET",    &ZG_sigMET);
  	ZGtree->SetBranchAddress("dPhiLepMET",&ZG_dPhiLepMET);
  	ZGtree->SetBranchAddress("dRPhoLep",  &ZG_dRPhoLep);
  	ZGtree->SetBranchAddress("HT",        &ZG_HT);
  	ZGtree->SetBranchAddress("nJet",      &ZG_nJet);
	ZGtree->SetBranchAddress("threeMass", &ZG_threeMass);
	ZGtree->SetBranchAddress("dilepMass", &ZG_dilepMass);
	ZGtree->SetBranchAddress("ISRJetPt",  &ZG_JetPt);
  	ZGtree->SetBranchAddress("mcPID",     &ZG_mcPID);
  	ZGtree->SetBranchAddress("mcEta",     &ZG_mcEta);
  	ZGtree->SetBranchAddress("mcPhi",     &ZG_mcPhi);
  	ZGtree->SetBranchAddress("mcPt",      &ZG_mcPt);
  	ZGtree->SetBranchAddress("mcMomPID",  &ZG_mcMomPID);
  	ZGtree->SetBranchAddress("mcGMomPID", &ZG_mcGMomPID);

  	for(unsigned ievt(0); ievt<ZGtree->GetEntries(); ++ievt){//loop on entries
		ZGtree->GetEntry(ievt);

		//if(doEB && fabs(ZG_phoEta) > 1.4442)continue;
		//if(!doEB && (fabs(ZG_phoEta)< 1.56 || fabs(ZG_phoEta) > 2.4))continue;
		double scalefactor = GetChannelSF(objectESF, channel, ZG_lepPt, ZG_lepEta, ZG_phoEt, ZG_phoEta);

		if(ZG_phoEt > 799)ZG_phoEt = 799;
		if(ZG_JetPt > 799)ZG_JetPt = 799;
		if(ZG_dRPhoLep < 0.8)continue;
		if(ZG_dilepMass < 80 || ZG_dilepMass >100)continue;
		double isrW = applyISR ? ISRWeightFromPt(ZG_JetPt, loadedISRWeights, loadedNormalization) : 1.0;
		double weight = ZG_MCweight*ZG_PUweight*scalefactor*isrW;

		bool   isTruePho(false);
		double mindR(0.3);
		unsigned phoIndex(0);
		for(unsigned iMC(0); iMC< ZG_mcPID->size(); iMC++){
			double dR1 = DeltaR((*ZG_mcEta)[iMC], (*ZG_mcPhi)[iMC], ZG_phoEta, ZG_phoPhi);
			double dE1 = fabs((*ZG_mcPt)[iMC] - ZG_phoEt)/ZG_phoEt;
			if(dR1 < mindR && dE1 < 0.2){mindR=dR1; phoIndex=iMC;}
		}
		if(mindR < 0.1){
			if((*ZG_mcPID)[phoIndex] == 22 && (fabs((*ZG_mcMomPID)[phoIndex]) == 23 || fabs((*ZG_mcMomPID)[phoIndex]) == 13 || fabs((*ZG_mcMomPID)[phoIndex])==999)){
				isTruePho=true;
			}
		}

		p_JetPt_ZG->Fill(ZG_JetPt, weight);	
		p_phoEt_ZG->Fill(ZG_phoEt, weight); 
		p_llmass_ZG->Fill(ZG_dilepMass, weight);
	}//loop on  events


	//************ Signal Tree **********************//
  	TChain *raretree = new TChain("ZTree");
	// using analysis_ISRMC.C
	raretree->Add(Form("/eos/uscms/store/user/tmishra/ISRweighting/%s/resTree_ISR_TTGJets_%d%s.root", whichChannel.c_str(), RunYear,whichVFP.c_str()));
	raretree->Add(Form("/eos/uscms/store/user/tmishra/ISRweighting/%s/resTree_ISR_TTJets_%d%s.root", whichChannel.c_str(), RunYear,whichVFP.c_str()));
	raretree->Add(Form("/eos/uscms/store/user/tmishra/ISRweighting/%s/resTree_ISR_WWG_%d%s.root", whichChannel.c_str(), RunYear,whichVFP.c_str()));
	raretree->Add(Form("/eos/uscms/store/user/tmishra/ISRweighting/%s/resTree_ISR_WZG_%d%s.root", whichChannel.c_str(),RunYear,whichVFP.c_str()));

	float rare_MCweight(0);
	float rare_PUweight(0);
  	float rare_phoEt(0);
  	float rare_phoEta(0);
  	float rare_phoPhi(0);
  	float rare_lepPt(0);
  	float rare_lepEta(0);
  	float rare_lepPhi(0);
  	float rare_sigMT(0);
  	float rare_sigMET(0);
  	float rare_dPhiLepMET(0);
  	float rare_dRPhoLep(0);
  	float rare_HT(0);
  	float rare_nJet(0);
	float rare_threeMass(0);
	float rare_dilepMass(0);
	float rare_JetPt(0);
  	std::vector<int>   *rare_mcPID=0;
  	std::vector<float> *rare_mcEta=0;
  	std::vector<float> *rare_mcPhi=0;
  	std::vector<float> *rare_mcPt=0;
  	std::vector<int>   *rare_mcMomPID=0;
  	std::vector<int>   *rare_mcGMomPID=0;

	raretree->SetBranchAddress("MCweight",  &rare_MCweight); 
	raretree->SetBranchAddress("PUweight",  &rare_PUweight);
  	raretree->SetBranchAddress("phoEt",     &rare_phoEt);
  	raretree->SetBranchAddress("phoEta",    &rare_phoEta);
  	raretree->SetBranchAddress("phoPhi",    &rare_phoPhi);
  	raretree->SetBranchAddress("lepPt",     &rare_lepPt);
  	raretree->SetBranchAddress("lepEta",    &rare_lepEta);
  	raretree->SetBranchAddress("lepPhi",    &rare_lepPhi);
  	raretree->SetBranchAddress("sigMT",     &rare_sigMT);
  	raretree->SetBranchAddress("sigMET",    &rare_sigMET);
  	raretree->SetBranchAddress("dPhiLepMET",&rare_dPhiLepMET);
  	raretree->SetBranchAddress("dRPhoLep",  &rare_dRPhoLep);
  	raretree->SetBranchAddress("HT",        &rare_HT);
  	raretree->SetBranchAddress("nJet",      &rare_nJet);
	raretree->SetBranchAddress("ISRJetPt",  &rare_JetPt);
	raretree->SetBranchAddress("threeMass", &rare_threeMass);
	raretree->SetBranchAddress("dilepMass", &rare_dilepMass);
  	raretree->SetBranchAddress("mcPID",    &rare_mcPID);
  	raretree->SetBranchAddress("mcEta",    &rare_mcEta);
  	raretree->SetBranchAddress("mcPhi",    &rare_mcPhi);
  	raretree->SetBranchAddress("mcPt",     &rare_mcPt);
  	raretree->SetBranchAddress("mcMomPID", &rare_mcMomPID);
  	raretree->SetBranchAddress("mcGMomPID",&rare_mcGMomPID);

  	for(unsigned ievt(0); ievt<raretree->GetEntries(); ++ievt){//loop on entries
		raretree->GetEntry(ievt);
		
		//if(doEB && fabs(rare_phoEta) > 1.4442)continue;
		//else if(!doEB && (fabs(rare_phoEta) < 1.56|| fabs(rare_phoEta) > 2.4))continue;
		double scalefactor = GetChannelSF(objectESF, channel, rare_lepPt, rare_lepEta, rare_phoEt, rare_phoEta);

		if(rare_phoEt > 799)rare_phoEt = 799;
		if(rare_JetPt > 799)rare_JetPt = 799;
		if(rare_dRPhoLep < 0.8)continue;
		if(rare_dilepMass < 80 || rare_dilepMass > 100)continue;
		// ISR weights are applied to V+gamma only (AN), not to the rare backgrounds
		double weight = rare_MCweight*rare_PUweight*scalefactor;

		p_phoEt_rare->Fill(rare_phoEt, weight); 
		p_JetPt_rare->Fill(rare_JetPt, weight);
		p_llmass_rare->Fill(rare_dilepMass, weight);

	}//loop on  events
	

	
	TCanvas *can_phoEt     = new TCanvas("can_phoEt",       "can_phoEt", 600,600); 
	TCanvas *can_JetPt        = new TCanvas("can_JetPt",          "can_JetPt",600,600); 
	TCanvas *can_JetPt_alter  = new TCanvas("can_JetPt_alter",    "can_JetPt_alter",600,600); 


	p_llmass_ZG->Add(p_llmass_rare);
        float scalefactor = 1.0;
        float scalefactorhighEt = 1.0;

	
	can_phoEt->cd();
	TPad *phoEt_pad1 = new TPad("phoEt_pad1", "phoEt_pad1", 0, 0.3, 1, 1.0);
	phoEt_pad1->SetBottomMargin(0);
        phoEt_pad1->Draw();
        phoEt_pad1->cd();
        gStyle->SetOptStat(0);
        phoEt_pad1->SetLogy();
        p_phoEt_data->GetXaxis()->SetTitle("ISR P_{T} (GeV)");
        p_phoEt_data->SetMinimum(0.05);

	p_phoEt_data->SetLineColor(kBlack);
	p_phoEt_data->SetMarkerStyle(20);
	p_phoEt_data->SetMarkerColor(kBlack);
	p_phoEt_data->Draw("P");
	p_phoEt_ZG->Add(p_phoEt_rare);
	p_phoEt_ZG->SetLineColor(6);
	p_phoEt_ZG->SetFillColor(6);
	p_phoEt_ZG->Scale(scalefactor);
	p_phoEt_ZG->Draw("hist same");
	p_phoEt_rare->SetLineColor(8);
	p_phoEt_rare->SetFillColor(8);
	p_phoEt_rare->Scale(scalefactor);
	p_phoEt_rare->Draw("hist same");
	p_phoEt_data->Draw("EP same");

	TLegend *leg_phoEt =  new TLegend(0.6,0.7,0.9,0.9);
	leg_phoEt->SetFillStyle(0);
	leg_phoEt->AddEntry(p_phoEt_data, "Data");
	leg_phoEt->AddEntry(p_phoEt_ZG,   "Z#gamma");
	leg_phoEt->AddEntry(p_phoEt_rare, "t#bar{t},t#bar{t}#gamma,WW#gamma,WZ#gamma");   
	leg_phoEt->Draw("same");
	
	can_phoEt->cd();
	TPad *phoEt_pad2 = new TPad("phoEt_pad2", "phoEt_pad2", 0, 0.05, 1, 0.3);
	phoEt_pad2->SetTopMargin(0);
	phoEt_pad2->SetBottomMargin(0.3);
	phoEt_pad2->Draw();
	phoEt_pad2->cd();
  	TLine *flatratio_phoEt = new TLine(0,1,800,1);
	TH1F *ratio_phoEt=(TH1F*)p_phoEt_data->Clone("transfer factor");
	ratio_phoEt->SetMarkerStyle(20);
	ratio_phoEt->SetLineColor(kBlack);
	ratio_phoEt->GetXaxis()->SetRangeUser(0,800);
	ratio_phoEt->GetYaxis()->SetRangeUser(0.95,1.05);
	ratio_phoEt->GetYaxis()->SetLabelSize(12);
	ratio_phoEt->SetMinimum(0.2);
	ratio_phoEt->SetMaximum(1.7);
	ratio_phoEt->Divide(p_phoEt_ZG);
	ratio_phoEt->SetTitle("");
	ratio_phoEt->GetYaxis()->SetTitle("Data/MC");
	ratio_phoEt->Draw();
	flatratio_phoEt->Draw("same");
//	can_phoEt->SaveAs(Form("/eos/uscms/store/user/tmishra/ISRweighting/PLOT_ISRweight_phoEt_%d%s.pdf",RunYear,whichVFP.c_str()));
	
	// ISR weight = (data - fakes) / Z+gamma MC per ISR pT bin, as in the 2016 analysis:
	// rare MC is not subtracted (taken before p_JetPt_ZG->Add(p_JetPt_rare) below).
	// Fakes are already removed by the per-event (1 - fakerate) data weight.
	// Normalization keeps the Z+gamma yield unchanged: sum_i w_i*N*ZG_i = ZG.
	const int nbinsISR = p_JetPt_data->GetNbinsX();
	std::vector<double> ISRwgt_norm(nbinsISR,1.0);
	for(int i=1; i<=nbinsISR; i++){
		double dataBin = p_JetPt_data->GetBinContent(i);
		double mcBin   = p_JetPt_ZG->GetBinContent(i);
		double ratio   = SafeRatioISR(dataBin, mcBin);
		double stat    = (mcBin>0 ? p_JetPt_data->GetBinError(i)/mcBin : 0.0);
		std::cout << "norm ratio " << i << " " << ratio << " stat " << stat
		          << "  (data " << dataBin << ", rare " << p_JetPt_rare->GetBinContent(i)
		          << ", ZG " << mcBin << ")" << std::endl;
		ISRwgt_norm[i-1] = ratio;
	}

	double normalization = SafeRatioISR(p_JetPt_ZG->Integral(), p_JetPt_data->Integral());

	if(!applyISR){
		if(WriteISRWeights(isrWgtFile, ISRwgt_norm, normalization))
			std::cout << "[INFO] wrote ISR weights to " << isrWgtFile << std::endl;
		else
			std::cout << "[WARN] failed to write ISR weights to " << isrWgtFile << std::endl;
	}

	gStyle->SetOptStat(0);
	can_JetPt->cd();
	TPad *JetPt_pad1 = new TPad("JetPt_pad1", "JetPt_pad1", 0, 0.3, 1, 1.0);
	JetPt_pad1->SetBottomMargin(0);
	JetPt_pad1->Draw();  
	JetPt_pad1->cd();  
	gStyle->SetOptStat(0);
	JetPt_pad1->SetLogy();
	p_JetPt_data->GetXaxis()->SetTitleOffset(0.9);
	p_JetPt_data->GetXaxis()->SetTitle("ISR P_{T} (GeV)");
	p_JetPt_data->SetMinimum(0.05);
	p_JetPt_data->SetLineColor(kBlack);
	p_JetPt_data->SetMarkerStyle(20);
	p_JetPt_data->SetMarkerColor(kBlack);
	p_JetPt_data->Draw("P");
	p_JetPt_ZG->Add(p_JetPt_rare);

	p_JetPt_ZG->SetLineColor(6);
	p_JetPt_ZG->SetFillStyle(1001);
	p_JetPt_ZG->SetFillColor(6);
	p_JetPt_ZG->Scale(scalefactor);
	p_JetPt_ZG->Draw("hist same");
	p_JetPt_rare->SetLineColor(8);
	p_JetPt_rare->SetFillStyle(1001);
	p_JetPt_rare->SetFillColor(8);
	p_JetPt_rare->Scale(scalefactor);
	p_JetPt_rare->Draw("hist same");
	p_JetPt_data->Draw("EP same");
        if(RunYear==2016 and preVFP == 1)       CMS_lumi( JetPt_pad1,1,1, 11 );
        else if(RunYear==2016 and preVFP == 0)  CMS_lumi( JetPt_pad1,2,1, 11 );
        else if(RunYear==2017)                  CMS_lumi( JetPt_pad1,3,1, 11 );
        else if(RunYear==2018)                  CMS_lumi( JetPt_pad1,4,1, 11 );
	
	if(applyISR == false)  cout<<"                 Normalization = "<<normalization<<";    }"<<endl;
	
	cout<<"ZG integral is "<<p_JetPt_ZG->Integral()<<endl;
	cout<<"Data integral is "<<p_JetPt_data->Integral()<<endl;
	cout<<"Ratio "<<p_JetPt_ZG->Integral()/p_JetPt_data->Integral()  <<endl;
	TLegend *leg_JetPt =  new TLegend(0.6,0.7,0.9,0.9);
	leg_JetPt->SetFillStyle(0);
	leg_JetPt->AddEntry(p_JetPt_data, "Data");
	leg_JetPt->AddEntry(p_JetPt_ZG,   "Z#gamma");
	leg_JetPt->AddEntry(p_JetPt_rare, "t#bar{t},t#bar{t}#gamma,WW#gamma,WZ#gamma");   
	leg_JetPt->Draw("same");

	can_JetPt->cd();
	TPad *JetPt_pad2 = new TPad("JetPt_pad2", "JetPt_pad2", 0, 0.05, 1, 0.3);
	JetPt_pad2->SetTopMargin(0);
	JetPt_pad2->SetBottomMargin(0.4);
	JetPt_pad2->Draw();
	JetPt_pad2->cd();
  	TLine *flatratio_JetPt = new TLine(0,1,800,1);
	TH1F *ratio_JetPt=(TH1F*)p_JetPt_data->Clone("transfer factor");
	ratio_JetPt->SetMarkerStyle(20);
	ratio_JetPt->SetLineColor(kBlack);
	ratio_JetPt->GetXaxis()->SetRangeUser(0,800);
	ratio_JetPt->GetYaxis()->SetRangeUser(0.5,1.35);
	ratio_JetPt->GetYaxis()->SetLabelSize(12);
	
	for(int i=1; i<=nbinsISR; ++i){
		double mcBin = p_JetPt_ZG->GetBinContent(i);
		if(mcBin<=0){
			ratio_JetPt->SetBinContent(i,1.0);
			ratio_JetPt->SetBinError(i,0.0);
		}else{
			ratio_JetPt->SetBinContent(i, p_JetPt_data->GetBinContent(i)/mcBin);
			ratio_JetPt->SetBinError(i, p_JetPt_data->GetBinError(i)/mcBin);
		}
	}
	ratio_JetPt->SetTitle("");
	ratio_JetPt->GetYaxis()->SetTitle("Data/MC");
	ratio_JetPt->Draw();
	flatratio_JetPt->Draw("same");

	TLine *ratioValue_JetPt = new TLine(0, p_JetPt_data->Integral() / p_JetPt_ZG->Integral(),
                                    800, p_JetPt_data->Integral() / p_JetPt_ZG->Integral());
	ratioValue_JetPt->SetLineColor(kRed);
	ratioValue_JetPt->Draw("same");

	std::ostringstream outPdf;
	outPdf << "/eos/uscms/store/user/tmishra/ISRweighting/PLOT_ISRweight_" << GetEraTag(RunYear, preVFP) << "_" <<whichChannel;
	if(applyISR) outPdf << "_ISRweighted";
	outPdf << ".png";
	can_JetPt->SaveAs(outPdf.str().c_str());

	delete ratio_JetPt;
	delete flatratio_JetPt;
	delete leg_JetPt;
	delete JetPt_pad2;
	delete JetPt_pad1;
	delete can_JetPt;
	delete can_JetPt_alter;
	delete can_phoEt;
	delete tree;
	delete ZGtree;
	delete raretree;
}
int main(int argc, char** argv){
	int runYear = atoi(argv[1]);
	bool preVFP = atoi(argv[2]) == 1;
	bool applyISR = atoi(argv[3]) == 1;
	int channel = atoi(argv[4]);
	plot_ISRweight(runYear, preVFP, applyISR, channel);
	return 0;
}
