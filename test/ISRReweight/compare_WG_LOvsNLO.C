// LO (madgraphMLM) vs NLO (amcatnloFXFX) W+gamma: ISR pT and related shapes
// at the W+gamma selection of VGamma/src/analysis_VGamma_Run2.C, for the reviewer
// question on applying Z+gamma(NLO)-derived ISR weights to LO W+gamma.
//
// Selection and ISR pT copied from analysis_VGamma_Run2.C (signal photon with
// pixel/GSF/FSR vetoes and trigger matching, signal lepton, dR(pho,lep) > 0.8,
// eg Z veto, MET filters; ISR pT = vector sum of signal jets not matched to the
// prompt W/gamma decay products), plus the analysis_VGBkg.C cuts: photon
// pT > 35 GeV, |eta| < 1.4442, truth-matched prompt photon.
// LO stitching as in the skim: WGToLNuG for reco photon Et <= 50, WGJet40 above,
// WGJet130 always; LO weight = xsec * k-factor(photon pT) / N.
// NLO weight = sign(genWeight); the signed event count is stored in h_count so
// that the merge step can normalise each sample (h_count: bin 1 = sum of signs).
// No PU or object scale factors: both samples are treated the same way.
// Triggers are off by default (kUseTrigger): the eg trigger bits of the NLO ggNtuples
// (newer ntuplizer) do not match the LO ones (HLTPho bit 14 fires for 28% vs 43% of
// the same offline selection), and trigger emulation is not a generator effect.
// ISRJetPtReco: ISR pT with jets within dR < 0.4 of the selected photon and lepton
// removed, as in the derivation (analysis_ISRMC.C), independent of the generator record.
//
// LO samples also fill, event by event, <var>_<ch>_<region>_LOw (x data/NLO-Z+gamma ISR weight w of the era,
// ISRweights_<era>_<ch>.txt) and <var>_<ch>_<region>_LOr (x LO->NLO factor r, as Background/analysis_VGBkg.C),
// for the closure of r in variables other than ISR pT.
//
// Output: one ROOT file per job with histograms <var>_<ch>_<region>[_LOw|_LOr].
//   ./compare_WG_LOvsNLO.exe <RunYear> <preVFP 0|1> <WGToLNuG|WGJet40|WGJet130|NLO> <firstEntry> <lastEntry> <output.root> [nloList]

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_jet.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"

// as in VGamma/src/analysis_VGamma_Run2.C
static bool nisrMatch(float jetEta, float jetPhi, std::vector<mcData>& genParticles){
	for(std::vector<mcData>::iterator itMC = genParticles.begin(); itMC != genParticles.end(); itMC++){
		bool isPromptFinal(false);
		int momid = abs(itMC->getMomPID());
		int pid   = abs(itMC->getPID());
		if(pid <= 5 && (momid==6 || momid==23 || momid==24 || momid==25))isPromptFinal = true;
		else if(pid == 11 && (momid<=6 || momid==23 || momid==24 || momid==22 || momid==15))isPromptFinal = true;
		else if(pid == 13 && (momid<=6 || momid==23 || momid==24 || momid==22 || momid==15))isPromptFinal = true;
		else if(pid == 15 && (momid<=6 || momid==23 || momid==24 || momid==22))isPromptFinal = true;
		else if(pid == 22 && (momid==11 || momid==13 || momid==15 || momid==23 || momid==24))isPromptFinal = true;
		if(!isPromptFinal)continue;
		if(DeltaR(jetEta, jetPhi, itMC->getEta(), itMC->getPhi()) < 0.3)return true;
	}
	return false;
}

static float getWGKfactor(float phoPt){
	if     (phoPt < 190.) return 1.40;
	else if(phoPt < 250.) return 1.37;
	else if(phoPt < 400.) return 1.29;
	else if(phoPt < 700.) return 1.26;
	else                  return 1.15;
}

// truth-matched prompt photon, as in Background/analysis_VGBkg.C
static bool isTruePhoton(recoPhoton& pho, std::vector<mcData>& genParticles){
	double mindR(0.3);
	std::vector<mcData>::iterator best = genParticles.end();
	for(std::vector<mcData>::iterator itMC = genParticles.begin(); itMC != genParticles.end(); itMC++){
		double dR = DeltaR(itMC->getEta(), itMC->getPhi(), pho.getEta(), pho.getPhi());
		double dE = fabs(itMC->getPt() - pho.getCalibEt())/pho.getCalibEt();
		if(dR < mindR && dE < 0.5){mindR = dR; best = itMC;}
	}
	if(best == genParticles.end() || mindR >= 0.2 || best->getPID() != 22)return false;
	int mom = abs(best->getMomPID());
	return (mom <= 6 || mom == 21 || mom == 999 || mom == 11 || mom == 13 || mom == 15 || mom == 23 || mom == 24);
}

const double kISRedges[7] = {0, 50, 100, 150, 200, 250, 300};

// LO -> NLO factors r, as LOtoNLOFactor() in Background/analysis_VGBkg.C (ch 0 = eg, 1 = mg), in ISRJetPtReco
static double LOtoNLOFactor(int ch, float isrPt){
	static const double eg[7] = {1.06025, 0.94482, 0.89026, 0.89977, 0.89691, 0.89908, 0.87877};
	static const double mg[7] = {1.06165, 0.93595, 0.89678, 0.87963, 0.85610, 0.86764, 0.85067};
	int i = 6;
	while(i > 0 && isrPt < kISRedges[i]) i--;
	return ch == 0 ? eg[i] : mg[i];
}

// ISR weights w (x normalization) of one era and channel, as in plot_WG_LOvsNLO.C
static std::vector<double> readISRWeights(const std::string& era, int ch){
	std::vector<double> w;
	double norm = 1;
	std::ifstream in(Form("/eos/uscms/store/user/tmishra/ISRweighting/ISRweights_%s_%s.txt", era.c_str(), ch == 0 ? "egamma" : "mgamma"));
	std::string line;
	while(std::getline(in, line)){
		if(line.rfind("# normalization", 0) == 0){ std::istringstream s(line.substr(15)); s >> norm; continue; }
		if(line.empty() || line[0] == '#')continue;
		std::istringstream s(line);
		double lo, hi, wt;
		if(s >> lo >> hi >> wt) w.push_back(wt*norm);
	}
	if(w.size() != 7){ std::cerr << "unexpected ISR weight file for " << era << " ch " << ch << std::endl; exit(1); }
	return w;
}

static double isrWeightAt(const std::vector<double>& w, float isrPt){
	int i = 6;
	while(i > 0 && isrPt < kISRedges[i]) i--;
	return w[i];
}

const bool kUseTrigger = false;

const int nCh = 2, nReg = 3;
const char* chName[nCh]  = {"eg", "mg"};
const char* regName[nReg] = {"presel", "CR", "SR"};   // CR: MET < 70; SR: MET > 120, MT > 100

struct HistSet {
	TH1D *isr, *isrFine, *isrReco, *nJet, *HT, *phoEt, *MET, *bosonPt;
	void book(const char* tag){
		const double isrBins[] = {0,50,100,150,200,250,300,800};
		isr     = new TH1D(Form("ISRJetPt_%s", tag),     ";ISR p_{T} (GeV)", 7, isrBins);
		isrFine = new TH1D(Form("ISRJetPtFine_%s", tag), ";ISR p_{T} (GeV)", 32, 0, 800);
		isrReco = new TH1D(Form("ISRJetPtReco_%s", tag), ";ISR p_{T}, reco-cleaned (GeV)", 7, isrBins);
		nJet    = new TH1D(Form("nJet_%s", tag),         ";N_{jet}", 10, -0.5, 9.5);
		HT      = new TH1D(Form("HT_%s", tag),           ";H_{T} (GeV)", 30, 0, 1500);
		phoEt   = new TH1D(Form("phoEt_%s", tag),        ";photon p_{T} (GeV)", 30, 0, 1200);
		MET     = new TH1D(Form("MET_%s", tag),          ";p_{T}^{miss} (GeV)", 30, 0, 900);
		bosonPt = new TH1D(Form("bosonPt_%s", tag),      ";gen W p_{T} (GeV)", 30, 0, 900);
		for(TH1D* h : {isr, isrFine, isrReco, nJet, HT, phoEt, MET, bosonPt})h->Sumw2();
	}
	void fill(float isrPt, float isrPtReco, float nj, float ht, float et, float met, float wpt, double w){
		// overflow into the last bin, as the >300 GeV ISR weight bin
		isr->Fill(std::min(isrPt, 799.f), w);
		isrFine->Fill(std::min(isrPt, 799.f), w);
		isrReco->Fill(std::min(isrPtReco, 799.f), w);
		nJet->Fill(std::min(nj, 9.f), w);
		HT->Fill(std::min(ht, 1499.f), w);
		phoEt->Fill(std::min(et, 1199.f), w);
		MET->Fill(std::min(met, 899.f), w);
		bosonPt->Fill(std::min(wpt, 899.f), w);
	}
};

void compare_WG_LOvsNLO(int RunYear, bool preVFP, const char* sample, long long first, long long last,
                        const char* outName, const char* nloList){
	const std::string vfp = (RunYear == 2016) ? (preVFP ? "preVFP" : "postVFP") : "";
	const std::string s(sample);
	const bool isNLO = (s == "NLO");
	int mcType = -1;
	if(s == "WGToLNuG")     mcType = MCType::WGJetInclusive;
	else if(s == "WGJet40") mcType = MCType::WGJet40;
	else if(s == "WGJet130")mcType = MCType::WGJet130;
	else if(!isNLO){ std::cerr << "unknown sample " << s << std::endl; exit(1); }

	TChain* es = new TChain("ggNtuplizer/EventTree");
	if(isNLO){
		std::ifstream in(nloList);
		std::string line;
		while(std::getline(in, line)) if(!line.empty()) es->Add(line.c_str());
	}
	else es->Add(Form("/eos/uscms/store/user/tmishra/InputFilesMC/%s/%s_%d%s.root", sample, sample, RunYear, vfp.c_str()));
	const long long nEvts = es->GetEntries();
	if(nEvts == 0){ std::cerr << "no input events" << std::endl; exit(1); }
	if(last < 0 || last > nEvts) last = nEvts;
	std::cout << sample << " " << RunYear << vfp << " entries " << nEvts << ", processing [" << first << ", " << last << ")" << std::endl;

	TFile* fOut = TFile::Open(outName, "RECREATE");
	// bin 1: sum of sign(genWeight) over all processed events (NLO) or number of events (LO)
	TH1D* h_count = new TH1D("h_count", "", 1, 0, 1);
	h_count->Sumw2();
	HistSet H[nCh][nReg], HW[nCh][nReg], HR[nCh][nReg];
	std::vector<double> isrW[nCh];
	for(int c = 0; c < nCh; c++) for(int r = 0; r < nReg; r++){
		H[c][r].book(Form("%s_%s", chName[c], regName[r]));
		if(isNLO)continue;
		HW[c][r].book(Form("%s_%s_LOw", chName[c], regName[r]));
		HR[c][r].book(Form("%s_%s_LOr", chName[c], regName[r]));
	}
	if(!isNLO) for(int c = 0; c < nCh; c++) isrW[c] = readISRWeights(Form("%d%s", RunYear, vfp.c_str()), c);

	rawData raw(es, MC);
	float genWeight(1);
	es->SetBranchAddress("genWeight", &genWeight);

	std::vector<mcData> MCData;
	std::vector<recoPhoton> Photon;
	std::vector<recoMuon> Muon;
	std::vector<recoEle> Ele;
	std::vector<recoJet> JetCollection;

	for(long long ievt = first; ievt < last; ++ievt){
		if((ievt - first) % 200000 == 0) std::cout << " -- event " << ievt << std::endl;
		raw.GetData(es, ievt);
		const double sign = isNLO ? (genWeight < 0 ? -1. : 1.) : 1.;
		h_count->Fill(0.5, sign);
		if(raw.nPho < 1)continue;

		MCData.clear(); Photon.clear(); Muon.clear(); Ele.clear(); JetCollection.clear();
		for(int i = 0; i < raw.nMC; i++)  MCData.push_back(mcData(raw, i));
		for(int i = 0; i < raw.nPho; i++) Photon.push_back(recoPhoton(raw, i));
		for(int i = 0; i < raw.nMu; i++)  Muon.push_back(recoMuon(raw, i));
		for(int i = 0; i < raw.nEle; i++) Ele.push_back(recoEle(raw, i));
		for(int i = 0; i < raw.nJet; i++) JetCollection.push_back(recoJet(raw, i));
		const float MET = raw.pfMET, METPhi = raw.pfMETPhi;

		int Nmedpho(0);
		for(std::vector<recoPhoton>::iterator itpho = Photon.begin(); itpho != Photon.end(); ++itpho){
			if(!itpho->isMedium())continue;
			if(itpho->getR9() < 0.5 || itpho->getR9() > 1.0)continue;
			if(fabs(itpho->getEta()) > 1.4442 || itpho->getCalibEt() < 40)continue;
			if(itpho->getSigma() < 0.005)continue;
			if(itpho->PixelSeed() != 0)continue;
			Nmedpho += 1;
		}
		if(Nmedpho >= 2)continue;

		float WGPt(-1);
		for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC != MCData.end(); itMC++)
			if(abs(itMC->getPID()) == 24 || itMC->getPID() == 23) WGPt = itMC->getPt();

		bool hasegPho(false), hasmgPho(false);
		std::vector<recoPhoton>::iterator egPho = Photon.begin(), mgPho = Photon.begin();
		for(std::vector<recoPhoton>::iterator itpho = Photon.begin(); itpho != Photon.end(); ++itpho){
			if(itpho->getR9() < 0.5)continue;
			if(!itpho->passSignalSelection())continue;
			bool PixelVeto = itpho->PixelSeed() == 0;
			bool GSFveto(true), FSRVeto(true);
			for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
				if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.02)GSFveto = false;
				if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3)FSRVeto = false;
			}
			for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
				if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getPt() > 2.0)FSRVeto = false;
			if(!(GSFveto && PixelVeto && FSRVeto))continue;
			bool passEGPhoTrig = !kUseTrigger || itpho->fireDoubleTrg(5) || itpho->fireDoubleTrg(6);
			bool passMGPhoTrig = !kUseTrigger || ((RunYear == 2016)
			                   ? (itpho->fireDoubleTrg(28) || itpho->fireDoubleTrg(29) || itpho->fireDoubleTrg(30))
			                   : (itpho->fireDoubleTrg(27) || itpho->fireDoubleTrg(36)));
			if(passEGPhoTrig && !hasegPho){ hasegPho = true; egPho = itpho; }
			if(passMGPhoTrig && !hasmgPho){ hasmgPho = true; mgPho = itpho; }
		}
		// LO stitching on reco photon Et, as in the skim
		if(hasegPho){
			if(egPho->getCalibEt() > 50  && mcType == MCType::WGJetInclusive)hasegPho = false;
			if(egPho->getCalibEt() <= 50 && mcType == MCType::WGJet40)hasegPho = false;
		}
		if(hasmgPho){
			if(mgPho->getCalibEt() > 50  && mcType == MCType::WGJetInclusive)hasmgPho = false;
			if(mgPho->getCalibEt() <= 50 && mcType == MCType::WGJet40)hasmgPho = false;
		}

		bool hasEle(false);
		std::vector<recoEle>::iterator signalEle = Ele.begin();
		if(hasegPho){
			for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
				if((itEle->isEB() && itEle->getR9() < 0.5) || (itEle->isEE() && itEle->getR9() < 0.8))continue;
				bool passEleTrig = !kUseTrigger || ((RunYear == 2016) ? (itEle->fireTrgs(21) || itEle->fireTrgs(22))
				                                   : (itEle->fireTrgs(43) || itEle->fireTrgs(44)));
				if(passEleTrig && itEle->passSignalSelection()){ hasEle = true; signalEle = itEle; break; }
			}
		}
		bool hasMu(false);
		std::vector<recoMuon>::iterator signalMu = Muon.begin();
		if(hasmgPho){
			for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
				bool passMuTrig = !kUseTrigger || ((RunYear == 2016)
				                ? (itMu->fireSingleTrg(2) || itMu->fireSingleTrg(21) || itMu->fireSingleTrg(22))
				                : (itMu->fireSingleTrg(0) || itMu->fireSingleTrg(32)));
				if(passMuTrig && itMu->passSignalSelection()){ hasMu = true; signalMu = itMu; break; }
			}
		}

		TLorentzVector JetVec(0,0,0,0);
		for(std::vector<recoJet>::iterator itJet = JetCollection.begin(); itJet != JetCollection.end(); ++itJet){
			if(!itJet->passSignalSelection())continue;
			if(!nisrMatch(itJet->getEta(), itJet->getPhi(), MCData))JetVec = JetVec + itJet->getP4();
		}
		const float ISRJetPt = JetVec.Pt();

		// channel 0 = eg, 1 = mg; mg only if the event is not already eg (as in the skim)
		bool hasDoubleEG(false);
		for(int c = 0; c < nCh; c++){
			recoPhoton* pho(0);
			float lepPt(0), lepEta(0), lepPhi(0);
			if(c == 0){
				if(!(hasegPho && hasEle && (!kUseTrigger || ((raw.HLTPho >> 14)&1) == 1)))continue;
				if(DeltaR(egPho->getEta(), egPho->getPhi(), signalEle->getEta(), signalEle->getPhi()) <= 0.8)continue;
				if(fabs((egPho->getP4() + signalEle->getP4()).M() - 91.188) <= 10.0)continue;
				pho = &(*egPho); lepPt = signalEle->getPt(); lepEta = signalEle->getEta(); lepPhi = signalEle->getPhi();
			}
			else{
				bool passMGHLT = !kUseTrigger || ((RunYear == 2016) ? ((((raw.HLTEleMuX >> 8)&1) != 0) || (((raw.HLTEleMuX >> 51)&1) != 0))
				                                  : ((((raw.HLTEleMuX >> 8)&1) != 0) || (((raw.HLTEleMuX >> 57)&1) != 0)));
				if(!(hasmgPho && hasMu && passMGHLT) || hasDoubleEG)continue;
				if(DeltaR(mgPho->getEta(), mgPho->getPhi(), signalMu->getEta(), signalMu->getPhi()) <= 0.8)continue;
				pho = &(*mgPho); lepPt = signalMu->getPt(); lepEta = signalMu->getEta(); lepPhi = signalMu->getPhi();
			}
			if(!raw.passMETFilter(RunYear, raw.metFilters))continue;
			if(c == 0) hasDoubleEG = true;

			const float phoEt = pho->getCalibEt();
			if(phoEt < 35 || fabs(pho->getEta()) > 1.4442)continue;
			if(!isTruePhoton(*pho, MCData))continue;

			float HT(0), nJet(0);
			TLorentzVector JetVecReco(0,0,0,0);
			for(std::vector<recoJet>::iterator itJet = JetCollection.begin(); itJet != JetCollection.end(); ++itJet){
				if(!itJet->passSignalSelection())continue;
				if(DeltaR(itJet->getEta(), itJet->getPhi(), lepEta, lepPhi) <= 0.4)continue;
				if(DeltaR(itJet->getEta(), itJet->getPhi(), pho->getEta(), pho->getPhi()) <= 0.4)continue;
				nJet += 1; HT += itJet->getPt();
				JetVecReco += itJet->getP4();
			}
			const float MT = sqrt(2*MET*lepPt*(1 - std::cos(DeltaPhi(lepPhi, METPhi))));

			const float ISRJetPtReco = JetVecReco.Pt();
			const double w = isNLO ? sign : getWGKfactor(phoEt)*MC_XS[mcType];
			const bool inReg[nReg] = {true, MET < 70, MET > 120 && MT > 100};
			for(int r = 0; r < nReg; r++){
				if(!inReg[r])continue;
				H[c][r].fill(ISRJetPt, ISRJetPtReco, nJet, HT, phoEt, MET, WGPt, w);
				if(isNLO)continue;
				// w and r are applied at ISRJetPtReco, as in the V+gamma skim and the background macros
				HW[c][r].fill(ISRJetPt, ISRJetPtReco, nJet, HT, phoEt, MET, WGPt, w*isrWeightAt(isrW[c], ISRJetPtReco));
				HR[c][r].fill(ISRJetPt, ISRJetPtReco, nJet, HT, phoEt, MET, WGPt, w*LOtoNLOFactor(c, ISRJetPtReco));
			}
		}
	}
	fOut->Write();
	fOut->Close();
	std::cout << "wrote " << outName << std::endl;
}

int main(int argc, char** argv){
	if(argc < 7){
		std::cerr << "Usage: " << argv[0] << " <RunYear> <preVFP 0|1> <WGToLNuG|WGJet40|WGJet130|NLO> <first> <last> <out.root> [nloList]\n";
		return 1;
	}
	compare_WG_LOvsNLO(atoi(argv[1]), atoi(argv[2]), argv[3], atoll(argv[4]), atoll(argv[5]), argv[6], argc > 7 ? argv[7] : "");
	return 0;
}
