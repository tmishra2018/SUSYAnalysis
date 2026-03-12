#include<string>
#include "../../include/analysis_commoncode.h"

#define NTOY 1000
bool useGaussFit=false;


namespace {

double fakerate_pt(double pt, double slope, double constant, double index)
{
	double x = TMath::Max(pt, 0.000001);
	return pow(slope*x + constant, index);
}

bool loadBinnedFakeRate(const std::string& fileName, std::vector<double>& bins, std::vector<double>& weights)
{
	bins.clear();
	weights.clear();

	std::ifstream infile(fileName.c_str());
	if(!infile.is_open()) return false;

	std::vector<double> denBins;
	std::vector<double> denVals;
	std::vector<double> numBins;
	std::vector<double> numVals;

	std::string bintype, numtype;
	double lowcut(0), signal1(0), error1(0), fitmean(0), fitrms(0);
	while(infile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms){
		if(numtype == "den"){
			denBins.push_back(lowcut);
			denVals.push_back(fitmean);
		}
		else if(numtype == "num"){
			numBins.push_back(lowcut);
			numVals.push_back(fitmean);
		}
	}
	infile.close();

	if(denBins.empty() || denBins.size() != numBins.size()) return false;

	for(unsigned i(0); i<denBins.size(); ++i){
		int matchIndex(-1);
		for(unsigned j(0); j<numBins.size(); ++j){
			if(fabs(numBins[j] - denBins[i]) < 1e-6){
				matchIndex = (int)j;
				break;
			}
		}
		if(matchIndex < 0) return false;

		bins.push_back(denBins[i]);
		double den = denVals[i];
		double num = numVals[matchIndex];
		weights.push_back((den != 0) ? num/den : 0);
	}
	return true;
}


double getBinnedWeight(double x, const std::vector<double>& bins, const std::vector<double>& weights)
{
	if(bins.empty() || weights.empty()) return 0;
	for(unsigned i(0); i<bins.size(); ++i){
		if(i+1 < bins.size()){
			if(x > bins[i] && x < bins[i+1]) return weights[i];
		}
		else if(x > bins[i]) return weights[i];
	}
	return weights[0];
}

}

void analysis_eleBkg_HT(){
	SetRunConfig();
	setTDRStyle();
 	gSystem->Load("../../lib/libAnaClasses.so");
   	int channelType = ichannel; // eg = 1; mg =2;
        std::string whichVFP;
        if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
        if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
        if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

	std::ostringstream elefake_config;
        elefake_config.str("");
        elefake_config << "/eos/uscms/store/user/tmishra/elefakepho/DATAResult"<<RunYear<<whichVFP<<"/EleFakeRate-Data-ByPt-EB.txt";
        std::ifstream elefake_file(elefake_config.str().c_str());
        double scalefactor(0);
        double ptslope(0);
        double ptconstant(0);
        double ptindex(0);
        std::string line;
        if(elefake_file.is_open()){
                while(std::getline(elefake_file, line)){
                        if(line.find("scalefactor") != std::string::npos){
                                std::size_t pos = line.find_last_of(' ');
                                if(pos != std::string::npos) scalefactor = atof(line.substr(pos+1).c_str());
                        }
                        else if(line.find("ptslope") != std::string::npos){
                                std::stringstream ss(line);
                                std::string label;
                                ss >> label >> ptslope;
                        }
                        else if(line.find("ptconstant") != std::string::npos){
                                std::stringstream ss(line);
                                std::string label;
                                ss >> label >> ptconstant;
                        }
                        else if(line.find("ptindex") != std::string::npos){
                                std::stringstream ss(line);
                                std::string label;
                                ss >> label >> ptindex;
                        }
                }
        }

        elefake_file.close();

	std::vector<double> etaBins, etaWeights, htBins, htWeights;
	bool loadedEta = loadBinnedFakeRate(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-ker-eta-60-120.txt",RunYear,whichVFP.c_str()), etaBins, etaWeights);
	bool loadedHT = loadBinnedFakeRate(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/EleFakeRate-DATA-Bw-ker-ht-60-120.txt",RunYear,whichVFP.c_str()), htBins, htWeights);

	std::vector<double> toy_scalefactor(NTOY, scalefactor);
	std::vector<double> toy_ptslope(NTOY, ptslope);
	std::vector<double> toy_ptconstant(NTOY, ptconstant);
	std::vector<double> toy_ptindex(NTOY, ptindex);
	std::ifstream elefake_toyfile(Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s/ToyFakeRate_Data_EB_HT-added.txt",RunYear,whichVFP.c_str()));
 	if(elefake_toyfile.is_open()){
		for(int i(0); i<NTOY; i++){
			double toyScale(0), toySlope(0), toyConst(0), toyIndex(0);
			if(!(elefake_toyfile >> toyScale >> toySlope >> toyConst >> toyIndex)) break;
			toy_scalefactor[i] = toyScale;
			toy_ptslope[i] = toySlope;
			toy_ptconstant[i] = toyConst;
			toy_ptindex[i] = toyIndex;
		}
 	}
 	elefake_toyfile.close();
 
 	//*********** histo list **********************//
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
		// different MET, MT and lep pT cuts
		double weight_pt = fakerate_pt(phoEt, ptslope, ptconstant, ptindex);
		double weight_eta = getBinnedWeight(fabs(phoEta), etaBins, etaWeights);
		double weight_ht = getBinnedWeight(HT, htBins, htWeights);
		double w_ele = scalefactor*weight_pt*weight_eta*weight_ht;
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
			double toy_pt = fakerate_pt(phoEt, toy_ptslope[it], toy_ptconstant[it], toy_ptindex[it]);
			double toy_ele = toy_scalefactor[it]*toy_pt*weight_eta*weight_ht;
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
				double toy_pt = fakerate_pt(phoEt, toy_ptslope[it], toy_ptconstant[it], toy_ptindex[it]);
				double toy_ele = toy_scalefactor[it]*toy_pt*weight_eta*weight_ht;
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
