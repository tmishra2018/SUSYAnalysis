#include "../../include/analysis_commoncode.h"

#define NTOY 1000
bool useGaussFit;
bool dopostfit=true;

void pred_eleBkg(){
	SetSignalConfig();
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);
	setTDRStyle();
  	gSystem->Load("../../lib/libAnaClasses.so");
  	int channelType = ichannel; // eg = 1; mg =2;
        std::string whichVFP;	
        if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
        else if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
        else whichVFP = "";
  	/************************************/
	/*	double normfactor = par[0]; */  
  	/*	double slope = par[1];				*/
  	/*	double constant = par[2];			*/
  	/*	double index = par[3];				*/
  	/*	double coeff = par[4]; 				*/
	/*	double vtx_constant = par[5];	*/
	/*	double vtx_slope = par[6];		*/
  	/**********************************/
	std::ifstream elefake_file;
	elefake_file.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/EleFakeRate-Data-ByPtVtx-EB.txt",RunYear,whichVFP.c_str()));
	double scalefactor(0);
	double ptslope(0);
	double ptconstant(0);
	double ptindex(0);
	double ptcoeff(0);
	double vtxconst(0);
	double vtxslope(0);
	std::string variabletype;
	double variablevalue;
	if(elefake_file.is_open()){
  	for(int i(0); i<7; i++){ 
			elefake_file >> variabletype >> variablevalue; 
			if(variabletype.find("scalefactor")!=std::string::npos)scalefactor = variablevalue;
			else if(variabletype.find("ptslope")!=std::string::npos)ptslope = variablevalue;
			else if(variabletype.find("ptconstant")!=std::string::npos)ptconstant = variablevalue;
			else if(variabletype.find("ptindex")!=std::string::npos)ptindex = variablevalue;
			else if(variabletype.find("ptcoeff")!=std::string::npos)ptcoeff = variablevalue;
			else if(variabletype.find("vtxconst")!=std::string::npos)vtxconst = variablevalue;
			else if(variabletype.find("vtxslope")!=std::string::npos)vtxslope = variablevalue;
	  }
	}
	elefake_file.close();
	// fakerate_func  defined in ../../include/analysis_fakes.h
	TF3 h_nominal_fakerate("h_nominal_fakerate", fakerate_func,10,1000,0,100,0,1.5,7);
	h_nominal_fakerate.SetParameters(scalefactor, ptslope, ptconstant, ptindex, ptcoeff, vtxconst, vtxslope);

	TF3 *h_toymc_fakerate[NTOY];
	std::ostringstream funcname;
	std::ifstream elefake_toyfile;
	elefake_toyfile.open(Form("/eos/uscms/store/user/tmishra/elefakepho/DataResult%d%s/ToyFakeRate_Data_EB.txt",RunYear,whichVFP.c_str()));
	if(elefake_toyfile.is_open()){
  	for(int i(0); i<NTOY; i++){ 
			elefake_toyfile >> scalefactor >> ptslope >> ptconstant >> ptindex >>  vtxconst >> vtxslope;
			funcname.str("");
			funcname << "h_toymc_fakerate_" << i;
			h_toymc_fakerate[i] = new TF3(funcname.str().c_str(), fakerate_func,10,1000,0,100,0,1.5,7);
			h_toymc_fakerate[i]->SetParameters(scalefactor, ptslope, ptconstant, ptindex, 1, vtxconst, vtxslope); 
	  }
	}
	elefake_toyfile.close();

	//*********** histo list **********************//
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
	//************ Proxy Tree **********************//
	TChain *proxytree = new TChain("proxyTree");
        if(channelType==1)proxytree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s.root",RunYear,whichVFP.c_str()));
        //if(channelType==1)proxytree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s_isFakeProxyOLD.root",RunYear,whichVFP.c_str()));
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
	float nJet(0);
	
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
	proxytree->SetBranchAddress("nJet",      &nJet);
 
	for (unsigned ievt(0); ievt<proxytree->GetEntries(); ++ievt){//loop on entries
		proxytree->GetEntry(ievt);
		p_PU->Fill(nVertex);
		/** cut flow *****/
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
		if(sigMT < lowMt)continue;
		if(highMt > 0 && sigMT > highMt)continue;
		if(lepPt < lowPt)continue;
		if(highPt > 0 && lepPt > highPt)continue;

		double w_ele = h_nominal_fakerate(phoEt, nVertex, fabs(phoEta));
		p_MET->Fill(sigMET, w_ele);

		for(unsigned it(0); it < NTOY; it++){
			double toy_ele = h_toymc_fakerate[it]->Eval(phoEt,nVertex,fabs(phoEta));
			toy_MET[it]->Fill(sigMET, toy_ele);
		}

		/** cut flow *****/
		if(sigMET < lowMET)continue;
		if(highMET > 0 && sigMET > highMET)continue;

		w_ele = h_nominal_fakerate(phoEt, nVertex, fabs(phoEta));
		p_PhoEt->Fill(phoEt,w_ele);
		p_PhoEta->Fill(phoEta, w_ele);
		p_Mt->Fill(sigMT, w_ele);
		p_HT->Fill(HT, w_ele);
		p_LepPt->Fill(lepPt, w_ele);
		p_LepEta->Fill(lepEta, w_ele);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), w_ele);
		p_nJet->Fill(nJet, w_ele);

		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET, HT, phoEt); 
		if(SigBinIndex >=0){
			h_elefakepho_norm->Fill( SigBinIndex, w_ele);
			h_elefakepho_controlsample->Fill( SigBinIndex );
			if(phoEt < 200){
				predict_lowEt += w_ele;
				proxy_lowEt   += 1;
			}
			else{
				predict_highEt += w_ele;
				proxy_highEt   += 1;
			}
		}	

		for(unsigned it(0); it < NTOY; it++){
			double toy_ele = h_toymc_fakerate[it]->Eval(phoEt,nVertex,fabs(phoEta));
			toy_PhoEt[it]->Fill(phoEt,toy_ele);
			toy_Mt[it]->Fill(sigMT, toy_ele);
			toy_HT[it]->Fill(HT, toy_ele);
			toy_LepPt[it]->Fill(lepPt, toy_ele);
			toy_dPhiEleMET[it]->Fill(fabs(dPhiLepMET), toy_ele);
			if( SigBinIndex >= 0){
				toy_eventcount[it]->Fill(SigBinIndex, toy_ele); 
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


