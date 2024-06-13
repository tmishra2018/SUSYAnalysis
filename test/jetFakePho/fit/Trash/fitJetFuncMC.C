#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

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

#define NTOY 10
#define NBIN 18
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
	//return (jetfakes_highedge + jetfakes_lowedge)/2.0*REBINSIZE;
	return (jetfakes_highedge - jetfakes_lowedge);
}

void fitJetFuncMC(int detType,int  RunYear,bool ISpreVFP){
	int channel = 5; // 1 = eg; 2 = mg; 3 = egloose; 4 = mgloose; 5 = MCEG;

	setTDRStyle();
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetErrorX(0.5);
	gStyle->SetTitleX(0.5);
	gROOT->SetBatch(kTRUE);

	std::string whichVFP;
        if(RunYear==2016 and ISpreVFP == true) whichVFP = "preVFP";
        if(RunYear==2016 and ISpreVFP == false) whichVFP = "postVFP";
        if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

	TChain *sigtree = new TChain("signalTree");
	if(channel == 1)
		sigtree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s.root",RunYear,whichVFP.c_str()));
	else if(channel ==2)
		sigtree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s.root",RunYear,whichVFP.c_str()));
	else if(channel ==5){ 
		sigtree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_DYJetsToLL_%d%s.root",RunYear,whichVFP.c_str()));
		sigtree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_WGToLNuG_%d%s.root",RunYear,whichVFP.c_str()));}

	TChain *controltree = new TChain("jetTree");
	if(channel == 1)
		controltree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s.root",RunYear,whichVFP.c_str()));
	else if(channel ==2)
		controltree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s.root",RunYear,whichVFP.c_str()));
	else if(channel ==5){
		controltree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_DYJetsToLL_%d%s.root",RunYear,whichVFP.c_str()));
		controltree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_WGToLNuG_%d%s.root",RunYear,whichVFP.c_str()));}


	std::stringstream fakerate_filename;
	fakerate_filename.str("");
	if(channel == 1) fakerate_filename << "/eos/uscms/store/user/tmishra/jetfakepho/txt"<<RunYear<<whichVFP<<"/JetFakeRate-DoubleEG-MC";  // some issue here
        if(channel == 2) fakerate_filename << "/eos/uscms/store/user/tmishra/jetfakepho/txt"<<RunYear<<whichVFP<<"/JetFakeRate-MuonEG-MC";
	if(detType == 1)fakerate_filename << "EB.txt";
	else if(detType == 2)fakerate_filename << "EE.txt";

	std::ifstream jetfake_file(fakerate_filename.str().c_str());
	double PtBin[NBIN];
	double fracHad[NBIN];
	double fracHadError[NBIN];
	int i(0);
	float pt_lower(0), pt_upper(0);
	float fakerate(0), error(0), systematic(0);
	float truefake; 
	if(jetfake_file.is_open()){
		for(int i(0); i < NBIN; i++){ 
			jetfake_file >> pt_lower >> pt_upper >> fakerate >> error >> systematic >> truefake;
			PtBin[i] = pt_lower;
			fracHad[i]= fakerate;
			fracHadError[i] = error;
		}
		jetfake_file.close(); 
	}

	std::ostringstream elefake_config;
	elefake_config.str("");
	if(detType == 1) elefake_config << "/eos/uscms/store/user/tmishra/elefakepho/DataResult"<<RunYear<<whichVFP<<"/EleFakeRate-Data-ByPtVtx-EB_" << RunYear <<whichVFP<<".txt";
	else if(detType == 2) elefake_config << "/eos/uscms/store/user/tmishra/elefakepho/DataResult"<<RunYear<<whichVFP<<"/EleFakeRate-Data-ByPtVtx-EE_" << RunYear <<whichVFP<<".txt";
	std::ifstream elefake_file(elefake_config.str().c_str());
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
	/**********************************/
	/*  double normfactor = par[0];   */  
	/*  double slope = par[1];        */
	/*  double constant = par[2];     */
	/*  double index = par[3];        */
	/*  double coeff = par[4];        */
	/*  double vtx_constant = par[5]; */
	/*  double vtx_slope = par[6];    */
	/**********************************/
	TF3 f3("f3", fakerate_func,10,1000,0,100,0,2.5,7);
	f3.SetParameters(scalefactor, ptslope, ptconstant, ptindex, 1.0, vtxconst, vtxslope);


	TH1D *p_controlPhoEt = new TH1D("p_controlPhoEt",";p_{T} (GeV);Events",10,35,135);
	TH1D *p_sigPhoEt  = new TH1D("p_sigPhoEt",";p_{T} (GeV);Events",10,35,135);
	TH1D *p_fakesPhoEt = new TH1D("p_fakesPhoEt",";p_{T} (GeV);Events",10,35,135);
	TH1D *p_fakesPhoEt_stat = new TH1D("p_fakesPhoEt_stat",";p_{T} (GeV);Events",10,35,135);
	TH1D *p_elebkgPhoEt = new TH1D("p_elebkgPhoEt",";p_{T} (GeV);Events",10,35,135);

	if(channel <= 4){
		if(detType == 1){
			sigtree->Draw("phoEt >> p_sigPhoEt", " phoEt >35 && sigMET < 70 && fabs(phoEta) < 1.4442");
			std::cout << "signal " << p_sigPhoEt->GetEntries() << std::endl;
			controltree->Draw("phoEt >> p_controlPhoEt", "phoEt > 35 && sigMET < 70 && fabs(phoEta) < 1.4442");
		}
		if(detType == 2){
			sigtree->Draw("phoEt >> p_sigPhoEt", " phoEt >35 && sigMET < 70 && fabs(phoEta) > 1.56 && fabs(phoEta) < 2.1");
			std::cout << "signal " << p_sigPhoEt->GetEntries() << std::endl;
			controltree->Draw("phoEt >> p_controlPhoEt", "phoEt > 35 && sigMET < 70 && fabs(phoEta) > 1.56 && fabs(phoEta) < 2.1");
		}
	}
	else{
		float crosssection(1);
	  float ntotalevent(1);
		float phoEt(0);
		float phoEta(0);
		float sigMET(0);
		sigtree->SetBranchAddress("crosssection",&crosssection);
		sigtree->SetBranchAddress("ntotalevent", &ntotalevent);
		sigtree->SetBranchAddress("phoEt",       &phoEt);
		sigtree->SetBranchAddress("phoEta",      &phoEta);
		sigtree->SetBranchAddress("sigMET",      &sigMET);
		for(unsigned ievt(0); ievt < sigtree->GetEntries(); ievt++){
			sigtree->GetEntry(ievt);
			if( phoEt < 35 || fabs(phoEta) > 1.4442 || sigMET > 70)continue;
			p_sigPhoEt->Fill( phoEt, crosssection/ntotalevent);
		}	
		float control_crosssection(1);
	  float control_ntotalevent(1);
		float control_phoEt(0);
		float control_phoEta(0);
		float control_sigMET(0);
		controltree->SetBranchAddress("crosssection",&control_crosssection);
		controltree->SetBranchAddress("ntotalevent", &control_ntotalevent);
		controltree->SetBranchAddress("phoEt",       &control_phoEt);
		controltree->SetBranchAddress("phoEta",      &control_phoEta);
		controltree->SetBranchAddress("sigMET",      &control_sigMET);
		for(unsigned ievt(0); ievt < controltree->GetEntries(); ievt++){
			controltree->GetEntry(ievt);
			if( control_phoEt < 35 || fabs(control_phoEta) > 1.4442 || control_sigMET > 70)continue;
			p_controlPhoEt->Fill( control_phoEt, control_crosssection/control_ntotalevent);
		}
	}
  
	//p_sigPhoEt->Sumw2();
	//p_controlPhoEt->Sumw2();


	//************ Proxy Tree **********************//
	TChain *proxytree = new TChain("proxyTree");
	if(channel == 1)proxytree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s.root",RunYear,whichVFP.c_str()));
	else if(channel ==2)proxytree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s.root",RunYear,whichVFP.c_str()));
	else if(channel ==5){
			      proxytree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_DYJetsToLL_%d%s.root",RunYear,whichVFP.c_str()));
			      proxytree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_WGToLNuG_%d%s.root",RunYear,whichVFP.c_str()));
			    }

	float proxy_crosssection(1);
	float proxy_ntotalevent(1);
	float proxyphoEt(0);
	float proxyphoEta(0);
	float proxylepPt(0);
	float proxysigMET(0);
	int   proxynVertex(0);
 
	if(channel ==5){ 
		proxytree->SetBranchAddress("crosssection",&proxy_crosssection);
		proxytree->SetBranchAddress("ntotalevent", &proxy_ntotalevent);
	}
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
		double w_ele = 1;
		w_ele = f3(proxyphoEt, proxynVertex, fabs(proxyphoEta)); 
		if(channel ==5){w_ele = w_ele*proxy_crosssection/proxy_ntotalevent;}
		p_elebkgPhoEt->Fill(proxyphoEt,w_ele);
	}
	p_elebkgPhoEt->Sumw2();
	p_sigPhoEt->Add(p_elebkgPhoEt, -1);

	p_sigPhoEt->Sumw2();
	for(int ibin(1); ibin < p_sigPhoEt->GetSize()-1; ibin++){
		double xvalue = p_sigPhoEt->GetBinCenter(ibin);
		double frac(0),fracerror(0);
		for(unsigned i(0); i< NBIN-1; i++) 
			if(xvalue >= PtBin[i] && xvalue < PtBin[i+1]){
				frac = fracHad[i];
				fracerror = fracHadError[i];
			}
		if(xvalue >= PtBin[NBIN-1]){ frac = fracHad[NBIN-1]; fracerror = fracHadError[NBIN-1];}
    
		double binvalue = p_sigPhoEt->GetBinContent(ibin)*frac;
		if(binvalue == 0)continue;
		double totalerror = p_sigPhoEt->GetBinContent(ibin)*fracerror;
		std::cout << "fake " << p_fakesPhoEt->GetBinCenter(ibin) << " fakerate=" << frac << " fakeraterror " << fracerror << std::endl;
		p_fakesPhoEt_stat->SetBinContent(ibin, p_sigPhoEt->GetBinContent(ibin)); 
		p_fakesPhoEt->SetBinContent(ibin, binvalue);
		p_fakesPhoEt->SetBinError(ibin, totalerror);
	}


	TCanvas *c_pt = new TCanvas("Photon_Pt", "Photon P_{T}",800,800);
	c_pt->cd();
	TPad *can_pad1 = new TPad("can_pad1", "pad1", 0, 0.05, 1, 1.0);
	can_pad1->SetBottomMargin(0.1);
	can_pad1->Draw();          
	can_pad1->cd();          

	TGraphAsymmErrors *new_controlPhoEt = new TGraphAsymmErrors(10);
	TGraphAsymmErrors *new_fakesPhoEt   = new TGraphAsymmErrors(10);
	for(unsigned ibin(1); ibin <= 10; ibin++){
		double new_control_value = p_controlPhoEt->GetBinContent(ibin); 
		double new_fakes_value = p_fakesPhoEt->GetBinContent(ibin);
		double new_control_error = p_controlPhoEt->GetBinError(ibin); 
		double new_fakes_error = p_fakesPhoEt->GetBinError(ibin);
		double staterror = sqrt(p_fakesPhoEt_stat->GetBinContent(ibin))/p_fakesPhoEt_stat->GetBinContent(ibin)*new_fakes_value;
		new_fakes_error = sqrt(new_fakes_error*new_fakes_error + staterror);
		new_controlPhoEt->SetPoint(ibin,p_controlPhoEt->GetBinCenter(ibin),  new_control_value);
		new_controlPhoEt->SetPointError(ibin,0, 0, min(new_control_error, new_control_value-0.00001), new_control_error);
		new_fakesPhoEt->SetPoint(ibin, p_fakesPhoEt->GetBinCenter(ibin),  new_fakes_value);
		new_fakesPhoEt->SetPointError(ibin,0, 0, min(new_fakes_error, new_fakes_value - 0.00001), new_fakes_error);
	}
	
	new_fakesPhoEt->GetXaxis()->SetTitle("p_{T} (GeV)");
	new_controlPhoEt->GetXaxis()->SetTitle("p_{T} (GeV)");
	new_controlPhoEt->GetXaxis()->SetTitleOffset(1);
	new_controlPhoEt->GetXaxis()->SetTitleSize(20);
	new_controlPhoEt->GetXaxis()->SetRangeUser(35,135);
	new_fakesPhoEt->GetXaxis()->SetRangeUser(35,135);
	gPad->SetLogy();
	TH1D *new_dummy=new TH1D("dummy",";p_{T} (GeV); Event/GeV",10,35,135);
	new_dummy->SetMinimum(10);
	new_dummy->SetMaximum(1000000);
	new_dummy->GetXaxis()->SetTitle("p_{T} (GeV)");
	new_dummy->GetXaxis()->SetTitleOffset(1.0);
	new_dummy->Draw();
	new_controlPhoEt->Draw("P same");
	new_controlPhoEt->SetLineColor(kBlack);
	new_controlPhoEt->SetMarkerColor(kBlack);
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
	new_controlPhoEt->Fit("expo");
	TF1 *inifitden = new_controlPhoEt->GetFunction("expo");
        inifitden->SetLineColor(kBlack);
	double iniLambda_den1 = inifitden->GetParameter(1);
	double iniCoeff_den   = exp(inifitden->GetParameter(0))/2;
	double iniLambda_den2 =  iniLambda_den1/2; 
 	TF1 *fitfunc_den= new TF1("fitfunc_den", tmpjetfake_func, 35, 135, 4);
	if(iniLambda_den2 > 0 && iniLambda_den2 < 1e6)fitfunc_den->SetParameters(iniCoeff_den, iniCoeff_den/10, iniLambda_den1, iniLambda_den2);
	else fitfunc_den->SetParameters(iniCoeff_den, iniCoeff_den/10, iniLambda_den1, iniLambda_den1);
	new_controlPhoEt->Fit("fitfunc_den","S");
	TF1 *fitden = new_controlPhoEt->GetFunction("fitfunc_den");
        fitden->SetLineColor(kBlack);
	ofstream myfile;
	if(detType == 1)myfile.open(Form("/eos/uscms/store/user/tmishra/jetfakepho/txt%d%s/JetFakeRate-transferfactor-MCEG-EB.txt",RunYear,whichVFP.c_str()), std::ios_base::app | std::ios_base::out);
        if(detType == 2)myfile.open(Form("/eos/uscms/store/user/tmishra/jetfakepho/txt%d%s/JetFakeRate-transferfactor-MCEG-EE.txt",RunYear,whichVFP.c_str()), std::ios_base::app | std::ios_base::out);

	TH1D *ratio = new TH1D("transfer fraction","",int(p_fakesPhoEt->GetXaxis()->GetNbins()/REBINSIZE),35,135);
	for(int ibin(1); ibin < int(p_fakesPhoEt->GetXaxis()->GetNbins()/REBINSIZE); ibin++){
		ratio->SetBinContent(ibin, new_fakesPhoEt->Eval(35+(ibin-1)*2 +1)/new_controlPhoEt->Eval(35+(ibin-1)*2 +1));
	}
 
	myfile << "den_coeff1 " << fitden->GetParameter(0) << std::endl;
	myfile << "den_coeff2 " << fitden->GetParameter(1) << std::endl;
	myfile << "den_lambd1 " << fitden->GetParameter(2) << std::endl;
	myfile << "den_lambd2 " << fitden->GetParameter(3) << std::endl;
   
	TFitResultPtr rden = new_controlPhoEt->Fit("fitfunc_den","S");
	TMatrixDSym covden = rden->GetCovarianceMatrix(); 
	rden->Print("V");     
	//fitfunc_den->Draw("same");
	TVectorD muden(4) ;
	muden(0) = rden->Parameter(0); 
	muden(1) = rden->Parameter(1);
	muden(2) = rden->Parameter(2);
	muden(3) = rden->Parameter(3);
	RooRealVar central_coeff1_den("central_coeff1_den","central_coeff1_den",muden(0)-rden->ParError(0), muden(0)+rden->ParError(0));
	RooRealVar central_coeff2_den("central_coeff2_den","central_coeff2_den",muden(1)-rden->ParError(1), muden(1)+rden->ParError(1));
	RooRealVar central_lambda1_den("central_lambda1_den","central_lambda1_den",muden(2)-rden->ParError(2),muden(2)+rden->ParError(2));
	RooRealVar central_lambda2_den("central_lambda2_den","central_lambda2_den",muden(3)-rden->ParError(3),muden(3)+rden->ParError(3));
	
	RooMultiVarGaussian mvgden("mvgden","mvgden",RooArgList(central_coeff1_den,central_coeff2_den,central_lambda1_den,central_lambda2_den),muden,covden);
	RooDataSet* toymcdataden = mvgden.generate(RooArgSet(central_coeff1_den,central_coeff2_den,central_lambda1_den,central_lambda2_den),NTOY);
	std::ostringstream modelnameden;
	TF1 *gen_den[NTOY];
	TH1D *den_upper = new TH1D("den_upper","den_upper",10,35,135);
	TH1D *den_lower = new TH1D("den_lower","den_lower",10,35,135);
	can_pad1->cd();          
	for(unsigned ibin(1); ibin <= 10; ibin++)den_lower->SetBinContent(ibin, fitfunc_den->Eval(den_upper->GetBinCenter(ibin)));
	for(int i(0); i<NTOY; i++){
		double data1 = toymcdataden->get(i)->getRealValue("central_coeff1_den");
		double data2 = toymcdataden->get(i)->getRealValue("central_coeff2_den");
		double data3 = toymcdataden->get(i)->getRealValue("central_lambda1_den");
		double data4 = toymcdataden->get(i)->getRealValue("central_lambda2_den");
		modelnameden.str("");
		modelnameden << "gen_den_" << i;
		gen_den[i] = new TF1(modelnameden.str().c_str(), tmpjetfake_func, 35, 135, 4);
		gen_den[i]->SetParameters(data1, data2, data3, data4);
		gen_den[i]->SetLineColorAlpha(kBlue, 0.35);
		//gen_den[i]->Draw("same");
	 	for(unsigned ibin(1); ibin <= 10; ibin++){
			double estimated = gen_den[i]->Eval(den_upper->GetBinCenter(ibin));
			if(den_upper->GetBinContent(ibin) < estimated)den_upper->SetBinContent(ibin, estimated);
			if(den_lower->GetBinContent(ibin) > estimated)den_lower->SetBinContent(ibin, estimated);
		}
	}
	
//	den_upper->Draw("L same");
//	den_lower->Draw("L same");
// *************************  Numerator ******************************************************************//
 
new_fakesPhoEt->Fit("expo");
TF1 *inifit = new_fakesPhoEt->GetFunction("expo");
   inifit->SetLineColor(kRed);
double iniLambda_num1 = inifit->GetParameter(1);
double iniCoeff_num  	= exp(inifit->GetParameter(0))/2;
double iniLambda_num2 = iniLambda_num1/2; 
TF1 *fitfunc_num= new TF1("fitfunc_num", tmpjetfake_func, 35, 135, 4);
fitfunc_num->SetParameters(iniCoeff_num, 20, iniLambda_num1, -0.02);
TVirtualFitter::SetMaxIterations(1000000);
TFitResultPtr r = new_fakesPhoEt->Fit("fitfunc_num","R S");
TF1 *fit = new_fakesPhoEt->GetFunction("fitfunc_num");
        fit->SetLineColor(kRed);
myfile << "num_coeff1 " << fit->GetParameter(0) << std::endl;
myfile << "num_coeff2 " << fit->GetParameter(1) << std::endl;
myfile << "num_lambd1 " << fit->GetParameter(2) << std::endl;
myfile << "num_lambd2 " << fit->GetParameter(3) << std::endl;
 
//std::ostringstream testname;
//TF1 *test_num[200][200];
//for(unsigned i(0); i < 200; i++){
//	for(unsigned j(0); j < 200; j++){
//		testname.str("");
//		testname << "test_" << i << "_" << j;
//		test_num[i][j] = new TF1(testname.str().c_str(), tmpjetfake_func, 35, 135, 4);
//		test_num[i][j]->SetParameters(iniCoeff_num, 20+i*5, iniLambda_num1, -0.001-0.001*j);
//		//TVirtualFitter::SetMaxIterations(1000000);
//	  int status = new_fakesPhoEt->Fit(testname.str().c_str(),"R");
//		TF1 *tmpf1 = new_fakesPhoEt->GetFunction(testname.str().c_str());
//		float diff = fabs(tmpf1->Eval(200) - new_fakesPhoEt->Eval(200))/new_fakesPhoEt->Eval(200);
//		if(tmpf1->GetParError(0)/tmpf1->GetParameter(0) < 1 && tmpf1->GetParError(1)/tmpf1->GetParameter(1) < 1 && tmpf1->GetParError(2)/tmpf1->GetParameter(2) < 1 && tmpf1->GetParError(3)/tmpf1->GetParameter(3) < 1	)std::cout << "good point " << i*5 << " " << -1+0.005*j << " status = " << status << " diff " << diff <<  std::endl;
//		delete tmpf1;
//	}
//}
 
	can_pad1->cd();          
	TMatrixDSym cov = r->GetCovarianceMatrix(); 
	r->Print("V");     
	//fitfunc_num->Draw("same");
	float nominalvalue_num[10];
	for(unsigned ibin(0); ibin < 10; ibin++){
		nominalvalue_num[ibin] =  fitfunc_num->Eval(35+ibin);
	}
	TVectorD mu(4) ;
	mu(0) = r->Parameter(0); 
	mu(1) = r->Parameter(1);
	mu(2) = r->Parameter(2);
	mu(3) = r->Parameter(3);
	RooRealVar central_coeff1_num("central_coeff1_num","central_coeff1_num",mu(0)-r->ParError(0), mu(0)+r->ParError(0));
	RooRealVar central_coeff2_num("central_coeff2_num","central_coeff2_num",mu(1)-r->ParError(1), mu(1)+r->ParError(1));
	RooRealVar central_lambda1_num("central_lambda1_num","central_lambda1_num",mu(2)-r->ParError(2),mu(2)+r->ParError(2));
	RooRealVar central_lambda2_num("central_lambda2_num","central_lambda2_num",mu(3)-r->ParError(3),mu(3)+r->ParError(3));
	
	RooMultiVarGaussian mvg("mvg","mvg",RooArgList(central_coeff1_num,central_coeff2_num,central_lambda1_num,central_lambda2_num),mu,cov);
	RooDataSet* toymcdata = mvg.generate(RooArgSet(central_coeff1_num,central_coeff2_num,central_lambda1_num,central_lambda2_num),NTOY);
	std::ostringstream modelname;
	TF1 *gen_num[NTOY];
	TH1D *num_upper = new TH1D("num_upper","num_upper",10,35,135);
	TH1D *num_lower = new TH1D("num_lower","num_lower",10,35,135);
	float toyptvalue[10][NTOY];
	float lowtoyptvalue[10];
	float hightoyptvalue[10];
	for(unsigned ii(0); ii < 10; ii++){
		lowtoyptvalue[ii] = 100000;
		hightoyptvalue[ii] = 0;
	}
	
	for(int i(0); i<NTOY; i++){
		double data1 = toymcdata->get(i)->getRealValue("central_coeff1_num");
		double data2 = toymcdata->get(i)->getRealValue("central_coeff2_num");
		double data3 = toymcdata->get(i)->getRealValue("central_lambda1_num");
		double data4 = toymcdata->get(i)->getRealValue("central_lambda2_num");
		if(data1 < 1e8 && data2 < 1e8 && data3 < 1e8 && data4 < 1e8){
			modelname.str("");
			modelname << "gen_num_" << i;
			gen_num[i] = new TF1(modelname.str().c_str(), tmpjetfake_func, 35, 135, 4);
			gen_num[i]->SetParameters(data1, data2, data3, data4);
			gen_num[i]->SetLineColorAlpha(kBlue, 0.35);
			//gen_num[i]->Draw("same");
			bool exception(false);
			//std::cout << std::endl;
			//for(unsigned ibin(1); ibin < new_fakesPhoEt->GetSize(); ibin++){
			//	std::cout << gen_num[i]->Eval(new_fakesPhoEt->GetBinCenter(ibin)) << " " << new_fakesPhoEt->GetBinContent(ibin) << std::endl; 
			//	if(gen_num[i]->Eval(new_fakesPhoEt->GetBinCenter(ibin)) - new_fakesPhoEt->GetBinContent(ibin) < -1*1.5*new_fakesPhoEt->GetBinError(ibin))exception=true;
			//}
			if(exception)continue; 
			for(unsigned ibin(0); ibin < 10; ibin++){
				double estimated = gen_num[i]->Eval(35+ibin);
				toyptvalue[ibin][i] = estimated;
				if(lowtoyptvalue[ibin] > estimated)lowtoyptvalue[ibin] = estimated;
				if(hightoyptvalue[ibin]< estimated)hightoyptvalue[ibin] = estimated;
			}
		}
	}
	
	TCanvas *cangaus = new TCanvas("cangaus","",600,600);
	cangaus->cd();
	for(unsigned ibin(0); ibin < 10; ibin++){
		num_upper->SetBinContent(ibin, hightoyptvalue[ibin]); 
		num_lower->SetBinContent(ibin, lowtoyptvalue[ibin]);
	}
	
	can_pad1->cd();          
//	num_upper->Draw("L same");
//	num_lower->Draw("L same");
	new_controlPhoEt->Draw("P same");
	new_fakesPhoEt->Draw("EP same");
	
	
	for(unsigned ii(1); ii <= 10; ii++){
		myfile << "numerror " << ii-1 << " " <<  num_upper->GetBinContent(ii) - num_lower->GetBinContent(ii) << std::endl;
	} 
	for(unsigned ii(1); ii <= 10; ii++){
		myfile << "denerror " << ii-1 << " " <<  den_upper->GetBinContent(ii) - den_lower->GetBinContent(ii) << std::endl;
	}
	
	std::cout << std::endl; 
	for(int ibin(1); ibin < ratio->GetSize(); ibin++){
		myfile << "ratio " << ibin << " = " << ratio->GetBinContent(ibin) << "  est=" << fitfunc_num->Eval(ratio->GetBinCenter(ibin))/fitfunc_den->Eval(ratio->GetBinCenter(ibin)) << std::endl;
	}
	myfile.close();
	
	if(channel == 1 and detType == 1)c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Plots%d%s/JetFakeRate_transfer_MCEG_EB.png",RunYear,whichVFP.c_str()));
        else if(channel == 1 and detType == 2)c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Plots%d%s/JetFakeRate_transfer_MCEG_EE.png",RunYear,whichVFP.c_str()));
        else if(channel == 2 and detType == 1)c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Plots%d%s/JetFakeRate_transfer_MCMG_EB.png",RunYear,whichVFP.c_str()));
        else if(channel == 2 and detType == 2)c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Plots%d%s/JetFakeRate_transfer_MCMG_EE.png",RunYear,whichVFP.c_str()));
}
