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
#include "TROOT.h"

#define NTOY 1000
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
// root -l -q "fitJetFunc.C(1,2016,1)"
void fitJetFunc_oldFunc(int ichannel,int  RunYear,bool ISpreVFP){
	int detType = 1;
	int channel = ichannel;
	gROOT->SetBatch(kTRUE);
	setTDRStyle();
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetErrorX(0.5);
	gStyle->SetTitleX(0.5);
	
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
	if(channel == 1) fakerate_filename << "/eos/uscms/store/user/tmishra/jetfakepho/txt"<<RunYear<<whichVFP<<"/JetFakeRate-DoubleEG-";  // some issue here
        if(channel == 2) fakerate_filename << "/eos/uscms/store/user/tmishra/jetfakepho/txt"<<RunYear<<whichVFP<<"/JetFakeRate-MuonEG-";

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

	std::ostringstream elefake_config;
	elefake_config.str("");
	if(detType == 1) elefake_config << "/eos/uscms/store/user/tmishra/elefakepho/DataResult"<<RunYear<<whichVFP<<"/EleFakeRate-Data-ByPtVtx-EB.txt";
	else if(detType == 2) elefake_config << "/eos/uscms/store/user/tmishra/elefakepho/DataResult"<<RunYear<<whichVFP<<"/EleFakeRate-Data-ByPtVtx-EE.txt";
	
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
	// ele fake file as input
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
        else cout<<"file does not exist..  "<<"\n";
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


	TH1D *p_controlPhoEt = new TH1D("p_controlPhoEt",";p_{T} (GeV);Events",265,35,300);
	TH1D *p_sigPhoEt  = new TH1D("p_sigPhoEt",";p_{T} (GeV);Events",265,35,300);
	TH1D *p_fakesPhoEt = new TH1D("p_fakesPhoEt",";p_{T} (GeV);Events",265,35,300);
	TH1D *p_fakesPhoEt_stat = new TH1D("p_fakesPhoEt_stat",";p_{T} (GeV);Events",265,35,300);
	TH1D *p_elebkgPhoEt = new TH1D("p_elebkgPhoEt",";p_{T} (GeV);Events",265,35,300);

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
		double w_ele = 1;
		w_ele = f3(proxyphoEt, proxynVertex, fabs(proxyphoEta));
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
	TPad *can_pad1 = new TPad("can_pad1", "pad1", 0, 0.05, 1, 1.0);
	can_pad1->SetBottomMargin(0.1);
	can_pad1->Draw();          
	can_pad1->cd();          

	TGraphAsymmErrors *new_controlPhoEt = new TGraphAsymmErrors(95);
	TGraphAsymmErrors *new_fakesPhoEt   = new TGraphAsymmErrors(95);
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
	// bins 76-94, bin width =  10, pt 110-300
	for(unsigned ibin(76); ibin <= 94; ibin++){
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
		double new_control_value = p_controlPhoEt->GetBinContent(266); 
		double new_fakes_value = p_fakesPhoEt->GetBinContent(266);
		double new_control_error = p_controlPhoEt->GetBinError(266); 
		double new_fakes_error = p_fakesPhoEt->GetBinError(266);
		new_controlPhoEt->SetPoint(95,300,  new_control_value);
		new_controlPhoEt->SetPointError(95,5, 5, min(new_control_error, new_control_value-0.001), new_control_error);
		new_fakesPhoEt->SetPoint(95, 300,  new_fakes_value);
		new_fakesPhoEt->SetPointError(95,5, 5, min(new_fakes_error, new_fakes_value - 0.001), new_fakes_error);
	}
	
	new_fakesPhoEt->GetXaxis()->SetTitle("p_{T} (GeV)");
	new_controlPhoEt->GetXaxis()->SetTitle("p_{T} (GeV)");
	new_controlPhoEt->GetXaxis()->SetTitleOffset(1);
	new_controlPhoEt->GetXaxis()->SetTitleSize(20);
	new_controlPhoEt->GetXaxis()->SetRangeUser(35,300);
	new_fakesPhoEt->GetXaxis()->SetRangeUser(35,300);
	gPad->SetLogy();
	TH1D *new_dummy=new TH1D("dummy",";p_{T} (GeV); Event/GeV",265,35,300);
	new_dummy->SetMinimum(0.001);
	new_dummy->SetMaximum(100000);
	new_dummy->GetXaxis()->SetTitle("p_{T} (GeV)");
	new_dummy->GetXaxis()->SetTitleOffset(1.0);
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
	new_controlPhoEt->Fit("expo");
	TF1 *inifitden = new_controlPhoEt->GetFunction("expo");
	double iniLambda_den1 = inifitden->GetParameter(1);
	double iniCoeff_den   = exp(inifitden->GetParameter(0))/2;
	double iniLambda_den2 = (log(new_controlPhoEt->Eval(111)- 2*iniCoeff_den*exp(iniLambda_den1*111.0)) - log(new_controlPhoEt->Eval(250)- 2*iniCoeff_den*exp(iniLambda_den1*250.0)) )/( 111.0-250.0);

	cout<<"iniLambda_den2  : "<<iniLambda_den2<<endl; 
	TF1 *fitfunc_den= new TF1("fitfunc_den", tmpjetfake_func, 35, 300, 4);
	fitfunc_den->SetParameters(iniCoeff_den, iniCoeff_den/10, iniLambda_den1, iniLambda_den1);
	new_controlPhoEt->Fit("fitfunc_den","S");
	TF1 *fitden = new_controlPhoEt->GetFunction("fitfunc_den");

	ofstream myfile;
	if(channel == 1 and detType == 1)myfile.open(Form("/eos/uscms/store/user/tmishra/jetfakepho/txt%d%s/JetFakeRate-transferfactor-DoubleEG-EB.txt",RunYear,whichVFP.c_str()), std::ios_base::ate | std::ios_base::out);
	if(channel == 1 and detType == 2)myfile.open(Form("/eos/uscms/store/user/tmishra/jetfakepho/txt%d%s/JetFakeRate-transferfactor-DoubleEG-EE.txt",RunYear,whichVFP.c_str()), std::ios_base::ate | std::ios_base::out);

	if(channel == 2 and detType == 1)myfile.open(Form("/eos/uscms/store/user/tmishra/jetfakepho/txt%d%s/JetFakeRate-transferfactor-MuonEG-EB_oldFunc.txt",RunYear,whichVFP.c_str()), std::ios_base::ate | std::ios_base::out);
	if(channel == 2 and detType == 2)myfile.open(Form("/eos/uscms/store/user/tmishra/jetfakepho/txt%d%s/JetFakeRate-transferfactor-MuonEG-EE.txt",RunYear,whichVFP.c_str()), std::ios_base::ate | std::ios_base::out);
	
	TH1D *ratio = new TH1D("transfer fraction","",int(p_fakesPhoEt->GetXaxis()->GetNbins()/REBINSIZE),35,300);
	for(int ibin(1); ibin < int(p_fakesPhoEt->GetXaxis()->GetNbins()/REBINSIZE); ibin++){
		ratio->SetBinContent(ibin, new_fakesPhoEt->Eval(35+(ibin-1)*2 +1)/new_controlPhoEt->Eval(35+(ibin-1)*2 +1));
	}
 
	myfile << "den_coeff1 " << fitden->GetParameter(0) << std::endl;
	myfile << "den_coeff2 " << fitden->GetParameter(1) << std::endl;
	myfile << "den_lambd1 " << fitden->GetParameter(2) << std::endl;
	myfile << "den_lambd2 " << fitden->GetParameter(3) << std::endl;

	cout << "den_coeff1 " << fitden->GetParameter(0) << std::endl;
	cout << "den_coeff2 " << fitden->GetParameter(1) << std::endl;
	cout << "den_lambd1 " << fitden->GetParameter(2) << std::endl;
	cout << "den_lambd2 " << fitden->GetParameter(3) << std::endl;
	cout << "chiSquare/ndf:  " << fitden->GetChisquare()/fitden->GetNDF() << std::endl;
	cout << "chiSquare:  " << fitden->GetChisquare() << std::endl;

	TFitResultPtr rden = new_controlPhoEt->Fit("fitfunc_den","S");
	// covariance matrix
	TMatrixDSym covden = rden->GetCovarianceMatrix(); 
	rden->Print("V");
	fitfunc_den->Draw("same");
	TVectorD muden(4) ;
	muden(0) = rden->Parameter(0); 
	muden(1) = rden->Parameter(1);
	muden(2) = rden->Parameter(2);
	muden(3) = rden->Parameter(3);
	RooRealVar central_coeff1_den("central_coeff1_den","central_coeff1_den",muden(0)-rden->ParError(0), muden(0)+rden->ParError(0));
	RooRealVar central_coeff2_den("central_coeff2_den","central_coeff2_den",muden(1)-rden->ParError(1), muden(1)+rden->ParError(1));
	RooRealVar central_lambda1_den("central_lambda1_den","central_lambda1_den",muden(2)-rden->ParError(2),muden(2)+rden->ParError(2));
	RooRealVar central_lambda2_den("central_lambda2_den","central_lambda2_den",muden(3)-rden->ParError(3),muden(3)+rden->ParError(3));
	// Multivariant Gaussian pdf, get NTOY numbers of paramters
	RooMultiVarGaussian mvgden("mvgden","mvgden",RooArgList(central_coeff1_den,central_coeff2_den,central_lambda1_den,central_lambda2_den),muden,covden);
	RooDataSet* toymcdataden = mvgden.generate(RooArgSet(central_coeff1_den,central_coeff2_den,central_lambda1_den,central_lambda2_den),NTOY);
	std::ostringstream modelnameden;
	TF1 *gen_den[NTOY];
	TH1D *den_upper = new TH1D("den_upper","den_upper",265,35,300);
	TH1D *den_lower = new TH1D("den_lower","den_lower",265,35,300);
	can_pad1->cd();          
	for(unsigned ibin(1); ibin <= 265; ibin++)den_lower->SetBinContent(ibin, fitfunc_den->Eval(den_upper->GetBinCenter(ibin)));
	for(int i(0); i<NTOY; i++){
		double data1 = toymcdataden->get(i)->getRealValue("central_coeff1_den");
		double data2 = toymcdataden->get(i)->getRealValue("central_coeff2_den");
		double data3 = toymcdataden->get(i)->getRealValue("central_lambda1_den");
		double data4 = toymcdataden->get(i)->getRealValue("central_lambda2_den");
		modelnameden.str("");
		modelnameden << "gen_den_" << i;
		gen_den[i] = new TF1(modelnameden.str().c_str(), tmpjetfake_func, 35, 300, 4);
		gen_den[i]->SetParameters(data1, data2, data3, data4);
		gen_den[i]->SetLineColorAlpha(kBlue, 0.35);
		//gen_den[i]->Draw("same");
	 	for(unsigned ibin(1); ibin <= 265; ibin++){
			double estimated = gen_den[i]->Eval(den_upper->GetBinCenter(ibin));
			if(den_upper->GetBinContent(ibin) < estimated)den_upper->SetBinContent(ibin, estimated);
			if(den_lower->GetBinContent(ibin) > estimated)den_lower->SetBinContent(ibin, estimated);
		}
	}
	
	den_upper->Draw("L same");
	den_lower->Draw("L same");
// *************************  Numerator ******************************************************************//
 
new_fakesPhoEt->Fit("expo");
TF1 *inifit = new_fakesPhoEt->GetFunction("expo");
double iniLambda_num1 = inifit->GetParameter(1);
double iniCoeff_num  	= exp(inifit->GetParameter(0))/2;
//double iniLambda_num2 = (log(new_fakesPhoEt->GetBinContent(38)- 2*iniCoeff_num*exp(iniLambda_num1*new_fakesPhoEt->GetBinCenter(38))) - log(new_fakesPhoEt->GetBinContent(56)- 2*iniCoeff_num*exp(iniLambda_num1*new_fakesPhoEt->GetBinCenter(56))) )/(new_fakesPhoEt->GetBinCenter(38)-new_fakesPhoEt->GetBinCenter(56));
double iniLambda_num2 = (log(new_fakesPhoEt->Eval(111)- 2*iniCoeff_num*exp(iniLambda_num1*111.0)) - log(new_fakesPhoEt->Eval(200)- 2*iniCoeff_num*exp(iniLambda_num1*200.0)) )/(111.0-200.0);



std::cout << "log(new_fakesPhoEt->Eval(111)- 2*iniCoeff_num*exp(iniLambda_num1*111.0)) = " << log(new_fakesPhoEt->Eval(111)- 2*iniCoeff_num*exp(iniLambda_num1*111.0)) << "  log(new_fakesPhoEt->Eval(200)- 2*iniCoeff_num*exp(iniLambda_num1*200.0)) = " << log(new_fakesPhoEt->Eval(200)- 2*iniCoeff_num*exp(iniLambda_num1*200.0))  << "  lambda2 = " << iniLambda_num2 << std::endl;

TF1 *fitfunc_num= new TF1("fitfunc_num", tmpjetfake_func, 35, 300, 4);
// Eqn 6, 7 AN
cout<<iniCoeff_num<<"\t"<<iniCoeff_num/10<<"\t"<<iniLambda_num1<<"\t"<<iniLambda_num2<<endl;
//if(channel == 1)fitfunc_num->SetParameters(iniCoeff_num, iniCoeff_num/10, iniLambda_num1, iniLambda_num2);
if(channel == 1)fitfunc_num->SetParameters(31140, 2199.48, -0.0862613, -0.0402501);
//if(channel == 1)fitfunc_num->SetParameters(14300, 300, -0.070, -0.04);
//else if(channel == 2)fitfunc_num->SetParameters(iniCoeff_num, iniCoeff_num/10, iniLambda_num1, iniLambda_num2);
else if(channel == 2) fitfunc_num->SetParameters(23900, 300, -0.082, -0.01);

TVirtualFitter::SetMaxIterations(1000000);
TFitResultPtr r = new_fakesPhoEt->Fit("fitfunc_num","R S");
TF1 *fit = new_fakesPhoEt->GetFunction("fitfunc_num");
myfile << "num_coeff1 " << fit->GetParameter(0) << std::endl;
myfile << "num_coeff2 " << fit->GetParameter(1) << std::endl;
myfile << "num_lambd1 " << fit->GetParameter(2) << std::endl;
myfile << "num_lambd2 " << fit->GetParameter(3) << std::endl;

cout << "num_coeff1 " << fit->GetParameter(0) << std::endl;
cout << "num_coeff2 " << fit->GetParameter(1) << std::endl;
cout << "num_lambd1 " << fit->GetParameter(2) << std::endl;
cout << "num_lambd2 " << fit->GetParameter(3) << std::endl;
cout << "chiSquare/ndf:  " << fit->GetChisquare()/fit->GetNDF() << std::endl;
cout << "chiSquare:  " << fit->GetChisquare() << std::endl;

//std::ostringstream testname;
//TF1 *test_num[200][200];
//for(unsigned i(0); i < 200; i++){
//	for(unsigned j(0); j < 200; j++){
//		testname.str("");
//		testname << "test_" << i << "_" << j;
//		test_num[i][j] = new TF1(testname.str().c_str(), tmpjetfake_func, 35, 300, 4);
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
	fitfunc_num->Draw("same");
	float nominalvalue_num[265];
	for(unsigned ibin(0); ibin < 265; ibin++){
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
	TH1D *num_upper = new TH1D("num_upper","num_upper",265,35,300);
	TH1D *num_lower = new TH1D("num_lower","num_lower",265,35,300);
	float toyptvalue[265][NTOY];
	float lowtoyptvalue[265];
	float hightoyptvalue[265];
	for(unsigned ii(0); ii < 265; ii++){
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
			gen_num[i] = new TF1(modelname.str().c_str(), tmpjetfake_func, 35, 300, 4);
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
			for(unsigned ibin(0); ibin < 265; ibin++){
				double estimated = gen_num[i]->Eval(35+ibin);
				toyptvalue[ibin][i] = estimated;
				if(lowtoyptvalue[ibin] > estimated)lowtoyptvalue[ibin] = estimated;
				if(hightoyptvalue[ibin]< estimated)hightoyptvalue[ibin] = estimated;
			}
		}
	}
	
	TCanvas *cangaus = new TCanvas("cangaus","",600,600);
	cangaus->cd();
	for(unsigned ibin(0); ibin < 265; ibin++){
		num_upper->SetBinContent(ibin, hightoyptvalue[ibin]); 
		num_lower->SetBinContent(ibin, lowtoyptvalue[ibin]);
	}
	
	can_pad1->cd();          
	num_upper->Draw("L same");
	num_lower->Draw("L same");
	new_controlPhoEt->Draw("EP same");
	new_fakesPhoEt->Draw("EP same");
	
        if(RunYear==2016 and ISpreVFP == 1)       CMS_lumi( can_pad1, 1, ichannel, 11 );
        else if(RunYear==2016 and ISpreVFP == 0)  CMS_lumi( can_pad1, 2, ichannel, 11 );
        else if(RunYear==2017)                    CMS_lumi( can_pad1, 3, ichannel, 11 );
        else if(RunYear==2018)                    CMS_lumi( can_pad1, 4, ichannel, 11 );
	
	
	
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
	       
	if(channel == 1 and detType == 1)c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Plots%d%s/DoubleEG/JetFakeRate_transfer_DoubleEG_EB_%d%s.png", RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));
	else if(channel == 1 and detType == 2)c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Plots%d%s/DoubleEG/JetFakeRate_transfer_DoubleEG_EE_%d%s.png", RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));
	else if(channel == 2 and detType == 1)c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Plots%d%s/MuonEG/JetFakeRate_transfer_MuonEG_EB_%d%s.png", RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));
	else if(channel == 2 and detType == 2)c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Plots%d%s/MuonEG/JetFakeRate_transfer_MuonEG_EE_%d%s.png", RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));
}
