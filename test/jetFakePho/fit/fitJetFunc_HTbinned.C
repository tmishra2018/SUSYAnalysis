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
		weights.push_back((denVals[i] != 0) ? numVals[matchIndex]/denVals[i] : 0);
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

void fitJetFunc_HTbinned(int ichannel,int  RunYear,bool ISpreVFP){
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
	if(channel == 1) fakerate_filename << "/eos/uscms/store/user/tmishra/jetfakepho/txt"<<RunYear<<whichVFP<<"/JetFakeRate-DoubleEG-";
        if(channel == 2) fakerate_filename << "/eos/uscms/store/user/tmishra/jetfakepho/txt"<<RunYear<<whichVFP<<"/JetFakeRate-MuonEG-";

	if(detType == 1)fakerate_filename << "EB_HTbinned.txt";
	else if(detType == 2)fakerate_filename << "EE_HTbinned.txt";
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

	// ele fake file as input
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
	else if(channel ==2)proxytree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/NEW/resTree_mgsignal_MuonEG_%d%s.root",RunYear,whichVFP.c_str()));

	float proxyphoEt(0);
	float proxyphoEta(0);
	float proxylepPt(0);
	float proxysigMET(0);
	float proxyHT(0);
 
	proxytree->SetBranchAddress("phoEt",     &proxyphoEt);
	proxytree->SetBranchAddress("phoEta",    &proxyphoEta);
	proxytree->SetBranchAddress("lepPt",     &proxylepPt);
	proxytree->SetBranchAddress("sigMET",    &proxysigMET);
	proxytree->SetBranchAddress("HT",        &proxyHT);
 
	for (unsigned ievt(0); ievt<proxytree->GetEntries(); ++ievt){//loop on entries
		proxytree->GetEntry(ievt);
		if(proxysigMET > 70)continue;
		if(detType == 1 && fabs(proxyphoEta) > 1.4442)continue;
		else if(detType == 2 && (fabs(proxyphoEta) < 1.56 || fabs(proxyphoEta) > 2.1))continue;
		if(proxyphoEt < 35)continue;
		double weight_pt = fakerate_pt(proxyphoEt, ptslope, ptconstant, ptindex);
		double weight_eta = getBinnedWeight(fabs(proxyphoEta), etaBins, etaWeights);
		double weight_ht = getBinnedWeight(proxyHT, htBins, htWeights);
		double w_ele = scalefactor*weight_pt*weight_eta*weight_ht;
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
    		fout = new TFile(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/fitJetFunc_Output_DoubleEG_%d%s.root", RunYear, whichVFP.c_str()), "RECREATE");
	else if (channel == 2)
    		fout = new TFile(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/fitJetFunc_Output_MuonEG_%d%s.root", RunYear, whichVFP.c_str()), "RECREATE");

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
	TF1 *fitden = new_controlPhoEt->GetFunction("fitfunc_den");

	ofstream myfile;
	if(channel == 1 and detType == 1)myfile.open(Form("/eos/uscms/store/user/tmishra/jetfakepho/txt%d%s/JetFakeRate-transferfactor-DoubleEG-EB_HTbinned.txt",RunYear,whichVFP.c_str()), std::ios_base::ate | std::ios_base::out);
	if(channel == 1 and detType == 2)myfile.open(Form("/eos/uscms/store/user/tmishra/jetfakepho/txt%d%s/JetFakeRate-transferfactor-DoubleEG-EE_HTbinned.txt",RunYear,whichVFP.c_str()), std::ios_base::ate | std::ios_base::out);

	if(channel == 2 and detType == 1)myfile.open(Form("/eos/uscms/store/user/tmishra/jetfakepho/txt%d%s/JetFakeRate-transferfactor-MuonEG-EB_HTbinned.txt",RunYear,whichVFP.c_str()), std::ios_base::ate | std::ios_base::out);
	if(channel == 2 and detType == 2)myfile.open(Form("/eos/uscms/store/user/tmishra/jetfakepho/txt%d%s/JetFakeRate-transferfactor-MuonEG-EE_HTbinned.txt",RunYear,whichVFP.c_str()), std::ios_base::ate | std::ios_base::out);
	
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
	for(unsigned ibin(1); ibin <= 565; ibin++)den_lower->SetBinContent(ibin, fitfunc_den->Eval(den_upper->GetBinCenter(ibin)));
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
			if(den_upper->GetBinContent(ibin) < estimated)den_upper->SetBinContent(ibin, estimated);
			if(den_lower->GetBinContent(ibin) > estimated)den_lower->SetBinContent(ibin, estimated);
		}
	}
	
	den_upper->Draw("L same");
	den_lower->Draw("L same");
// *************************  Numerator ******************************************************************//

TF1 *fitfunc_num = new TF1("fitfunc_num", exp2c_func, 20, 600, 6);

fitfunc_num->SetParameters(51521, 1.13974, 0.0560658, 0.0200371, 125.136, 0.0970124);
fitfunc_num->SetParLimits(0, 0, 1e9);      // c1: positive
fitfunc_num->SetParLimits(1, 0, 1e5);      // c2: positive
fitfunc_num->SetParLimits(2, 0.0001, 1);   // lamda1: decay rate
fitfunc_num->SetParLimits(3, 0.0001, 1);   // lamda2: decay rate
fitfunc_num->SetParLimits(4, 0, 300);      // shift
fitfunc_num->SetParLimits(5, -100, 100);   // offset

if (ichannel == 1 and RunYear == 2016 and ISpreVFP == 1)	
	fitfunc_num->FixParameter(5, 0);

if (ichannel == 2 and RunYear == 2016 and ISpreVFP == 0){	
	fitfunc_num->FixParameter(2, 0.066);    
	fitfunc_num->FixParameter(3, 0.02);  
	fitfunc_num->FixParameter(4, 25.0);    
	fitfunc_num->FixParameter(5, 0);    
}
//if (ichannel == 2 and RunYear == 2016 and ISpreVFP == 1){
//	fitfunc_num->SetParLimits(0, 0, 1e6);
//	fitfunc_num->SetParLimits(1, 0, 1e3);
//	fitfunc_num->SetParLimits(2, 0.0001, 0.5);
//	fitfunc_num->SetParLimits(3, 0.0001, 0.5);
//	fitfunc_num->SetParLimits(4, 0, 100);
//	fitfunc_num->SetParLimits(5, -10, 10);
//}

if (ichannel == 1 and RunYear == 2016)
    fitfunc_num->SetParameters(800, 0.001, 0.056, 0.02, 0.0, 0.0);
if (ichannel == 1 && RunYear == 2016 && ISpreVFP == true) {
    fitfunc_num->SetParLimits(0, 0, 1e8);
    fitfunc_num->SetParLimits(1, 0, 1e6);
    fitfunc_num->SetParLimits(2, 0.001, 0.5);
    fitfunc_num->SetParLimits(3, 0.0001, 0.1);
    fitfunc_num->SetParLimits(4, -300, -165); // Force shift+200 < 35 GeV
    fitfunc_num->FixParameter(5, 0);          // offset=0
}
/*
if (ichannel == 1 && RunYear == 2016 && ISpreVFP == false) {
    fitfunc_num->SetParLimits(0, 0, 1e8);
    fitfunc_num->SetParLimits(1, 0, 1e5);
    fitfunc_num->SetParLimits(2, 0.001, 0.5);
    fitfunc_num->SetParLimits(3, 0.0001, 0.05);
    fitfunc_num->SetParLimits(4, -300, -165);
    fitfunc_num->FixParameter(5, 0);
}
*/

cout << "[DEBUG] Initial Parameters Set: " << endl;
for (int i = 0; i < 6; ++i) {
    cout << "  p" << i << " = " << fitfunc_num->GetParameter(i) << endl;
}

cout << "[DEBUG] Histogram integral (total entries): " << new_fakesPhoEt->Integral() << endl;
if (new_fakesPhoEt->Integral() == 0) {
    cout << "[ERROR] Histogram has zero total content. Fit will fail." << endl;
}

new_fakesPhoEt->Fit("fitfunc_num", "S"); // silent fit, retrieve result below
TF1 *fitnum = new_fakesPhoEt->GetFunction("fitfunc_num");

TVirtualFitter::SetMaxIterations(1000000);

TFitResultPtr r = new_fakesPhoEt->Fit("fitfunc_num", "R S");

fitfunc_num->SetNpx(1000);

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
		lowtoyptvalue[ii] = 100000;
		hightoyptvalue[ii] = 0;
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
	//			if (estimated < 0){
        //				std::cout << "[WARNING] Negative toy value: " << estimated << " at bin " << ibin << ", toy " << i << std::endl;
        //				estimated = 0;
    	//			}
				toyptvalue[ibin][i] = estimated;
				if(lowtoyptvalue[ibin] > estimated)lowtoyptvalue[ibin] = estimated;
				if(hightoyptvalue[ibin]< estimated)hightoyptvalue[ibin] = estimated;
			}

		}
		if (data3 < 0 || data4 < 0) 
			    std::cout << "[DEBUG] λ1 or λ2 negative: λ1 = " << data3 << ", λ2 = " << data4 << std::endl;
		if (data5 + 200 > 600) 
    			std::cout << "[DEBUG] Shift too large: shift+200 = " << (data5 + 200) << std::endl;
	}
	
	TCanvas *cangaus = new TCanvas("cangaus","",600,600);
	cangaus->cd();
	for(unsigned ibin(0); ibin < 565; ibin++){
		num_upper->SetBinContent(ibin, hightoyptvalue[ibin]); 
		num_lower->SetBinContent(ibin, lowtoyptvalue[ibin]);
		cout<< ibin << " "<< hightoyptvalue[ibin] << " " << lowtoyptvalue[ibin]<<endl;
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
	       
	if(channel == 1 and detType == 1)c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Plots%d%s/DoubleEG/JetFakeRate_transfer_DoubleEG_EB_%d%s_HTbinned.png", RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));
	else if(channel == 1 and detType == 2)c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Plots%d%s/DoubleEG/JetFakeRate_transfer_DoubleEG_EE_%d%s_HTbinned.png", RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));
	else if(channel == 2 and detType == 1)c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Plots%d%s/MuonEG/JetFakeRate_transfer_MuonEG_EB_%d%s_HTbinned.png", RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));
	else if(channel == 2 and detType == 2)c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/jetfakepho/Plots%d%s/MuonEG/JetFakeRate_transfer_MuonEG_EE_%d%s_HTbinned.png", RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));
}
