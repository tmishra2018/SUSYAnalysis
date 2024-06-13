#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<ctime>
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TF3.h"
#include "TH2F.h"
#include "TMath.h"
#include "TSystem.h"
#include "TFractionFitter.h"
#include "TLatex.h"
#include "TStyle.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooAbsReal.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooKeysPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "RooNumIntConfig.h"
#include "RooFFTConvPdf.h"
#include "RooPolynomial.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooFitResult.h"

#include "../../../include/tdrstyle.C"
#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_fakes.h"

bool useMC = true;
bool doIterate = true;

int  eventType = 2; // 1 = eg, 2 = mg, 3 = MCeg, 4 = MCmg
// to do in control region
float METLOWCUT = 0;
float METHIGHCUT = 70;

bool ISpreVFP = false;
int  RunYear = 2016;
int FitJetFake(float lowercut, float uppercut, int detType, float loweta, float higheta){

	setTDRStyle();   
	time_t now = time(0);
	int nBinTotal = 20;
	// No of bins in HadIsoCharge distribution
	if(lowercut > 150)nBinTotal = 10;

	char lowername[5];
	sprintf(lowername, "%d", (int)lowercut);
	// Low pT cut
	char uppername[5];
	if(uppercut < 1000){
	  sprintf(uppername, "%d", (int)uppercut);
	// High pT cut
	}
	else sprintf(uppername, "Inf");

	std::cout << "start fitting " << std::endl;
	gSystem->Load("/uscms/home/tmishra/work/CMSSW_10_2_22/src/SUSYAnalysis/lib/libAnaClasses.so");
	ofstream myfile;
	std::string  outputType, outputName, outputT;
	// Channel looking for
	std::string  outputDet;
	// Barrel or EndCap
	std::ostringstream  myfilename;  myfilename.str("");
	std::ostringstream  Hist1Dname;  Hist1Dname.str("");
	std::ostringstream  Hist2Dname;  Hist2Dname.str("");
	switch(eventType){
		case 1: outputT = "DoubleEG";	outputType = "DoubleEG-"; outputName = "eg-"; break;
		case 2: outputT = "MuonEG";   outputType = "MuonEG-"; outputName = "mg-"; break;
		case 3: outputType = "MCEG-"; outputName = "MCeg-"; break;
		case 4: outputType = "MCMG-"; outputName = "MCmg-"; break;
	}
	if(detType == 1)outputDet = "EB";
	else if(detType == 2)outputDet = "EE";

	std::string whichVFP;
	if(RunYear==2016 and ISpreVFP == true) whichVFP = "preVFP";
	if(RunYear==2016 and ISpreVFP == false) whichVFP = "postVFP";
	if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

        // ------------------------------------------------------------------------------------------------------------

	myfilename << "/eos/uscms/store/user/tmishra/jetfakepho/txt"<<RunYear<< whichVFP<<"/JetFakeRate-" << outputType << outputDet << ".txt";
	Hist1Dname << "/eos/uscms/store/user/tmishra/jetfakepho/Plots"<<RunYear<< whichVFP<<"/"<<outputT<<"/frac-" << lowername << "-" << uppername << "-" << outputType << outputDet << ".png";
	Hist2Dname << "/eos/uscms/store/user/tmishra/jetfakepho/Plots"<<RunYear<< whichVFP<<"/"<<outputT<<"/can2D-" << lowername << "-" << uppername << "-" << outputType << outputDet << ".png";
	myfile.open(myfilename.str().c_str(), std::ios_base::app | std::ios_base::out);

	ofstream logfile;
	logfile.open(Form("/eos/uscms/store/user/tmishra/jetfakepho/logs/fittingquality_%d%s.log",RunYear,whichVFP.c_str()));
	logfile << outputType << " " << lowername << " " << uppername << std::endl;
	
        // ------------------------------------------------------------------------------------------------------------
        // tree from Data
	TChain *datatree = new TChain("hadronTree");
	
	switch(eventType){
		//case 1: datatree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_egsignal_DoubleEG_ReMiniAOD_FullEcal.root"); break;
		case 1: datatree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s.root",RunYear,whichVFP.c_str())); break;
		case 2: datatree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s.root",RunYear,whichVFP.c_str())); break;
		case 3: datatree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_DYJetsToLL_%d.root",RunYear));
			datatree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_WJetsToLNu_%d.root",RunYear));
			datatree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_QCD_DoubleEM_%d.root",RunYear));
	
						break;
	}	

	std::ostringstream mcTreename; 
	mcTreename.str("");
	if(eventType == 1 || eventType == 3)mcTreename << "egTree";
	else if(eventType == 2 || eventType == 4)mcTreename << "mgTree";
	TChain *mctree = new TChain(mcTreename.str().c_str());
        mctree->Add(Form("/eos/uscms/store/user/tmishra/jetfakepho/files/plot_hadron_GJet_%d%s.root",RunYear,whichVFP.c_str()));
        //mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/plot_hadron_GJet.root");



        // ------------------------------------------------------------------------------------------------------------
	
	/*************************************/
	/*      double normfactor = par[0];  */  
	/*      double slope = par[1];       */
	/*      double constant = par[2];    */
	/*      double index = par[3];       */
	/*      double coeff = par[4];       */
	/*      double vtx_constant = par[5];*/
	/*      double vtx_slope = par[6];   */
	/*************************************/
	std::ostringstream elefake_config;
	elefake_config.str("");
	elefake_config << "/eos/uscms/store/user/tmishra/elefakepho/DataResult"<<RunYear<<whichVFP<<"/EleFakeRate-Data-ByPtVtx-" << outputDet <<".txt";
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
        else cout<<"file does not exist..  "<<"\n";
	elefake_file.close();
	TF3 f_elefake("f_elefake", fakerate_func,10,1000,0,100,0,2.5,7);
        // Electron faking as photon rate as parametrised function of Pt, nVtx
	f_elefake.SetParameters(scalefactor, ptslope, ptconstant, ptindex, ptcoeff, vtxconst, vtxslope);

        // ------------------------------------------------------------------------------------------------------------

	double   StandardCut = 0;
	double   StandardIso = 0;

//	Double_t SigmaCutLower_EB[]={0.0103,0.0104,0.0105,0.0106,0.0107,0.0108,0.0109,0.0110,0.0111,0.0112};
//      Double_t SigmaCutUpper_EB[]={0.0140,0.0145,0.0150,0.0155,0.0160,0.0165,0.0170,0.0175,0.0180,0.0185};
//      Double_t SigmaCutLower_EE[]={0.03013,0.0302,0.0303,0.0304,0.0305,0.0306,0.0307,0.0308,0.0309,0.031};
//      Double_t SigmaCutUpper_EE[]={0.035,0.036,0.037,0.038,0.039,0.04};

	Double_t SigmaCutLower_EB[]={0.0106,0.0107,0.0108,0.0109,0.0110,0.0111,0.0112,0.0113,0.0114,0.0115};
	Double_t SigmaCutUpper_EB[]={0.0140,0.0145,0.0150,0.0155,0.0160,0.0165,0.0170,0.0175,0.0180,0.0185};
	Double_t SigmaCutLower_EE[]={0.0272,0.0273,0.0274,0.0275,0.0276,0.0277,0.0278,0.0279,0.028,0.0281};
	Double_t SigmaCutUpper_EE[]={0.035,0.036,0.037,0.038,0.039,0.04};
	std::vector<double> SigmaCutLower;
	std::vector<double> SigmaCutUpper;
	SigmaCutLower.clear();
	SigmaCutUpper.clear();
//cut based photon ID 94X-V2 ID for Run2
	if(detType == 1){
	//	StandardCut = 0.0103;
	//	StandardIso = 1.295;
		StandardCut = 0.0106;
		StandardIso = 1.694;
		for(unsigned i(0); i<sizeof(SigmaCutLower_EB)/sizeof(Double_t); i++)SigmaCutLower.push_back(SigmaCutLower_EB[i]);
		for(unsigned i(0); i<sizeof(SigmaCutUpper_EB)/sizeof(Double_t); i++)SigmaCutUpper.push_back(SigmaCutUpper_EB[i]);
	}
	else if(detType == 2){
	//	StandardCut = 0.03013;
	//	StandardIso = 1.011;
		StandardCut = 0.0272;
		StandardIso = 2.089;
		for(unsigned i(0); i<sizeof(SigmaCutLower_EE)/sizeof(Double_t); i++)SigmaCutLower.push_back(SigmaCutLower_EE[i]);
		for(unsigned i(0); i<sizeof(SigmaCutUpper_EE)/sizeof(Double_t); i++)SigmaCutUpper.push_back(SigmaCutUpper_EE[i]);
	}
	unsigned nUpper = SigmaCutUpper.size();
	unsigned nLower = SigmaCutLower.size();
	std::ostringstream hname;
	hname.str(""); hname <<"ChIsoTar-pt-" << lowername << "-" << uppername;
	TH1D* h_target= new TH1D(hname.str().c_str(),";Iso_{h^{#pm}} (GeV);",nBinTotal, 0.0, 20);
	// fitting target
	hname.str(""); hname <<"mcChIso-pt-" <<  lowername << "-" << uppername;
	TH1D* mc_sig = new TH1D(hname.str().c_str(), ";Iso_{h^#pm} (GeV);",nBinTotal, 0.0, 20);
	// true photon
	hname.str(""); hname <<"mcChIso-fake-pt-" << lowername << "-" << uppername;
	TH1D* mc_fake = new TH1D(hname.str().c_str(),";Iso_{h^#pm} (GeV);",nBinTotal, 0.0, 20);
	TH1D* h_bg[nLower][nUpper];
	// jet fake photon
	
	TH1*  mc_predict[nLower][nUpper];
	// Final jet fake photon
	TH1D* mc_sbcontamination[nLower][nUpper];
	// side band contamination
	float templateCorrFactor[nLower][nUpper];
	TObjArray *templateContainer[nLower][nUpper];
	TFractionFitter* fitter[nLower][nUpper];
	// Fitter
	TH1D* result[nLower][nUpper];
	TCanvas *can[nLower][nUpper]; 
	hname.str(""); hname <<"can2D" << lowername << "-" << uppername;
	TCanvas *can2D = new TCanvas(hname.str().c_str(), hname.str().c_str(),1200,1200);
        can2D->SetBottomMargin(0.15);
        can2D->SetTopMargin(0.1);
        can2D->SetLeftMargin(0.1);
        can2D->SetRightMargin(0.1);
	// Loop over sigma cut lower and upper to name hostograms
	for(unsigned iUpper(0); iUpper<nUpper; iUpper++){
		for(unsigned iLower(0); iLower<nLower; iLower++){
			templateContainer[iLower][iUpper] = new TObjArray(2);
			hname.str("");
			hname << "can" << lowername << "-" << uppername << "_" << SigmaCutLower[iLower] << "_" << SigmaCutUpper[iUpper] ;
			can[iLower][iUpper]=new TCanvas(hname.str().c_str(),hname.str().c_str(),600,600);
			can[iLower][iUpper]->SetBottomMargin(0.15);
			hname.str(""); hname <<"ChIsoHad-pt-" <<  lowername << "-" << uppername << "_" << SigmaCutLower[iLower] << "_" << SigmaCutUpper[iUpper];
			h_bg[iLower][iUpper] = new TH1D(hname.str().c_str(), hname.str().c_str(), nBinTotal, 0.0, 20);
			hname.str(""); hname <<"mcChIso-sideband-pt-" <<  lowername << "-" << uppername << "_" << SigmaCutLower[iLower] << "_" << SigmaCutUpper[iUpper];
			mc_sbcontamination[iLower][iUpper] = new TH1D(hname.str().c_str(), hname.str().c_str(), nBinTotal, 0.0, 20);
		}
	}

	TFile *outputfile;
	outputfile = TFile::Open(Form("/uscms_data/d3/tmishra/jetfakepho/JetFakeRate-%s%d%s.root",outputName.c_str(),RunYear,whichVFP.c_str()),"RECREATE");
	outputfile->cd();
	hname.str("");
	hname << "fracHad2D_" << lowername << "-" << uppername;
	TH2F* fracHad2D = new TH2F(hname.str().c_str(),"  Hadron fraction;Lower #sigma_{i#etai#eta} Threshold(GeV);Upper #sigma_{i#etai#eta} Threshold(GeV)",9,SigmaCutLower[0], SigmaCutLower[SigmaCutLower.size()-1],10,0.014,0.019);
	// Fraction of hadrons for photon
        fracHad2D->GetYaxis()->SetTitleOffset(2.0);
        fracHad2D->GetXaxis()->SetTitleOffset(1.5);
        fracHad2D->SetTitleOffset(1.5);
	hname.str("");
	hname << "fracHad1D_" << lowername << "-" << uppername;
	TH1D* fracHad1D = new TH1D(hname.str().c_str(),hname.str().c_str(),500,0,1.0);
	// Fraction of hadrons for photon


//************ Signal(true photon) Tree **********************//
	float mc_phoEt(0);
	float mc_phoEta(0); 
	float mc_phoPhi(0); 
	float mc_phoSigma(0);
	float mc_phoChIso(0);
	std::vector<int> *mc_mcPID=0;
	std::vector<float> *mc_mcEta=0;
	std::vector<float> *mc_mcPhi=0;
	std::vector<float> *mc_mcPt=0;
	std::vector<int> *mc_mcMomPID=0; 

	mctree->SetBranchAddress("phoEt",     &mc_phoEt);
	mctree->SetBranchAddress("phoEta",    &mc_phoEta);
	mctree->SetBranchAddress("phoPhi",    &mc_phoPhi);
	mctree->SetBranchAddress("phoSigma",  &mc_phoSigma);
	mctree->SetBranchAddress("phoChIso",  &mc_phoChIso);
	if(useMC){
		mctree->SetBranchAddress("mcPID",   &mc_mcPID);
		mctree->SetBranchAddress("mcEta",   &mc_mcEta);
		mctree->SetBranchAddress("mcPhi",   &mc_mcPhi);
		mctree->SetBranchAddress("mcPt",    &mc_mcPt);
		mctree->SetBranchAddress("mcMomPID",&mc_mcMomPID);
	}

//************ data Tree **********************//

	float crosssection(1);
	float ntotalevent(1);
	float phoEt(0);
	float phoEta(0);
	float sigMET(0);
	float phoSigma(0);
	float phoChIso(0);
	std::vector<float> *eleproxyEt=0;
	std::vector<float> *eleproxyEta=0;
	std::vector<float> *eleproxyPhi=0;
	std::vector<float> *eleproxySigma=0;
	std::vector<float> *eleproxyChIso=0;
	std::vector<int>   *eleproxynVertex=0;
	if(eventType == 3 || eventType == 4){
		datatree->SetBranchAddress("crosssection",&crosssection);
		datatree->SetBranchAddress("ntotalevent", &ntotalevent);
	}
	datatree->SetBranchAddress("phoEt",     &phoEt);
	datatree->SetBranchAddress("phoEta",    &phoEta);
	datatree->SetBranchAddress("sigMET",    &sigMET);
	datatree->SetBranchAddress("phoSigma",  &phoSigma);
	datatree->SetBranchAddress("phoChIso",  &phoChIso);
	datatree->SetBranchAddress("eleproxyEt",&eleproxyEt);
	datatree->SetBranchAddress("eleproxyEta",&eleproxyEta);
	datatree->SetBranchAddress("eleproxyPhi",&eleproxyPhi);
	datatree->SetBranchAddress("eleproxySigma", &eleproxySigma);
	datatree->SetBranchAddress("eleproxyChIso", &eleproxyChIso);
	datatree->SetBranchAddress("eleproxynVertex",&eleproxynVertex);

        // ------------------------------------------------------------------------------------------------------------
        // Fake ratio = N_{fake} / N_{proxy}

	double hadfrac(0);
	double fittingError(0);
	double systematicError(0), lowestfrac(1), highestfrac(0);
	double mcfakerate(0);
	double nden(0), totalden(0), nsignalregion(0);
	double nMCtarget(0),  nMCfaketarget(0);
	double nnum[nLower][nUpper], totalnum[nLower][nUpper];
	double nsideband[nLower][nUpper], nsbPassIso[nLower][nUpper];
	for(unsigned ii(0); ii < nLower; ii++){
	// for N_{fake} estimated using template fit on Charge had Isolation
		for(unsigned jj(0); jj < nUpper; jj++){
			nnum[ii][jj] = 0;
			totalnum[ii][jj] = 0;
			nsideband[ii][jj] = 0;
			nsbPassIso[ii][jj] = 0;
		}
	}

	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);

		if(mc_phoEt < lowercut || mc_phoEt > uppercut)continue;
		if(detType == 1 && fabs(mc_phoEta) > 1.4442)continue;
		else if(detType == 2 && (fabs(mc_phoEta) < loweta || fabs(mc_phoEta) > higheta) )continue;

		bool isFake(true);
		unsigned mcIndex(0);
		unsigned mcPhoIndex(0);
		float mindR(0.7);
		float phodR(0.7);
		bool hasMatch(false);
		// photon PID = 22, checks if MC match with photon, using GJets MC
		for(unsigned ii(0); ii < mc_mcPID->size(); ii++){
			float dR = DeltaR(mc_phoEta, mc_phoPhi, (*mc_mcEta)[ii], (*mc_mcPhi)[ii]);
			float dE = fabs(mc_phoEt - (*mc_mcPt)[ii])/mc_phoEt;
			if((*mc_mcPID)[ii] == 22){
				phodR = DeltaR(mc_phoEta, mc_phoPhi, (*mc_mcEta)[ii], (*mc_mcPhi)[ii]);
				mcPhoIndex = ii;
			}
			if(dR < mindR && dE < 0.7){mindR = dR; mcIndex = ii; hasMatch = true;}
		}
		if(phodR < 0.1)mcIndex = mcPhoIndex;
		if(hasMatch)isFake = isHad(fabs((*mc_mcPID)[mcIndex]), fabs((*mc_mcMomPID)[mcIndex]));
		// checking motherID if photon has a matching

		bool isTrueTemplate = (mc_phoSigma <= StandardCut && !isFake);
		// true photon condition
		bool isInTargetRegion = (mc_phoSigma <= StandardCut && mc_phoChIso < StandardIso);
		// fitting target region

		if(mc_phoSigma <= StandardCut && isFake)mc_fake->Fill(mc_phoChIso);
		if(isTrueTemplate){
			mc_sig->Fill(mc_phoChIso);
			// SIGNAL photon
			nsignalregion += 1;
		}
		if(isInTargetRegion){
			nMCtarget += 1;
			if(isFake)nMCfaketarget += 1; 
		}
		for(unsigned iUpper(0); iUpper<nUpper; iUpper++){
			for(unsigned iLower(0); iLower<nLower; iLower++){
				if(mc_phoSigma > SigmaCutLower[iLower] && mc_phoSigma < SigmaCutUpper[iUpper] && !isFake){
					mc_sbcontamination[iLower][iUpper]->Fill(mc_phoChIso);
					// true photon in sideband
					nsideband[iLower][iUpper] += 1;
					if(mc_phoChIso < StandardIso)nsbPassIso[iLower][iUpper] += 1;
					//  also passing Iso cut
				}
			}
		}					
	}

	std::cout << "MC tree has been filled " <<mctree->GetEntries()<<std::endl;

	if(nMCtarget > 0)mcfakerate = 1.0*nMCfaketarget/nMCtarget;

        // -----------------------------------------------------------------------------------------------------------

	for(unsigned ievt(0); ievt < datatree->GetEntries(); ievt++){
		datatree->GetEntry(ievt);
		if(sigMET < METLOWCUT || sigMET > METHIGHCUT)continue;
		// pTmiss < 70 GeV CR

		double LumiWeight = 1;
		if(eventType == 3 || eventType == 4)	LumiWeight = 35.8*crosssection*1000/ntotalevent;
		//if(eventType == 3 || eventType == 4)	LumiWeight = getEvtWeight(RunYear,crosssection,ntotalevent);

                // counting electron fake photon, storing totalden, totalnum, nden, nnum

		for(unsigned iele(0); iele < eleproxyEt->size(); iele++){
			if(detType == 1 && fabs( (*eleproxyEta)[iele]) > 1.4442)continue;
			else if(detType == 2 && ( fabs( (*eleproxyEta)[iele]) < loweta || fabs( (*eleproxyEta)[iele]) > higheta) )continue;
			double w_ele = -1.0*f_elefake( (*eleproxyEt)[iele], (*eleproxynVertex)[iele], fabs( (*eleproxyEta)[iele]));
			w_ele = w_ele*LumiWeight;
			if( (*eleproxyEt)[iele] < lowercut || (*eleproxyEt)[iele]  > uppercut)continue;

			if( (*eleproxySigma)[iele] <= StandardCut){
				h_target->Fill( (*eleproxyChIso)[iele], w_ele);
				// Fitting TARGET
				totalden+= w_ele;
				if( (*eleproxyChIso)[iele] < StandardIso)nden+= w_ele;
			}
			else{
			// Numerator is stored in bins of sigmaIetaIeta
				for(unsigned iUpper(0); iUpper<nUpper; iUpper++){
					for(unsigned iLower(0); iLower<nLower; iLower++){
						if( (*eleproxySigma)[iele] > SigmaCutLower[iLower] && (*eleproxySigma)[iele] < SigmaCutUpper[iUpper]){
							h_bg[iLower][iUpper]->Fill( (*eleproxyChIso)[iele], w_ele);
							// Fake photon for SIetaIeta range
							totalnum[iLower][iUpper] += w_ele;
							if( (*eleproxyChIso)[iele] < StandardIso)nnum[iLower][iUpper] += w_ele;
						}
					}
				}
			}

		}

		if(phoEt < lowercut || phoEt > uppercut )continue;
		if(sigMET < METLOWCUT || sigMET > METHIGHCUT)continue;
		if(detType == 1 && fabs(phoEta) > 1.4442)continue;
		else if(detType == 2 && ( fabs(phoEta) < loweta || fabs(phoEta) > higheta) )continue;

                //counting photon candidate, storing totalden,totalnum,nden,nnum

		if(phoSigma <= StandardCut){
			h_target->Fill(phoChIso,LumiWeight);
			// Fitting TARGET
			totalden+=1*LumiWeight;
			if(phoChIso < StandardIso)nden+=1*LumiWeight;
		}
		else{
			// Numerator is stored in bins of sigmaIetaIeta
			for(unsigned iUpper(0); iUpper<nUpper; iUpper++){
				for(unsigned iLower(0); iLower<nLower; iLower++){
					if(phoSigma > SigmaCutLower[iLower] && phoSigma < SigmaCutUpper[iUpper]){
						h_bg[iLower][iUpper]->Fill(phoChIso, LumiWeight);
						// Fake photon for SIetaIeta range
						totalnum[iLower][iUpper] +=1*LumiWeight;
						if(phoChIso < StandardIso)nnum[iLower][iUpper] +=1*LumiWeight;
					}
				}
			}
		}

	}
	//h_target->Sumw2();
	//for(unsigned iUpper(0); iUpper<nUpper; iUpper++)
	//	for(unsigned iLower(0); iLower<nLower; iLower++)
	//		h_bg[iLower][iUpper]->Sumw2();

	std::cout << "data tree has been filled " <<datatree->GetEntries()<< std::endl;

	TCanvas *districan = new TCanvas("districan","",1200,600);
	districan->Divide(2);
	districan->cd(1);
	h_target->Draw();
	mc_sig->Scale(h_target->Integral(1,nBinTotal)/mc_sig->Integral(1,nBinTotal));
	mc_sig->SetLineColor(kRed);
	mc_sig->Draw("hist same");	
	districan->cd(2);
	h_bg[0][0]->Draw();

	std::cout << "target " << h_target->Integral(1,nBinTotal) << std::endl;

	for(unsigned iUpper(0); iUpper<nUpper; iUpper++){
		for(unsigned iLower(0); iLower<nLower; iLower++){
			double iteratorUncertainty(1), lastUpdatePurity(1);
			double frac_passSigma = nsbPassIso[iLower][iUpper]/(1.0*mc_sbcontamination[iLower][iUpper]->GetEntries());
			if(mc_sbcontamination[iLower][iUpper]->GetEntries()==0)frac_passSigma = 0;
			double phoPurity(0);
			templateCorrFactor[iLower][iUpper] = 0.0;
	int itertime = 0;
                        // Iterative method to remove the true photon sigmaIetaIeta sideband contaminating hadronic template
			if(doIterate){
				while(iteratorUncertainty > 0.01){

					itertime += 1; logfile << "iter " << itertime << std::endl;
                                        // RAW hadronic template - true photon in sideband
					TH1D *hist_new=(TH1D*)h_bg[iLower][iUpper]->Clone();
					hist_new->Add(mc_sbcontamination[iLower][iUpper], -1*templateCorrFactor[iLower][iUpper]/mc_sbcontamination[iLower][iUpper]->GetEntries());
					hist_new->Sumw2();
					templateContainer[iLower][iUpper]->Add(mc_sig);
					templateContainer[iLower][iUpper]->Add(hist_new);
                                        // Fitter with target, signal and background(fake)
					fitter[iLower][iUpper]= new TFractionFitter(h_target, templateContainer[iLower][iUpper], "V");
					fitter[iLower][iUpper]->Constrain(1,0.0,1.0);
					Int_t status = fitter[iLower][iUpper]->Fit();
					double tmpSig, tmpErrorSig;
					if(status == 0){
						fitter[iLower][iUpper]->GetResult(0, tmpSig,  tmpErrorSig);
						phoPurity = tmpSig;
					}
                                        // f_i - f_{i+1}
					iteratorUncertainty = fabs(lastUpdatePurity - phoPurity)/phoPurity;
					lastUpdatePurity = phoPurity;
                                        //  N_{true-in-sb} = N_{sig} * (1-f_i) * (p_{sb} / p_{sig})
					templateCorrFactor[iLower][iUpper] = totalden*phoPurity*nsideband[iLower][iUpper]/nsignalregion;
					templateContainer[iLower][iUpper]->Clear();
					if( fabs(SigmaCutLower[iLower]- StandardCut) < 0.00005 && iUpper == 0)logfile << "  Standard "  << "iteratorUncertainty " << iteratorUncertainty << " phoPurity " << phoPurity << "  templateCorrFactor " << templateCorrFactor[iLower][iUpper] << std::endl;
					else logfile << "  " << SigmaCutLower[iLower] << " " << SigmaCutLower[iUpper] << " " << "iteratorUncertainty " << iteratorUncertainty << " phoPurity " << phoPurity << "  templateCorrFactor " << templateCorrFactor[iLower][iUpper] << std::endl;
					hist_new->Delete();
				}
			}
                        //  No iteration method
			else{
				iteratorUncertainty = 0;
				frac_passSigma = 0;
				templateCorrFactor[iLower][iUpper] = 0;
				templateContainer[iLower][iUpper]->Add(mc_sig);
				templateContainer[iLower][iUpper]->Add(h_bg[iLower][iUpper]);
				fitter[iLower][iUpper]= new TFractionFitter(h_target, templateContainer[iLower][iUpper], "V");
				fitter[iLower][iUpper]->Constrain(1,0.0,1.0);
				Int_t status = fitter[iLower][iUpper]->Fit();
			}


			// different loop
			if(iteratorUncertainty <= 0.01 || !doIterate){
				double den = nden/totalden;  
				double denerror = (1.0/sqrt(nden)+ 1.0/sqrt(totalden))*den;
				double num = (nnum[iLower][iUpper] - templateCorrFactor[iLower][iUpper]*frac_passSigma)/(totalnum[iLower][iUpper] - templateCorrFactor[iLower][iUpper]); 
				double numerror = (1.0/sqrt(nnum[iLower][iUpper]) + 1.0/sqrt(totalnum[iLower][iUpper]))*num;
				double fracBkg, errorBkg;
				double fracSig, errorSig;
                                // histogram corresponding to the total MC prediction (which is not the same as a simple weighted sum of the input MC distributions)
                                // can be obtained by GetPlot()
                                // result[][]  = signal + bkg, as fitter is defined
				result[iLower][iUpper] = (TH1D*)fitter[iLower][iUpper]->GetPlot();
				result[iLower][iUpper]->SetMinimum(10);
                                // For fits that converged, parameter values and errors can be obtained through GetResult()
				fitter[iLower][iUpper]->GetResult(0, fracSig, errorSig);
				fitter[iLower][iUpper]->GetResult(1, fracBkg, errorBkg);
				can[iLower][iUpper]->cd();
				gStyle->SetOptStat(0);
				hname.str("");
				hname << lowername << "< pt <" << uppername; 
				h_target->SetTitle(hname.str().c_str());
				h_target->SetMinimum(h_target->GetBinContent(nBinTotal));
				result[iLower][iUpper]->SetMinimum(100);
                                gPad->SetTopMargin(0.1);
				gPad->SetLogy();
				h_target->GetXaxis()->SetTitleOffset(0.8);
				h_target->SetMarkerStyle(20);
				h_target->Draw("Ep");
				result[iLower][iUpper]->SetLineColor(kBlue);
				result[iLower][iUpper]->SetFillStyle(1001);
				result[iLower][iUpper]->SetFillColor(kBlue);
				result[iLower][iUpper]->SetMarkerColor(kBlue);
				result[iLower][iUpper]->Draw("same");
                                // GetMCPrediction (parm) Return the adjusted MC template (Aji) for template (parm).
				mc_predict[iLower][iUpper] = fitter[iLower][iUpper]->GetMCPrediction(1);
				mc_predict[iLower][iUpper]->SetMinimum(100);
                                // normalised to total, scaled to fracBkg
				mc_predict[iLower][iUpper]->Scale(1.0/mc_predict[iLower][iUpper]->Integral(1,21)*result[iLower][iUpper]->Integral(1,21)*fracBkg);
				mc_predict[iLower][iUpper]->SetLineColor(kGreen);
				mc_predict[iLower][iUpper]->SetFillStyle(1001);
				mc_predict[iLower][iUpper]->SetFillColor(kGreen);
				mc_predict[iLower][iUpper]->SetMarkerColor(kGreen);
				mc_predict[iLower][iUpper]->Draw("hist same");

		////		mc_fake->SetLineColor(kRed);
		////		mc_fake->Draw("same");

				h_target->Draw("Ep same");
				TLatex* latex = new TLatex();
		//		hname.str("");
		//		hname << "Bkg% = " << fracBkg << " #pm " << errorBkg << std::endl;
		//		latex->DrawLatex(2,0.9*h_target->GetMaximum(),hname.str().c_str()); 
				hname.str("");
                                // NDF = the number of points used in the fit minus the number of templates.
                                float chindf = fitter[iLower][iUpper]->GetChisquare()/fitter[iLower][iUpper]->GetNDF();
				hname << "#chi^{2}/ndof = " << chindf;
				latex->DrawLatex(5,0.15*h_target->GetMaximum(),hname.str().c_str());
		//		hname.str("");
		//		hname << "Hadfrac="<< num*fracBkg/den << " #pm " << (numerror/num + denerror/den + errorBkg/fracBkg)*num*fracBkg/den;
		//		latex->DrawLatex(2,0.15*h_target->GetMaximum(),hname.str().c_str());
		// for sigma 0.0106-0.0140, 1D histograms
				if(SigmaCutLower[iLower] == StandardCut && iUpper == 0){
					// imp! num=nnum/totalnum and en=nden/totalden
					hadfrac =  num*fracBkg/den;
					fittingError = (numerror/num + denerror/den + errorBkg/fracBkg)*num*fracBkg/den ;
					TLegend *leg = new TLegend(0.6,0.75,0.87,0.88);
					leg->AddEntry(h_target, "observed");
					leg->AddEntry(result[iLower][iUpper], "true photons");
					leg->AddEntry(mc_predict[iLower][iUpper], "hadrons");
					leg->Draw("same");
					can[iLower][iUpper]->SaveAs(Hist1Dname.str().c_str());
					
					mc_sig->Write();
					h_bg[iLower][iUpper]->Write();	
					h_target->Write();
					//result[iLower][iUpper]->Write();
					//mc_predict[iLower][iUpper]->Write();
				}
				if(fracBkg >0 && fracBkg < 1)fracHad2D->Fill(SigmaCutLower[iLower]+0.00005, SigmaCutUpper[iUpper]+0.00005, num*fracBkg/den);
				if(fracBkg >0 && fracBkg < 1)fracHad1D->Fill(num*fracBkg/den);
				if(lowestfrac > num*fracBkg/den)lowestfrac = num*fracBkg/den;
				if(highestfrac < num*fracBkg/den)highestfrac = num*fracBkg/den; 
			}
		}
	}

	if(hadfrac < 0.001){
		hadfrac = fracHad1D->GetMean();	
         }

	can2D->cd();
        fracHad2D->SetTitle("hadron fraction");
        can2D->SetTitle("hadron fraction");
	fracHad2D->Draw("colz");
	can2D->SaveAs(Hist2Dname.str().c_str());
	systematicError = highestfrac - lowestfrac;
	myfile << " " << lowername << " " << uppername << " " << hadfrac << " "; 
	// why this large error in 3 bins ???
	//if(RunYear == 2016 and ISpreVFP == true and lowercut == 150 and uppercut == 180) myfile << sqrt(fittingError*fittingError)<< " ";
	myfile << sqrt(fittingError*fittingError + systematicError*systematicError) << " ";
	myfile << systematicError << "	" << " ";
	myfile << mcfakerate << std::endl;
	char* dt = ctime(&now);
	if(uppercut > 500){
		myfile<< std::endl;
		myfile<< "ibin pt_lower pt_upper fakerate error errorsystematic" << std::endl;
		myfile << dt << std::endl;
		myfile << "SigmaCutLower = ";
		for(unsigned l(0); l < nLower; l++)myfile << SigmaCutLower[l] << " ";
		myfile << std::endl;
		for(unsigned u(0); u < nUpper; u++)myfile << SigmaCutUpper[u] << " ";
		myfile << std::endl;
	} 
	myfile.close();
	logfile.close();
//	outputfile->Write();
	return 1;
}
