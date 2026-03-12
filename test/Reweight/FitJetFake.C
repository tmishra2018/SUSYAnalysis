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

#include "../../include/tdrstyle.C"
#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_jet.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_fakes.h"

bool useMC = true;
bool doIterate = true;


int FitJetFake(float lowercut, float uppercut, int detType, int RunYear, bool preVFP) {
	setTDRStyle();   
	time_t now = time(0);

	std::cout << "start fitting " << std::endl;
	gSystem->Load("../../lib/libAnaClasses.so");
	
	std::string whichVFP;
        if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
  	else if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
  	else if(RunYear==2017 or  RunYear == 2018) whichVFP = "";
	
	std::string whichRegion;
	if(detType==1) whichRegion = "EB";
	else if(detType==2) whichRegion = "EE";
	
	ofstream myfile;
	myfile.open(Form("/eos/uscms/store/user/tmishra/jetfakepho/files/JetFakeRate-ISR-%s_%d%s.txt",whichRegion.c_str(),RunYear,whichVFP.c_str()), std::ios_base::app | std::ios_base::out);

	std::ostringstream datasetname;
	datasetname.str("");
	datasetname << "/eos/uscms/store/user/tmishra/jetfakepho/files/plot_hadron_ISR_"<<RunYear<<whichVFP<<".root";
	// made with analysis_mgHadron.C
	char lowername[4];
	sprintf(lowername, "%d", (int)lowercut);
	char uppername[4];
	if(uppercut < 1000){
	  sprintf(uppername, "%d", (int)uppercut);
	}
	else sprintf(uppername, "Inf");

	double   StandardCut = 0;
	double   StandardIso = 0;
	Double_t SigmaCutLower_EB[]={0.0106,0.0107,0.0108,0.0109,0.0110,0.0111,0.0112,0.0113,0.0114,0.0115};
        Double_t SigmaCutUpper_EB[]={0.0140,0.0145,0.0150,0.0155,0.0160,0.0165,0.0170,0.0175,0.0180,0.0185};
        Double_t SigmaCutLower_EE[]={0.0272,0.0273,0.0274,0.0275,0.0276,0.0277,0.0278,0.0279,0.028,0.0281};
        Double_t SigmaCutUpper_EE[]={0.035,0.036,0.037,0.038,0.039,0.04};

	std::vector<double> SigmaCutLower;
	std::vector<double> SigmaCutUpper;
	SigmaCutLower.clear();
	SigmaCutUpper.clear();
	if(detType == 1){
		StandardCut = 0.0106;
		StandardIso = 1.694;
		for(unsigned i(0); i<sizeof(SigmaCutLower_EB)/sizeof(Double_t); i++)SigmaCutLower.push_back(SigmaCutLower_EB[i]);
		for(unsigned i(0); i<sizeof(SigmaCutUpper_EB)/sizeof(Double_t); i++)SigmaCutUpper.push_back(SigmaCutUpper_EB[i]);
	}
	else if(detType == 2){
		StandardCut = 0.0272;
                StandardIso = 2.089;
		for(unsigned i(0); i<sizeof(SigmaCutLower_EE)/sizeof(Double_t); i++)SigmaCutLower.push_back(SigmaCutLower_EE[i]);
		for(unsigned i(0); i<sizeof(SigmaCutUpper_EE)/sizeof(Double_t); i++)SigmaCutUpper.push_back(SigmaCutUpper_EE[i]);
	}
	unsigned nUpper = SigmaCutUpper.size(); 
	unsigned nLower = SigmaCutLower.size();
	std::ostringstream hname;
	hname.str(""); hname <<"ChIsoTar-pt-" << lowername << "-" << uppername;
	TH1D* h_target= new TH1D(hname.str().c_str(),";Iso_{h^{#pm}} (GeV);",20, 0.0, 20); 
	hname.str(""); hname <<"mcChIso-pt-" <<  lowername << "-" << uppername;
	TH1D* mc_sig = new TH1D(hname.str().c_str(), ";Iso_{h^#pm} (GeV);",20, 0.0, 20);
	hname.str(""); hname <<"mcChIso-fake-pt-" << lowername << "-" << uppername;
	TH1D* mc_fake = new TH1D(hname.str().c_str(),";Iso_{h^#pm} (GeV);",20, 0.0, 20);
	TH1D* h_bg[nLower][nUpper];
	
	TH1*  mc_predict[nLower][nUpper];
	TH1D* mc_sbcontamination[nLower][nUpper];
	float templateCorrFactor[nLower][nUpper];
	TObjArray *templateContainer[nLower][nUpper];
	TFractionFitter* fitter[nLower][nUpper];
	TH1D* result[nLower][nUpper];
	TCanvas *can[nLower][nUpper]; 
	hname.str(""); hname <<"can2D" << lowername << "-" << uppername;
	TCanvas *can2D = new TCanvas(hname.str().c_str(), hname.str().c_str(),600,600);
	for(unsigned iUpper(0); iUpper<nUpper; iUpper++){
		for(unsigned iLower(0); iLower<nLower; iLower++){
			templateContainer[iLower][iUpper] = new TObjArray(2);
			hname.str("");
			hname << "can" << lowername << "-" << uppername << "_" << SigmaCutLower[iLower] << "_" << SigmaCutUpper[iUpper] ;
			can[iLower][iUpper]=new TCanvas(hname.str().c_str(),hname.str().c_str(),600,600);
			hname.str(""); hname <<"ChIsoHad-pt-" <<  lowername << "-" << uppername << "_" << SigmaCutLower[iLower] << "_" << SigmaCutUpper[iUpper];
			h_bg[iLower][iUpper] = new TH1D(hname.str().c_str(), hname.str().c_str(), 20, 0.0, 20);
			hname.str(""); hname <<"mcChIso-sideband-pt-" <<  lowername << "-" << uppername << "_" << SigmaCutLower[iLower] << "_" << SigmaCutUpper[iUpper];
			mc_sbcontamination[iLower][iUpper] = new TH1D(hname.str().c_str(), hname.str().c_str(),20, 0.0, 20);   
		}
	}

	hname.str("");
	hname << "fracHad2D_" << lowername << "-" << uppername;
	TH2F* fracHad2D = new TH2F(hname.str().c_str(),"  Hadron fraction;Lower #sigma_{i#etai#eta} Threshold(GeV);Upper #sigma_{i#etai#eta} Threshold(GeV)",9,SigmaCutLower[0], SigmaCutLower[SigmaCutLower.size()-1],10,0.014,0.019);
	hname.str("");
	hname << "fracHad1D_" << lowername << "-" << uppername;
	TH1D* fracHad1D = new TH1D(hname.str().c_str(),hname.str().c_str(),500,0,1.0);


	//************ Signal Tree **********************//
	TChain *mctree = new TChain("egTree");
	mctree->Add(Form("/eos/uscms/store/user/tmishra/jetfakepho/files/plot_hadron_mgGJet_%d%s.root",RunYear,whichVFP.c_str()));
	// made with ../jetFakePho/src/analysis_mgGJet.C
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

	//TChain *egtree = new TChain("hadronTree");
	TChain *egtree = new TChain("signalTree");
	egtree->Add(datasetname.str().c_str());
	float phoEt(0);
	float phoEta(0);
	float phoSigma(0);
	float phoChIso(0);
	float llmass(0);
	egtree->SetBranchAddress("phoEt",     &phoEt);
	egtree->SetBranchAddress("phoEta",    &phoEta);
	egtree->SetBranchAddress("phoSigma",  &phoSigma);
	egtree->SetBranchAddress("phoChIso",  &phoChIso);
	egtree->SetBranchAddress("dilepMass", &llmass);


	double hadfrac(0);
	double fittingError(0);
	double systematicError(0), lowestfrac(1), highestfrac(0);
	double mcfakerate(0);
	double nden(0), totalden(0), nsignalregion(0);
	double nMCtarget(0),  nMCfaketarget(0);
	double nnum[nLower][nUpper], totalnum[nLower][nUpper];
	double nsideband[nLower][nUpper], nsbPassIso[nLower][nUpper];
	for(unsigned ii(0); ii < nLower; ii++){
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
		else if(detType == 2 && (fabs(mc_phoEta) < 1.56 || fabs(mc_phoEta) > 2.4) )continue;

		bool isFake(true);
		unsigned mcIndex(0);
		unsigned mcPhoIndex(0);
		float mindR(0.7);
		float phodR(0.7);
		bool hasMatch(false);
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

		bool isTrueTemplate = (mc_phoSigma <= StandardCut && !isFake);
		bool isInTargetRegion = (mc_phoSigma <= StandardCut && mc_phoChIso < StandardIso);

		if(mc_phoSigma <= StandardCut && isFake)mc_fake->Fill(mc_phoChIso);
		if(isTrueTemplate){
			mc_sig->Fill(mc_phoChIso);
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
					nsideband[iLower][iUpper] += 1;
					if(mc_phoChIso < StandardIso)nsbPassIso[iLower][iUpper] += 1;
				}
			}
		}					
	}
	std::cout << "MC tree has been filled" << std::endl;

	if(nMCtarget > 0)mcfakerate = 1.0*nMCfaketarget/nMCtarget;

	for(unsigned ievt(0); ievt < egtree->GetEntries(); ievt++){
		egtree->GetEntry(ievt);
		
		if(llmass < 80 || llmass > 100)continue;

		if(phoEt < lowercut || phoEt > uppercut )continue;
		if(detType == 1 && fabs(phoEta) > 1.4442)continue;
		else if(detType == 2 && ( fabs(phoEta) < 1.56 || fabs(phoEta) > 2.4) )continue;

		if(phoSigma <= StandardCut){
			h_target->Fill(phoChIso);
			totalden+=1;
			if(phoChIso < StandardIso)nden+=1;
		}
		else{
			for(unsigned iUpper(0); iUpper<nUpper; iUpper++){
				for(unsigned iLower(0); iLower<nLower; iLower++){
					if(phoSigma > SigmaCutLower[iLower] && phoSigma < SigmaCutUpper[iUpper]){
						h_bg[iLower][iUpper]->Fill(phoChIso);
						totalnum[iLower][iUpper] +=1;
						if(phoChIso < StandardIso)nnum[iLower][iUpper] +=1;
					}
				}
			}
		}

	}

	std::cout << "data tree has been filled" << std::endl;

	TCanvas *districan = new TCanvas("districan","",1200,600);
	districan->Divide(2);
	districan->cd(1);
	h_target->Draw();
	mc_sig->Scale(h_target->Integral(1,20)/mc_sig->Integral(1,20));
	mc_sig->SetLineColor(kRed);
	mc_sig->Draw("hist same");	
	districan->cd(2);
	h_bg[0][0]->Draw();

	std::cout << "target " << h_target->Integral(1,20) << std::endl;
	std::cout << "MC " << mc_sig->Integral(1,20) << std::endl;

	for(unsigned iUpper(0); iUpper<nUpper; iUpper++){
		for(unsigned iLower(0); iLower<nLower; iLower++){
			double iteratorUncertainty(1), lastUpdatePurity(1);
			double frac_passSigma = nsbPassIso[iLower][iUpper]/(1.0*mc_sbcontamination[iLower][iUpper]->GetEntries());
			if(mc_sbcontamination[iLower][iUpper]->GetEntries()==0)frac_passSigma = 0;
			double phoPurity(0);
			templateCorrFactor[iLower][iUpper] = 0.0;
			if(doIterate){
				while(iteratorUncertainty > 0.01){

					TH1D *hist_new=(TH1D*)h_bg[iLower][iUpper]->Clone();
					hist_new->Add(mc_sbcontamination[iLower][iUpper], -1*templateCorrFactor[iLower][iUpper]/mc_sbcontamination[iLower][iUpper]->GetEntries());
					hist_new->Sumw2();
					templateContainer[iLower][iUpper]->Add(mc_sig);
					templateContainer[iLower][iUpper]->Add(hist_new);
					fitter[iLower][iUpper]= new TFractionFitter(h_target, templateContainer[iLower][iUpper], "V");
					fitter[iLower][iUpper]->Constrain(1,0.0,1.0);
					Int_t status = fitter[iLower][iUpper]->Fit();
					double tmpSig, tmpErrorSig;
					if(status == 0){
						fitter[iLower][iUpper]->GetResult(0, tmpSig,  tmpErrorSig);
						phoPurity = tmpSig;
					}
					iteratorUncertainty = fabs(lastUpdatePurity - phoPurity)/phoPurity;
					lastUpdatePurity = phoPurity;
					templateCorrFactor[iLower][iUpper] = totalden*phoPurity*nsideband[iLower][iUpper]/nsignalregion;
					templateContainer[iLower][iUpper]->Clear();
					hist_new->Delete();
				}
			}
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
			if(iteratorUncertainty <= 0.01 || !doIterate){
				double den = nden/totalden; 
				double denerror = (1.0/sqrt(nden)+ 1.0/sqrt(totalden))*den;
				double num = (nnum[iLower][iUpper] - templateCorrFactor[iLower][iUpper]*frac_passSigma)/(totalnum[iLower][iUpper] - templateCorrFactor[iLower][iUpper]); 
				double numerror = (1.0/sqrt(nnum[iLower][iUpper]) + 1.0/sqrt(totalnum[iLower][iUpper]))*num;
				double fracBkg, errorBkg;
				double fracSig, errorSig;
				result[iLower][iUpper] = (TH1D*)fitter[iLower][iUpper]->GetPlot();
				result[iLower][iUpper]->SetMinimum(10);
				fitter[iLower][iUpper]->GetResult(0, fracSig, errorSig);
				fitter[iLower][iUpper]->GetResult(1, fracBkg, errorBkg);
				can[iLower][iUpper]->cd();
				gStyle->SetOptStat(0);
				hname.str("");
				hname << lowername << "< pt <" << uppername; 
				h_target->SetTitle(hname.str().c_str());
				h_target->SetMinimum(h_target->GetBinContent(20));
				result[iLower][iUpper]->SetMinimum(100);
				//gPad->SetLogy();
				h_target->SetMarkerStyle(20);
				h_target->Draw("Ep");
				result[iLower][iUpper]->SetLineColor(kBlue);
				result[iLower][iUpper]->SetFillStyle(1001);
				result[iLower][iUpper]->SetFillColor(kBlue);
				result[iLower][iUpper]->SetMarkerColor(kBlue);
				result[iLower][iUpper]->Draw("same");

				mc_predict[iLower][iUpper] = fitter[iLower][iUpper]->GetMCPrediction(1);
				mc_predict[iLower][iUpper]->SetMinimum(100);
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
				hname << "#chi^{2}/ndof = " << fitter[iLower][iUpper]->GetChisquare() << "/" << fitter[iLower][iUpper]->GetNDF(); 
				latex->DrawLatex(5,0.15*h_target->GetMaximum(),hname.str().c_str());
		//		hname.str("");
		//		hname << "Hadfrac="<< num*fracBkg/den << " #pm " << (numerror/num + denerror/den + errorBkg/fracBkg)*num*fracBkg/den;
		//		latex->DrawLatex(2,0.15*h_target->GetMaximum(),hname.str().c_str());
		 
				if(SigmaCutLower[iLower] == StandardCut && iUpper == 0){
					hadfrac =  num*fracBkg/den;
					fittingError = (numerror/num + denerror/den + errorBkg/fracBkg)*num*fracBkg/den ;
					TLegend *leg = new TLegend(0.6,0.75,0.87,0.9);
					leg->AddEntry(h_target, "observed");
					leg->AddEntry(result[iLower][iUpper], "true photons");
					leg->AddEntry(mc_predict[iLower][iUpper], "hadrons");
					leg->Draw("same");
					hname.str("");
					hname << "/eos/uscms/store/user/tmishra/ISRweighting/JetFakePlots/frac-" << lowername << "-" << uppername << "-ISR-FullEcal";
					if(detType == 1)hname << "-EB.pdf";
					else if(detType == 2)hname << "-EE.pdf";
					can[iLower][iUpper]->SaveAs(hname.str().c_str());
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

	hname.str(""); 
	hname << "/eos/uscms/store/user/tmishra/ISRweighting/JetFakePlots/can2D-" << lowername << "-" << uppername << "-ISR-FullEcal";
	if(detType == 1)hname << "-EB.pdf";
	else if(detType == 2)hname << "-EE.pdf";
	
	can2D->cd();
	fracHad2D->Draw("colz");
	can2D->SaveAs(hname.str().c_str());
	systematicError = highestfrac - lowestfrac;
	myfile << " " << lowername << " " << uppername << " " << hadfrac << " "; 
	myfile << sqrt(fittingError*fittingError + systematicError*systematicError) << " ";
	myfile << systematicError << "	" << " ";
	myfile << mcfakerate;
        if(detType == 1)myfile << " EB" <<  std::endl;
        else if(detType == 2)myfile << " EE" <<  std::endl;
	char* dt = ctime(&now);
	//if(uppercut > 500){
	//	myfile<< std::endl;
	//	myfile<< "ibin pt_lower pt_upper fakerate error errorsystematic" << std::endl;
	//	myfile << dt << std::endl;
	//	myfile << "SigmaCutLower = ";
	//	for(unsigned l(0); l < nLower; l++)myfile << SigmaCutLower[l] << " ";
	//	myfile << std::endl;
	//	for(unsigned u(0); u < nUpper; u++)myfile << SigmaCutUpper[u] << " ";
	//	myfile << std::endl;
	//	myfile << datasetname.str().c_str() << std::endl;
	//} 
	myfile.close();
	return 1;
}
