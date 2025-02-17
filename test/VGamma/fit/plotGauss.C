// Distribution of the scale factors derived from 1000 toy MC experiments. plot 29, AN
// this is to get the error on scalefactor
#include "TPaveText.h"
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1F.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TH2F.h"
#include "TMath.h"
#include "TSystem.h"
#include "TFractionFitter.h"
#include "TLatex.h"
#include "TGraphErrors.h"

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
#include "RooNumIntConfig.h"
#include "TROOT.h"
#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/tdrstyle.C"
//int ichannel = 1;
//int RunYear = 2016;
//bool isQCD=true;
//bool preVFP = true;

std::string whichVFP;

void plotGauss(int RunYear, int preVFP, int isQCD, int ichannel){
	if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
	else if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
        else whichVFP = "";

	gStyle->SetOptStat(0);
	gROOT->SetBatch(kTRUE);
	std::ifstream vgammascalefile;
	if(ichannel == 1) vgammascalefile.open(Form("/eos/uscms/store/user/tmishra/VGamma/VGamma_scalefactor_eg_%d%s.txt",RunYear,whichVFP.c_str()));
	if(ichannel == 2) vgammascalefile.open(Form("/eos/uscms/store/user/tmishra/VGamma/VGamma_scalefactor_mg_%d%s.txt",RunYear,whichVFP.c_str()));
	float fakescale, vgammascale, fakescaleerror, vgammascaleerror;
	float leplow, lephigh;
	
	float lowRange,highRange;
	for(unsigned i(0);  i < 1; i++){
		vgammascalefile >> leplow >> lephigh >> fakescale >> fakescaleerror >> vgammascale >> vgammascaleerror;
			if(isQCD==1){
				lowRange = 0.7*fakescale;
				highRange = 1.3*fakescale;}
			else{
				lowRange = 0.7*vgammascale;
				highRange = 1.3*vgammascale;}
	}
	cout<<"lowRange = "<<lowRange <<" highRange :"<<highRange<<endl;

	TH1F *p_frac_0;
	p_frac_0 = new TH1F("p_frac_0","",100,lowRange,highRange);

	float highest(0), lowest(1);
	float fittingerror(0), systematicerror(0), totalerror(0);
	
	double lowbound(2), highbound(0), norm(0);
	for(unsigned i(0);  i < 1000; i++){
		vgammascalefile >> leplow >> lephigh >> fakescale >> fakescaleerror >> vgammascale >> vgammascaleerror;
		if(isQCD==1){
			if(i == 0)norm = fakescale;
			if(i == 0)fittingerror = fakescaleerror;
			p_frac_0->Fill(fakescale);
			if(fakescale < lowbound)lowbound = fakescale;
			if(fakescale > highbound)highbound = fakescale;}
		else{
			if(i == 0)norm = vgammascale;   // for VGamma scale
			if(i == 0)fittingerror = vgammascaleerror;
			p_frac_0->Fill(vgammascale);
			if(vgammascale < lowbound)lowbound = vgammascale;
			if(vgammascale > highbound)highbound = vgammascale;}
	}
	float h_mean = p_frac_0->GetMean();
	float h_std = p_frac_0->GetMeanError();
	cout << "Histogram mean: " << h_mean << ", std error: " << h_std << endl;

	float fitRangeLow = std::max(lowRange, static_cast<float>(h_mean - 0.1));
	float fitRangeHigh = std::min(highRange, static_cast<float>(h_mean + 0.1));
	
	//float fitRangeLow = h_mean - 2 * h_std;
	//float fitRangeHigh = h_mean + 2 * h_std;

	if (p_frac_0->GetEntries() == 0) {
    		cout << "Histogram is empty. Skipping fit." << endl;
    		return;
	}

	if (fitRangeLow >= fitRangeHigh) {
    		cout << "Invalid fit range: [" << fitRangeLow << ", " << fitRangeHigh << "]" << endl;
    		return;
	}

	p_frac_0->Fit("gaus", "", "", fitRangeLow, fitRangeHigh);
	TF1* MyFit1 = p_frac_0->GetFunction("gaus");
	if (!MyFit1) {
    		cout << "Fit failed. Skipping further steps." << endl;
    		return;
	}

	float mu = MyFit1->GetParameter(1);
	float sigma = MyFit1->GetParameter(2);
	cout << "Fit results - Mu: " << mu << ", Sigma: " << sigma << endl;
	p_frac_0->Fit("gaus", "", "", mu - 0.1, mu + 0.1);
	//p_frac_0->Fit("gaus", "", "", mu - sigma, mu + sigma);

	MyFit1 = p_frac_0->GetFunction("gaus");
	float Mu = MyFit1->GetParameter(1);
        float Sigma = MyFit1->GetParameter(2);
	
	//systematicerror = 3*p_frac_0->GetFunction("gaus")->GetParameter(2);
	systematicerror = (highbound-norm) > (norm - lowbound)? (highbound-norm):(norm - lowbound); 
	if(p_frac_0->GetFunction("gaus")->GetParameter(1) + systematicerror > highest)highest=p_frac_0->GetFunction("gaus")->GetParameter(1) + systematicerror;
	if(p_frac_0->GetFunction("gaus")->GetParameter(1) - systematicerror < lowest)lowest=p_frac_0->GetFunction("gaus")->GetParameter(1) - systematicerror;
	totalerror = sqrt(systematicerror*systematicerror + fittingerror*fittingerror);


	std::string whichVFP, bkg, channel;
	if(ichannel == 1)	channel="eg";
	if(ichannel == 2)	channel="mg";

	if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
        else if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
        else whichVFP = "";

	if(isQCD)	bkg="QCD";
	else		bkg="VGamma";

	ofstream myfile;
	myfile.open("scaleFactors.txt", std::ios_base::app | std::ios_base::out);
	myfile << "highbound = " << highbound << std::endl;
	myfile << "lowbound = " << lowbound << std::endl<<endl;

	myfile <<"double factor_"<<channel<<bkg<<"_"<<RunYear<<whichVFP<<"("<<Mu<<");"<<endl;
	myfile <<"double factorerror_"<<channel<<bkg<<"_"<<RunYear<<whichVFP<<"("<<totalerror<<");"<<endl<<endl;
	myfile.close();

	float peakHeight = MyFit1->GetParameter(0);
    	float histMax = p_frac_0->GetMaximum();
   	float yMax = std::max(peakHeight, histMax) * 1.2;
	p_frac_0->SetMaximum(yMax);
	cout<< "Maximum:  "<< yMax<< endl;
		
	double chi2 = MyFit1->GetChisquare();
    	int ndf = MyFit1->GetNDF();
    	double chi2ndf = (ndf > 0) ? chi2 / ndf : 0; // Avoid division by zero
    	cout << "Chi2/ndf = " << chi2 << "/" << ndf << " = " << chi2ndf << endl;

	p_frac_0->SetMaximum(yMax);
	TCanvas *canscale = new TCanvas("canscale", "", 600, static_cast<int>(0.7 * 600)); // 0.65 corresponds to (1.0 - 0.35)


	p_frac_0->Draw();	

	TLatex title;
    	title.SetTextSize(0.1); 
    	title.SetTextAlign(22); 
    	title.SetTextFont(42); 

    	if (isQCD)
        	title.DrawLatexNDC(0.5, 0.95, "a_{fake}");
    	else
        	title.DrawLatexNDC(0.5, 0.95, "a_{V#gamma}");
	TLatex chi2Text;
    	chi2Text.SetTextSize(0.06);
    	chi2Text.SetNDC();          
    	chi2Text.DrawLatex(0.6, 0.8, Form("#chi^{2}/ndf = %.2f", chi2ndf));


	if(RunYear==2016 and preVFP == 1)       CMS_lumi(canscale, 1,ichannel, 11);
        else if(RunYear==2016 and preVFP == 0)  CMS_lumi(canscale, 2,ichannel, 11);
        else if(RunYear==2017)                  CMS_lumi(canscale, 3,ichannel, 11);
        else if(RunYear==2018)                  CMS_lumi(canscale, 4,ichannel, 11);

	if(isQCD==1){
		if(ichannel == 1) canscale->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d%s/QCDScale_eg_%d%s.png",RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));	
		if(ichannel == 2) canscale->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d%s/QCDScale_mg_%d%s.png",RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));	}
	else{
		if(ichannel == 1) canscale->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d%s/VGammaScale_eg_%d%s.png",RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));	
		if(ichannel == 2) canscale->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d%s/VGammaScale_mg_%d%s.png",RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));	}
}
int main(int argc, char** argv)
{
    if(argc < 3)
      cout << "You have to provide two arguments!!\n";
    plotGauss(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
    return 0;
}
