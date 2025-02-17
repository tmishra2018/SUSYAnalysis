#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1D.h"
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
#include "TRandom3.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TROOT.h"
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

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/tdrstyle.C"

//int RunYear = 2017;
//bool preVFP = false;
std::string whichVFP;

int Fitfractionmg(int ih,int metlow, int methigh, int leplow, int lephigh, int isocut, int RunYear, bool preVFP){
	gROOT->SetBatch(kTRUE);
	gStyle->SetOptStat(0);
	setTDRStyle();
	gStyle->SetErrorX(0);
	gStyle->SetTitleXOffset(2.5);
  	RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setRealValue("maxSteps",50000); 
	std::ostringstream histname;
	if(RunYear==2016 and preVFP == true) whichVFP = "preVFP";
	else if(RunYear==2016 and preVFP == false) whichVFP = "postVFP";
        else whichVFP = "";

	ofstream myfile;
	myfile.open(Form("/eos/uscms/store/user/tmishra/VGamma/VGamma_scalefactor_mg_%d%s.txt",RunYear,whichVFP.c_str()), std::ios_base::app | std::ios_base::out);
	TString filepath = "/eos/uscms/store/user/tmishra/Background/";
	std::ostringstream filename;
	filename.str("");
	// // data in the control region
	filename << filepath << "controlTree_mg_signal_" << "met" << metlow << "_" << methigh << "_pt" << leplow << "_" << lephigh << "_" << RunYear << whichVFP <<".root";
	TFile *file_sig = TFile::Open(filename.str().c_str());
	filename.str("");
	filename << filepath << "controlTree_mg_eleBkg_" << "met" << metlow << "_" << methigh << "_pt" << leplow << "_" << lephigh << "_" << RunYear << whichVFP <<".root";
	TFile *file_ele = TFile::Open(filename.str().c_str());
	filename.str("");
	filename << filepath << "controlTree_mg_jetbkg_" << "met" << metlow << "_" << methigh << "_pt" << leplow << "_" << lephigh << "_" << RunYear << whichVFP <<".root";
	TFile *file_jet = TFile::Open(filename.str().c_str());
	filename.str("");
	filename << filepath << "controlTree_mg_rareBkg_" << "met" << metlow << "_" << methigh << "_pt" << leplow << "_" << lephigh << "_" << RunYear << whichVFP <<".root";
	TFile *file_rare = TFile::Open(filename.str().c_str());
	filename.str("");
	filename << filepath << "controlTree_mg_qcd_" << "met" << metlow << "_" << methigh << "_pt" << leplow << "_" << lephigh <<  "_" << RunYear << whichVFP <<".root";
	TFile *file_qcd = TFile::Open(filename.str().c_str());
	filename.str("");
	filename << filepath << "controlTree_mg_VGBkg_" << "met" << metlow << "_" << methigh << "_pt" << leplow << "_" << lephigh << "_" << RunYear << whichVFP <<".root";
	TFile *file_VG = TFile::Open(filename.str().c_str());

	TH1D  *p_sig = (TH1D*)file_sig->Get("p_dPhiEleMET");
	TH1D  *p_rare = (TH1D*)file_rare->Get("p_dPhiEleMET");
	histname.str("");
	if(ih == 0)histname << "p_dPhiEleMET";
	else histname << "toy_eledPhiEleMET_" << ih;
	TH1D  *p_ele = (TH1D*)file_ele->Get(histname.str().c_str());
	TH1D  *p_jet = (TH1D*)file_jet->Get(histname.str().c_str());

	TH1D *p_rawsig = (TH1D*)p_sig->Clone("p_rawsig");
	TH1D *p_target = (TH1D*)p_sig->Clone("fit_target");
	//  CHECK !!!!
	//  JetFakePhoton is thrice for menglei compare to me in the Control region
	cout<< "\tRare : "<< p_rare->Integral()/p_target->Integral() <<"\t efakePho : "<<p_ele->Integral()/p_target->Integral() <<"\t jetfakePho : "<<p_jet->Integral()/p_target->Integral() << "\n";

	//  Target distribution is deltaPhi shape of the data in the CR, with jet fake photon, ele fake photon, rare EWK backgrounds subtracted.
	p_target->Add(p_rare, -1);
	p_target->Add(p_ele, -1);
	p_target->Add(p_jet, -1);
	p_target->Sumw2();
	//fake lepton template
  TH1D *p_proxy = (TH1D*)file_qcd->Get("p_dPhiEleMET");
  TH1D *p_MC;
	if(ih < 1000){
		// VGamma template
		p_MC = (TH1D*)file_VG->Get("p_dPhiEleMET"); 
		gRandom = new TRandom3(0);
		gRandom->SetSeed(0);
		double radomMC = gRandom->Gaus(0, 1);
		if(ih == 0)radomMC = 0;
		for(int ibin(1); ibin < p_MC->GetSize(); ibin++){
			radomMC =  gRandom->Gaus(0, 1);
			// VGamma get error from random number times error
			if(ih != 0)p_MC->SetBinContent(ibin, p_MC->GetBinContent(ibin)+radomMC*p_MC->GetBinError(ibin));
		}
	}
	else{
		histname.str("");
		histname << "toy_VGdPhiEleMET_" << ih-500;
		p_MC = (TH1D*)file_VG->Get(histname.str().c_str());
	}
	p_MC->Sumw2();

  RooRealVar dphi("dphi","",0,3.2);
	// RooDataHist for target and two fitting templates
  RooDataHist* h_target = new RooDataHist("h_target","h_target",dphi,p_target);
  RooDataHist* h_proxy  = new RooDataHist("h_proxy", "h_proxy", dphi,p_proxy );
  RooDataHist* h_MC     = new RooDataHist("h_MC",    "h_MC",    dphi,p_MC);
	// RooHistPdf for fake lepton and VGamma
  RooHistPdf*  pdf_proxy= new RooHistPdf("pdf_proxy","pdf_proxy",dphi, *h_proxy);
  RooHistPdf*  pdf_MC   = new RooHistPdf("pdf_MC",   "pdf_MC",   dphi, *h_MC);
  RooRealVar nFake("nFake","nFake", 0.5*p_target->Integral(1,32), 0, p_target->Integral(1,32)*2);
  RooRealVar nMC("nMC","nMC",0.5*p_target->Integral(1,32),0, p_target->Integral(1,32)*2);
  RooRealVar fakefrac("fakefrac","fakefrac",0.2,0,1.0);
  RooAddPdf model2("model2","",RooArgList(*pdf_proxy, *pdf_MC),RooArgList(nFake,nMC));
	// // combined pdf of fake lepton and VGamma
  RooAddPdf model("model","",RooArgList(*pdf_proxy, *pdf_MC),fakefrac);

  RooPlot* frame = dphi.frame(RooFit::Title("#Delta#phi(l,MET) (40 GeV < MET < 70 GeV)"));
  TCanvas *can=new TCanvas("can","",600,600);
  can->cd(1);
	// fitting dPhi of two templates to data
	model.fitTo(*h_target,RooFit::SumW2Error(kTRUE),RooFit::Save());
  RooFitResult* result = model.fitTo(*h_target,RooFit::SumW2Error(kTRUE),RooFit::Save());
  h_target->plotOn(frame);
  model.plotOn(frame,
		   RooFit::FillColor(kAzure-9));
  model.plotOn(frame, RooFit::Components(*pdf_proxy),
               RooFit::LineStyle(kDashed),
               RooFit::LineColor(kRed),
               RooFit::Normalization(1.0));
  model.plotOn(frame, RooFit::Components(*pdf_MC),
               RooFit::LineStyle(kDashed),
               RooFit::LineColor(kViolet),
               RooFit::Normalization(1.0));
  frame->Draw();
	std::ostringstream figurename;
  figurename.str("");
  figurename << "/eos/uscms/store/user/tmishra/VGamma/"<<RunYear<<whichVFP <<"/fit_dPhi_mg_met" << metlow << "_" << methigh << "_pt" <<  leplow << "_" << lephigh <<  "_iso" << isocut << ".png";
  //can->SaveAs(figurename.str().c_str());

	TH1D *p_combine = new TH1D("p_combine","",32,0,3.2);
	TH1D *p_qcd  = new TH1D("p_qcd","",32,0,3.2);
	TH1D *p_VGAMMA = new TH1D("p_VGAMMA","",32,0,3.2);
	// normalised template histograms
	TH1D *norm_MC = (TH1D*)p_MC->Clone("norm_MC");
	TH1D *norm_proxy = (TH1D*)p_proxy->Clone("norm_proxy");
	norm_MC->Scale(1.0/norm_MC->Integral(1, norm_MC->GetSize()) );
	norm_proxy->Scale(1.0/norm_proxy->Integral(1, norm_proxy->GetSize()));

	// unit normalised templates scaled by fakeFrac (from template fit) and combined
	p_combine->Add(norm_proxy, fakefrac.getVal());
	p_combine->Add(norm_MC, 1-fakefrac.getVal());
	// p_combine should of Integral 1, It is scaled to data template
	p_combine->Scale(p_target->Integral(1, p_target->GetSize()));
	for(int ibin(1); ibin <= 32; ibin++)p_combine->SetBinError(ibin, fakefrac.getError()*p_combine->GetBinContent(ibin));
	// unit normalised templates scaled to a fraction from fitting
	p_qcd->Add(norm_proxy, fakefrac.getVal()*p_target->Integral(1, p_target->GetSize()));
	p_VGAMMA->Add(norm_MC, (1-fakefrac.getVal())*p_target->Integral(1, p_target->GetSize()));
	TH1D *fitratio = new TH1D("fitratio",";#Delta#phi(l, E_{T}^{miss}) (radians);Fit/Data", 32,0,3.2);
	TGraphErrors *fitratio_error = new TGraphErrors(32);
	TGraphErrors *p_combine_error = new TGraphErrors(32);
  double staterror = fakefrac.getError()/fakefrac.getVal();
	for(int ibin(1); ibin <= 32; ibin++)fitratio->SetBinContent(ibin, p_combine->GetBinContent(ibin)/p_target->GetBinContent(ibin));
	for(int ibin(1); ibin <= 32; ibin++)fitratio->SetBinError(ibin, p_target->GetBinError(ibin)/p_target->GetBinContent(ibin));
	for(int ibin(0); ibin < 32; ibin++)fitratio_error->SetPoint(ibin, 0.05+ibin*0.1, 1);
	for(int ibin(0); ibin < 32; ibin++)fitratio_error->SetPointError(ibin, 0.05, staterror);
	for(int ibin(0); ibin < 32; ibin++)p_combine_error->SetPoint(ibin, -0.05+ibin*0.1, p_combine->GetBinContent(ibin));
	for(int ibin(0); ibin < 32; ibin++)p_combine_error->SetPointError(ibin, 0.05, staterror*p_combine->GetBinContent(ibin));

  TCanvas *canres=new TCanvas("canres","",600,600);
	setCanvas(canres); 
  canres->cd();
	TPad *canpt_pad1 = new TPad("canpt_pad1", "pad1", 0, 0.35, 1, 1.0);
	setTopPad(canpt_pad1); 
	canpt_pad1->Draw();          
	canpt_pad1->cd();  
	p_target->GetYaxis()->SetTitleOffset(1.4);
	p_target->GetYaxis()->SetTitle("Events / 0.1 radians");
	p_target->SetTitle(""); 
	p_target->SetMaximum(1.5*p_target->GetBinContent(p_target->GetMaximumBin()));
	p_target->SetMinimum(1);
	p_target->SetMarkerStyle(20);
	p_target->Draw("EP");
	p_combine->SetLineColor(kAzure-9);
	p_combine->SetLineWidth(0);
	p_combine->SetFillColor(kAzure-9);
	p_combine->Draw("hist same");
	p_combine_error->SetFillColor(kBlack);
	p_combine_error->SetFillStyle(3345);
	p_combine_error->Draw("E2");
	p_target->Draw("EP same");
	p_qcd->SetLineColor(kRed);
	p_qcd->SetLineWidth(3);
	p_qcd->Draw("hist same");
	p_VGAMMA->SetLineColor(kGreen+2);
	p_VGAMMA->SetLineStyle(2);
	p_VGAMMA->SetLineWidth(5);
	p_VGAMMA->Draw("hist same");
	
	TLegend *leg=new TLegend(0.45,0.65,0.9,0.9);
	leg->SetNColumns(2);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	p_combine->SetMarkerSize(0);
	p_VGAMMA->SetMarkerSize(0);
	p_qcd->SetMarkerSize(0);
	p_combine_error->SetMarkerSize(0);
	p_combine_error->SetLineWidth(0);
	leg->AddEntry(p_target, "Data", "ep");
	leg->AddEntry(p_combine,"Total fit");
	leg->AddEntry(p_VGAMMA, "V#gamma","l");
	leg->AddEntry(p_qcd, "Misid. #mu","l");
	leg->AddEntry(p_combine_error, "Fit uncertainty");
	leg->Draw("same");
 	gPad->RedrawAxis();
	if(RunYear==2016 and preVFP == 1)       CMS_lumi( canpt_pad1,1,2, 11 );
        else if(RunYear==2016 and preVFP == 0)  CMS_lumi( canpt_pad1,2,2, 11 );
        else if(RunYear==2017)                  CMS_lumi( canpt_pad1,3,2, 11 );
        else if(RunYear==2018)                  CMS_lumi( canpt_pad1,4,2, 11 );

  TLatex chantex;
  chantex.SetNDC();
  chantex.SetTextFont(42);
  chantex.SetTextSize(0.07);
  chantex.DrawLatex(0.2,0.7," #mu#gamma");
 	gPad->RedrawAxis();

	canres->cd();   
	TPad *canpt_pad2 = new TPad("canpt_pad2", "pad2", 0, 0, 1, 0.35);
	canpt_pad2->SetBottomMargin(0.3);
	canpt_pad2->Draw();
	canpt_pad2->cd(); 	
	TH1D *dummy_ptratio = new TH1D("dummy_ptratio",";|#Delta#phi(#mu, #vec{p}_{T}^{ miss})| (radians);#frac{Fit}{Data}",32,0,3.2);
	dummy_ptratio->SetMaximum(1.3);
	dummy_ptratio->SetMinimum(0.7);
	dummy_ptratio->GetYaxis()->SetNdivisions(504);
	dummy_ptratio->Draw();
  TLine *flatratio = new TLine(0,1,3.2,1);
	flatratio->Draw("same");
	fitratio->SetMarkerStyle(20);
	fitratio->Draw("EP same");
	fitratio_error->SetFillColor(kBlack);
	fitratio_error->SetFillStyle(3345);
	fitratio_error->Draw("E2 same");	

	if(ih == 0 and leplow==0 and lephigh==1000) canres->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/%d%s/fit_dPhi_mg_%d%s.png",RunYear,whichVFP.c_str(),RunYear,whichVFP.c_str()));	
	//else{
	//	std::ostringstream savename;
	//	savename.str("");
	//	savename << "/eos/uscms/store/user/tmishra/VGamma/fit_dPhi_mg_sys" << ih << ".png";
	//	canres->SaveAs(savename.str().c_str());
	//}

  double nMCtotal(0);
  for(unsigned i(1); i<=32; i++)nMCtotal+= p_MC->GetBinContent(i);
  std::cout << "factor for fake: " << fakefrac.getVal()*p_target->Integral(1, p_target->GetSize()) << "/" << p_proxy->Integral(1, p_proxy->GetSize()) << std::endl;
  std::cout << "factor for MC: " << (1-fakefrac.getVal())*p_target->Integral(1, p_target->GetSize()) << "/" << p_MC->Integral(1, p_MC->GetSize())  << std::endl;
  std::cout << "fakefrac=" << fakefrac.getVal()<< std::endl;
	
	float fakescale = fakefrac.getVal()*p_target->Integral(1, p_target->GetSize())/p_proxy->Integral(1, p_proxy->GetSize());
	float fakescaleerror = fakefrac.getError()*p_target->Integral(1, p_target->GetSize())/p_proxy->Integral(1, p_proxy->GetSize());
	float vgammascale = (1-fakefrac.getVal())*p_target->Integral(1, p_target->GetSize())/p_MC->Integral(1, p_MC->GetSize());
	float vgammascaleerror = fakefrac.getError()*p_target->Integral(1, p_target->GetSize())/p_MC->Integral(1, p_MC->GetSize()); 
	myfile << leplow << " " << lephigh << " " << fakescale << " " << fakescaleerror  << " " << vgammascale << " " << vgammascaleerror << std::endl; 
	cout << leplow << " " << lephigh << " " << fakescale << " " << fakescaleerror  << " " << vgammascale << " " << vgammascaleerror << std::endl;
	myfile.close();

  return 1;
}
