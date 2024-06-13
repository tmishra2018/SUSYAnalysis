#include "TH2D.h"
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
#include "RooMultiVarGaussian.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "TGraphErrors.h"
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

float analysis_PU(unsigned i, unsigned j, TH2D *p_lowPU_t5wg_pass, TH2D *p_lowPU_t5wg_all, TH2D *p_highPU_t5wg_pass, TH2D *p_highPU_t5wg_all, TH1D *p_PU_data){ 
	
	std::ostringstream histname;
	histname.str("");
	histname << "T5WG_" << int(p_lowPU_t5wg_all->GetXaxis()->GetBinCenter(i)) << "_" << int(p_lowPU_t5wg_all->GetYaxis()->GetBinCenter(j));
  float pass_L = p_lowPU_t5wg_pass->GetBinContent(i,j);
  float fail_L = p_lowPU_t5wg_all->GetBinContent(i,j) - pass_L;
  float tot_L  = pass_L + fail_L;
  float pass_H = p_highPU_t5wg_pass->GetBinContent(i,j);
  float fail_H = p_highPU_t5wg_all->GetBinContent(i,j) - pass_H;
  float tot_H  = pass_H + fail_H;
  float pass_A = pass_L + pass_H;
  float fail_A = fail_L + fail_H;
  float tot_A = pass_A + fail_A;
  
  float errpass_L = sqrt(pass_L);
  float errfail_L = sqrt(fail_L);
  float errpass_H = sqrt(pass_H);
  float errfail_H = sqrt(fail_H);
  
  
  double eff_L = (pass_L/tot_L)/(pass_A/tot_A);
  double eff_H = (pass_H/tot_H)/(pass_A/tot_A);
  double par_L_toPL = (fail_A*(fail_L + 2*pass_L + pass_H) + fail_H*(fail_A - pass_A*pass_A))/( tot_L*tot_L*pass_A*pass_A);
  double par_L_toFL = -1*pass_L*tot_H/( tot_L*tot_L*pass_A);
  double par_L_toPH = -1*pass_L*fail_A/( tot_L*pass_A*pass_A);
  double par_L_toFH = pass_L/( tot_L*pass_A);
  
  double par_H_toPH = (fail_A*(fail_H + 2*pass_H + pass_L) + fail_L*(fail_A - pass_A*pass_A))/( tot_H*tot_H*pass_A*pass_A);
  double par_H_toFH = -1*pass_H*tot_L/( tot_H*tot_H*pass_A);
  double par_H_toPL = -1*pass_H*fail_A/( tot_H*pass_A*pass_A);
  double par_H_toFL = pass_H/( tot_H*pass_A);
  
  double unc_L = sqrt( par_L_toPL*par_L_toPL*errpass_L*errpass_L + par_L_toFL*par_L_toFL*errfail_L*errfail_L + par_L_toPH*par_L_toPH*errpass_H*errpass_H + par_L_toFH*par_L_toFH*errfail_H*errfail_H);	
  double unc_H = sqrt( par_H_toPH*par_H_toPH*errpass_H*errpass_H + par_H_toFH*par_H_toFH*errfail_H*errfail_H + par_H_toPL*par_H_toPL*errpass_L*errpass_L + par_H_toFL*par_H_toFL*errfail_L*errfail_L);	
  
  TGraphErrors *p_eff = new TGraphErrors(2);
	//std::cout << eff_L << " " << unc_L << " " << eff_H << " " << unc_H << std::endl;
  p_eff->SetPoint(1,10,eff_L);
  p_eff->SetPointError(1, 0, unc_L);
  p_eff->SetPoint(2,30,eff_H);
  p_eff->SetPointError(2,0,unc_H);
  p_eff->SetMarkerStyle(20);
	TF1 *func = new TF1("linear","pol1",0,100);
	func->SetParameter(1, (eff_H-eff_L)/20);
  TFitResultPtr rPU = p_eff->Fit("linear","S");
  TH1F *p_effband = new TH1F(histname.str().c_str(),"",100,0,100); 
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(p_effband);
  p_effband->SetFillColor(kYellow);
  p_effband->Draw("E3 same");

	TVectorD pufit(2) ;
	pufit(0) = rPU->Parameter(0); 
	pufit(1) = rPU->Parameter(1);
	TMatrixDSym covden = rPU->GetCovarianceMatrix(); 
	RooRealVar offset("offset","offset",pufit(0)-rPU->ParError(0), pufit(0)+rPU->ParError(0));
	RooRealVar slope("slope","slope",pufit(1)-rPU->ParError(1), pufit(1)+rPU->ParError(1));
	RooMultiVarGaussian mvgpu("mvgpu","mvgpu",RooArgList(offset,slope),pufit,covden);
	RooDataSet* toymc = mvgpu.generate(RooArgSet(offset,slope),10000);
	std::ostringstream modelnameden;
	TH1D *band_upper = new TH1D("band_upper","band_upper",100,0,100);
	TH1D *band_lower = new TH1D("band_lower","band_lower",100,0,100);
	for(unsigned ibin(1); ibin <= 100; ibin++)band_lower->SetBinContent(ibin,2);
	for(int i(0); i<10000; i++){
		double data1 = toymc->get(i)->getRealValue("offset");
		double data2 = toymc->get(i)->getRealValue("slope");
		for(unsigned ibin(1); ibin <= 100; ibin++){
			double estimated = data1 + data2*ibin; 
			if(band_upper->GetBinContent(ibin) < estimated)band_upper->SetBinContent(ibin, estimated);
			if(band_lower->GetBinContent(ibin) > estimated)band_lower->SetBinContent(ibin, estimated);
		}
	}
	
	band_upper->Draw("L same");
	band_lower->Draw("L same");
  
  double c_center(0), c_up(0), c_down(0);
  for(unsigned i(1); i <=100; i++){
  	c_center += p_effband->GetBinContent(i)*p_PU_data->GetBinContent(i);
  	c_up += (band_upper->GetBinContent(i))*p_PU_data->GetBinContent(i);
  }
  
	TCanvas *can = new TCanvas(histname.str().c_str(),"",600,600);
	can->cd();
	p_eff->Draw("AP");
	p_effband->Draw("E3 same");
	p_eff->Draw("P same");
  p_eff->Fit("linear");
	band_upper->Draw("L same");
	band_lower->Draw("L same");
	histname.str("");
	histname << "T5WG_" << int(p_lowPU_t5wg_all->GetXaxis()->GetBinCenter(i)) << "_" << int(p_lowPU_t5wg_all->GetYaxis()->GetBinCenter(j)) << ".pdf";
	//can->SaveAs(histname.str().c_str());	
	//std::cout << i << " " << j  << " " << c_center << " " << c_up << " " << (c_up-c_center)/c_center <<  std::endl;
  delete band_upper;
  delete band_lower;
	return (c_up-c_center)/c_center;	
}


float analysis_PU(unsigned i, TH1D *p_lowPU_t5wg_pass, TH1D *p_lowPU_t5wg_all, TH1D *p_highPU_t5wg_pass, TH1D *p_highPU_t5wg_all, TH1D *p_PU_data){ 
	
	std::ostringstream histname;
	histname.str("");
	histname << "TChiWG_" << int(p_lowPU_t5wg_all->GetXaxis()->GetBinCenter(i));
  float pass_L = p_lowPU_t5wg_pass->GetBinContent(i);
  float fail_L = p_lowPU_t5wg_all->GetBinContent(i) - pass_L;
  float tot_L  = pass_L + fail_L;
  float pass_H = p_highPU_t5wg_pass->GetBinContent(i);
  float fail_H = p_highPU_t5wg_all->GetBinContent(i) - pass_H;
  float tot_H  = pass_H + fail_H;
  float pass_A = pass_L + pass_H;
  float fail_A = fail_L + fail_H;
  float tot_A = pass_A + fail_A;
  
  float errpass_L = sqrt(pass_L);
  float errfail_L = sqrt(fail_L);
  float errpass_H = sqrt(pass_H);
  float errfail_H = sqrt(fail_H);
  
  
  double eff_L = (pass_L/tot_L)/(pass_A/tot_A);
  double eff_H = (pass_H/tot_H)/(pass_A/tot_A);
  double par_L_toPL = (fail_A*(fail_L + 2*pass_L + pass_H) + fail_H*(fail_A - pass_A*pass_A))/( tot_L*tot_L*pass_A*pass_A);
  double par_L_toFL = -1*pass_L*tot_H/( tot_L*tot_L*pass_A);
  double par_L_toPH = -1*pass_L*fail_A/( tot_L*pass_A*pass_A);
  double par_L_toFH = pass_L/( tot_L*pass_A);
  
  double par_H_toPH = (fail_A*(fail_H + 2*pass_H + pass_L) + fail_L*(fail_A - pass_A*pass_A))/( tot_H*tot_H*pass_A*pass_A);
  double par_H_toFH = -1*pass_H*tot_L/( tot_H*tot_H*pass_A);
  double par_H_toPL = -1*pass_H*fail_A/( tot_H*pass_A*pass_A);
  double par_H_toFL = pass_H/( tot_H*pass_A);
  
  double unc_L = sqrt( par_L_toPL*par_L_toPL*errpass_L*errpass_L + par_L_toFL*par_L_toFL*errfail_L*errfail_L + par_L_toPH*par_L_toPH*errpass_H*errpass_H + par_L_toFH*par_L_toFH*errfail_H*errfail_H);	
  double unc_H = sqrt( par_H_toPH*par_H_toPH*errpass_H*errpass_H + par_H_toFH*par_H_toFH*errfail_H*errfail_H + par_H_toPL*par_H_toPL*errpass_L*errpass_L + par_H_toFL*par_H_toFL*errfail_L*errfail_L);	
  
  TGraphErrors *p_eff = new TGraphErrors(2);
  p_eff->SetPoint(1,10,eff_L);
  p_eff->SetPointError(1, 0, unc_L);
  p_eff->SetPoint(2,30,eff_H);
  p_eff->SetPointError(2,0,unc_H);
  p_eff->SetMarkerStyle(20);
	TF1 *func = new TF1("linear","pol1",0,100);
	func->SetParameter(1, (eff_H-eff_L)/20);
  TFitResultPtr rPU = p_eff->Fit("linear","S");
  TH1F *p_effband = new TH1F(histname.str().c_str(),"",100,0,100); 
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(p_effband);
  p_effband->SetFillColor(kYellow);
  p_effband->Draw("E3 same");

	TVectorD pufit(2) ;
	pufit(0) = rPU->Parameter(0); 
	pufit(1) = rPU->Parameter(1);
	TMatrixDSym covden = rPU->GetCovarianceMatrix(); 
	RooRealVar offset("offset","offset",pufit(0)-rPU->ParError(0), pufit(0)+rPU->ParError(0));
	RooRealVar slope("slope","slope",pufit(1)-rPU->ParError(1), pufit(1)+rPU->ParError(1));
	RooMultiVarGaussian mvgpu("mvgpu","mvgpu",RooArgList(offset,slope),pufit,covden);
	RooDataSet* toymc = mvgpu.generate(RooArgSet(offset,slope),10000);
	std::ostringstream modelnameden;
	TH1D *band_upper = new TH1D("band_upper","band_upper",100,0,100);
	TH1D *band_lower = new TH1D("band_lower","band_lower",100,0,100);
	for(unsigned ibin(1); ibin <= 100; ibin++)band_lower->SetBinContent(ibin,2);
	for(int i(0); i<10000; i++){
		double data1 = toymc->get(i)->getRealValue("offset");
		double data2 = toymc->get(i)->getRealValue("slope");
		for(unsigned ibin(1); ibin <= 100; ibin++){
			double estimated = data1 + data2*ibin; 
			if(band_upper->GetBinContent(ibin) < estimated)band_upper->SetBinContent(ibin, estimated);
			if(band_lower->GetBinContent(ibin) > estimated)band_lower->SetBinContent(ibin, estimated);
		}
	}
	
	band_upper->Draw("L same");
	band_lower->Draw("L same");
  
  double c_center(0), c_up(0), c_down(0);
  for(unsigned i(1); i <=100; i++){
  	c_center += p_effband->GetBinContent(i)*p_PU_data->GetBinContent(i);
  	c_up += (band_upper->GetBinContent(i))*p_PU_data->GetBinContent(i);
  }
  
	TCanvas *can = new TCanvas(histname.str().c_str(),"",600,600);
	can->cd();
	p_eff->Draw("AP");
	p_effband->Draw("E3 same");
	p_eff->Draw("P same");
  p_eff->Fit("linear");
	band_upper->Draw("L same");
	band_lower->Draw("L same");
	histname.str("");
	histname << "T5WG_" << int(p_lowPU_t5wg_all->GetXaxis()->GetBinCenter(i));
  delete band_upper;
  delete band_lower;
	return (c_up-c_center)/c_center;	
}
