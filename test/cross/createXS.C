#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TF3.h"
#include "TH1F.h"
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

void createXS(){
	gROOT->SetBatch(kTRUE);
	std::ifstream T5WG_file("ULCrossSectionT5WG.txt");
	std::ifstream T6WG_file("ULCrossSectionT6WG.txt");
	std::ifstream TChiWG_file("newCrossSectionTChiWG.txt");

  	TFile *outputfile = TFile::Open("ULsusyCrossSection.root","RECREATE");
  	outputfile->cd();
	TH1D *p_XS_T5WG   = new TH1D("p_gluinoxSec","p_gluinoxSec",500,197.5,2697.5);
	TH1D *p_XS_T6WG   = new TH1D("p_squarkxSec","p_squarkxSec",500,197.5,2697.5);
	TH1D *p_XS_TChiWG = new TH1D("p_charginoSec","p_charginoSec",77,87.5,2012.5);

	double susymass(0);
	double xsvalue(0);
	double xserror(0);
	if(T5WG_file.is_open()){
  	for(int i(0); i<561; i++){ 
			T5WG_file >> susymass >> xsvalue >> xserror;
			int bin = p_XS_T5WG->FindBin(susymass);
			p_XS_T5WG->SetBinContent(bin, xsvalue);
			p_XS_T5WG->SetBinError(  bin, xserror);	
	  }
	}
	T5WG_file.close();

	if(T6WG_file.is_open()){
  	for(int i(0); i<500; i++){ 
			T6WG_file >> susymass >> xsvalue >> xserror;
			int bin = p_XS_T6WG->FindBin(susymass);
			p_XS_T6WG->SetBinContent(bin, xsvalue);
			p_XS_T6WG->SetBinError(  bin, xserror);	
	  }
	}
	T6WG_file.close();

	if(TChiWG_file.is_open()){
		for(int i(0); i<77; i++){
			TChiWG_file >> susymass >> xsvalue >> xserror;
			int bin = p_XS_TChiWG->FindBin(susymass);
			p_XS_TChiWG->SetBinContent(bin, xsvalue);
			p_XS_TChiWG->SetBinError(  bin, xserror);
		}
	}
	TChiWG_file.close();

	outputfile->Write();
	outputfile->Close();
}

