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
#include "TLatex.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "../../../include/tdrstyle.C"
void plot_sigmap(int NBIN){//main  

	std::ostringstream susyname;
	gStyle->SetOptStat(0);
	setTDRStyle();
  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

	TFile *susy_sig = TFile::Open("signalTree_T5WG.root");
	TH2D *h_SUSY[36];
	for(unsigned i(0); i < 36; i++){
		susyname.str("");
		susyname << "h_chan" << i+1 << "_rate_nom";
		h_SUSY[i] = (TH2D*)susy_sig->Get(susyname.str().c_str());
	}
	TH2D *p_channel_map  = new TH2D("p_channel_map","", 27, 775.0, 2125.0, 80, 12.5, 2012.5);

	TFile *tchiwg_sig = TFile::Open("signalTree_TChiWG.root");
	TH1D  *h_tchiwg[36];
	for(unsigned i(0); i < 36; i++){
		susyname.str("");
		susyname << "h_chan" << i+1 << "_rate_nom";
		h_tchiwg[i] = (TH1D*)tchiwg_sig->Get(susyname.str().c_str());
	}
	TH1D *sens_map[36];
	for(unsigned i(0); i < 36; i++){
    susyname.str("");
    susyname << "sens_map_" << i+1;
		sens_map[i] = new TH1D(susyname.str().c_str(), susyname.str().c_str(), 40,287.5,1287.5);
	}
	TCanvas *can[36];
	for(unsigned i(0); i < 36; i++){
		susyname.str("");
		susyname << "sens_chan_" << i+1;
		can[i] = new TCanvas(susyname.str().c_str(), susyname.str().c_str(), 600,600);
	} 

	TFile *egfile_sig = TFile::Open("signalTree_egamma_signal.root");
	TFile *mgfile_sig = TFile::Open("signalTree_mg_signal.root");
	TH1D  *eg_sig     = (TH1D*)egfile_sig->Get("p_eventcount");
	TH1D  *mg_sig     = (TH1D*)mgfile_sig->Get("p_eventcount");
	TH1D  *h_sig      = new TH1D("h_sig","",NBIN*2,0,NBIN*2);

		h_sig->SetBinContent(1,	324.214);
		h_sig->SetBinContent(2,	467.461);
		h_sig->SetBinContent(3,	 98.865);
		h_sig->SetBinContent(4,	 27.271);
		h_sig->SetBinContent(5,	 61.307);
		h_sig->SetBinContent(6,	 45.150);
		h_sig->SetBinContent(7,	  1.236);
		h_sig->SetBinContent(8,	  2.914);
		h_sig->SetBinContent(9,	  5.311);
		h_sig->SetBinContent(10,		  6.538);
		h_sig->SetBinContent(11,		 21.444);
		h_sig->SetBinContent(12,		 15.580);
		h_sig->SetBinContent(13,		  5.012);
		h_sig->SetBinContent(14,		  8.368);
		h_sig->SetBinContent(15,		  5.459);
		h_sig->SetBinContent(16,		  0.766);
		h_sig->SetBinContent(17,		  0.634);
		h_sig->SetBinContent(18,		  0.539);
		h_sig->SetBinContent(19,		163.626);
		h_sig->SetBinContent(20,		256.955);
		h_sig->SetBinContent(21,		 78.950);
		h_sig->SetBinContent(22,		 17.309);
		h_sig->SetBinContent(23,		 49.956); 	
		h_sig->SetBinContent(24,		 28.164); 	
		h_sig->SetBinContent(25,		  1.302); 	
		h_sig->SetBinContent(26,		  1.160); 	
		h_sig->SetBinContent(27,		  2.804); 	
		h_sig->SetBinContent(28,		  5.927); 	
		h_sig->SetBinContent(29,		 21.614); 	
		h_sig->SetBinContent(30,		 11.677); 	
		h_sig->SetBinContent(31,		  4.135);
		h_sig->SetBinContent(32,		  8.876);
		h_sig->SetBinContent(33,		  5.123);
		h_sig->SetBinContent(34,		  0.394);
		h_sig->SetBinContent(35,		  0.532);
		h_sig->SetBinContent(36,		  0.848);


	for(unsigned ix(0); ix < 27; ix++){
		std::cout << "process " << ix << std::endl;
		for(unsigned iy(0); iy < 80; iy++){
			unsigned xbin = ix+1;
			unsigned ybin = iy+1;
			double stob = 0.00001;
			double stob2 = 0.00001;
			int    bestchan = -1;
			int    secondchan = -1;
			for(unsigned ih(0); ih < 36; ih++){
				if( h_SUSY[ih]->GetBinContent(xbin, ybin)/sqrt(h_sig->GetBinContent(ih+1)) > stob){
					stob = h_SUSY[ih]->GetBinContent(xbin, ybin)/sqrt(h_sig->GetBinContent(ih+1));
					bestchan = ih;
				}
				if( h_SUSY[ih]->GetBinContent(xbin, ybin)/sqrt(h_sig->GetBinContent(ih+1)) > stob2 && h_SUSY[ih]->GetBinContent(xbin, ybin)/sqrt(h_sig->GetBinContent(ih+1)) < stob){
					stob2 = h_SUSY[ih]->GetBinContent(xbin, ybin)/sqrt(h_sig->GetBinContent(ih+1));
					secondchan = ih;
				}
			}
			std::cout << bestchan << " " << secondchan << std::endl;
			p_channel_map->SetBinContent(xbin, ybin, bestchan);
		}
	}


	for(unsigned i(1); i <= 40 ;i++){
		for(unsigned ih(0); ih < 36; ih++){
			sens_map[ih]->SetBinContent(i, h_tchiwg[ih]->GetBinContent(i)/sqrt( h_sig->GetBinContent(ih+1)));
		}
	}

	p_channel_map->Draw("colz");
	for(unsigned i(0); i < 36; i++){
		susyname.str("");
		susyname << "./sens/sens_chan_" << i+1 << ".pdf";
		can[i]->cd();
		sens_map[i]->Draw();
		can[i]->SaveAs( susyname.str().c_str());
	} 
	

}
