#include "../../../include/analysis_commoncode.h"
#include "TProfile2D.h"
#define NBIN 18
#define NCH  4

int DIM;
int NBIN_2DX;
int NBIN_2DY;
double XBIN_low;
double XBIN_high;
double YBIN_low;
double YBIN_high;

TFile *outputfile;
TH2D *p_SUSYMass;
TH2D *p_SUSYGEN[NCH];
TH2D *p_SUSYISR;
TH2D *p_SUSYMASS_pdf[9];
TH2D *p_lowPU_pass;
TH2D *p_lowPU_all;
TH2D *p_highPU_pass;
TH2D *p_highPU_all;

TH2D *susy_chan_rate_nom[NCH][NBIN*2]; 
TH2D *susy_chan_rate_control[NCH][NBIN*2]; 
TH2D *susy_chan_rate_jesUp[NBIN*2]; 
TH2D *susy_chan_rate_jesDown[NBIN*2];    
TH2D *susy_chan_rate_isr[NBIN*2];    
TH2D *susy_chan_rate_phoesfUp[NBIN*2];       
TH2D *susy_chan_rate_lepesfUp[NBIN*2];       
TH2D *susy_chan_rate_pdf[NBIN*2][9]; 
TH2D *susy_chan_rate_genMET[NBIN*2];
                                
TH2D *susy_syserr_PU;
TH2D *susy_chan_syserr_jes[NBIN*2];      
TH2D *susy_chan_syserr_phoesf[NBIN*2];     
TH2D *susy_chan_syserr_lepesf[NBIN*2];     
TH2D *susy_chan_syserr_scale[NBIN*2];   
TH2D *susy_chan_syserr_eleshape[NBIN*2]; 
TH2D *susy_chan_syserr_jetshape[NBIN*2];
TH2D *susy_chan_syserr_qcdshape[NBIN*2];
TH2D *susy_chan_syserr_xs[NBIN*2];
TH2D *susy_chan_syserr_lumi[NBIN*2];     
TH2D *susy_chan_syserr_isr[NBIN*2];     
TH2D *susy_chan_syserr_pdf[NBIN*2];
TH2D *susy_chan_syserr_met[NBIN*2];

void init_histo(TString modelname, TString filename, int dim, int nxbin, double xbinlow, double xbinhigh, int nybin, double ybinlow, double ybinhigh){
	
	DIM = dim;	
	NBIN_2DX = nxbin;
	NBIN_2DY = nybin;
	XBIN_low = xbinlow;
	XBIN_high = xbinhigh;
	YBIN_low = ybinlow;
	YBIN_high = ybinhigh;
	
	std::ostringstream outputfilename;
	outputfilename.str("");
	outputfilename << filename;
	outputfile = TFile::Open(outputfilename.str().c_str(),"RECREATE");
	outputfile->cd();
	
	p_SUSYMass         = new TH2D("SUSYMass","",NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
	p_SUSYISR          = new TH2D("SUSYISR","",NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
	p_lowPU_pass       = new TH2D("p_lowPU_pass","", NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
	p_lowPU_all        = new TH2D("p_lowPU_all", "", NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
	p_highPU_pass      = new TH2D("p_highPU_pass","",NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
	p_highPU_all       = new TH2D("p_highPU_all","", NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
	
	std::ostringstream histname;
	
	for(unsigned i(0); i < 9; i++){
		histname.str("");
		histname << "p_SUSYMASS_pdf" << i;
		p_SUSYMASS_pdf[i] = new TH2D(histname.str().c_str(), histname.str().c_str(), NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
	}
	for(unsigned ich(0); ich < NCH; ich++){

  	histname.str("");
		if(ich == 0)histname << "ALL_SUSYGEN" << std::endl;
		else if(ich == 1)histname << "WG_SUSYGEN" << std::endl;
		else if(ich == 2)histname << "GG_SUSYGEN" << std::endl;
		else if(ich == 3)histname << "WW_SUSYGEN" << std::endl;
		p_SUSYGEN[ich]= new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
	}

	susy_syserr_PU 	 = new TH2D("susy_syserr_PU","susy_syserr_PU", NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);

	for(unsigned ich(0); ich < NCH; ich++){

		TString channame;
		if(ich == 0)channame = "ALL_";
		else if(ich == 1)channame = "WG_";
		else if(ich == 2)channame = "GG_";
		else if(ich == 3)channame = "WW_";

  	for(unsigned i(0); i < NBIN*2; i++){
  		histname.str("");
  		histname << channame << "chan" << i+1 << "_rate_nom";
  		susy_chan_rate_nom[ich][i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
  		histname.str("");
  		histname << channame << "chan" << i+1 << "_rate_control";
			susy_chan_rate_control[ich][i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		}
	}


  for(unsigned i(0); i < NBIN*2; i++){
	  for(unsigned is(0); is<9; is++){
	  	histname.str("");
	  	histname << "susy_chan" << i+1 << "_rate_pdf" << is;
	  	susy_chan_rate_pdf[i][is] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high); 
	  }		

  	histname.str("");
  	histname << "susy_chan" << i+1 << "_rate_jesUp";
  	susy_chan_rate_jesUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
  	histname.str("");
  	histname << "susy_chan" << i+1 << "_rate_jesDown";
  	susy_chan_rate_jesDown[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
  	histname.str("");
  	histname << "susy_chan" << i+1 << "_rate_isr";
  	susy_chan_rate_isr[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
  	histname.str("");
  	histname << "susy_chan" << i+1 << "_rate_phoesfUp";
  	susy_chan_rate_phoesfUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
  	histname.str("");
  	histname << "susy_chan" << i+1 << "_rate_lepesfUp";
  	susy_chan_rate_lepesfUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
  	histname.str("");
  	histname << "susy_chan" << i+1 << "_rate_genMET";
  	susy_chan_rate_genMET[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
  															
  	histname.str("");
  	histname << "susy_chan" << i+1 << "_syserr_jes";
  	susy_chan_syserr_jes[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
  	histname.str("");
  	histname << "susy_chan" << i+1 << "_syserr_phoesf";
  	susy_chan_syserr_phoesf[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
  	histname.str("");
  	histname << "susy_chan" << i+1 << "_syserr_lepesf";
  	susy_chan_syserr_lepesf[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
  	histname.str("");
  	histname << "susy_chan" << i+1 << "_syserr_scale";
  	susy_chan_syserr_scale[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
  	histname.str("");
  	histname << "susy_chan" << i+1 << "_syserr_eleshape";
  	susy_chan_syserr_eleshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
  	histname.str("");
  	histname << "susy_chan" << i+1 << "_syserr_jetshape";
  	susy_chan_syserr_jetshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
  	histname.str("");
  	histname << "susy_chan" << i+1 << "_syserr_qcdshape";
  	susy_chan_syserr_qcdshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
  	histname.str("");
  	histname << "susy_chan" << i+1 << "_syserr_xs";
  	susy_chan_syserr_xs[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
  	histname.str("");
  	histname << "susy_chan" << i+1 << "_syserr_lumi";
  	susy_chan_syserr_lumi[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
  	histname.str("");
  	histname << "susy_chan" << i+1 << "_syserr_isr";
  	susy_chan_syserr_isr[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
  	histname.str("");
  	histname << "susy_chan" << i+1 << "_syserr_pdf";
  	susy_chan_syserr_pdf[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high); 
  	histname.str("");
  	histname << "susy_chan" << i+1 << "_syserr_met";
  	susy_chan_syserr_met[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high); 
  }

}


void write_histo(){		
	outputfile->Write();
	outputfile->Close();
}









//	p_SUSYselect       = new TH2D("p_SUSYselect","",NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high); 
//  TH2D *p_SUSYselect;
//  
//  TH1D *p_MET_signal_mg;
//  TH1D *p_HT_signal_mg;
//  TH1D *p_PhoEt_signal_mg;
//  TH1D *p_MET_signal_eg;
//  TH1D *p_HT_signal_eg;
//  TH1D *p_PhoEt_signal_eg;
//  
//  		
//  	p_MET_signal_mg = new TH1D("p_MET_signal_1700_1000_mg","", nSigMETBins, sigMETBins);
//  	p_HT_signal_mg  = new TH1D("p_HT_signal_1700_1000_mg","",  nSigHTBins, sigHTBins);
//  	p_PhoEt_signal_mg = new TH1D("p_PhoEt_signal_1700_1000_mg","",nSigEtBins, sigEtBins);
//  	p_MET_signal_eg = new TH1D("p_MET_signal_1700_1000_eg","", nSigMETBins, sigMETBins);
//  	p_HT_signal_eg  = new TH1D("p_HT_signal_1700_1000_eg","",  nSigHTBins, sigHTBins);
//  	p_PhoEt_signal_eg = new TH1D("p_PhoEt_signal_1700_1000_eg","",nSigEtBins, sigEtBins);
//  
//  }
