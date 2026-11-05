#include "../../../include/analysis_commoncode.h"
#include "TProfile2D.h"
#include "analysis_PU.C"

void analysis_TChi(){//main  

	SetSignalConfig();
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);
	esfScaleFactor  objectESF;

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
	
	std::ifstream C1N1_file("../../cross/CrossSectionTChiWG.txt");
	std::ifstream C1C1_file("../../cross/CrossSectionC1C1.txt");

	std::map<int, float> xsmap;

	int susymass1(0),susymass2(0);
	double xsvalue1(0),xsvalue2(0);
	double xserror1(0),xserror2(0);
	if(C1N1_file.is_open() && C1C1_file.is_open()){
  	for(int i(0); i<60; i++){ 
			C1N1_file >> susymass1 >> xsvalue1 >> xserror1;
			C1C1_file >> susymass2 >> xsvalue2 >> xserror2;
			if(susymass1 != susymass2){
				std::cout << "wrong XS value" << std::endl;
				abort();
			}
			//xsmap.insert(make_pair(susymass1, xsvalue1+xsvalue2));
			xsmap.insert(make_pair(susymass1, xsvalue1));
	  }
	}
	C1N1_file.close();
	C1C1_file.close();

//	for(map<pair<int, int>, float>::const_iterator it = xsmap.begin(); it != xsmap.end(); ++it)
//    std::cout << it->first.first << " " << it->first.second << " " << it->second << " " << xsmap[make_pair(it->first.first,it->first.second)] << "\n";

	TChain *datachain = new TChain("signalTree");
	datachain->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_mgsignal_MuonEG_FullEcal.root");
	TH1D *p_PU_data = new TH1D("p_PU_data",";N_{vtx};",100,0,100); 
  datachain->Draw("nVertex >> p_PU_data");
	p_PU_data->Scale(1.0/p_PU_data->Integral(1,101));

	std::ostringstream histname;
	//**************   M1M2  ***************************//
  TFile *file_m1m2 = TFile::Open("/uscms_data/d3/mengleis/FullStatusOct/resTree_TChiNg_test.root");
  TTree *tree_m1m2 = (TTree*)file_m1m2->Get("SUSYtree");
	float Mass1(0);
  float Mass2(0);
	int   nVertex(0);
	int   branch(0);
	tree_m1m2->SetBranchAddress("Mass1",      &Mass1);
  tree_m1m2->SetBranchAddress("Mass2",      &Mass2);
  tree_m1m2->SetBranchAddress("nVertex",    &nVertex);
	tree_m1m2->SetBranchAddress("branch",     &branch);

	TFile *outputfile_m1m2 = TFile::Open("signalTree_TChiNg_C1N1.root","RECREATE");
	outputfile_m1m2->cd();

	TH2D *p_M1M2MASS         = new TH2D("SUSYMass","", 48,287.5,1487.5,10,0,10);
	TH2D *p_M1M2accept       = new TH2D("p_M1M2accept","",48,287.5,1487.5,10,0,10);
	TH2D *p_M1M2_XS          = new TH2D("p_M1M2_XS","", 48,287.5,1487.5,10,0,10);
	TH2D *p_lowPU_m1m2_pass  = new TH2D("p_lowPU_m1m2_pass","", 48,287.5,1487.5,10,0,10);
	TH2D *p_lowPU_m1m2_all   = new TH2D("p_lowPU_m1m2_all", "", 48,287.5,1487.5,10,0,10);
	TH2D *p_highPU_m1m2_pass = new TH2D("p_highPU_m1m2_pass","",48,287.5,1487.5,10,0,10);
	TH2D *p_highPU_m1m2_all  = new TH2D("p_highPU_m1m2_all","", 48,287.5,1487.5,10,0,10);

	for(unsigned ievt(0); ievt < tree_m1m2->GetEntries(); ievt++){
		tree_m1m2->GetEntry(ievt);
    if (ievt%1000000==0) std::cout << " -- Processing event " << ievt << std::endl;
    	
		p_M1M2MASS->Fill(Mass1,branch);
		if(nVertex < 20)p_lowPU_m1m2_all->Fill(Mass1,branch, 1);
		else p_highPU_m1m2_all->Fill(Mass1,branch, 1);
	}

	TH2D *m1m2_h_chan_rate_nom[NBIN*2]; 
	for(unsigned i(0); i < NBIN*2; i++){
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_nom";
		m1m2_h_chan_rate_nom[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),48,287.5,1487.5,10,0,10);
	}

	TH2D *m1m2_h_chan_rate_jesUp[NBIN*2]; 
	TH2D *m1m2_h_chan_rate_jesDown[NBIN*2];    
	TH2D *m1m2_h_chan_rate_jerUp[NBIN*2];    
	TH2D *m1m2_h_chan_rate_jerDown[NBIN*2];    
	TH2D *m1m2_h_chan_rate_xsUp[NBIN*2];       
	TH2D *m1m2_h_chan_rate_esfUp[NBIN*2];       
                                  
	TH2D *m1m2_h_syserr_PU = new TH2D("m1m2_h_syserr_PU","m1m2_h_syserr_PU", 48,287.5,1487.5,10,0,10); 
	TH2D *m1m2_h_chan_syserr_jes[NBIN*2];      
	TH2D *m1m2_h_chan_syserr_jer[NBIN*2];     
	TH2D *m1m2_h_chan_syserr_esf[NBIN*2];     
	TH2D *m1m2_h_chan_syserr_scale[NBIN*2];   
	TH2D *m1m2_h_chan_syserr_eleshape[NBIN*2]; 
	TH2D *m1m2_h_chan_syserr_jetshape[NBIN*2];
	TH2D *m1m2_h_chan_syserr_qcdshape[NBIN*2];
	TH2D *m1m2_h_chan_syserr_xs[NBIN*2];
	TH2D *m1m2_h_chan_syserr_lumi[NBIN*2];     
	TH2D *m1m2_h_chan_syserr_isr[NBIN*2];     
	TH2D *m1m2_h_chan_staterr[NBIN*2];     

	for(unsigned i(0); i < NBIN*2; i++){
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jesUp";
		m1m2_h_chan_rate_jesUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),48,287.5,1487.5,10,0,10);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jesDown";
		m1m2_h_chan_rate_jesDown[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),48,287.5,1487.5,10,0,10);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jerUp";
		m1m2_h_chan_rate_jerUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),48,287.5,1487.5,10,0,10);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jerDown";
		m1m2_h_chan_rate_jerDown[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),48,287.5,1487.5,10,0,10);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_xsUp";
		m1m2_h_chan_rate_xsUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),48,287.5,1487.5,10,0,10);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_esfUp";
		m1m2_h_chan_rate_esfUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),48,287.5,1487.5,10,0,10);
															
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jes";
		m1m2_h_chan_syserr_jes[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),48,287.5,1487.5,10,0,10);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jer";
		m1m2_h_chan_syserr_jer[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),48,287.5,1487.5,10,0,10);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_esf";
		m1m2_h_chan_syserr_esf[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),48,287.5,1487.5,10,0,10);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_scale";
		m1m2_h_chan_syserr_scale[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),48,287.5,1487.5,10,0,10);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_eleshape";
		m1m2_h_chan_syserr_eleshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),48,287.5,1487.5,10,0,10);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jetshape";
		m1m2_h_chan_syserr_jetshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),48,287.5,1487.5,10,0,10);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_qcdshape";
		m1m2_h_chan_syserr_qcdshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),48,287.5,1487.5,10,0,10);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_xs";
		m1m2_h_chan_syserr_xs[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),48,287.5,1487.5,10,0,10);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_lumi";
		m1m2_h_chan_syserr_lumi[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),48,287.5,1487.5,10,0,10);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_isr";
		m1m2_h_chan_syserr_isr[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),48,287.5,1487.5,10,0,10);
		histname.str("");
		histname << "h_chan" << i+1 << "_staterr";
		m1m2_h_chan_staterr[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),48,287.5,1487.5,10,0,10);
	}
	
	TH1D *p_m1m2_MET_signal_mg = new TH1D("p_m1m2_MET_signal_1700_1000_mg","", nSigMETBins, sigMETBins);
	TH1D *p_m1m2_HT_signal_mg  = new TH1D("p_m1m2_HT_signal_1700_1000_mg","",  nSigHTBins, sigHTBins);
	TH1D *p_m1m2_PhoEt_signal_mg = new TH1D("p_m1m2_PhoEt_signal_1700_1000_mg","",nSigEtBins, sigEtBins);
	TH1D *p_m1m2_MET_signal_eg = new TH1D("p_m1m2_MET_signal_1700_1000_eg","", nSigMETBins, sigMETBins);
	TH1D *p_m1m2_HT_signal_eg  = new TH1D("p_m1m2_HT_signal_1700_1000_eg","",  nSigHTBins, sigHTBins);
	TH1D *p_m1m2_PhoEt_signal_eg = new TH1D("p_m1m2_PhoEt_signal_1700_1000_eg","",nSigEtBins, sigEtBins);
		
  TChain *mgtree_m1m2;
  mgtree_m1m2 = new TChain("mgTree","mgTree");
  mgtree_m1m2->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_TChiNg_test.root");
  float phoEt_m1m2_mg(0);
  float phoEta_m1m2_mg(0);
  float lepPt_m1m2_mg(0);
  float lepEta_m1m2_mg(0);
  float sigMT_m1m2_mg(0);
  float sigMET_m1m2_mg(0);
  float HT_m1m2_mg(0);
	int   nVertex_m1m2_mg(0);
	float sigMETJESup_m1m2_mg(0);
	float sigMETJESdo_m1m2_mg(0);
	float sigMETJERup_m1m2_mg(0);
	float sigMETJERdo_m1m2_mg(0);
	float sigMTJESup_m1m2_mg(0);
	float sigMTJESdo_m1m2_mg(0);
	float sigMTJERup_m1m2_mg(0);
	float sigMTJERdo_m1m2_mg(0);
	float HTJESup_m1m2_mg(0);
	float HTJESdo_m1m2_mg(0);
	float Mass1_m1m2_mg(0); 
  float Mass2_m1m2_mg(0);
	int   branch_mg(0);
  mgtree_m1m2->SetBranchAddress("phoEt",      &phoEt_m1m2_mg);
  mgtree_m1m2->SetBranchAddress("phoEta",     &phoEta_m1m2_mg);
  mgtree_m1m2->SetBranchAddress("lepPt",      &lepPt_m1m2_mg);
  mgtree_m1m2->SetBranchAddress("lepEta",     &lepEta_m1m2_mg);
  mgtree_m1m2->SetBranchAddress("sigMT",      &sigMT_m1m2_mg);
  mgtree_m1m2->SetBranchAddress("sigMET",     &sigMET_m1m2_mg);
  mgtree_m1m2->SetBranchAddress("HT",         &HT_m1m2_mg);
	mgtree_m1m2->SetBranchAddress("nVertex",    &nVertex_m1m2_mg);
	mgtree_m1m2->SetBranchAddress("Mass1",      &Mass1_m1m2_mg);
  mgtree_m1m2->SetBranchAddress("Mass2",      &Mass2_m1m2_mg);
	mgtree_m1m2->SetBranchAddress("branch",     &branch_mg);
	mgtree_m1m2->SetBranchAddress("sigMETJESup",&sigMETJESup_m1m2_mg);
	mgtree_m1m2->SetBranchAddress("sigMETJESdo",&sigMETJESdo_m1m2_mg);
	mgtree_m1m2->SetBranchAddress("sigMETJERup",&sigMETJERup_m1m2_mg);
	mgtree_m1m2->SetBranchAddress("sigMETJERdo",&sigMETJERdo_m1m2_mg);
	mgtree_m1m2->SetBranchAddress("sigMTJESup", &sigMTJESup_m1m2_mg);
	mgtree_m1m2->SetBranchAddress("sigMTJESdo", &sigMTJESdo_m1m2_mg);
	mgtree_m1m2->SetBranchAddress("sigMTJERup", &sigMTJERup_m1m2_mg);
	mgtree_m1m2->SetBranchAddress("sigMTJERdo", &sigMTJERdo_m1m2_mg);
	mgtree_m1m2->SetBranchAddress("HTJESup",    &HTJESup_m1m2_mg);
	mgtree_m1m2->SetBranchAddress("HTJESdo",    &HTJESdo_m1m2_mg);

	for(unsigned ievt(0); ievt < mgtree_m1m2->GetEntries(); ievt++){
		mgtree_m1m2->GetEntry(ievt);
	
		/** cut flow *****/
		if(phoEt_m1m2_mg < 35 || lepPt_m1m2_mg < 25)continue;
		if(fabs(phoEta_m1m2_mg) > 1.4442 || fabs(lepEta_m1m2_mg) > 2.5)continue;

		double scalefactor(0);
		double scalefactorup(0);
		scalefactor = objectESF.getFastMuonESF(lepPt_m1m2_mg,lepEta_m1m2_mg)*objectESF.getPhotonESF(phoEt_m1m2_mg, phoEta_m1m2_mg)*objectESF.getFastMuonEGTRGESF(phoEt_m1m2_mg, lepPt_m1m2_mg);
		double s_mu_error = objectESF.getFastMuonESFError(lepPt_m1m2_mg,lepEta_m1m2_mg)*objectESF.getPhotonESF(phoEt_m1m2_mg, phoEta_m1m2_mg)*objectESF.getFastMuonEGTRGESF(phoEt_m1m2_mg, lepPt_m1m2_mg);
    double s_pho_error = objectESF.getPhotonESFError(phoEt_m1m2_mg, phoEta_m1m2_mg)*objectESF.getFastMuonESF(lepPt_m1m2_mg,lepEta_m1m2_mg)*objectESF.getFastMuonEGTRGESF(phoEt_m1m2_mg, lepPt_m1m2_mg);
    double s_trg_error = objectESF.getFastMuonEGTRGESFError(phoEt_m1m2_mg, lepPt_m1m2_mg)*objectESF.getFastMuonESF(lepPt_m1m2_mg,lepEta_m1m2_mg)*objectESF.getPhotonESF(phoEt_m1m2_mg, phoEta_m1m2_mg);
		double s_error = sqrt(pow(s_mu_error,2) + pow(s_pho_error,2) + pow(s_trg_error,2));
		scalefactorup = scalefactor + s_error; 

		if( fabs(Mass1_m1m2_mg - 900) < 10){ 
			if(sigMT_m1m2_mg > 100){
				p_m1m2_MET_signal_mg->Fill(sigMET_m1m2_mg);
				if(sigMET_m1m2_mg > 120)p_m1m2_HT_signal_mg->Fill(HT_m1m2_mg);
				if(sigMET_m1m2_mg > 120)p_m1m2_PhoEt_signal_mg->Fill(phoEt_m1m2_mg);
			}
		}

		if(sigMET_m1m2_mg < 120 || sigMT_m1m2_mg < 100)continue;

		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET_m1m2_mg, HT_m1m2_mg, phoEt_m1m2_mg);
		if(SigBinIndex >=0){
			m1m2_h_chan_rate_nom[SigBinIndex]->Fill( Mass1_m1m2_mg,branch_mg, scalefactor);
			m1m2_h_chan_rate_xsUp[SigBinIndex]->Fill( Mass1_m1m2_mg,branch_mg, scalefactor); 
			m1m2_h_chan_rate_esfUp[SigBinIndex]->Fill( Mass1_m1m2_mg,branch_mg, scalefactorup); 


			if(nVertex_m1m2_mg < 20)p_lowPU_m1m2_pass->Fill( Mass1_m1m2_mg,branch_mg, 1);
			else p_highPU_m1m2_pass->Fill( Mass1_m1m2_mg,branch_mg, 1);
		}
		int jesupBinIndex(-1);	
		jesupBinIndex = Bin.findSignalBin(sigMETJESup_m1m2_mg, HTJESup_m1m2_mg, phoEt_m1m2_mg); 
		if(jesupBinIndex >=0){
			m1m2_h_chan_rate_jesUp[jesupBinIndex]->Fill( Mass1_m1m2_mg,branch_mg, scalefactor); 
		}
		int jesdoBinIndex(-1);
		jesdoBinIndex = Bin.findSignalBin(sigMETJESdo_m1m2_mg, HTJESdo_m1m2_mg, phoEt_m1m2_mg);
		if(jesdoBinIndex >=0){
			m1m2_h_chan_rate_jesDown[jesdoBinIndex]->Fill( Mass1_m1m2_mg,branch_mg, scalefactor);   
		}	
		int jerupBinIndex(-1);
		jerupBinIndex = Bin.findSignalBin(sigMETJERup_m1m2_mg, HT_m1m2_mg, phoEt_m1m2_mg); 
		if( jerupBinIndex >=0){
			m1m2_h_chan_rate_jerUp[jerupBinIndex]->Fill( Mass1_m1m2_mg,branch_mg, scalefactor);
		}
		int jerdoBinIndex(-1);
		jerdoBinIndex = Bin.findSignalBin(sigMETJERdo_m1m2_mg, HT_m1m2_mg, phoEt_m1m2_mg);
		if(jerdoBinIndex >= 0){
			m1m2_h_chan_rate_jerDown[jerdoBinIndex]->Fill( Mass1_m1m2_mg,branch_mg, scalefactor);
		}  
	}

  TChain *egtree_m1m2;
  egtree_m1m2 = new TChain("egTree","egTree");
  egtree_m1m2->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_TChiNg_test.root");
  float phoEt_m1m2_eg(0);
  float phoEta_m1m2_eg(0);
  float lepPt_m1m2_eg(0);
  float lepEta_m1m2_eg(0);
  float sigMT_m1m2_eg(0);
  float sigMET_m1m2_eg(0);
  float HT_m1m2_eg(0);
  int   nVertex_m1m2_eg(0); 
	float sigMETJESup_m1m2_eg(0);
	float sigMETJESdo_m1m2_eg(0);
	float sigMETJERup_m1m2_eg(0);
	float sigMETJERdo_m1m2_eg(0);
	float sigMTJESup_m1m2_eg(0);
	float sigMTJESdo_m1m2_eg(0);
	float sigMTJERup_m1m2_eg(0);
	float sigMTJERdo_m1m2_eg(0);
	float HTJESup_m1m2_eg(0);
	float HTJESdo_m1m2_eg(0);
	float Mass1_m1m2_eg(0);
	float Mass2_m1m2_eg(0);
	int   branch_eg(0);
  egtree_m1m2->SetBranchAddress("phoEt",      &phoEt_m1m2_eg);
  egtree_m1m2->SetBranchAddress("phoEta",     &phoEta_m1m2_eg);
  egtree_m1m2->SetBranchAddress("lepPt",      &lepPt_m1m2_eg);
  egtree_m1m2->SetBranchAddress("lepEta",     &lepEta_m1m2_eg);
  egtree_m1m2->SetBranchAddress("sigMT",      &sigMT_m1m2_eg);
  egtree_m1m2->SetBranchAddress("sigMET",     &sigMET_m1m2_eg);
	egtree_m1m2->SetBranchAddress("nVertex",    &nVertex_m1m2_eg);
  egtree_m1m2->SetBranchAddress("HT",         &HT_m1m2_eg);
	egtree_m1m2->SetBranchAddress("Mass1",      &Mass1_m1m2_eg);
	egtree_m1m2->SetBranchAddress("Mass2",      &Mass2_m1m2_eg);
	egtree_m1m2->SetBranchAddress("branch",     &branch_eg);
	egtree_m1m2->SetBranchAddress("sigMETJESup",&sigMETJESup_m1m2_eg);
	egtree_m1m2->SetBranchAddress("sigMETJESdo",&sigMETJESdo_m1m2_eg);
	egtree_m1m2->SetBranchAddress("sigMETJERup",&sigMETJERup_m1m2_eg);
	egtree_m1m2->SetBranchAddress("sigMETJERdo",&sigMETJERdo_m1m2_eg);
	egtree_m1m2->SetBranchAddress("sigMTJESup", &sigMTJESup_m1m2_eg);
	egtree_m1m2->SetBranchAddress("sigMTJESdo", &sigMTJESdo_m1m2_eg);
	egtree_m1m2->SetBranchAddress("sigMTJERup", &sigMTJERup_m1m2_eg);
	egtree_m1m2->SetBranchAddress("sigMTJERdo", &sigMTJERdo_m1m2_eg);
	egtree_m1m2->SetBranchAddress("HTJESup",    &HTJESup_m1m2_eg);
	egtree_m1m2->SetBranchAddress("HTJESdo",    &HTJESdo_m1m2_eg);

	for(unsigned ievt(0); ievt < egtree_m1m2->GetEntries(); ievt++){
		egtree_m1m2->GetEntry(ievt);

		/** cut flow *****/
		if(phoEt_m1m2_eg < 35 || lepPt_m1m2_eg < 25)continue;
		if(fabs(phoEta_m1m2_eg) > 1.4442 || fabs(lepEta_m1m2_eg) > 2.5)continue;

		double scalefactor(1);
		double scalefactorup(0);
		scalefactor = objectESF.getFastElectronESF(lepPt_m1m2_eg,lepEta_m1m2_eg)*objectESF.getPhotonESF(phoEt_m1m2_eg,phoEta_m1m2_eg)*objectESF.getFastegPhotonTRGESF(phoEt_m1m2_eg,phoEta_m1m2_eg)*objectESF.getFastElectronTRGESF(lepPt_m1m2_eg,lepEta_m1m2_eg);
		double s_ele_error = objectESF.getFastElectronESFError(lepPt_m1m2_eg,lepEta_m1m2_eg)*objectESF.getPhotonESF(phoEt_m1m2_eg,phoEta_m1m2_eg)*objectESF.getFastegPhotonTRGESF(phoEt_m1m2_eg,phoEta_m1m2_eg)*objectESF.getFastElectronTRGESF(lepPt_m1m2_eg,lepEta_m1m2_eg);
		double s_pho_error = objectESF.getPhotonESFError(phoEt_m1m2_eg,phoEta_m1m2_eg)*objectESF.getFastElectronESF(lepPt_m1m2_eg,lepEta_m1m2_eg)*objectESF.getFastegPhotonTRGESF(phoEt_m1m2_eg,phoEta_m1m2_eg)*objectESF.getFastElectronTRGESF(lepPt_m1m2_eg,lepEta_m1m2_eg);
		double s_eletrg_error = objectESF.getFastElectronTRGESFError(lepPt_m1m2_eg,lepEta_m1m2_eg)*objectESF.getFastElectronESF(lepPt_m1m2_eg,lepEta_m1m2_eg)*objectESF.getPhotonESF(phoEt_m1m2_eg,phoEta_m1m2_eg)*objectESF.getFastegPhotonTRGESF(phoEt_m1m2_eg,phoEta_m1m2_eg);
		double s_photrg_error = objectESF.getFastegPhotonTRGESFError(phoEt_m1m2_eg,phoEta_m1m2_eg)*objectESF.getFastElectronESF(lepPt_m1m2_eg,lepEta_m1m2_eg)*objectESF.getPhotonESF(phoEt_m1m2_eg,phoEta_m1m2_eg)*objectESF.getFastElectronTRGESF(lepPt_m1m2_eg,lepEta_m1m2_eg);
		double s_error = sqrt(pow(s_ele_error,2) + pow(s_pho_error,2) + pow(s_eletrg_error,2) + pow(s_photrg_error, 2)); 
		scalefactorup = scalefactor + s_error; 

		if( fabs(Mass1_m1m2_eg - 900) < 10 ){
			if(sigMT_m1m2_eg > 100){
				p_m1m2_MET_signal_eg->Fill(sigMET_m1m2_eg);
				if(sigMET_m1m2_eg > 120)p_m1m2_HT_signal_eg->Fill(HT_m1m2_eg);
				if(sigMET_m1m2_eg > 120)p_m1m2_PhoEt_signal_eg->Fill(phoEt_m1m2_eg);
			}
		}

		if(sigMET_m1m2_eg < 120 || sigMT_m1m2_eg < 100)continue;

		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET_m1m2_eg, HT_m1m2_eg, phoEt_m1m2_eg) + NBIN;
		if(SigBinIndex >=0){
			m1m2_h_chan_rate_nom[SigBinIndex]->Fill( Mass1_m1m2_eg,branch_eg, scalefactor);
			m1m2_h_chan_rate_xsUp[SigBinIndex]->Fill( Mass1_m1m2_eg,branch_eg,  scalefactor); 
			m1m2_h_chan_rate_esfUp[SigBinIndex]->Fill( Mass1_m1m2_eg,branch_eg,  scalefactorup); 

	    p_M1M2accept->Fill(Mass1_m1m2_eg,branch_eg);
			if(nVertex_m1m2_eg < 20)p_lowPU_m1m2_pass->Fill( Mass1_m1m2_eg, branch_eg, 1);
			else p_highPU_m1m2_pass->Fill( Mass1_m1m2_eg,branch_eg,  1);
		}
		int jesupBinIndex(-1);	
		jesupBinIndex = Bin.findSignalBin(sigMETJESup_m1m2_eg, HTJESup_m1m2_eg, phoEt_m1m2_eg) + NBIN;
		if(jesupBinIndex >=0){
			m1m2_h_chan_rate_jesUp[jesupBinIndex]->Fill( Mass1_m1m2_eg, branch_eg, scalefactor); 
		}
		int jesdoBinIndex(-1);
		jesdoBinIndex = Bin.findSignalBin(sigMETJESdo_m1m2_eg, HTJESdo_m1m2_eg, phoEt_m1m2_eg) + NBIN;
		if(jesdoBinIndex >=0){
			m1m2_h_chan_rate_jesDown[jesdoBinIndex]->Fill( Mass1_m1m2_eg, branch_eg, scalefactor);   
		}	
		int jerupBinIndex(-1);
		jerupBinIndex = Bin.findSignalBin(sigMETJERup_m1m2_eg, HT_m1m2_eg, phoEt_m1m2_eg) + NBIN;
		if( jerupBinIndex >=0){
			m1m2_h_chan_rate_jerUp[jerupBinIndex]->Fill( Mass1_m1m2_eg, branch_eg, scalefactor);
		}
		int jerdoBinIndex(-1);
		jerdoBinIndex = Bin.findSignalBin(sigMETJERdo_m1m2_eg, HT_m1m2_eg, phoEt_m1m2_eg) + NBIN;	
		if(jerdoBinIndex >= 0){
			m1m2_h_chan_rate_jerDown[jerdoBinIndex]->Fill( Mass1_m1m2_eg, branch_eg,  scalefactor);
		}  
	}

	for(unsigned i(1); i < p_M1M2accept->GetXaxis()->GetNbins() + 1; i++){
		for(unsigned j(1); j < p_M1M2accept->GetYaxis()->GetNbins() + 1; j++){
			float nSelect = p_M1M2accept->GetBinContent(i,j);
			if(p_M1M2MASS->GetBinContent(i,j) > 1)p_M1M2accept->SetBinContent(i,j, nSelect/p_M1M2MASS->GetBinContent(i,j));
		}
	}

	for(unsigned ih(0); ih < NBIN*2; ih++){
		for(unsigned i(1); i < m1m2_h_chan_rate_nom[ih]->GetXaxis()->GetNbins() + 1; i++){
			for(unsigned j(1); j < m1m2_h_chan_rate_nom[ih]->GetYaxis()->GetNbins() + 1; j++){
				if(p_M1M2MASS->GetBinContent(i,j) < 100){
					p_M1M2accept->SetBinContent(i,j,-1);
					p_M1M2MASS->SetBinContent(i,j,-1);
				}
	
				if(p_M1M2MASS->GetBinContent(i,j) <= 0){
					m1m2_h_chan_rate_nom[ih]->SetBinContent(i,j,-1); 
					m1m2_h_chan_rate_jesUp[ih]->SetBinContent(i,j,-1);
					m1m2_h_chan_rate_jesDown[ih]->SetBinContent(i,j,-1);
					m1m2_h_chan_rate_jerUp[ih]->SetBinContent(i,j,-1);
					m1m2_h_chan_rate_jerDown[ih]->SetBinContent(i,j,-1);
					m1m2_h_chan_rate_xsUp[ih]->SetBinContent(i,j, -1);
							
					m1m2_h_syserr_PU->SetBinContent(i,j, -1);
					m1m2_h_chan_syserr_jes[ih]->SetBinContent(i,j, -1); 
					m1m2_h_chan_syserr_jer[ih]->SetBinContent(i,j, -1); 
					m1m2_h_chan_syserr_esf[ih]->SetBinContent(i,j, -1); 
					m1m2_h_chan_syserr_scale[ih]->SetBinContent(i,j, -1); 
					m1m2_h_chan_syserr_eleshape[ih]->SetBinContent(i,j, -1); 
					m1m2_h_chan_syserr_jetshape[ih]->SetBinContent(i,j, -1); 
					m1m2_h_chan_syserr_qcdshape[ih]->SetBinContent(i,j, -1); 
					m1m2_h_chan_syserr_xs[ih]->SetBinContent(i,j, -1); 
					m1m2_h_chan_syserr_lumi[ih]->SetBinContent(i,j, -1); 
					m1m2_h_chan_syserr_isr[ih]->SetBinContent(i,j, -1);
					m1m2_h_chan_staterr[ih]->SetBinContent(i,j, -1); 
				}

				else{
					float noe = p_M1M2MASS->GetBinContent(i,j);
					int sparticleMass1 = int(p_M1M2MASS->GetXaxis()->GetBinCenter(i));
					float crosssection = xsmap[sparticleMass1]; 
					float crosssectionUp = xsmap[sparticleMass1]; 
					//std::cout << sparticleMass1 << " " << sparticleMass2 << " " << crosssection << std::endl;
					p_M1M2_XS->SetBinContent(i,j, crosssection);

					m1m2_h_chan_staterr[ih]->SetBinContent(i,j, sqrt(m1m2_h_chan_rate_nom[ih]->GetBinContent(i,j))/m1m2_h_chan_rate_nom[ih]->GetBinContent(i));
					m1m2_h_chan_rate_nom[ih]->SetBinError(i,j, sqrt(m1m2_h_chan_rate_nom[ih]->GetBinContent(i,j))*35.9*crosssection/noe);
					m1m2_h_chan_rate_nom[ih]->SetBinContent(i,j, m1m2_h_chan_rate_nom[ih]->GetBinContent(i,j)*35.9*crosssection/noe); 
					m1m2_h_chan_rate_jesUp[ih]->SetBinContent(i,j, m1m2_h_chan_rate_jesUp[ih]->GetBinContent(i,j)*35.9*crosssection/noe); 
					m1m2_h_chan_rate_jesDown[ih]->SetBinContent(i,j, m1m2_h_chan_rate_jesDown[ih]->GetBinContent(i,j)*35.9*crosssection/noe);
					m1m2_h_chan_rate_jerUp[ih]->SetBinContent(i,j, m1m2_h_chan_rate_jerUp[ih]->GetBinContent(i,j)*35.9*crosssection/noe);
					m1m2_h_chan_rate_jerDown[ih]->SetBinContent(i,j, m1m2_h_chan_rate_jerDown[ih]->GetBinContent(i,j)*35.9*crosssection/noe);
					m1m2_h_chan_rate_xsUp[ih]->SetBinContent(i,j, m1m2_h_chan_rate_xsUp[ih]->GetBinContent(i,j)*35.9*crosssectionUp/noe);
					m1m2_h_chan_rate_esfUp[ih]->SetBinContent(i,j, m1m2_h_chan_rate_esfUp[ih]->GetBinContent(i,j)*35.9*crosssection/noe);
									
					m1m2_h_chan_syserr_jes[ih]->SetBinContent(i,j, max( fabs(m1m2_h_chan_rate_jesUp[ih]->GetBinContent(i,j)-m1m2_h_chan_rate_nom[ih]->GetBinContent(i,j)), fabs(m1m2_h_chan_rate_jesDown[ih]->GetBinContent(i,j)-m1m2_h_chan_rate_nom[ih]->GetBinContent(i,j)) ) ); 
					m1m2_h_chan_syserr_jer[ih]->SetBinContent(i,j, max( fabs(m1m2_h_chan_rate_jerUp[ih]->GetBinContent(i,j)-m1m2_h_chan_rate_nom[ih]->GetBinContent(i,j)), fabs(m1m2_h_chan_rate_jerDown[ih]->GetBinContent(i,j)-m1m2_h_chan_rate_nom[ih]->GetBinContent(i,j)) ) ); 
					m1m2_h_chan_syserr_esf[ih]->SetBinContent(i,j, fabs( m1m2_h_chan_rate_esfUp[ih]->GetBinContent(i,j)-m1m2_h_chan_rate_nom[ih]->GetBinContent(i,j)) );
					m1m2_h_chan_syserr_scale[ih]->SetBinContent(i,j, -1);
					m1m2_h_chan_syserr_eleshape[ih]->SetBinContent(i,j, -1);
					m1m2_h_chan_syserr_jetshape[ih]->SetBinContent(i,j, -1);
					m1m2_h_chan_syserr_qcdshape[ih]->SetBinContent(i,j, -1);
					m1m2_h_chan_syserr_xs[ih]->SetBinContent(i,j, -1); 
					m1m2_h_chan_syserr_lumi[ih]->SetBinContent(i,j, 0.026*m1m2_h_chan_rate_nom[ih]->GetBinContent(i,j));   
					m1m2_h_chan_syserr_isr[ih]->SetBinContent(i,j, -1);
					if(ih==0)m1m2_h_syserr_PU->SetBinContent(i,j, analysis_PU(i,j, p_lowPU_m1m2_pass, p_lowPU_m1m2_all, p_highPU_m1m2_pass, p_highPU_m1m2_all, p_PU_data)); 
				} 
			}
		} 
	} 

////		for(unsigned i(1); i < p_M1M2MASS->GetXaxis()->GetNbins() + 1; i++){
////			for(unsigned j(1); j < p_M1M2MASS->GetYaxis()->GetNbins() + 1; j++){
////				float sparticleMass = p_M1M2MASS->GetXaxis()->GetBinCenter(i);
////				if( fabs(sparticleMass - 1700 ) < 10 && fabs( p_M1M2MASS->GetYaxis()->GetBinCenter(j) - 1000) < 5){
////					float noe = p_M1M2MASS->GetBinContent(i,j);
////					float crosssection = 0.5*p_crosssection_m1m2->GetBinContent( p_crosssection_m1m2->FindBin(sparticleMass) );
////					p_m1m2_MET_signal_mg->Scale(35.8*crosssection/noe);
////					p_m1m2_HT_signal_mg->Scale(35.8*crosssection/noe);
////					p_m1m2_PhoEt_signal_mg->Scale(35.8*crosssection/noe);
////					p_m1m2_MET_signal_eg->Scale(35.8*crosssection/noe);
////					p_m1m2_HT_signal_eg->Scale(35.8*crosssection/noe);
////					p_m1m2_PhoEt_signal_eg->Scale(35.8*crosssection/noe);
////			}
////		}
////	}

	outputfile_m1m2->Write();
	outputfile_m1m2->Close();

}


