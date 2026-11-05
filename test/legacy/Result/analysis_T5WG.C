#include "analysis_PU.C"
#include "analysis_SUSY.h"

void analysis_T5WG(){//main  

	SetSignalConfig();
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);
	esfScaleFactor  objectESF;

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

	TFile xSecFile("../../cross/susyCrossSection.root");
	TH1D *p_crosssection_susy   = (TH1D*)xSecFile.Get("p_gluinoxSec");

	TChain *datachain = new TChain("signalTree");
	datachain->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_mgsignal_MuonEG_FullEcal.root");
	TH1D *p_PU_data = new TH1D("p_PU_data",";N_{vtx};",100,0,100); 
  datachain->Draw("nVertex >> p_PU_data");
	p_PU_data->Scale(1.0/p_PU_data->Integral(1,101));

	init_histo("T5WG","test_T5WG_1.root",2,27, 775.0, 2125.0, 420, 2.5, 2102.5);

	std::ostringstream histname;
	//**************   T5WG  ***************************//
	TFile *file_susy = TFile::Open("/uscms_data/d3/mengleis/FullStatusOct/resTree_T5WG_debug2.root");
  TTree *tree_susy = (TTree*)file_susy->Get("SUSYtree");
	float Mgluino_susy(0);
  float Mchargino_susy(0);
  float Mneutralino_susy(0);
	float Mass1_susy(0);
	float Mass2_susy(0);
	float ISRPt(0);
	int   nVertex(0);
	std::vector<float> *ScaleSystWeight = 0;
	tree_susy->SetBranchAddress("MsGsQ",      &Mgluino_susy);  
  tree_susy->SetBranchAddress("Mchargino",  &Mchargino_susy);
  tree_susy->SetBranchAddress("Mneutralino",&Mneutralino_susy);
  tree_susy->SetBranchAddress("Mass1",      &Mass1_susy);
  tree_susy->SetBranchAddress("Mass2",      &Mass2_susy);
  tree_susy->SetBranchAddress("ISRPt",      &ISRPt);
  tree_susy->SetBranchAddress("nVertex",    &nVertex);
	tree_susy->SetBranchAddress("ScaleSystWeight",&ScaleSystWeight);

	for(unsigned ievt(0); ievt < tree_susy->GetEntries(); ievt++){
		tree_susy->GetEntry(ievt);
    if (ievt%1000000==0) std::cout << " -- Processing event " << ievt << std::endl;

		p_SUSYMass->Fill(Mass1_susy, Mass2_susy);
		for(unsigned i(0); i < 9; i++)p_SUSYMASS_pdf[i]->Fill(Mass1_susy, Mass2_susy, (*ScaleSystWeight)[i]);
		if(nVertex < 20)p_lowPU_susy_all->Fill(Mass1_susy, Mass2_susy, 1);
		else p_highPU_susy_all->Fill(Mass1_susy, Mass2_susy, 1);

		double reweightF = 1;
    if(ISRPt < 50)reweightF = 1.015;
    else if(ISRPt >= 50 && ISRPt < 100)reweightF  = 1.110;
    else if(ISRPt >= 100 && ISRPt < 150)reweightF = 0.845;
    else if(ISRPt >= 150 && ISRPt < 200)reweightF = 0.715;
    else if(ISRPt >= 200 && ISRPt < 250)reweightF = 0.730;
    else if(ISRPt >= 250 && ISRPt < 300)reweightF = 0.732;
    else if(ISRPt >= 300)reweightF =  0.642;
		p_SUSYISR->Fill(Mass1_susy, Mass2_susy, reweightF);
	}

  TChain *mgtree_susy;
  mgtree_susy = new TChain("mgTree","mgTree");
	mgtree_susy->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_T5WG_debug2.root");
  float phoEt_susy_mg(0);
  float phoEta_susy_mg(0);
  float lepPt_susy_mg(0);
  float lepEta_susy_mg(0);
  float sigMT_susy_mg(0);
  float sigMET_susy_mg(0);
  float genMET_susy_mg(0);
  float HT_susy_mg(0);
	int   nVertex_susy_mg(0);
	float sigMETJESup_susy_mg(0);
	float sigMETJESdo_susy_mg(0);
	float sigMTJESup_susy_mg(0);
	float sigMTJESdo_susy_mg(0);
	float HTJESup_susy_mg(0);
	float HTJESdo_susy_mg(0);
	float gluinoMass_susy_mg(0);
  float charginoMass_susy_mg(0);
  float neutralinoMass_susy_mg(0);
  float Mass1_susy_mg(0);
  float Mass2_susy_mg(0);
	float ISRPt_mg(0);
	std::vector<float> *ScaleSystWeight_mg = 0;
  mgtree_susy->SetBranchAddress("phoEt",      &phoEt_susy_mg);
  mgtree_susy->SetBranchAddress("phoEta",     &phoEta_susy_mg);
  mgtree_susy->SetBranchAddress("lepPt",      &lepPt_susy_mg);
  mgtree_susy->SetBranchAddress("lepEta",     &lepEta_susy_mg);
  mgtree_susy->SetBranchAddress("sigMT",      &sigMT_susy_mg);
  mgtree_susy->SetBranchAddress("sigMET",     &sigMET_susy_mg);
  mgtree_susy->SetBranchAddress("genMET",     &genMET_susy_mg);
  mgtree_susy->SetBranchAddress("HT",         &HT_susy_mg);
	mgtree_susy->SetBranchAddress("nVertex",    &nVertex_susy_mg);
	mgtree_susy->SetBranchAddress("MsGsQ",      &gluinoMass_susy_mg);
  mgtree_susy->SetBranchAddress("Mchargino",  &charginoMass_susy_mg);
  mgtree_susy->SetBranchAddress("Mneutralino",&neutralinoMass_susy_mg);
	mgtree_susy->SetBranchAddress("Mass1",      &Mass1_susy_mg);
	mgtree_susy->SetBranchAddress("Mass2",      &Mass2_susy_mg);
  mgtree_susy->SetBranchAddress("ISRPt",      &ISRPt_mg);
	mgtree_susy->SetBranchAddress("sigMETJESup",&sigMETJESup_susy_mg);
	mgtree_susy->SetBranchAddress("sigMETJESdo",&sigMETJESdo_susy_mg);
	mgtree_susy->SetBranchAddress("sigMTJESup", &sigMTJESup_susy_mg);
	mgtree_susy->SetBranchAddress("sigMTJESdo", &sigMTJESdo_susy_mg);
	mgtree_susy->SetBranchAddress("HTJESup",    &HTJESup_susy_mg);
	mgtree_susy->SetBranchAddress("HTJESdo",    &HTJESdo_susy_mg);
	mgtree_susy->SetBranchAddress("ScaleSystWeight", &ScaleSystWeight_mg);

	for(unsigned ievt(0); ievt < mgtree_susy->GetEntries(); ievt++){
		mgtree_susy->GetEntry(ievt);
	
		/** cut flow *****/
		if(phoEt_susy_mg < 35 || lepPt_susy_mg < 25)continue;
		if(fabs(phoEta_susy_mg) > 1.4442 || fabs(lepEta_susy_mg) > 2.5)continue;

    p_SUSYselect->Fill(Mass1_susy_mg, Mass2_susy_mg);

		double scalefactor(0);
		double scalefactorup(0);
		scalefactor = objectESF.getFastMuonESF(lepPt_susy_mg,lepEta_susy_mg)*objectESF.getPhotonESF(phoEt_susy_mg, phoEta_susy_mg)*objectESF.getFastMuonEGTRGESF(phoEt_susy_mg, lepPt_susy_mg);
		double s_mu_error = objectESF.getFastMuonESFError(lepPt_susy_mg,lepEta_susy_mg)*objectESF.getPhotonESF(phoEt_susy_mg, phoEta_susy_mg)*objectESF.getFastMuonEGTRGESF(phoEt_susy_mg, lepPt_susy_mg);
    double s_pho_error = objectESF.getPhotonESFError(phoEt_susy_mg, phoEta_susy_mg)*objectESF.getFastMuonESF(lepPt_susy_mg,lepEta_susy_mg)*objectESF.getFastMuonEGTRGESF(phoEt_susy_mg, lepPt_susy_mg);
    double s_trg_error = objectESF.getFastMuonEGTRGESFError(phoEt_susy_mg, lepPt_susy_mg)*objectESF.getFastMuonESF(lepPt_susy_mg,lepEta_susy_mg)*objectESF.getPhotonESF(phoEt_susy_mg, phoEta_susy_mg);
		double s_error = sqrt(pow(s_mu_error,2) + pow(s_pho_error,2) + pow(s_trg_error,2));
		scalefactorup = scalefactor + s_error; 

		double reweightF = 1;
    if(ISRPt_mg < 50)reweightF = 1.015;
    else if(ISRPt_mg >= 50 && ISRPt_mg < 100)reweightF  = 1.110;
    else if(ISRPt_mg >= 100 && ISRPt_mg < 150)reweightF = 0.845;
    else if(ISRPt_mg >= 150 && ISRPt_mg < 200)reweightF = 0.715;
    else if(ISRPt_mg >= 200 && ISRPt_mg < 250)reweightF = 0.730;
    else if(ISRPt_mg >= 250 && ISRPt_mg < 300)reweightF = 0.732;
    else if(ISRPt_mg >= 300)reweightF =  0.642;

		if( fabs(Mass1_susy_mg - 1700) < 10 && fabs(Mass2_susy_mg - 1000) < 5){
			if(sigMT_susy_mg > 100){
				p_susy_MET_signal_mg->Fill(sigMET_susy_mg);
				if(sigMET_susy_mg > 120)p_susy_HT_signal_mg->Fill(HT_susy_mg);
				if(sigMET_susy_mg > 120)p_susy_PhoEt_signal_mg->Fill(phoEt_susy_mg);
			}
		}

		if(sigMT_susy_mg > 100){
			int pfMETBinIndex(-1);
			pfMETBinIndex = Bin.findSignalBin(sigMET_susy_mg, HT_susy_mg, phoEt_susy_mg);
			if(pfMETBinIndex >= 0)pfMETBinIndex += NBIN;
			int genMETBinIndex(-1);
			genMETBinIndex = Bin.findSignalBin(genMET_susy_mg, HT_susy_mg, phoEt_susy_mg);
			if(genMETBinIndex >= 0)genMETBinIndex += NBIN;
			if(pfMETBinIndex >=0)susy_all_chan_rate_pfMET[pfMETBinIndex]->Fill(Mass1_susy_mg, Mass2_susy_mg);
			if(genMETBinIndex >=0)susy_all_chan_rate_genMET[genMETBinIndex]->Fill(Mass1_susy_mg, Mass2_susy_mg);
			if(charginoMass_susy_mg > 0 && neutralinoMass_susy_mg > 0){
				if(pfMETBinIndex >=0)susy_wg_chan_rate_pfMET[pfMETBinIndex]->Fill(Mass1_susy_mg,Mass2_susy_mg);
				if(genMETBinIndex >=0)susy_wg_chan_rate_genMET[genMETBinIndex]->Fill(Mass1_susy_mg,Mass2_susy_mg);
			}
			else if(charginoMass_susy_mg <= 0 && neutralinoMass_susy_mg > 0){
				if(pfMETBinIndex >=0)susy_gg_chan_rate_pfMET[pfMETBinIndex]->Fill(Mass1_susy_mg,Mass2_susy_mg); 
				if(genMETBinIndex >=0)susy_gg_chan_rate_genMET[genMETBinIndex]->Fill(Mass1_susy_mg,Mass2_susy_mg);
			}
		}	

		if(sigMET_susy_mg < 120 || sigMT_susy_mg < 100)continue;
		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET_susy_mg, HT_susy_mg, phoEt_susy_mg);
		if(SigBinIndex >= 0)SigBinIndex += NBIN;
		int jesupBinIndex(-1);	
		jesupBinIndex = Bin.findSignalBin(sigMETJESup_susy_mg, HTJESup_susy_mg, phoEt_susy_mg);
		if(jesupBinIndex >= 0)jesupBinIndex += NBIN; 
		int jesdoBinIndex(-1);
		jesdoBinIndex = Bin.findSignalBin(sigMETJESdo_susy_mg, HTJESdo_susy_mg, phoEt_susy_mg);
		if(jesdoBinIndex >= 0)jesdoBinIndex += NBIN;

		if(Mass1_susy_mg > 0 && Mass2_susy_mg > 0){ 
			if(SigBinIndex >=0){
				susy_all_chan_rate_nom[SigBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor);
				susy_all_chan_rate_esfUp[SigBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg,scalefactorup); 
				susy_all_chan_rate_isr[SigBinIndex]->Fill(Mass1_susy_mg, Mass2_susy_mg, scalefactor*reweightF);	
	
				if(nVertex_susy_mg < 20)p_lowPU_susy_pass->Fill( Mass1_susy_mg, Mass2_susy_mg, 1);
				else p_highPU_susy_pass->Fill( Mass1_susy_mg, Mass2_susy_mg, 1);

			  for(unsigned is(0); is < 9; is++)
					susy_all_chan_rate_pdf[SigBinIndex][is]->Fill(Mass1_susy_mg, Mass2_susy_mg, scalefactor*(*ScaleSystWeight_mg)[is]);
			}
			if(jesupBinIndex >=0){
				susy_all_chan_rate_jesUp[jesupBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor); 
			}
			if(jesdoBinIndex >=0){
				susy_all_chan_rate_jesDown[jesdoBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor);   
			}	
		}

		if(charginoMass_susy_mg > 0 && neutralinoMass_susy_mg > 0){
			if(SigBinIndex >=0){
				susy_wg_chan_rate_nom[SigBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor);
				susy_wg_chan_rate_esfUp[SigBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg,scalefactorup); 
				susy_wg_chan_rate_isr[SigBinIndex]->Fill(Mass1_susy_mg, Mass2_susy_mg, scalefactor*reweightF);	
	
				if(nVertex_susy_mg < 20)p_lowPU_susy_pass->Fill( Mass1_susy_mg, Mass2_susy_mg, 1);
				else p_highPU_susy_pass->Fill( Mass1_susy_mg, Mass2_susy_mg, 1);

			  for(unsigned is(0); is < 9; is++)
					susy_wg_chan_rate_pdf[SigBinIndex][is]->Fill(Mass1_susy_mg, Mass2_susy_mg, scalefactor*(*ScaleSystWeight_mg)[is]);
			}
			if(jesupBinIndex >=0){
				susy_wg_chan_rate_jesUp[jesupBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor); 
			}
			if(jesdoBinIndex >=0){
				susy_wg_chan_rate_jesDown[jesdoBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor);   
			}	
		} 
		else if(charginoMass_susy_mg <= 0 && neutralinoMass_susy_mg > 0){
			if(SigBinIndex >=0){
				susy_gg_chan_rate_nom[SigBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor);
				susy_gg_chan_rate_esfUp[SigBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg,scalefactorup); 
				susy_gg_chan_rate_isr[SigBinIndex]->Fill(Mass1_susy_mg, Mass2_susy_mg, scalefactor*reweightF);	
	
				if(nVertex_susy_mg < 20)p_lowPU_susy_pass->Fill( Mass1_susy_mg, Mass2_susy_mg, 1);
				else p_highPU_susy_pass->Fill( Mass1_susy_mg, Mass2_susy_mg, 1);

			  for(unsigned is(0); is < 9; is++)
					susy_gg_chan_rate_pdf[SigBinIndex][is]->Fill(Mass1_susy_mg, Mass2_susy_mg, scalefactor*(*ScaleSystWeight_mg)[is]);
			}
			if(jesupBinIndex >=0){
				susy_gg_chan_rate_jesUp[jesupBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor); 
			}
			if(jesdoBinIndex >=0){
				susy_gg_chan_rate_jesDown[jesdoBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor);   
			}	
		} 
	}

  TChain *egtree_susy;
  egtree_susy = new TChain("egTree","egTree");
	egtree_susy->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_T5WG_debug2.root");
  float phoEt_susy_eg(0);
  float phoEta_susy_eg(0);
  float lepPt_susy_eg(0);
  float lepEta_susy_eg(0);
  float sigMT_susy_eg(0);
  float sigMET_susy_eg(0);
  float HT_susy_eg(0);
  int   nVertex_susy_eg(0); 
	float sigMETJESup_susy_eg(0);
	float sigMETJESdo_susy_eg(0);
	float sigMTJESup_susy_eg(0);
	float sigMTJESdo_susy_eg(0);
	float HTJESup_susy_eg(0);
	float HTJESdo_susy_eg(0);
	float gluinoMass_susy_eg(0);
  float charginoMass_susy_eg(0);
  float neutralinoMass_susy_eg(0);
	float Mass1_susy_eg(0);
  float Mass2_susy_eg(0);
	float ISRPt_eg(0);
  float genMET_susy_eg(0);
	std::vector<float> *ScaleSystWeight_eg = 0;
  egtree_susy->SetBranchAddress("phoEt",      &phoEt_susy_eg);
  egtree_susy->SetBranchAddress("phoEta",     &phoEta_susy_eg);
  egtree_susy->SetBranchAddress("lepPt",      &lepPt_susy_eg);
  egtree_susy->SetBranchAddress("lepEta",     &lepEta_susy_eg);
  egtree_susy->SetBranchAddress("sigMT",      &sigMT_susy_eg);
  egtree_susy->SetBranchAddress("sigMET",     &sigMET_susy_eg);
  egtree_susy->SetBranchAddress("genMET",     &genMET_susy_eg);
	egtree_susy->SetBranchAddress("nVertex",    &nVertex_susy_eg);
  egtree_susy->SetBranchAddress("HT",         &HT_susy_eg);
	egtree_susy->SetBranchAddress("MsGsQ",      &gluinoMass_susy_eg);
  egtree_susy->SetBranchAddress("Mchargino",  &charginoMass_susy_eg);
  egtree_susy->SetBranchAddress("Mneutralino",&neutralinoMass_susy_eg);
	egtree_susy->SetBranchAddress("Mass1",      &Mass1_susy_eg);
	egtree_susy->SetBranchAddress("Mass2",      &Mass2_susy_eg);
  egtree_susy->SetBranchAddress("ISRPt",      &ISRPt_eg);
	egtree_susy->SetBranchAddress("sigMETJESup",&sigMETJESup_susy_eg);
	egtree_susy->SetBranchAddress("sigMETJESdo",&sigMETJESdo_susy_eg);
	egtree_susy->SetBranchAddress("sigMTJESup", &sigMTJESup_susy_eg);
	egtree_susy->SetBranchAddress("sigMTJESdo", &sigMTJESdo_susy_eg);
	egtree_susy->SetBranchAddress("HTJESup",    &HTJESup_susy_eg);
	egtree_susy->SetBranchAddress("HTJESdo",    &HTJESdo_susy_eg);
	egtree_susy->SetBranchAddress("ScaleSystWeight",&ScaleSystWeight_eg);	

	for(unsigned ievt(0); ievt < egtree_susy->GetEntries(); ievt++){
		egtree_susy->GetEntry(ievt);

		/** cut flow *****/
		if(phoEt_susy_eg < 35 || lepPt_susy_eg < 25)continue;
		if(fabs(phoEta_susy_eg) > 1.4442 || fabs(lepEta_susy_eg) > 2.5)continue;

		double scalefactor(1);
		double scalefactorup(0);
		scalefactor = objectESF.getFastElectronESF(lepPt_susy_eg,lepEta_susy_eg)*objectESF.getPhotonESF(phoEt_susy_eg,phoEta_susy_eg)*objectESF.getFastegPhotonTRGESF(phoEt_susy_eg,phoEta_susy_eg)*objectESF.getFastElectronTRGESF(lepPt_susy_eg,lepEta_susy_eg);
		double s_ele_error = objectESF.getFastElectronESFError(lepPt_susy_eg,lepEta_susy_eg)*objectESF.getPhotonESF(phoEt_susy_eg,phoEta_susy_eg)*objectESF.getFastegPhotonTRGESF(phoEt_susy_eg,phoEta_susy_eg)*objectESF.getFastElectronTRGESF(lepPt_susy_eg,lepEta_susy_eg);
		double s_pho_error = objectESF.getPhotonESFError(phoEt_susy_eg,phoEta_susy_eg)*objectESF.getFastElectronESF(lepPt_susy_eg,lepEta_susy_eg)*objectESF.getFastegPhotonTRGESF(phoEt_susy_eg,phoEta_susy_eg)*objectESF.getFastElectronTRGESF(lepPt_susy_eg,lepEta_susy_eg);
		double s_eletrg_error = objectESF.getFastElectronTRGESFError(lepPt_susy_eg,lepEta_susy_eg)*objectESF.getFastElectronESF(lepPt_susy_eg,lepEta_susy_eg)*objectESF.getPhotonESF(phoEt_susy_eg,phoEta_susy_eg)*objectESF.getFastegPhotonTRGESF(phoEt_susy_eg,phoEta_susy_eg);
		double s_photrg_error = objectESF.getFastegPhotonTRGESFError(phoEt_susy_eg,phoEta_susy_eg)*objectESF.getFastElectronESF(lepPt_susy_eg,lepEta_susy_eg)*objectESF.getPhotonESF(phoEt_susy_eg,phoEta_susy_eg)*objectESF.getFastElectronTRGESF(lepPt_susy_eg,lepEta_susy_eg);
		double s_error = sqrt(pow(s_ele_error,2) + pow(s_pho_error,2) + pow(s_eletrg_error,2) + pow(s_photrg_error, 2)); 
		scalefactorup = scalefactor + s_error; 

		double reweightF = 1;
    if(ISRPt_eg < 50)reweightF = 1.015;
    else if(ISRPt_eg >= 50 && ISRPt_eg < 100)reweightF  = 1.110;
    else if(ISRPt_eg >= 100 && ISRPt_eg < 150)reweightF = 0.845;
    else if(ISRPt_eg >= 150 && ISRPt_eg < 200)reweightF = 0.715;
    else if(ISRPt_eg >= 200 && ISRPt_eg < 250)reweightF = 0.730;
    else if(ISRPt_eg >= 250 && ISRPt_eg < 300)reweightF = 0.732;
    else if(ISRPt_eg >= 300)reweightF =  0.642;

		if( fabs(Mass1_susy_eg - 1700) < 10 && fabs(Mass2_susy_eg - 1000) < 5){
			if(sigMT_susy_eg > 100){
				p_susy_MET_signal_eg->Fill(sigMET_susy_eg);
				if(sigMET_susy_eg > 120)p_susy_HT_signal_eg->Fill(HT_susy_eg);
				if(sigMET_susy_eg > 120)p_susy_PhoEt_signal_eg->Fill(phoEt_susy_eg);
			}
		}

		if(sigMT_susy_eg > 100){
			int pfMETBinIndex(-1);
			pfMETBinIndex = Bin.findSignalBin(sigMET_susy_eg, HT_susy_eg, phoEt_susy_eg);
			int genMETBinIndex(-1);
			genMETBinIndex = Bin.findSignalBin(genMET_susy_eg, HT_susy_eg, phoEt_susy_eg);
			if(pfMETBinIndex >=0)susy_all_chan_rate_pfMET[pfMETBinIndex]->Fill(Mass1_susy_eg, Mass2_susy_eg);
			if(genMETBinIndex >=0)susy_all_chan_rate_genMET[genMETBinIndex]->Fill(Mass1_susy_eg, Mass2_susy_eg);
			if(charginoMass_susy_eg > 0 && neutralinoMass_susy_eg > 0){
				if(pfMETBinIndex >=0)susy_wg_chan_rate_pfMET[pfMETBinIndex]->Fill(Mass1_susy_eg,Mass2_susy_eg);
				if(genMETBinIndex >=0)susy_wg_chan_rate_genMET[genMETBinIndex]->Fill(Mass1_susy_eg,Mass2_susy_eg);
			}
			else if(charginoMass_susy_eg <= 0 && neutralinoMass_susy_eg > 0){
				if(pfMETBinIndex >=0)susy_gg_chan_rate_pfMET[pfMETBinIndex]->Fill(Mass1_susy_eg,Mass2_susy_eg); 
				if(genMETBinIndex >=0)susy_gg_chan_rate_genMET[genMETBinIndex]->Fill(Mass1_susy_eg,Mass2_susy_eg);
			}
		}	

		if(sigMET_susy_eg < 120 || sigMT_susy_eg < 100)continue;
		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET_susy_eg, HT_susy_eg, phoEt_susy_eg);
		int jesupBinIndex(-1);	
		jesupBinIndex = Bin.findSignalBin(sigMETJESup_susy_eg, HTJESup_susy_eg, phoEt_susy_eg);
		int jesdoBinIndex(-1);
		jesdoBinIndex = Bin.findSignalBin(sigMETJESdo_susy_eg, HTJESdo_susy_eg, phoEt_susy_eg);

		if(Mass1_susy_eg > 0 && Mass2_susy_eg > 0){ 
			if(SigBinIndex >=0){
				susy_all_chan_rate_nom[SigBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor);
				susy_all_chan_rate_esfUp[SigBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg,scalefactorup); 
				susy_all_chan_rate_isr[SigBinIndex]->Fill(Mass1_susy_eg, Mass2_susy_eg, scalefactor*reweightF);	
	
				if(nVertex_susy_eg < 20)p_lowPU_susy_pass->Fill( Mass1_susy_eg, Mass2_susy_eg, 1);
				else p_highPU_susy_pass->Fill( Mass1_susy_eg, Mass2_susy_eg, 1);

			  for(unsigned is(0); is < 9; is++)
					susy_all_chan_rate_pdf[SigBinIndex][is]->Fill(Mass1_susy_eg, Mass2_susy_eg, scalefactor*(*ScaleSystWeight_eg)[is]);
			}
			if(jesupBinIndex >=0){
				susy_all_chan_rate_jesUp[jesupBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor); 
			}
			if(jesdoBinIndex >=0){
				susy_all_chan_rate_jesDown[jesdoBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor);   
			}	
		}

		if(charginoMass_susy_eg > 0 && neutralinoMass_susy_eg > 0){
			if(SigBinIndex >=0){
				susy_wg_chan_rate_nom[SigBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor);
				susy_wg_chan_rate_esfUp[SigBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg,scalefactorup); 
				susy_wg_chan_rate_isr[SigBinIndex]->Fill(Mass1_susy_eg, Mass2_susy_eg, scalefactor*reweightF);	
	
				if(nVertex_susy_eg < 20)p_lowPU_susy_pass->Fill( Mass1_susy_eg, Mass2_susy_eg, 1);
				else p_highPU_susy_pass->Fill( Mass1_susy_eg, Mass2_susy_eg, 1);

			  for(unsigned is(0); is < 9; is++)
					susy_wg_chan_rate_pdf[SigBinIndex][is]->Fill(Mass1_susy_eg, Mass2_susy_eg, scalefactor*(*ScaleSystWeight_eg)[is]);
			}
			if(jesupBinIndex >=0){
				susy_wg_chan_rate_jesUp[jesupBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor); 
			}
			if(jesdoBinIndex >=0){
				susy_wg_chan_rate_jesDown[jesdoBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor);   
			}	
		} 
		else if(charginoMass_susy_eg <= 0 && neutralinoMass_susy_eg > 0){
			if(SigBinIndex >=0){
				susy_gg_chan_rate_nom[SigBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor);
				susy_gg_chan_rate_esfUp[SigBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg,scalefactorup); 
				susy_gg_chan_rate_isr[SigBinIndex]->Fill(Mass1_susy_eg, Mass2_susy_eg, scalefactor*reweightF);	
	
				if(nVertex_susy_eg < 20)p_lowPU_susy_pass->Fill( Mass1_susy_eg, Mass2_susy_eg, 1);
				else p_highPU_susy_pass->Fill( Mass1_susy_eg, Mass2_susy_eg, 1);

			  for(unsigned is(0); is < 9; is++)
					susy_gg_chan_rate_pdf[SigBinIndex][is]->Fill(Mass1_susy_eg, Mass2_susy_eg, scalefactor*(*ScaleSystWeight_eg)[is]);
			}
			if(jesupBinIndex >=0){
				susy_gg_chan_rate_jesUp[jesupBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor); 
			}
			if(jesdoBinIndex >=0){
				susy_gg_chan_rate_jesDown[jesdoBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor);   
			}	
		} 
	}


	for(unsigned ih(0); ih < NBIN*2; ih++){
		for(unsigned i(1); i < susy_wg_chan_rate_nom[ih]->GetXaxis()->GetNbins() + 1; i++){
			for(unsigned j(1); j < susy_wg_chan_rate_nom[ih]->GetYaxis()->GetNbins() + 1; j++){
				if(p_SUSYMass->GetBinContent(i,j) < 1000){
					susy_all_chan_syserr_met[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_met[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_met[ih]->SetBinContent(i,j, -1); 
				}
				else{
					int ibin = (ih%9)/3;
					ibin = ibin*3;
					double pf_met = susy_all_chan_rate_pfMET[ibin]->GetBinContent(i,j)   + susy_all_chan_rate_pfMET[ibin+9]->GetBinContent(i,j) 
													+ susy_all_chan_rate_pfMET[ibin+1]->GetBinContent(i,j) + susy_all_chan_rate_pfMET[ibin+10]->GetBinContent(i,j)
													+ susy_all_chan_rate_pfMET[ibin+2]->GetBinContent(i,j) + susy_all_chan_rate_pfMET[ibin+11]->GetBinContent(i,j)
													+ susy_all_chan_rate_pfMET[ibin+18]->GetBinContent(i,j)   + susy_all_chan_rate_pfMET[ibin+18+9]->GetBinContent(i,j) 
													+ susy_all_chan_rate_pfMET[ibin+18+1]->GetBinContent(i,j) + susy_all_chan_rate_pfMET[ibin+18+10]->GetBinContent(i,j)
													+ susy_all_chan_rate_pfMET[ibin+18+2]->GetBinContent(i,j) + susy_all_chan_rate_pfMET[ibin+18+11]->GetBinContent(i,j);
					double gen_met = susy_all_chan_rate_genMET[ibin]->GetBinContent(i,j)   + susy_all_chan_rate_genMET[ibin+9]->GetBinContent(i,j) 
													+ susy_all_chan_rate_genMET[ibin+1]->GetBinContent(i,j) + susy_all_chan_rate_genMET[ibin+10]->GetBinContent(i,j)
													+ susy_all_chan_rate_genMET[ibin+2]->GetBinContent(i,j) + susy_all_chan_rate_genMET[ibin+11]->GetBinContent(i,j)
													+ susy_all_chan_rate_genMET[ibin+18]->GetBinContent(i,j)   + susy_all_chan_rate_genMET[ibin+18+9]->GetBinContent(i,j) 
													+ susy_all_chan_rate_genMET[ibin+18+1]->GetBinContent(i,j) + susy_all_chan_rate_genMET[ibin+18+10]->GetBinContent(i,j)
													+ susy_all_chan_rate_genMET[ibin+18+2]->GetBinContent(i,j) + susy_all_chan_rate_genMET[ibin+18+11]->GetBinContent(i,j);
					
					if(susy_all_chan_rate_nom[ih]->GetBinContent(i,j) > 0)
						susy_all_chan_syserr_met[ih]->SetBinContent(i,j, fabs((gen_met-pf_met)/2)/pf_met); 
			  	else susy_all_chan_syserr_met[ih]->SetBinContent(i,j, -1);
					if(susy_wg_chan_rate_nom[ih]->GetBinContent(i,j) > 0)
						susy_wg_chan_syserr_met[ih]->SetBinContent(i,j, fabs((gen_met-pf_met)/2)/pf_met);
			  	else susy_wg_chan_syserr_met[ih]->SetBinContent(i,j, -1);
					if(susy_gg_chan_rate_nom[ih]->GetBinContent(i,j) > 0)
						susy_gg_chan_syserr_met[ih]->SetBinContent(i,j, fabs((gen_met-pf_met)/2)/pf_met);
			  	else susy_gg_chan_syserr_met[ih]->SetBinContent(i,j, -1);
				}
			} 
		}
	} 

	for(unsigned ih(0); ih < NBIN*2; ih++){
		for(unsigned i(1); i < susy_wg_chan_rate_nom[ih]->GetXaxis()->GetNbins() + 1; i++){
			for(unsigned j(1); j < susy_wg_chan_rate_nom[ih]->GetYaxis()->GetNbins() + 1; j++){
				if(p_SUSYMass->GetBinContent(i,j) < 1000){
					p_SUSYselect->SetBinContent(i,j,-1);
					p_SUSYMass->SetBinContent(i,j,-1);
				}
	
				if(p_SUSYMass->GetBinContent(i,j) <= 0){
					susy_all_chan_rate_nom[ih]->SetBinContent(i,j, -1); 
							
					susy_all_syserr_PU->SetBinContent(i,j, -1);
					susy_all_chan_syserr_jes[ih]->SetBinContent(i,j, -1); 
					susy_all_chan_syserr_esf[ih]->SetBinContent(i,j, -1); 
					susy_all_chan_syserr_scale[ih]->SetBinContent(i,j, -1); 
					susy_all_chan_syserr_eleshape[ih]->SetBinContent(i,j, -1); 
					susy_all_chan_syserr_jetshape[ih]->SetBinContent(i,j, -1); 
					susy_all_chan_syserr_qcdshape[ih]->SetBinContent(i,j, -1); 
					susy_all_chan_syserr_xs[ih]->SetBinContent(i,j, -1); 
					susy_all_chan_syserr_lumi[ih]->SetBinContent(i,j, -1); 
					susy_all_chan_syserr_isr[ih]->SetBinContent(i,j, -1); 
					susy_all_chan_syserr_pdf[ih]->SetBinContent(i,j, -1); 
					susy_all_chan_syserr_met[ih]->SetBinContent(i,j, -1); 

					susy_wg_chan_rate_nom[ih]->SetBinContent(i,j, -1); 
							
					susy_wg_syserr_PU->SetBinContent(i,j, -1);
					susy_wg_chan_syserr_jes[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_esf[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_scale[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_eleshape[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_jetshape[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_qcdshape[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_xs[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_lumi[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_isr[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_pdf[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_met[ih]->SetBinContent(i,j, -1); 

					susy_gg_chan_rate_nom[ih]->SetBinContent(i,j, -1); 
							
					susy_gg_syserr_PU->SetBinContent(i,j, -1);
					susy_gg_chan_syserr_jes[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_esf[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_scale[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_eleshape[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_jetshape[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_qcdshape[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_xs[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_lumi[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_isr[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_pdf[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_met[ih]->SetBinContent(i,j, -1); 
				}
				else{
					float noe = p_SUSYMass->GetBinContent(i,j);
					float noeisr = p_SUSYISR->GetBinContent(i,j);
					float sparticleMass = p_SUSYMass->GetXaxis()->GetBinCenter(i);
					float crosssection = p_crosssection_susy->GetBinContent( p_crosssection_susy->FindBin(sparticleMass) );

					if(ih%3 == 0){
						double nominal(0);
						for(unsigned binHT(0); binHT < 3; binHT++)
							nominal += susy_all_chan_rate_nom[ih+binHT]->GetBinContent(i,j)/noe;
			  		double pdfdiff(0); 
			  		for(unsigned is(0); is < 9; is++){
							double noepdf = p_SUSYMASS_pdf[is]->GetBinContent(i,j);
							double rescale(0);
							for(unsigned binHT(0); binHT < 3; binHT++)
								rescale += susy_all_chan_rate_pdf[ih+binHT][is]->GetBinContent(i,j)/noepdf;
							pdfdiff = pdfdiff > fabs(rescale - nominal)? pdfdiff : fabs(rescale - nominal);
			  		}
						for(unsigned binHT(0); binHT < 3; binHT++){
							if(susy_all_chan_rate_nom[ih+binHT]->GetBinContent(i,j) > 0)
								susy_all_chan_syserr_pdf[ih+binHT]->SetBinContent(i,j, pdfdiff/nominal); 
			  			else susy_all_chan_syserr_pdf[ih+binHT]->SetBinContent(i,j, -1);
							if(susy_wg_chan_rate_nom[ih+binHT]->GetBinContent(i,j) > 0)
								susy_wg_chan_syserr_pdf[ih+binHT]->SetBinContent(i,j, pdfdiff/nominal); 
			  			else susy_wg_chan_syserr_pdf[ih+binHT]->SetBinContent(i,j, -1);
							if(susy_gg_chan_rate_nom[ih+binHT]->GetBinContent(i,j) > 0)
								susy_gg_chan_syserr_pdf[ih+binHT]->SetBinContent(i,j, pdfdiff/nominal); 
			  			else susy_gg_chan_syserr_pdf[ih+binHT]->SetBinContent(i,j, -1);
						}
			
						double isrdiff(0);
						for(unsigned binHT(0); binHT < 3; binHT++)
							isrdiff += susy_all_chan_rate_isr[ih+binHT]->GetBinContent(i,j)/noeisr;
						isrdiff = fabs(isrdiff - nominal);
						for(unsigned binHT(0); binHT < 3; binHT++){
							if(susy_all_chan_rate_nom[ih+binHT]->GetBinContent(i,j) > 0)
								susy_all_chan_syserr_isr[ih+binHT]->SetBinContent(i,j, isrdiff/nominal); 
			  			else susy_all_chan_syserr_isr[ih+binHT]->SetBinContent(i,j, -1);
							if(susy_wg_chan_rate_nom[ih+binHT]->GetBinContent(i,j) > 0)
								susy_wg_chan_syserr_isr[ih+binHT]->SetBinContent(i,j, isrdiff/nominal); 
			  			else susy_wg_chan_syserr_isr[ih+binHT]->SetBinContent(i,j, -1);
							if(susy_gg_chan_rate_nom[ih+binHT]->GetBinContent(i,j) > 0)
								susy_gg_chan_syserr_isr[ih+binHT]->SetBinContent(i,j, isrdiff/nominal); 
			  			else susy_gg_chan_syserr_isr[ih+binHT]->SetBinContent(i,j, -1);
						} 
					}
					
					susy_all_chan_rate_nom[ih]->SetBinError(i,j, sqrt(susy_all_chan_rate_nom[ih]->GetBinContent(i,j))*35.9*crosssection/noe);
					susy_all_chan_rate_nom[ih]->SetBinContent(i,j, susy_all_chan_rate_nom[ih]->GetBinContent(i,j)*35.9*crosssection/noe); 
					susy_all_chan_rate_jesUp[ih]->SetBinContent(i,j, susy_all_chan_rate_jesUp[ih]->GetBinContent(i,j)*35.9*crosssection/noe); 
					susy_all_chan_rate_jesDown[ih]->SetBinContent(i,j, susy_all_chan_rate_jesDown[ih]->GetBinContent(i,j)*35.9*crosssection/noe);
					susy_all_chan_rate_esfUp[ih]->SetBinContent(i,j, susy_all_chan_rate_esfUp[ih]->GetBinContent(i,j)*35.9*crosssection/noe);
									
					susy_all_chan_syserr_jes[ih]->SetBinContent(i,j, max( fabs(susy_all_chan_rate_jesUp[ih]->GetBinContent(i,j)-susy_all_chan_rate_nom[ih]->GetBinContent(i,j)), fabs(susy_all_chan_rate_jesDown[ih]->GetBinContent(i,j)-susy_all_chan_rate_nom[ih]->GetBinContent(i,j)) ) ); 
					susy_all_chan_syserr_esf[ih]->SetBinContent(i,j, fabs( susy_all_chan_rate_esfUp[ih]->GetBinContent(i,j)-susy_all_chan_rate_nom[ih]->GetBinContent(i,j)) );
					susy_all_chan_syserr_scale[ih]->SetBinContent(i,j, -1);
					susy_all_chan_syserr_eleshape[ih]->SetBinContent(i,j, -1);
					susy_all_chan_syserr_jetshape[ih]->SetBinContent(i,j, -1);
					susy_all_chan_syserr_qcdshape[ih]->SetBinContent(i,j, -1);
					susy_all_chan_syserr_xs[ih]->SetBinContent(i,j, -1); 
					susy_all_chan_syserr_lumi[ih]->SetBinContent(i,j, 0.025*susy_all_chan_rate_nom[ih]->GetBinContent(i,j));   
					if(ih==0)susy_all_syserr_PU->SetBinContent(i,j, analysis_PU(i,j, p_lowPU_susy_pass, p_lowPU_susy_all, p_highPU_susy_pass, p_highPU_susy_all, p_PU_data));  


					susy_wg_chan_rate_nom[ih]->SetBinError(i,j, sqrt(susy_wg_chan_rate_nom[ih]->GetBinContent(i,j))*35.9*crosssection/noe);
					susy_wg_chan_rate_nom[ih]->SetBinContent(i,j, susy_wg_chan_rate_nom[ih]->GetBinContent(i,j)*35.9*crosssection/noe); 
					susy_wg_chan_rate_jesUp[ih]->SetBinContent(i,j, susy_wg_chan_rate_jesUp[ih]->GetBinContent(i,j)*35.9*crosssection/noe); 
					susy_wg_chan_rate_jesDown[ih]->SetBinContent(i,j, susy_wg_chan_rate_jesDown[ih]->GetBinContent(i,j)*35.9*crosssection/noe);
					susy_wg_chan_rate_esfUp[ih]->SetBinContent(i,j, susy_wg_chan_rate_esfUp[ih]->GetBinContent(i,j)*35.9*crosssection/noe);
									
					susy_wg_chan_syserr_jes[ih]->SetBinContent(i,j, max( fabs(susy_wg_chan_rate_jesUp[ih]->GetBinContent(i,j)-susy_wg_chan_rate_nom[ih]->GetBinContent(i,j)), fabs(susy_wg_chan_rate_jesDown[ih]->GetBinContent(i,j)-susy_wg_chan_rate_nom[ih]->GetBinContent(i,j)) ) ); 
					susy_wg_chan_syserr_esf[ih]->SetBinContent(i,j, fabs( susy_wg_chan_rate_esfUp[ih]->GetBinContent(i,j)-susy_wg_chan_rate_nom[ih]->GetBinContent(i,j)) );
					susy_wg_chan_syserr_scale[ih]->SetBinContent(i,j, -1);
					susy_wg_chan_syserr_eleshape[ih]->SetBinContent(i,j, -1);
					susy_wg_chan_syserr_jetshape[ih]->SetBinContent(i,j, -1);
					susy_wg_chan_syserr_qcdshape[ih]->SetBinContent(i,j, -1);
					susy_wg_chan_syserr_xs[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_lumi[ih]->SetBinContent(i,j, 0.025*susy_wg_chan_rate_nom[ih]->GetBinContent(i,j));   
					if(ih==0)susy_wg_syserr_PU->SetBinContent(i,j, analysis_PU(i,j, p_lowPU_susy_pass, p_lowPU_susy_all, p_highPU_susy_pass, p_highPU_susy_all, p_PU_data));  

					susy_gg_chan_rate_nom[ih]->SetBinError(i,j, sqrt(susy_gg_chan_rate_nom[ih]->GetBinContent(i,j))*35.9*crosssection/noe);
					susy_gg_chan_rate_nom[ih]->SetBinContent(i,j, susy_gg_chan_rate_nom[ih]->GetBinContent(i,j)*35.9*crosssection/noe); 
					susy_gg_chan_rate_jesUp[ih]->SetBinContent(i,j, susy_gg_chan_rate_jesUp[ih]->GetBinContent(i,j)*35.9*crosssection/noe); 
					susy_gg_chan_rate_jesDown[ih]->SetBinContent(i,j, susy_gg_chan_rate_jesDown[ih]->GetBinContent(i,j)*35.9*crosssection/noe);
					susy_gg_chan_rate_esfUp[ih]->SetBinContent(i,j, susy_gg_chan_rate_esfUp[ih]->GetBinContent(i,j)*35.9*crosssection/noe);
									
					susy_gg_chan_syserr_jes[ih]->SetBinContent(i,j, max( fabs(susy_gg_chan_rate_jesUp[ih]->GetBinContent(i,j)-susy_gg_chan_rate_nom[ih]->GetBinContent(i,j)), fabs(susy_gg_chan_rate_jesDown[ih]->GetBinContent(i,j)-susy_gg_chan_rate_nom[ih]->GetBinContent(i,j)) ) ); 
					susy_gg_chan_syserr_esf[ih]->SetBinContent(i,j, fabs( susy_gg_chan_rate_esfUp[ih]->GetBinContent(i,j)-susy_gg_chan_rate_nom[ih]->GetBinContent(i,j)) );
					susy_gg_chan_syserr_scale[ih]->SetBinContent(i,j, -1);
					susy_gg_chan_syserr_eleshape[ih]->SetBinContent(i,j, -1);
					susy_gg_chan_syserr_jetshape[ih]->SetBinContent(i,j, -1);
					susy_gg_chan_syserr_qcdshape[ih]->SetBinContent(i,j, -1);
					susy_gg_chan_syserr_xs[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_lumi[ih]->SetBinContent(i,j, 0.026*susy_gg_chan_rate_nom[ih]->GetBinContent(i,j));   
					if(ih==0)susy_gg_syserr_PU->SetBinContent(i,j, analysis_PU(i,j, p_lowPU_susy_pass, p_lowPU_susy_all, p_highPU_susy_pass, p_highPU_susy_all, p_PU_data));  
				}
			}
		} 
	} 

		for(unsigned i(1); i < p_SUSYMass->GetXaxis()->GetNbins() + 1; i++){
			for(unsigned j(1); j < p_SUSYMass->GetYaxis()->GetNbins() + 1; j++){
				float sparticleMass = p_SUSYMass->GetXaxis()->GetBinCenter(i);
				if( fabs(sparticleMass - 1700 ) < 10 && fabs( p_SUSYMass->GetYaxis()->GetBinCenter(j) - 1000) < 5){
					float noe = p_SUSYMass->GetBinContent(i,j);
					float crosssection = 0.5*p_crosssection_susy->GetBinContent( p_crosssection_susy->FindBin(sparticleMass) );
					p_susy_MET_signal_mg->Scale(35.8*crosssection/noe);
					p_susy_HT_signal_mg->Scale(35.8*crosssection/noe);
					p_susy_PhoEt_signal_mg->Scale(35.8*crosssection/noe);
					p_susy_MET_signal_eg->Scale(35.8*crosssection/noe);
					p_susy_HT_signal_eg->Scale(35.8*crosssection/noe);
					p_susy_PhoEt_signal_eg->Scale(35.8*crosssection/noe);
			}
		}
	}

	double largeunc(0),smallunc(2);
	for(unsigned ih(0); ih < NBIN*2; ih++){
		for(unsigned i(1); i < susy_wg_chan_rate_nom[ih]->GetXaxis()->GetNbins() + 1; i++){
			for(unsigned j(1); j < susy_wg_chan_rate_nom[ih]->GetYaxis()->GetNbins() + 1; j++){
				if(p_SUSYMass->GetBinContent(i,j) > 1000){
					if(susy_all_chan_syserr_pdf[ih]->GetBinContent(i,j) > 0.1)std::cout << " !!" << std::endl;
					std::cout << "pdf ch " << ih << " mass " << susy_wg_chan_rate_nom[ih]->GetXaxis()->GetBinCenter(i) << " " << susy_wg_chan_rate_nom[ih]->GetYaxis()->GetBinCenter(j) << " " << susy_all_chan_syserr_pdf[ih]->GetBinContent(i,j) << std::endl;
					if(susy_all_chan_syserr_pdf[ih]->GetBinContent(i,j) > largeunc)largeunc=susy_all_chan_syserr_pdf[ih]->GetBinContent(i,j);
					if(susy_all_chan_syserr_pdf[ih]->GetBinContent(i,j) >= 0 &&susy_all_chan_syserr_pdf[ih]->GetBinContent(i,j) < smallunc)smallunc = largeunc=susy_all_chan_syserr_pdf[ih]->GetBinContent(i,j);
				} 
			}
		} 
	} 
	std::cout << "large " << largeunc << std::endl;
	std::cout << "small " << smallunc << std::endl;

	double isrlargeunc(0);
	double isrsmallunc(2);
	for(unsigned ih(0); ih < NBIN*2; ih++){
		for(unsigned i(1); i < susy_wg_chan_rate_nom[ih]->GetXaxis()->GetNbins() + 1; i++){
			for(unsigned j(1); j < susy_wg_chan_rate_nom[ih]->GetYaxis()->GetNbins() + 1; j++){
				if(p_SUSYMass->GetBinContent(i,j) > 1000){
					if(susy_all_chan_rate_nom[ih]->GetBinContent(i,j) <=0)std::cout << "isr ch " << ih << " mass " << susy_wg_chan_rate_nom[ih]->GetXaxis()->GetBinCenter(i) << " " << susy_wg_chan_rate_nom[ih]->GetYaxis()->GetBinCenter(j) << " -" << std::endl;
					else{
					if(susy_all_chan_syserr_isr[ih]->GetBinContent(i,j) > 0.1)std::cout << " !!" << std::endl;
					std::cout << "isr ch " << ih << " mass " << susy_wg_chan_rate_nom[ih]->GetXaxis()->GetBinCenter(i) << " " << susy_wg_chan_rate_nom[ih]->GetYaxis()->GetBinCenter(j) << " " << susy_all_chan_syserr_isr[ih]->GetBinContent(i,j) << std::endl;
					if(susy_all_chan_syserr_isr[ih]->GetBinContent(i,j) > isrlargeunc)isrlargeunc=susy_all_chan_syserr_isr[ih]->GetBinContent(i,j);
					if(susy_all_chan_syserr_isr[ih]->GetBinContent(i,j) >= 0 && susy_all_chan_syserr_isr[ih]->GetBinContent(i,j) < isrsmallunc)isrsmallunc = susy_all_chan_syserr_isr[ih]->GetBinContent(i,j);
					}
				} 
			}
		} 
	} 
	std::cout << "isrlarge " << isrlargeunc << std::endl;
	std::cout << "isrsmall " << isrsmallunc << std::endl;

	double metlargeunc(0);
	double metsmallunc(2);
	for(unsigned ih(0); ih < NBIN*2; ih++){
		for(unsigned i(1); i < susy_wg_chan_rate_nom[ih]->GetXaxis()->GetNbins() + 1; i++){
			for(unsigned j(1); j < susy_wg_chan_rate_nom[ih]->GetYaxis()->GetNbins() + 1; j++){
				if(p_SUSYMass->GetBinContent(i,j) > 1000){
					if(susy_all_chan_rate_nom[ih]->GetBinContent(i,j) <=0)std::cout << "met ch " << ih << " mass " << susy_wg_chan_rate_nom[ih]->GetXaxis()->GetBinCenter(i) << " " << susy_wg_chan_rate_nom[ih]->GetYaxis()->GetBinCenter(j) << " -" << std::endl;
					else{
					if(susy_all_chan_syserr_met[ih]->GetBinContent(i,j) > 0.1)std::cout << " !!" << std::endl;
					std::cout << "met ch " << ih << " mass " << susy_wg_chan_rate_nom[ih]->GetXaxis()->GetBinCenter(i) << " " << susy_wg_chan_rate_nom[ih]->GetYaxis()->GetBinCenter(j) << " " << susy_all_chan_syserr_met[ih]->GetBinContent(i,j) << std::endl;
					if(susy_all_chan_syserr_met[ih]->GetBinContent(i,j) > metlargeunc)metlargeunc=susy_all_chan_syserr_met[ih]->GetBinContent(i,j);
					if(susy_all_chan_syserr_met[ih]->GetBinContent(i,j) >= 0 && susy_all_chan_syserr_met[ih]->GetBinContent(i,j) < metsmallunc)metsmallunc = susy_all_chan_syserr_met[ih]->GetBinContent(i,j);
					}
				} 
			}
		} 
	} 
	std::cout << "metlarge " << metlargeunc << std::endl;
	std::cout << "metsmall " << metsmallunc << std::endl;

	write_histo();
  xSecFile.Close();
}


