#include "analysis_PU.C"
#include "analysis_SUSY.h"

#define NPDF 9

void fill_T6WG(){//main  

	SetSignalConfig();
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);
	esfScaleFactor  objectESF;

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

	TFile xSecFile("../cross/susyCrossSection.root");
	TH1D *p_crosssection = (TH1D*)xSecFile.Get("p_squarkxSec");

	TChain *datachain = new TChain("signalTree");
	datachain->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_mgsignal_MuonEG_FullEcal.root");
	TH1D *p_PU_data = new TH1D("p_PU_data",";N_{vtx};",100,0,100); 
  datachain->Draw("nVertex >> p_PU_data");
	p_PU_data->Scale(1.0/p_PU_data->Integral(1,101));

	init_histo("T5WG","test_T6WG.root",2,27, 775.0, 2125.0, 420, 2.5, 2102.5);

	std::ostringstream histname;
	//**************   T5WG  ***************************//
	TFile *file = TFile::Open("/uscms_data/d3/mengleis/FullStatusOct/resTree_T6WG_debug.root");
  TTree *tree = (TTree*)file->Get("SUSYtree");
	float Mass1(0);
	float Mass2(0);
	float Mgluino(0);
  float Mchargino(0);
  float Mneutralino(0);
	float ISRPt(0);
	int   nVertex(0);
	std::vector<float> *ScaleSystWeight = 0;
  tree->SetBranchAddress("Mass1",      &Mass1);
  tree->SetBranchAddress("Mass2",      &Mass2);
	tree->SetBranchAddress("MsGsQ",      &Mgluino);  
  tree->SetBranchAddress("Mchargino",  &Mchargino);
  tree->SetBranchAddress("Mneutralino",&Mneutralino);
  tree->SetBranchAddress("ISRPt",      &ISRPt);
  tree->SetBranchAddress("nVertex",    &nVertex);
	tree->SetBranchAddress("ScaleSystWeight",&ScaleSystWeight);

	for(unsigned ievt(0); ievt < tree->GetEntries(); ievt++){
		tree->GetEntry(ievt);

		p_SUSYMass->Fill(Mass1, Mass2);
		p_SUSYGEN[0]->Fill(Mass1, Mass2);	
		if(Mchargino> 0 && Mneutralino > 0)p_SUSYGEN[1]->Fill(Mass1, Mass2); 
		else if(Mchargino <= 0 && Mneutralino > 0)p_SUSYGEN[2]->Fill(Mass1, Mass2);
		else if(Mchargino > 0 && Mneutralino <= 0)p_SUSYGEN[3]->Fill(Mass1, Mass2);
		for(unsigned i(0); i < NPDF; i++)p_SUSYMASS_pdf[i]->Fill(Mass1, Mass2, (*ScaleSystWeight)[i]);

		if(nVertex < 20)p_lowPU_all->Fill(Mass1, Mass2, 1);
		else p_highPU_all->Fill(Mass1, Mass2, 1);

		double reweightF = isrW(ISRPt);
		p_SUSYISR->Fill(Mass1, Mass2, reweightF);

		//if(Mgluino <= 0 || (Mchargino <= 0 && Mneutralino <= 0)){
		if((Mchargino <= 0 && Mneutralino <= 0)){
			std::cout << "no sG mass" << std::endl;
			exit (EXIT_FAILURE);
		}
	}

  TTree *mgtree = (TTree*)file->Get("mgTree");
  float Mass1_mg(0);
  float Mass2_mg(0);
	float Mgluino_mg(0);
  float Mchargino_mg(0);
  float Mneutralino_mg(0);
  float phoEt_mg(0);
  float phoEta_mg(0);
  float lepPt_mg(0);
  float lepEta_mg(0);
  float sigMT_mg(0);
  float sigMET_mg(0);
  float genMET_mg(0);
  float HT_mg(0);
	float ISRPt_mg(0);
	int   nVertex_mg(0);
	float sigMETJESup_mg(0);
	float sigMETJESdo_mg(0);
	float sigMTJESup_mg(0);
	float sigMTJESdo_mg(0);
	float HTJESup_mg(0);
	float HTJESdo_mg(0);
	std::vector<float> *ScaleSystWeight_mg = 0;
	mgtree->SetBranchAddress("Mass1",      &Mass1_mg);
	mgtree->SetBranchAddress("Mass2",      &Mass2_mg);
	mgtree->SetBranchAddress("MsGsQ",      &Mgluino_mg);
  mgtree->SetBranchAddress("Mchargino",  &Mchargino_mg);
  mgtree->SetBranchAddress("Mneutralino",&Mneutralino_mg);
  mgtree->SetBranchAddress("phoEt",      &phoEt_mg);
  mgtree->SetBranchAddress("phoEta",     &phoEta_mg);
  mgtree->SetBranchAddress("lepPt",      &lepPt_mg);
  mgtree->SetBranchAddress("lepEta",     &lepEta_mg);
  mgtree->SetBranchAddress("sigMT",      &sigMT_mg);
  mgtree->SetBranchAddress("sigMET",     &sigMET_mg);
  mgtree->SetBranchAddress("genMET",     &genMET_mg);
  mgtree->SetBranchAddress("HT",         &HT_mg);
	mgtree->SetBranchAddress("nVertex",    &nVertex_mg);
  mgtree->SetBranchAddress("ISRPt",      &ISRPt_mg);
	mgtree->SetBranchAddress("sigMETJESup",&sigMETJESup_mg);
	mgtree->SetBranchAddress("sigMETJESdo",&sigMETJESdo_mg);
	mgtree->SetBranchAddress("sigMTJESup", &sigMTJESup_mg);
	mgtree->SetBranchAddress("sigMTJESdo", &sigMTJESdo_mg);
	mgtree->SetBranchAddress("HTJESup",    &HTJESup_mg);
	mgtree->SetBranchAddress("HTJESdo",    &HTJESdo_mg);
	mgtree->SetBranchAddress("ScaleSystWeight", &ScaleSystWeight_mg);

	for(unsigned ievt(0); ievt < mgtree->GetEntries(); ievt++){
		mgtree->GetEntry(ievt);

		/** cut flow *****/
		if(phoEt_mg < 35 || lepPt_mg < 25)continue;
		if(fabs(phoEta_mg) > 1.4442 || fabs(lepEta_mg) > 2.5)continue;

		double phoESF(0),   lepESF(0);	
		double phoESFUp(0), lepESFUp(0);	
		phoESF = objectESF.getPhotonESF(phoEt_mg, phoEta_mg)*objectESF.getFastMuonEGTRGESF(phoEt_mg, lepPt_mg);
		lepESF = objectESF.getFastMuonESF(lepPt_mg,lepEta_mg);
		double s_mu_error  =objectESF.getFastMuonESFError(lepPt_mg,lepEta_mg);
    double s_pho_error =objectESF.getPhotonESFError(phoEt_mg, phoEta_mg)*objectESF.getFastMuonEGTRGESF(phoEt_mg, lepPt_mg);
    double s_trg_error =objectESF.getFastMuonEGTRGESFError(phoEt_mg, lepPt_mg)*objectESF.getPhotonESF(phoEt_mg, phoEta_mg);
		phoESFUp = phoESF + sqrt(pow(s_pho_error,2) + pow(s_trg_error,2));
		lepESFUp = lepESF + s_mu_error;
		double scalefactor = phoESF*lepESF;

		double reweightF = isrW(ISRPt_mg);

		if(sigMT_mg > 100){
			int genMETBinIndex(-1);
			genMETBinIndex = Bin.findSignalBin(genMET_mg, HT_mg, phoEt_mg);
			if(genMETBinIndex >= 0)genMETBinIndex += NBIN;
			if(genMETBinIndex >=0)susy_chan_rate_genMET[genMETBinIndex]->Fill(Mass1_mg, Mass2_mg);
		}	

		if(sigMET_mg < 120 || sigMT_mg < 100)continue;
		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET_mg, HT_mg, phoEt_mg);
		if(SigBinIndex >= 0)SigBinIndex += NBIN;
		int jesupBinIndex(-1);	
		jesupBinIndex = Bin.findSignalBin(sigMETJESup_mg, HTJESup_mg, phoEt_mg);
		if(jesupBinIndex >= 0)jesupBinIndex += NBIN; 
		int jesdoBinIndex(-1);
		jesdoBinIndex = Bin.findSignalBin(sigMETJESdo_mg, HTJESdo_mg, phoEt_mg);
		if(jesdoBinIndex >= 0)jesdoBinIndex += NBIN;

		if(Mass1_mg > 0 && Mass2_mg > 0){ 

			if(SigBinIndex >=0){
				susy_chan_rate_nom[0][SigBinIndex]->Fill( Mass1_mg, Mass2_mg, scalefactor);
				susy_chan_rate_control[0][SigBinIndex]->Fill( Mass1_mg, Mass2_mg);
	
				if(Mchargino_mg > 0 && Mneutralino_mg > 0){
					susy_chan_rate_nom[1][SigBinIndex]->Fill( Mass1_mg, Mass2_mg, scalefactor);
					susy_chan_rate_control[1][SigBinIndex]->Fill( Mass1_mg, Mass2_mg);
				}
				else if(Mchargino_mg <= 0 && Mneutralino_mg > 0){
					susy_chan_rate_nom[2][SigBinIndex]->Fill( Mass1_mg, Mass2_mg, scalefactor);
					susy_chan_rate_control[2][SigBinIndex]->Fill( Mass1_mg, Mass2_mg);
				}
				else if(Mchargino_mg > 0 && Mneutralino_mg <= 0){
					susy_chan_rate_nom[3][SigBinIndex]->Fill( Mass1_mg, Mass2_mg, scalefactor);
					susy_chan_rate_control[3][SigBinIndex]->Fill( Mass1_mg, Mass2_mg);
				}

				susy_chan_rate_phoesfUp[SigBinIndex]->Fill( Mass1_mg, Mass2_mg, phoESFUp*lepESF); 
				susy_chan_rate_lepesfUp[SigBinIndex]->Fill( Mass1_mg, Mass2_mg, phoESF*lepESFUp); 
				susy_chan_rate_isr[SigBinIndex]->Fill(Mass1_mg, Mass2_mg, scalefactor*reweightF);	
	
				if(nVertex_mg < 20)p_lowPU_pass->Fill( Mass1_mg, Mass2_mg, 1);
				else p_highPU_pass->Fill( Mass1_mg, Mass2_mg, 1);

			  for(unsigned is(0); is < 9; is++)
					susy_chan_rate_pdf[SigBinIndex][is]->Fill(Mass1_mg, Mass2_mg, scalefactor*(*ScaleSystWeight_mg)[is]);
			}
			if(jesupBinIndex >=0){
				susy_chan_rate_jesUp[jesupBinIndex]->Fill( Mass1_mg, Mass2_mg, scalefactor); 
			}
			if(jesdoBinIndex >=0){
				susy_chan_rate_jesDown[jesdoBinIndex]->Fill( Mass1_mg, Mass2_mg, scalefactor);   
			}	
		}
	}



  TTree *egtree = (TTree*)file->Get("egTree");
  float Mass1_eg(0);
  float Mass2_eg(0);
	float Mgluino_eg(0);
  float Mchargino_eg(0);
  float Mneutralino_eg(0);
  float phoEt_eg(0);
  float phoEta_eg(0);
  float lepPt_eg(0);
  float lepEta_eg(0);
  float sigMT_eg(0);
  float sigMET_eg(0);
  float genMET_eg(0);
  float HT_eg(0);
	float ISRPt_eg(0);
	int   nVertex_eg(0);
	float sigMETJESup_eg(0);
	float sigMETJESdo_eg(0);
	float sigMTJESup_eg(0);
	float sigMTJESdo_eg(0);
	float HTJESup_eg(0);
	float HTJESdo_eg(0);
	std::vector<float> *ScaleSystWeight_eg = 0;
	egtree->SetBranchAddress("Mass1",      &Mass1_eg);
	egtree->SetBranchAddress("Mass2",      &Mass2_eg);
	egtree->SetBranchAddress("MsGsQ",      &Mgluino_eg);
  egtree->SetBranchAddress("Mchargino",  &Mchargino_eg);
  egtree->SetBranchAddress("Mneutralino",&Mneutralino_eg);
  egtree->SetBranchAddress("phoEt",      &phoEt_eg);
  egtree->SetBranchAddress("phoEta",     &phoEta_eg);
  egtree->SetBranchAddress("lepPt",      &lepPt_eg);
  egtree->SetBranchAddress("lepEta",     &lepEta_eg);
  egtree->SetBranchAddress("sigMT",      &sigMT_eg);
  egtree->SetBranchAddress("sigMET",     &sigMET_eg);
  egtree->SetBranchAddress("genMET",     &genMET_eg);
  egtree->SetBranchAddress("HT",         &HT_eg);
	egtree->SetBranchAddress("nVertex",    &nVertex_eg);
  egtree->SetBranchAddress("ISRPt",      &ISRPt_eg);
	egtree->SetBranchAddress("sigMETJESup",&sigMETJESup_eg);
	egtree->SetBranchAddress("sigMETJESdo",&sigMETJESdo_eg);
	egtree->SetBranchAddress("sigMTJESup", &sigMTJESup_eg);
	egtree->SetBranchAddress("sigMTJESdo", &sigMTJESdo_eg);
	egtree->SetBranchAddress("HTJESup",    &HTJESup_eg);
	egtree->SetBranchAddress("HTJESdo",    &HTJESdo_eg);
	egtree->SetBranchAddress("ScaleSystWeight", &ScaleSystWeight_eg);

	for(unsigned ievt(0); ievt < egtree->GetEntries(); ievt++){
		egtree->GetEntry(ievt);

		/** cut flow *****/
		if(phoEt_eg < 35 || lepPt_eg < 25)continue;
		if(fabs(phoEta_eg) > 1.4442 || fabs(lepEta_eg) > 2.5)continue;

		double phoESF(0),   lepESF(0);	
		double phoESFUp(0), lepESFUp(0);	
		phoESF = objectESF.getPhotonESF(phoEt_eg,phoEta_eg)*objectESF.getFastegPhotonTRGESF(phoEt_eg,phoEta_eg); 
		lepESF = objectESF.getFastElectronESF(lepPt_eg,lepEta_eg)*objectESF.getFastElectronTRGESF(lepPt_eg,lepEta_eg);
		double s_ele_error = objectESF.getFastElectronESFError(lepPt_eg,lepEta_eg)*objectESF.getFastElectronTRGESF(lepPt_eg,lepEta_eg);
		double s_eletrg_error = objectESF.getFastElectronTRGESFError(lepPt_eg,lepEta_eg)*objectESF.getFastElectronESF(lepPt_eg,lepEta_eg);
		double s_pho_error = objectESF.getPhotonESFError(phoEt_eg,phoEta_eg)*objectESF.getFastegPhotonTRGESF(phoEt_eg,phoEta_eg);
		double s_photrg_error = objectESF.getFastegPhotonTRGESFError(phoEt_eg,phoEta_eg)*objectESF.getPhotonESF(phoEt_eg,phoEta_eg);
		phoESFUp = phoESF + sqrt( pow(s_pho_error,2) + pow(s_photrg_error,2));
		lepESFUp = lepESF + sqrt( pow(s_ele_error,2) + pow(s_eletrg_error,2));
		double scalefactor = phoESF*lepESF;

		double reweightF = isrW(ISRPt_eg);

		if(sigMT_eg > 100){
			int genMETBinIndex(-1);
			genMETBinIndex = Bin.findSignalBin(genMET_eg, HT_eg, phoEt_eg);
			if(genMETBinIndex >=0)susy_chan_rate_genMET[genMETBinIndex]->Fill(Mass1_eg, Mass2_eg);
		}	

		if(sigMET_eg < 120 || sigMT_eg < 100)continue;
		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET_eg, HT_eg, phoEt_eg);
		int jesupBinIndex(-1);	
		jesupBinIndex = Bin.findSignalBin(sigMETJESup_eg, HTJESup_eg, phoEt_eg);
		int jesdoBinIndex(-1);
		jesdoBinIndex = Bin.findSignalBin(sigMETJESdo_eg, HTJESdo_eg, phoEt_eg);

		if(Mass1_eg > 0 && Mass2_eg > 0){ 

			if(SigBinIndex >=0){
				susy_chan_rate_nom[0][SigBinIndex]->Fill( Mass1_eg, Mass2_eg, scalefactor);
				susy_chan_rate_control[0][SigBinIndex]->Fill( Mass1_eg, Mass2_eg);
	
				if(Mchargino_eg > 0 && Mneutralino_eg > 0){
					susy_chan_rate_nom[1][SigBinIndex]->Fill( Mass1_eg, Mass2_eg, scalefactor);
					susy_chan_rate_control[1][SigBinIndex]->Fill( Mass1_eg, Mass2_eg);
				}
				else if(Mchargino_eg <= 0 && Mneutralino_eg > 0){
					susy_chan_rate_nom[2][SigBinIndex]->Fill( Mass1_eg, Mass2_eg, scalefactor);
					susy_chan_rate_control[2][SigBinIndex]->Fill( Mass1_eg, Mass2_eg);
				}
				else if(Mchargino_eg > 0 && Mneutralino_eg <= 0){
					susy_chan_rate_nom[3][SigBinIndex]->Fill( Mass1_eg, Mass2_eg, scalefactor);
					susy_chan_rate_control[3][SigBinIndex]->Fill( Mass1_eg, Mass2_eg);
				}

				susy_chan_rate_phoesfUp[SigBinIndex]->Fill( Mass1_eg, Mass2_eg, phoESFUp*lepESF); 
				susy_chan_rate_lepesfUp[SigBinIndex]->Fill( Mass1_eg, Mass2_eg, phoESF*lepESFUp); 
				susy_chan_rate_isr[SigBinIndex]->Fill(Mass1_eg, Mass2_eg, scalefactor*reweightF);	
	
				if(nVertex_eg < 20)p_lowPU_pass->Fill( Mass1_eg, Mass2_eg, 1);
				else p_highPU_pass->Fill( Mass1_eg, Mass2_eg, 1);

			  for(unsigned is(0); is < 9; is++)
					susy_chan_rate_pdf[SigBinIndex][is]->Fill(Mass1_eg, Mass2_eg, scalefactor*(*ScaleSystWeight_eg)[is]);
			}
			if(jesupBinIndex >=0){
				susy_chan_rate_jesUp[jesupBinIndex]->Fill( Mass1_eg, Mass2_eg, scalefactor); 
			}
			if(jesdoBinIndex >=0){
				susy_chan_rate_jesDown[jesdoBinIndex]->Fill( Mass1_eg, Mass2_eg, scalefactor);   
			}	
		}
	}

//************ Normalize Hists **************** //
	std::vector< std::pair< int, int> > metIndex;
	std::vector< std::pair< int, int> > htIndex;
	metIndex.clear();
	htIndex.clear();
	for(int i(0); i < 36; i++)
		metIndex.push_back(make_pair(i, ((i%9)/3)));
	for(int i(0); i < 36; i++)
		std::cout << "index " << metIndex[i].first << " " << metIndex[i].second << std::endl;
	for(int i(0); i < 36; i++)
		htIndex.push_back(make_pair(i, i/3));

	for(unsigned i(1); i < p_SUSYMass->GetXaxis()->GetNbins() + 1; i++){
		for(unsigned j(1); j < p_SUSYMass->GetYaxis()->GetNbins() + 1; j++){
			if(p_SUSYMass->GetBinContent(i,j) < 1000){
				p_SUSYMass->SetBinContent(i,j, -1);
				continue;
			}

			susy_syserr_PU->SetBinContent(i,j, analysis_PU(i,j, p_lowPU_pass, p_lowPU_all, p_highPU_pass, p_highPU_all, p_PU_data)); 
			/*******  process the MET uncertainty: add bin 1-3, 4-6, 7-9, etc **********/
			for(unsigned ih(0); ih < NBIN*2; ih++){
				// MET unc //
				double gen_met(0);
				double pf_met(0);
				for(int iM(0); iM < NBIN*2; iM++){
					if(metIndex[iM].second != metIndex[ih].second)continue;
					gen_met += susy_chan_rate_genMET[iM]->GetBinContent(i,j);
					pf_met  += susy_chan_rate_nom[0][iM]->GetBinContent(i,j);
				}
				double sys_genMET(-1);
				if(pf_met >0)sys_genMET = fabs((gen_met-pf_met)/2)/pf_met; 
				susy_chan_syserr_met[ih]->SetBinContent(i,j,sys_genMET);

				// Scale unc //
				double nom_pdf(0);
				double shift_pdf[NPDF] = {0};
				double shift_isr(0);
				for(int iH(0); iH < NBIN*2; iH++){
					if(htIndex[iH].second != htIndex[ih].second)continue;
					nom_pdf += susy_chan_rate_nom[0][iH]->GetBinContent(i,j)/p_SUSYMass->GetBinContent(i,j);
					shift_isr += susy_chan_rate_isr[iH]->GetBinContent(i,j)/p_SUSYISR->GetBinContent(i,j);
					for(unsigned is(0); is < NPDF; is++){
						shift_pdf[is] += susy_chan_rate_pdf[iH][is]->GetBinContent(i,j)/p_SUSYMASS_pdf[is]->GetBinContent(i,j);
					}
				}
				for(unsigned is(0); is < NPDF; is++)shift_pdf[is] = fabs(shift_pdf[is] - nom_pdf); 
				double sys_pdf(-1);
				if(nom_pdf > 0)sys_pdf = (*std::max_element(shift_pdf, shift_pdf+NPDF))/nom_pdf;
				susy_chan_syserr_pdf[ih]->SetBinContent(i,j,sys_pdf);
				double sys_isr(-1);
				if(nom_pdf > 0)sys_isr = fabs(shift_isr - nom_pdf)/nom_pdf;
				susy_chan_syserr_isr[ih]->SetBinContent(i,j, sys_isr);	
			}

			/******  process other uncertainties **************/
			for(unsigned ih(0); ih < NBIN*2; ih++){
				double nom_yield = susy_chan_rate_nom[0][ih]->GetBinContent(i,j);
				if(nom_yield <= 0)continue;
	
				float noe = p_SUSYMass->GetBinContent(i,j);
				float sparticleMass = p_SUSYMass->GetXaxis()->GetBinCenter(i);
				float crosssection = p_crosssection->GetBinContent( p_crosssection->FindBin(sparticleMass) );
		
				double err_jesUp = fabs( susy_chan_rate_jesUp[ih]->GetBinContent(i,j) - nom_yield); 
				double err_jesDo = fabs( susy_chan_rate_jesDown[ih]->GetBinContent(i,j) - nom_yield);
				susy_chan_syserr_jes[ih]->SetBinContent(i,j, max( err_jesUp, err_jesDo)/nom_yield );
				susy_chan_syserr_phoesf[ih]->SetBinContent(i,j, fabs(susy_chan_rate_phoesfUp[ih]->GetBinContent(i,j) - nom_yield)/nom_yield);
				susy_chan_syserr_lepesf[ih]->SetBinContent(i,j, fabs(susy_chan_rate_lepesfUp[ih]->GetBinContent(i,j) - nom_yield)/nom_yield);
				susy_chan_syserr_scale[ih]->SetBinContent(i,j, -1);
				susy_chan_syserr_eleshape[ih]->SetBinContent(i,j, -1);
				susy_chan_syserr_jetshape[ih]->SetBinContent(i,j, -1);
				susy_chan_syserr_qcdshape[ih]->SetBinContent(i,j, -1);			
				susy_chan_syserr_xs[ih]->SetBinContent(i,j, -1);
				susy_chan_syserr_lumi[ih]->SetBinContent(i,j, 0.025);
			
				for(unsigned ich(0); ich < 4; ich++){	
					double noe_ch = p_SUSYGEN[ich]->GetBinContent(i,j);	
					susy_chan_rate_nom[ich][ih]->SetBinContent(i,j, susy_chan_rate_nom[ich][ih]->GetBinContent(i,j)*35.9*crosssection/noe_ch);
				}
			}
		}
	} 

	double low_jes(0), high_jes(0);
	double low_phoesf(0), high_phoesf(0);
	double low_lepesf(0), high_lepesf(0);
	double low_isr(0), high_isr(0);
	double low_pdf(0), high_pdf(0);
	double low_met(0), high_met(0);
	for(unsigned i(1); i < p_SUSYMass->GetXaxis()->GetNbins() + 1; i++){
		for(unsigned j(1); j < p_SUSYMass->GetYaxis()->GetNbins() + 1; j++){
			if(p_SUSYMass->GetBinContent(i,j) < 0)
				continue;
				
			int totalevt(0);
			double ave_jes(0), ave_phoesf(0), ave_lepesf(0), ave_isr(0), ave_pdf(0), ave_met(0); 	
			for(unsigned ih(0); ih < NBIN*2; ih++){
				totalevt += susy_chan_rate_control[0][ih]->GetBinContent(i,j);
				ave_jes  += susy_chan_syserr_jes[ih]->GetBinContent(i,j)*susy_chan_rate_control[0][ih]->GetBinContent(i,j);
				ave_phoesf += susy_chan_syserr_phoesf[ih]->GetBinContent(i,j)*susy_chan_rate_control[0][ih]->GetBinContent(i,j);
				ave_lepesf += susy_chan_syserr_lepesf[ih]->GetBinContent(i,j)*susy_chan_rate_control[0][ih]->GetBinContent(i,j);
				ave_isr += susy_chan_syserr_isr[ih]->GetBinContent(i,j)*susy_chan_rate_control[0][ih]->GetBinContent(i,j);
				ave_pdf += susy_chan_syserr_pdf[ih]->GetBinContent(i,j)*susy_chan_rate_control[0][ih]->GetBinContent(i,j);
				ave_met += susy_chan_syserr_met[ih]->GetBinContent(i,j)*susy_chan_rate_control[0][ih]->GetBinContent(i,j);
				if(susy_chan_rate_control[0][ih]->GetBinContent(i,j) > 2){
					if(susy_chan_syserr_jes[ih]->GetBinContent(i,j) < low_jes)low_jes = susy_chan_syserr_jes[ih]->GetBinContent(i,j);	
					if(susy_chan_syserr_jes[ih]->GetBinContent(i,j) > high_jes)high_jes = susy_chan_syserr_jes[ih]->GetBinContent(i,j);
					if(susy_chan_syserr_phoesf[ih]->GetBinContent(i,j) < low_phoesf)low_phoesf = susy_chan_syserr_phoesf[ih]->GetBinContent(i,j);	
					if(susy_chan_syserr_phoesf[ih]->GetBinContent(i,j) > high_phoesf)high_phoesf = susy_chan_syserr_phoesf[ih]->GetBinContent(i,j);
					if(susy_chan_syserr_lepesf[ih]->GetBinContent(i,j) < low_lepesf)low_lepesf = susy_chan_syserr_lepesf[ih]->GetBinContent(i,j);	
					if(susy_chan_syserr_lepesf[ih]->GetBinContent(i,j) > high_lepesf)high_lepesf = susy_chan_syserr_lepesf[ih]->GetBinContent(i,j);
					if(susy_chan_syserr_isr[ih]->GetBinContent(i,j) < low_isr)low_isr = susy_chan_syserr_isr[ih]->GetBinContent(i,j);	
					if(susy_chan_syserr_isr[ih]->GetBinContent(i,j) > high_isr)high_isr = susy_chan_syserr_isr[ih]->GetBinContent(i,j);
					if(susy_chan_syserr_pdf[ih]->GetBinContent(i,j) < low_pdf)low_pdf = susy_chan_syserr_pdf[ih]->GetBinContent(i,j);	
					if(susy_chan_syserr_pdf[ih]->GetBinContent(i,j) > high_pdf)high_pdf = susy_chan_syserr_pdf[ih]->GetBinContent(i,j);
					if(susy_chan_syserr_met[ih]->GetBinContent(i,j) < low_met)low_met = susy_chan_syserr_met[ih]->GetBinContent(i,j);	
					if(susy_chan_syserr_met[ih]->GetBinContent(i,j) > high_met)high_met = susy_chan_syserr_met[ih]->GetBinContent(i,j);
				}
			}

			for(unsigned ih(0); ih < NBIN*2; ih++){
				if(susy_chan_rate_control[0][ih]->GetBinContent(i,j) == 0){
					susy_chan_syserr_jes[ih]->SetBinContent(i,j, ave_jes/totalevt);
					susy_chan_syserr_phoesf[ih]->SetBinContent(i,j, ave_phoesf/totalevt);
					susy_chan_syserr_lepesf[ih]->SetBinContent(i,j, ave_lepesf/totalevt);
					susy_chan_syserr_isr[ih]->SetBinContent(i,j, ave_isr/totalevt);
					susy_chan_syserr_pdf[ih]->SetBinContent(i,j, ave_pdf/totalevt);
					susy_chan_syserr_met[ih]->SetBinContent(i,j, ave_met/totalevt);
				}
			}
		}
	}	 
	std::cout << "JES " << low_jes << "-" << high_jes << std::endl;
	std::cout << "phoESF " << low_phoesf << "-" << high_phoesf<< std::endl;
	std::cout << "lepESF " << low_lepesf << "-" << high_lepesf<< std::endl;
	std::cout << "isrESF " << low_isr << "-" << high_isr<< std::endl;
	std::cout << "pdfESF " << low_pdf << "-" << high_pdf<< std::endl;
	std::cout << "metESF " << low_met << "-" << high_met<< std::endl;

	write_histo();
  xSecFile.Close();
}


