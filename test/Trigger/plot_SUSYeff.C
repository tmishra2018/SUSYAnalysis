#include "../legacy/Result/analysis_PU.C"
#include "../legacy/Result/analysis_SUSY.h"
void plot_SUSYeff(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

	TFile *e_trg_file = TFile::Open("../../sf/diphoton_eleg.root","read");
	TFile *g_trg_file = TFile::Open("../../sf/diphoton_pholeg.root","read");
	TFile *mg_trg_file= TFile::Open("../../sf/muonphoton_trigger.root","read");	
	
	TProfile *e_trg_hist = (TProfile*)e_trg_file->Get("Traileff_MC");
	TProfile *g_trg_hist = (TProfile*)g_trg_file->Get("Leadeff_MC");

	TH2F *p_mgEff = (TH2F*)mg_trg_file->Get("p_mgEff");
	TH2F *p_mgESF = (TH2F*)mg_trg_file->Get("p_mgESF");

	long total_mg(0);
	double total_effmg(0);
	long total_eg(0);
	double total_effeg(0);

  TChain *mgtree_susy;
  mgtree_susy = new TChain("mgTree","mgTree");
  mgtree_susy->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_T5WG_string.root");
  float phoEt_susy_mg(0);
  float phoEta_susy_mg(0);
  float lepPt_susy_mg(0);
  float lepEta_susy_mg(0);
  float sigMT_susy_mg(0);
  float sigMET_susy_mg(0);
  float HT_susy_mg(0);
	int   nVertex_susy_mg(0);
	float gluinoMass_susy_mg(0);
  float charginoMass_susy_mg(0);
  float neutralinoMass_susy_mg(0);
  float Mass1_susy_mg(0);
  float Mass2_susy_mg(0);
  mgtree_susy->SetBranchAddress("phoEt",      &phoEt_susy_mg);
  mgtree_susy->SetBranchAddress("phoEta",     &phoEta_susy_mg);
  mgtree_susy->SetBranchAddress("lepPt",      &lepPt_susy_mg);
  mgtree_susy->SetBranchAddress("lepEta",     &lepEta_susy_mg);
  mgtree_susy->SetBranchAddress("sigMT",      &sigMT_susy_mg);
  mgtree_susy->SetBranchAddress("sigMET",     &sigMET_susy_mg);
  mgtree_susy->SetBranchAddress("HT",         &HT_susy_mg);
	mgtree_susy->SetBranchAddress("nVertex",    &nVertex_susy_mg);
	mgtree_susy->SetBranchAddress("MsGsQ",      &gluinoMass_susy_mg);
  mgtree_susy->SetBranchAddress("Mchargino",  &charginoMass_susy_mg);
  mgtree_susy->SetBranchAddress("Mneutralino",&neutralinoMass_susy_mg);
	mgtree_susy->SetBranchAddress("Mass1",      &Mass1_susy_mg);
	mgtree_susy->SetBranchAddress("Mass2",      &Mass2_susy_mg);

	for(unsigned ievt(0); ievt < mgtree_susy->GetEntries(); ievt++){
		mgtree_susy->GetEntry(ievt);
	
		/** cut flow *****/
		if(phoEt_susy_mg < 50 || lepPt_susy_mg < 50)continue;
		if(fabs(phoEta_susy_mg) > 1.4442 || fabs(lepEta_susy_mg) > 2.4)continue;
		if(sigMET_susy_mg < 120 || sigMT_susy_mg < 100)continue;

		if(charginoMass_susy_mg > 0 && neutralinoMass_susy_mg > 0){
			if(phoEt_susy_mg >= 200)phoEt_susy_mg = 199;
			if(lepPt_susy_mg >= 200)lepPt_susy_mg = 199;
			double mgESF = p_mgESF->GetBinContent( p_mgESF->GetXaxis()->FindBin( phoEt_susy_mg ), p_mgESF->GetYaxis()->FindBin( lepPt_susy_mg));
			double mgEFF = p_mgEff->GetBinContent( p_mgEff->GetXaxis()->FindBin( phoEt_susy_mg ), p_mgEff->GetYaxis()->FindBin( lepPt_susy_mg));
			total_effmg += mgEFF/mgESF; 
			total_mg += 1;
		} 
	}

  TChain *egtree_susy;
  egtree_susy = new TChain("egTree","egTree");
  egtree_susy->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_T5WG_string.root");
  float phoEt_susy_eg(0);
  float phoEta_susy_eg(0);
  float lepPt_susy_eg(0);
  float lepEta_susy_eg(0);
  float sigMT_susy_eg(0);
  float sigMET_susy_eg(0);
  float HT_susy_eg(0);
  int   nVertex_susy_eg(0); 
	float gluinoMass_susy_eg(0);
  float charginoMass_susy_eg(0);
  float neutralinoMass_susy_eg(0);
	float Mass1_susy_eg(0);
  float Mass2_susy_eg(0);
  egtree_susy->SetBranchAddress("phoEt",      &phoEt_susy_eg);
  egtree_susy->SetBranchAddress("phoEta",     &phoEta_susy_eg);
  egtree_susy->SetBranchAddress("lepPt",      &lepPt_susy_eg);
  egtree_susy->SetBranchAddress("lepEta",     &lepEta_susy_eg);
  egtree_susy->SetBranchAddress("sigMT",      &sigMT_susy_eg);
  egtree_susy->SetBranchAddress("sigMET",     &sigMET_susy_eg);
	egtree_susy->SetBranchAddress("nVertex",    &nVertex_susy_eg);
  egtree_susy->SetBranchAddress("HT",         &HT_susy_eg);
	egtree_susy->SetBranchAddress("MsGsQ",      &gluinoMass_susy_eg);
  egtree_susy->SetBranchAddress("Mchargino",  &charginoMass_susy_eg);
  egtree_susy->SetBranchAddress("Mneutralino",&neutralinoMass_susy_eg);
	egtree_susy->SetBranchAddress("Mass1",      &Mass1_susy_eg);
	egtree_susy->SetBranchAddress("Mass2",      &Mass2_susy_eg);

	for(unsigned ievt(0); ievt < egtree_susy->GetEntries(); ievt++){
		egtree_susy->GetEntry(ievt);

		/** cut flow *****/
		if(phoEt_susy_eg < 50 || lepPt_susy_eg < 50)continue;
		if(fabs(phoEta_susy_eg) > 1.4442 || fabs(lepEta_susy_eg) > 2.5)continue;
		if(sigMET_susy_eg < 120 || sigMT_susy_eg < 100)continue;
		if(charginoMass_susy_eg > 0 && neutralinoMass_susy_eg > 0){
			if(phoEt_susy_eg >= 200)phoEt_susy_eg = 199;
			double phoTRGESF = g_trg_hist->GetBinContent(g_trg_hist->GetXaxis()->FindBin(phoEt_susy_eg), g_trg_hist->GetYaxis()->FindBin(phoEta_susy_eg));
			if(lepPt_susy_eg >= 200)lepPt_susy_eg = 199;
			double eleTRGESF = e_trg_hist->GetBinContent(e_trg_hist->GetXaxis()->FindBin(lepPt_susy_eg), e_trg_hist->GetYaxis()->FindBin(fabs(lepEta_susy_eg)));
			total_effeg += phoTRGESF*eleTRGESF;
			total_eg += 1;
		}
	} 
	std::cout <<"eff (eg): " << total_effeg/total_eg << std::endl;
	std::cout <<"eff (mg): " << total_effmg/total_mg << std::endl;
}


