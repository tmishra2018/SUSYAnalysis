void CheckVar(){
	gROOT->SetBatch(1);
	TFile file1("/eos/uscms/store/user/tmishra/jetfakepho/files/plot_hadron_GJet_2016preVFP.root");
	TTree *mctree = (TTree*)file1.Get("egTree");
        float mc_phoEt(0);
        float mc_phoEta(0);
        float mc_phoPhi(0);
        float mc_phoSigma(0);
        float mc_phoChIso(0);
        mctree->SetBranchAddress("phoEt",     &mc_phoEt);
        mctree->SetBranchAddress("phoEta",    &mc_phoEta);
        mctree->SetBranchAddress("phoPhi",    &mc_phoPhi);
        mctree->SetBranchAddress("phoSigma",  &mc_phoSigma);
        mctree->SetBranchAddress("phoChIso",  &mc_phoChIso);

	TCanvas *can = new TCanvas("can","can",1200,600);
	TH1D *hPhoEt;
	mctree->Draw("phoChIso>>hPhoEt");

//      for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
//             mctree->GetEntry(ievt);
	TFile file2("/uscms_data/d3/mengleis/FullStatusOct/plot_hadron_GJet.root");
	TTree *mctree_Mengl = (TTree*)file2.Get("egTree");

        float mc_phoEt_Mengl(0);
        float mc_phoEta_Mengl(0);
        float mc_phoPhi_Mengl(0);
        float mc_phoSigma_Mengl(0);
        float mc_phoChIso_Mengl(0);
        mctree_Mengl->SetBranchAddress("phoEt",     &mc_phoEt_Mengl);
        mctree_Mengl->SetBranchAddress("phoEta",    &mc_phoEta_Mengl);
        mctree_Mengl->SetBranchAddress("phoPhi",    &mc_phoPhi_Mengl);
        mctree_Mengl->SetBranchAddress("phoSigma",  &mc_phoSigma_Mengl);
        mctree_Mengl->SetBranchAddress("phoChIso",  &mc_phoChIso_Mengl);
	TCanvas *can2 = new TCanvas("can2","can2",1200,600);
	TH1D *hPhoEt_Mengl;
	mctree_Mengl->Draw("phoChIso>>hPhoEt_Mengl");

	can->SaveAs("Check.png");
	can2->SaveAs("Check2.png");
	
}
