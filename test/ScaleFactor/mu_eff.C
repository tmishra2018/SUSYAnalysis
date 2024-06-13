#include "../../include/analysis_commoncode.h"
#include "TH2D.h"
#include "TCanvas.h"
void mu_eff(){
	
	setTDRStyle();
	gROOT->SetBatch(kTRUE);
  	gSystem->Load("../../lib/libAnaClasses.so");
	esfScaleFactor  objectESF;
	
  	int channelType = 2; 

	//*********** histo list **********************//
	TH2D *p_eff = new TH2D("p_eff","Efficiency;Muon p_{T} (GeV);#eta",5,0,5,4,0,4);
	double lep_eff[5][4];
	double mc_count[5][4];
	for(unsigned i(0); i< 5;i++)
		for(unsigned j(0); j<4;j++){
			lep_eff[i][j] = 0;
			mc_count[i][j] = 0;
		}

// ********  MC *************************//
	std::ostringstream chainname;
	chainname.str("");
	if(channelType == 1)chainname << "egTree";
	else if(channelType == 2)chainname << "mgTree";
  TChain *mctree = new TChain(chainname.str().c_str(), chainname.str().c_str());
  mctree->Add("/eos/uscms/store/user/tmishra/Trigger/Efficiency_TT.root");
  float lepPt(0);
  float lepEta(0);
  float lepPhi(0);
  int   mcPID=0;
  float mcEta=0;
  float mcPhi=0;
  float mcPt=0;
  int   mcMomPID=0;
//
  mctree->SetBranchAddress("lepPt",     &lepPt);
  mctree->SetBranchAddress("lepEta",    &lepEta);
  mctree->SetBranchAddress("lepPhi",    &lepPhi);
  mctree->SetBranchAddress("mcPID",     &mcPID);
  mctree->SetBranchAddress("mcEta",     &mcEta);
  mctree->SetBranchAddress("mcPhi",     &mcPhi);
  mctree->SetBranchAddress("mcPt",      &mcPt);
  mctree->SetBranchAddress("mcMomPID",  &mcMomPID);

	int nMatch = 0;
	int nEle = 0;
	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);
		double scalefactor(0);
		double s_error(0);

		int ptbin(-1);
		int etabin(-1);
		if(mcPt > 25 && mcPt < 30)
			ptbin = 0;
		else if(mcPt >= 30 && mcPt < 40)
			ptbin = 1;
		else if(mcPt >= 40 && mcPt <50)
			ptbin = 2;
		else if(mcPt  >= 50 && mcPt < 100)
			ptbin = 3;
		else if(mcPt > 100)
			ptbin = 4;
		else continue;
	
		if(fabs(mcEta) < 0.9)
			etabin = 0;
		else if(fabs(mcEta) > 0.9 && fabs(mcEta) < 1.2)
			etabin = 1;
		else if(fabs(mcEta) > 1.2 && fabs(mcEta) < 2.1)
			etabin = 2;
		else if(fabs(mcEta) > 2.1 && fabs(mcEta) < 2.4)
			etabin = 3;

		if(fabs(mcEta) > 2.4)continue;

		if(channelType == 1){
			scalefactor = objectESF.getElectronESF(mcPt,mcEta);
			s_error = objectESF.getElectronESFError(mcPt,mcEta);
		}
		if(channelType == 2){
			scalefactor = objectESF.getMuonESF(mcPt,mcEta);
			s_error = objectESF.getMuonESFError(mcPt,mcEta);
		}	
		
		double dR = DeltaR(mcEta, mcPhi, lepEta, lepPhi);
		double dE = fabs(mcPt - lepPt)/mcPt;

		mc_count[ptbin][etabin] += 1;
		if(dR < 0.5 && dE < 0.5){
			lep_eff[ptbin][etabin] += 1;
		}
		
	}

	double ptlist[] = {25,35,45,70,130};
	double etalist[] = {0.45, 1.1, 1.8, 2.2};
	for(unsigned i(0); i<5; i++){
		for(unsigned j(0); j<4; j++){
			double pt = ptlist[i];
			double eta = etalist[j];
			double scalefactor(0);
			double s_error(0);
			if(channelType == 1){
				scalefactor = objectESF.getElectronESF(pt,eta);
				s_error = objectESF.getElectronESFError(pt,eta);
			}
			if(channelType == 2){
				scalefactor = objectESF.getMuonESF(pt,eta);
				s_error = objectESF.getMuonESFError(pt,eta);
			}	

			double num = lep_eff[i][j];
			double den = mc_count[i][j];
			double stat_e = sqrt((num*num*den + den*den*num))/den/den;
			lep_eff[i][j] = double(lep_eff[i][j])/mc_count[i][j];
			lep_eff[i][j] = lep_eff[i][j]*scalefactor;	
			p_eff->SetBinContent(i+1, j+1, lep_eff[i][j]);
			p_eff->SetBinError(i+1, j+1, sqrt(lep_eff[i][j]*lep_eff[i][j]*s_error*s_error + stat_e*stat_e*scalefactor*scalefactor) );
			std::cout << i << " " << j << " " << lep_eff[i][j] << " " << num << " " << den << " " << stat_e << " " << s_error<< std::endl;
		}
	}
	TCanvas *can=new TCanvas("can","",1200,800);
  gStyle->SetPadLeftMargin(0.10);
  can->SetLeftMargin(0.10);
  can->SetRightMargin(0.1);
	can->SetTopMargin(0.08);
	can->SetBottomMargin(0.15);
	can->cd();
	p_eff->SetTitleOffset(0.01);
	p_eff->GetXaxis()->SetBinLabel(1,"25-30");
	p_eff->GetXaxis()->SetBinLabel(2,"30-40");
	p_eff->GetXaxis()->SetBinLabel(3,"40-50");
	p_eff->GetXaxis()->SetBinLabel(4,"50-100");
	p_eff->GetXaxis()->SetBinLabel(5,"> 100");
	p_eff->GetXaxis()->SetTitleOffset(1.2);

	p_eff->GetYaxis()->SetBinLabel(1,"0-0.9");
	p_eff->GetYaxis()->SetBinLabel(2,"0.9-1.2");
	p_eff->GetYaxis()->SetBinLabel(3,"1.2-2.1");
	p_eff->GetYaxis()->SetBinLabel(4,"2.1-2.4");
	
	gStyle->SetPaintTextFormat("4.3f");
	p_eff->GetZaxis()->SetRangeUser(0,1);
	p_eff->Draw("colz text E");
	can->SaveAs("/eos/uscms/store/user/tmishra/Trigger/eff_mu_17012.pdf");

	TFile *outputfile = new TFile("/eos/uscms/store/user/tmishra/Trigger/eff_mu_17012.root","RECREATE");
	outputfile->cd();
	p_eff->Write();
	outputfile->Close();	
}
