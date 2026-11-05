// input from analysis_Mixing.C  look at SUSYAnalysis/test/legacy/macro_old/plot_Mixing.C (older copy)
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<ctime>

#include "TFile.h"
#include "TTree.h"
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
#include "TProfile2D.h"

float DeltaR(float eta1,float phi1,float eta2,float phi2)
{
	float deltaPhi = TMath::Abs(phi1-phi2);
	float deltaEta = eta1-eta2;
	if(deltaPhi > TMath::Pi())
	deltaPhi = TMath::TwoPi() - deltaPhi;
		return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}
void plot_Mixing(){//main 
  gROOT->SetBatch(kTRUE);
  double scalefactor1 = 35.87*1000*489.0/6103732;
	double scalefactor2 = 35.87*1000*17.01/5077584.0;
	double scalefactor3 = 35.87*1000*0.87/2354481;
	//double scalefactor3 = 35.87*1000*0.87/1644983;

//  double scalefactor2 = 35.87*1000*117.8/14372399;
//	double scalefactor3 = 35.87*1000*0.143/489423;

  TH1D *mugamma_phoEt_1 = new TH1D("mugamma_phoEt_1","",40,0,400); 
  TH1D *mugamma_phoEt_2 = new TH1D("mugamma_phoEt_2","",40,0,400); 
  TH1D *mugamma_phoEt_3 = new TH1D("mugamma_phoEt_3","",40,0,400); 
	TH1D *egamma_phoEt_1 =  new TH1D("egamma_phoEt_1","", 40,0,400); 
	TH1D *egamma_phoEt_2 =  new TH1D("egamma_phoEt_2","", 40,0,400); 
	TH1D *egamma_phoEt_3 =  new TH1D("egamma_phoEt_3","", 40,0,400); 

  TH1D *mugamma_lepPt_1 = new TH1D("mugamma_lepPt_1","",40,0,400); 
  TH1D *mugamma_lepPt_2 = new TH1D("mugamma_lepPt_2","",40,0,400); 
  TH1D *mugamma_lepPt_3 = new TH1D("mugamma_lepPt_3","",40,0,400); 
  TH1D *mugamma_sigMET_1 = new TH1D("mugamma_sigMET_1","",40,0,400); 
  TH1D *mugamma_sigMET_2 = new TH1D("mugamma_sigMET_2","",40,0,400); 
  TH1D *mugamma_sigMET_3 = new TH1D("mugamma_sigMET_3","",40,0,400); 
  TH1D *mugamma_dR_1 = new TH1D("mugamma_dR_1","",32,0,3.2); 
  TH1D *mugamma_dR_2 = new TH1D("mugamma_dR_2","",32,0,3.2); 
  TH1D *mugamma_dR_3 = new TH1D("mugamma_dR_3","",32,0,3.2); 

  TChain *mgtree = new TChain("mgTree");
	mgtree->Add("/eos/uscms/store/user/tmishra/egMC/mixing_WGToLNuG_2016_TH1D.root");
  float mg_phoEt=0;
  float mg_phoEta=0;
  float mg_phoPhi=0;
  float mg_lepPt(0);
  float mg_sigMET(0);
  float mg_dRPhoLep(0);

  mgtree->SetBranchAddress("phoEt",    &mg_phoEt);
  mgtree->SetBranchAddress("phoEta",   &mg_phoEta);
  mgtree->SetBranchAddress("phoPhi",   &mg_phoPhi);
  mgtree->SetBranchAddress("lepPt",    &mg_lepPt);
  mgtree->SetBranchAddress("sigMET",   &mg_sigMET);
  mgtree->SetBranchAddress("dRPhoLep", &mg_dRPhoLep);

	for(unsigned ievt(0); ievt < mgtree->GetEntries(); ievt++){
		mgtree->GetEntry(ievt);
		//if(fabs(mg_phoEta) > 1.4442)continue;
		mugamma_phoEt_1->Fill(mg_phoEt, scalefactor1);		
		
		if(mg_phoEt > 140 && mg_phoEt < 300){
			mugamma_lepPt_1->Fill(mg_lepPt, scalefactor1);
			mugamma_sigMET_1->Fill(mg_sigMET, scalefactor1);
			mugamma_dR_1->Fill(mg_dRPhoLep, scalefactor1);
		}
	}


  TChain *mg40tree = new TChain("mgTree");
	mg40tree->Add("/eos/uscms/store/user/tmishra/egMC/mixing_WGJet40_2016_TH1D.root");
  float mg40_phoEt=0;
  float mg40_phoEta=0;
  float mg40_phoPhi=0;
  float mg40_lepPt(0);
  float mg40_sigMET(0);
  float mg40_dRPhoLep(0);

  mg40tree->SetBranchAddress("phoEt",    &mg40_phoEt);
  mg40tree->SetBranchAddress("phoEta",   &mg40_phoEta);
  mg40tree->SetBranchAddress("phoPhi",   &mg40_phoPhi);
  mg40tree->SetBranchAddress("lepPt",    &mg40_lepPt);
  mg40tree->SetBranchAddress("sigMET",   &mg40_sigMET);
  mg40tree->SetBranchAddress("dRPhoLep", &mg40_dRPhoLep);

	for(unsigned ievt(0); ievt < mg40tree->GetEntries(); ievt++){
		mg40tree->GetEntry(ievt);
		//if(fabs(mg40_phoEta) > 1.4442)continue;
		mugamma_phoEt_2->Fill(mg40_phoEt, scalefactor2);			
		if(mg40_phoEt > 50 && mg40_phoEt < 130){
			mugamma_lepPt_2->Fill(mg40_lepPt, scalefactor2);
			mugamma_sigMET_2->Fill(mg40_sigMET, scalefactor2);
			mugamma_dR_2->Fill(mg40_dRPhoLep, scalefactor2);
		}
	}

  TChain *mg130tree = new TChain("mgTree");
	mg130tree->Add("/eos/uscms/store/user/tmishra/egMC/mixing_WGJet130_2016_TH1D.root");
  float mg130_phoEt=0;
  float mg130_phoEta=0;
  float mg130_phoPhi=0;
  float mg130_lepPt(0);
  float mg130_sigMET(0);
  float mg130_dRPhoLep(0);

  mg130tree->SetBranchAddress("phoEt",    &mg130_phoEt);
  mg130tree->SetBranchAddress("phoEta",   &mg130_phoEta);
  mg130tree->SetBranchAddress("phoPhi",   &mg130_phoPhi);
  mg130tree->SetBranchAddress("lepPt",    &mg130_lepPt);
  mg130tree->SetBranchAddress("sigMET",   &mg130_sigMET);
  mg130tree->SetBranchAddress("dRPhoLep", &mg130_dRPhoLep);

	for(unsigned ievt(0); ievt < mg130tree->GetEntries(); ievt++){
		mg130tree->GetEntry(ievt);
		//if(fabs(mg130_phoEta) > 1.4442)continue;
		mugamma_phoEt_3->Fill(mg130_phoEt, scalefactor3);			
		if(mg130_phoEt > 140 && mg130_phoEt < 300){
			mugamma_lepPt_3->Fill(mg130_lepPt, scalefactor3);
			mugamma_sigMET_3->Fill(mg130_sigMET, scalefactor3);
			mugamma_dR_3->Fill(mg130_dRPhoLep, scalefactor3);
		}
	}

  TChain *egtree = new TChain("egTree");
	egtree->Add("/eos/uscms/store/user/tmishra/egMC/mixing_WGToLNuG_2016_TH1D.root");
  float eg_phoEt=0;
  float eg_phoEta=0;
  float eg_phoPhi=0;

  egtree->SetBranchAddress("phoEt",    &eg_phoEt);
  egtree->SetBranchAddress("phoEta",   &eg_phoEta);
  egtree->SetBranchAddress("phoPhi",   &eg_phoPhi);

	for(unsigned ievt(0); ievt < egtree->GetEntries(); ievt++){
		egtree->GetEntry(ievt);
		//if(fabs(eg_phoEta) > 1.4442)continue;
		egamma_phoEt_1->Fill(eg_phoEt, scalefactor1);			
	}


  TChain *eg40tree = new TChain("egTree");
	eg40tree->Add("/eos/uscms/store/user/tmishra/egMC/mixing_WGJet40_2016_TH1D.root");
  float eg40_phoEt=0;
  float eg40_phoEta=0;
  float eg40_phoPhi=0;

  eg40tree->SetBranchAddress("phoEt",    &eg40_phoEt);
  eg40tree->SetBranchAddress("phoEta",   &eg40_phoEta);
  eg40tree->SetBranchAddress("phoPhi",   &eg40_phoPhi);

	for(unsigned ievt(0); ievt < eg40tree->GetEntries(); ievt++){
		eg40tree->GetEntry(ievt);
		//if(fabs(eg40_phoEta) > 1.4442)continue;
		egamma_phoEt_2->Fill(eg40_phoEt, scalefactor2);			
	}

  TChain *eg130tree = new TChain("egTree");
	eg130tree->Add("/eos/uscms/store/user/tmishra/egMC/mixing_WGJet130_2016_TH1D.root");
  float eg130_phoEt=0;
  float eg130_phoEta=0;
  float eg130_phoPhi=0;

  eg130tree->SetBranchAddress("phoEt",    &eg130_phoEt);
  eg130tree->SetBranchAddress("phoEta",   &eg130_phoEta);
  eg130tree->SetBranchAddress("phoPhi",   &eg130_phoPhi);

	for(unsigned ievt(0); ievt < eg130tree->GetEntries(); ievt++){
		eg130tree->GetEntry(ievt);
		//if(fabs(eg130_phoEta) > 1.4442)continue;
		egamma_phoEt_3->Fill(eg130_phoEt, scalefactor3);			
	}

	egamma_phoEt_3->Add(egamma_phoEt_2);
	mugamma_phoEt_3->Add(mugamma_phoEt_2);
  double ratio1 = egamma_phoEt_1->Integral(6,8)/egamma_phoEt_2->Integral(6,8);
  double ratio2 = mugamma_phoEt_1->Integral(6,8)/mugamma_phoEt_2->Integral(6,8);

  double ratio1_70 = egamma_phoEt_1->Integral(7,10)/egamma_phoEt_2->Integral(7,10);
  double ratio2_70 = mugamma_phoEt_1->Integral(7,10)/mugamma_phoEt_2->Integral(7,10);

  double ratio1_130 = egamma_phoEt_1->Integral(15,30)/egamma_phoEt_3->Integral(15,30);
  double ratio2_130 = mugamma_phoEt_1->Integral(15,30)/mugamma_phoEt_3->Integral(15,30);

  double ratio3_130 = egamma_phoEt_1->Integral(16,40)/egamma_phoEt_3->Integral(16,40);
  double ratio4_130 = mugamma_phoEt_1->Integral(16,40)/mugamma_phoEt_3->Integral(16,40);

	std::cout << "ratio3_130=" <<ratio3_130 <<std::endl;
	std::cout << "ratio4_130=" <<ratio4_130 <<std::endl;

  TLatex* latex = new TLatex();
	latex->SetTextAlign(13);
	latex->SetTextSize(0.035);
	std::ostringstream textname;

	
  mugamma_phoEt_1->Rebin(5);
  mugamma_phoEt_2->Rebin(5);
  mugamma_phoEt_3->Rebin(5);



	TCanvas *can1 = new TCanvas("can1","",600,600);
  can1->cd();
	gPad->SetLogy();
  egamma_phoEt_1->SetTitle("e+gamma channel photon");
	egamma_phoEt_1->Draw();
	egamma_phoEt_2->SetLineColor(kRed);
	egamma_phoEt_2->Scale(ratio1);
	egamma_phoEt_2->Draw("same");
	egamma_phoEt_3->SetLineColor(kGreen);
	egamma_phoEt_3->Scale(ratio1_130);
	egamma_phoEt_3->Draw("same");
  textname.str("");
  textname << "ratio(50~70) =" << ratio1;
  latex->DrawLatex(200,0.01*egamma_phoEt_1->GetMaximum(), textname.str().c_str() );
  textname.str("");
  textname << "ratio(70~100) =" << ratio1_70;
  latex->DrawLatex(200,0.005*egamma_phoEt_1->GetMaximum(), textname.str().c_str() );
  textname.str("");
  textname << "WG130 =" << ratio1_130; 
  latex->DrawLatex(200,0.001*egamma_phoEt_1->GetMaximum(), textname.str().c_str() );
	can1->SaveAs("/eos/uscms/store/user/tmishra/VGamma/ZGMixingScale_eg.pdf");

	TCanvas *can2 = new TCanvas("can2","",600,600);
  can2->cd();
	gPad->SetLogy();
  mugamma_phoEt_1->SetTitle("mu+gamma channel photon");
	mugamma_phoEt_1->Draw();
	mugamma_phoEt_2->SetLineColor(kRed);
	mugamma_phoEt_2->Scale(ratio2);
	mugamma_phoEt_2->Draw("same");
	mugamma_phoEt_3->SetLineColor(kGreen);
	mugamma_phoEt_3->Scale(ratio2_130);
	mugamma_phoEt_3->Draw("same");
  textname.str("");
  textname << "ratio(50~70) =" << ratio2;
  latex->DrawLatex(200,0.01*mugamma_phoEt_1->GetMaximum(), textname.str().c_str() );
  textname.str("");
  textname << "ratio(70~100) =" << ratio2_70;
  latex->DrawLatex(200,0.005*mugamma_phoEt_1->GetMaximum(), textname.str().c_str() );
  textname.str("");
  textname << "WG130 =" << ratio2_130; 
  latex->DrawLatex(200,0.001*mugamma_phoEt_1->GetMaximum(), textname.str().c_str() );
	can2->SaveAs("/eos/uscms/store/user/tmishra/VGamma/ZGMixingScale_mg.pdf");


	TCanvas *can3 = new TCanvas("can3","",600,600);
  can3->cd();
	gPad->SetLogy();
  mugamma_lepPt_1->SetTitle("e+gamma channel photon");
	mugamma_lepPt_1->Draw();
	mugamma_lepPt_3->SetLineColor(kRed);
	mugamma_lepPt_3->Scale(ratio2_130);
	mugamma_lepPt_3->Draw("same");
	can3->SaveAs("/eos/uscms/store/user/tmishra/VGamma/lepPt_mg.pdf");

	TCanvas *can4 = new TCanvas("can4","",600,600);
  can4->cd();
	gPad->SetLogy();
  mugamma_sigMET_1->SetTitle("e+gamma channel photon");
	mugamma_sigMET_1->Draw();
	mugamma_sigMET_3->SetLineColor(kRed);
	mugamma_sigMET_3->Scale(ratio2_130);
	mugamma_sigMET_3->Draw("same");
	can4->SaveAs("/eos/uscms/store/user/tmishra/VGamma/MET_mg.pdf");

	TCanvas *can5 = new TCanvas("can5","",600,600);
  can5->cd();
	gPad->SetLogy();
  mugamma_dR_1->SetTitle("e+gamma channel photon");
	mugamma_dR_1->Draw();
	mugamma_dR_3->SetLineColor(kRed);
	mugamma_dR_3->Scale(ratio2_130);
	mugamma_dR_3->Draw("same");
	can5->SaveAs("/eos/uscms/store/user/tmishra/VGamma/MET_mg.pdf");
}


