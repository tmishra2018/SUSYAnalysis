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
#include "TLatex.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "../../include/tdrstyle.C"

int channel = 2; // 1 = eg, 2 = mg

void plot_bkg(){//main  

  gSystem->Load("../../lib/libAnaClasses.so");
	setTDRStyle();

	std::ostringstream signame;  signame.str("");
	std::ostringstream elename;  elename.str("");
	std::ostringstream jetname;  jetname.str("");
	std::ostringstream qcdname;  qcdname.str("");
	std::ostringstream VGname;   VGname.str("");
	std::ostringstream rarename; rarename.str(""); 

	if(channel == 1){
			signame <<  "../Result/vetoDiPhoData/signalTree_egamma_signal.root";
			elename <<  "postTree_egamma_eleBkg.root";  
			jetname <<  "postTree_egamma_jetbkg.root";
			qcdname <<  "postTree_egamma_qcd.root";
			VGname  <<  "postTree_egamma_VGBkg.root";
			rarename << "postTree_egamma_rareBkg.root";
	}
	else if(channel == 2){
			signame << "../Result/vetoDiPhoData/signalTree_mg_signal.root";
			elename << "postTree_mg_eleBkg.root";  
			jetname << "postTree_mg_jetbkg.root";
			qcdname << "postTree_mg_qcd.root";
			VGname  << "postTree_mg_VGBkg.root";
			rarename <<"postTree_mg_rareBkg.root";
	}

	TFile *file_sig = TFile::Open(signame.str().c_str());
	TFile *file_ele = TFile::Open(elename.str().c_str());
	TFile *file_jet = TFile::Open(jetname.str().c_str());
	TFile *file_qcd = TFile::Open(qcdname.str().c_str());
	TFile *file_VG  = TFile::Open(VGname.str().c_str());
	TFile *file_rare= TFile::Open(rarename.str().c_str()); 

	Double_t bkgEtBins[]={35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400,500,800};
	int nBkgEtBins= sizeof(bkgEtBins)/sizeof(bkgEtBins[0]) -1;
	Double_t bkgPtBins[]={25,30,35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400,500,800};
	int nBkgPtBins= sizeof(bkgPtBins)/sizeof(bkgPtBins[0])-1;
	Double_t bkgMETBins[]={0,40,60,80,100,120,140,160,180,210,240,280,320,400,600,1000};
	int nBkgMETBins= sizeof(bkgMETBins)/sizeof(bkgMETBins[0]) -1;
	Double_t bkgHTBins[]={0,40,60,80,100,120,140,160,180,200,225,250,275,300,340,380,420,500,600,1000};
	int nBkgHTBins= sizeof(bkgHTBins)/sizeof(bkgHTBins[0]) -1;

	TH1F *p_allPhoEt = (TH1F*)file_sig->Get("p_PhoEt");
	TH1F *p_allMET   = (TH1F*)file_sig->Get("p_MET");
	TH1F *p_allHT  = (TH1F*)file_sig->Get("p_HT");

	p_allPhoEt->SetMarkerColor(1);
	p_allMET->SetMarkerColor(1);
	p_allHT->SetMarkerColor(1);
	p_allPhoEt->SetTitle("");
	p_allMET->SetTitle("");
	p_allHT->SetTitle("");
	p_allPhoEt->GetYaxis()->SetTitle("Events / bin");
	p_allMET->GetYaxis()->SetTitle("Events / bin");
	p_allHT->GetYaxis()->SetTitle("Events / bin");
	
	TH1F *postfit_elePhoEt = (TH1F*)file_ele->Get("postfit_PhoEt");
	TH1F *postfit_eleMET   = (TH1F*)file_ele->Get("postfit_MET");
	TH1F *postfit_eleHT  = (TH1F*)file_ele->Get("postfit_HT");

	TH1F *postfit_jetPhoEt = (TH1F*)file_jet->Get("postfit_PhoEt");
	TH1F *postfit_jetMET   = (TH1F*)file_jet->Get("postfit_MET");
	TH1F *postfit_jetHT  = (TH1F*)file_jet->Get("postfit_HT");

	TH1F *postfit_qcdPhoEt = (TH1F*)file_qcd->Get("postfit_PhoEt");
	TH1F *postfit_qcdMET   = (TH1F*)file_qcd->Get("postfit_MET");
	TH1F *postfit_qcdHT  = (TH1F*)file_qcd->Get("postfit_HT");

	TH1F *postfit_VGPhoEt = (TH1F*)file_VG->Get("postfit_PhoEt");
	TH1F *postfit_VGMET   = (TH1F*)file_VG->Get("postfit_MET");
	TH1F *postfit_VGHT  = (TH1F*)file_VG->Get("postfit_HT");

	TH1F *postfit_rarePhoEt = (TH1F*)file_rare->Get("postfit_PhoEt");
	TH1F *postfit_rareMET   = (TH1F*)file_rare->Get("postfit_MET");
	TH1F *postfit_rareHT  = (TH1F*)file_rare->Get("postfit_HT");

  int binnumber;
  binnumber = p_allPhoEt->GetSize()-2;  p_allPhoEt->SetBinContent(binnumber,  p_allPhoEt->GetBinContent(binnumber) +   p_allPhoEt->GetBinContent(binnumber+1) );   
  binnumber = p_allMET->GetSize()-2;    p_allMET->SetBinContent(binnumber,    p_allMET->GetBinContent(binnumber) +     p_allMET->GetBinContent(binnumber+1) );   
  binnumber = p_allHT->GetSize()-2;     p_allHT->SetBinContent(binnumber,     p_allHT->GetBinContent(binnumber) +      p_allHT->GetBinContent(binnumber+1) );    
                                                                                                                           
  binnumber = postfit_elePhoEt->GetSize()-2;  postfit_elePhoEt->SetBinContent(binnumber,  postfit_elePhoEt->GetBinContent(binnumber) +   postfit_elePhoEt->GetBinContent(binnumber+1) ); 
  binnumber = postfit_eleMET->GetSize()-2;    postfit_eleMET->SetBinContent(binnumber,    postfit_eleMET->GetBinContent(binnumber) +     postfit_eleMET->GetBinContent(binnumber+1) );   
  binnumber = postfit_eleHT->GetSize()-2;     postfit_eleHT->SetBinContent(binnumber,     postfit_eleHT->GetBinContent(binnumber) +      postfit_eleHT->GetBinContent(binnumber+1) );    
                                                                                                                                                                            
  binnumber = postfit_jetPhoEt->GetSize()-2;  postfit_jetPhoEt->SetBinContent(binnumber,  postfit_jetPhoEt->GetBinContent(binnumber) +   postfit_jetPhoEt->GetBinContent(binnumber+1) );             
  binnumber = postfit_jetMET->GetSize()-2;    postfit_jetMET->SetBinContent(binnumber,    postfit_jetMET->GetBinContent(binnumber) +     postfit_jetMET->GetBinContent(binnumber+1) );               
  binnumber = postfit_jetHT->GetSize()-2;     postfit_jetHT->SetBinContent(binnumber,     postfit_jetHT->GetBinContent(binnumber) +      postfit_jetHT->GetBinContent(binnumber+1) );  
                                                                                                                                                                            
  binnumber = postfit_qcdPhoEt->GetSize()-2;  postfit_qcdPhoEt->SetBinContent(binnumber,  postfit_qcdPhoEt->GetBinContent(binnumber) +   postfit_qcdPhoEt->GetBinContent(binnumber+1) );  
  binnumber = postfit_qcdMET->GetSize()-2;    postfit_qcdMET->SetBinContent(binnumber,    postfit_qcdMET->GetBinContent(binnumber) +     postfit_qcdMET->GetBinContent(binnumber+1) );    
  binnumber = postfit_qcdHT->GetSize()-2;     postfit_qcdHT->SetBinContent(binnumber,     postfit_qcdHT->GetBinContent(binnumber) +      postfit_qcdHT->GetBinContent(binnumber+1) );  
                                                                                                                                                                            
  binnumber = postfit_VGPhoEt->GetSize()-2;   postfit_VGPhoEt->SetBinContent(binnumber,   postfit_VGPhoEt->GetBinContent(binnumber) +    postfit_VGPhoEt->GetBinContent(binnumber+1) );   
  binnumber = postfit_VGMET->GetSize()-2;     postfit_VGMET->SetBinContent(binnumber,     postfit_VGMET->GetBinContent(binnumber) +      postfit_VGMET->GetBinContent(binnumber+1) );     
  binnumber = postfit_VGHT->GetSize()-2;      postfit_VGHT->SetBinContent(binnumber,      postfit_VGHT->GetBinContent(binnumber) +       postfit_VGHT->GetBinContent(binnumber+1) );   
                                                                                                                                                                            
  binnumber = postfit_rarePhoEt->GetSize()-2; postfit_rarePhoEt->SetBinContent(binnumber, postfit_rarePhoEt->GetBinContent(binnumber) +  postfit_rarePhoEt->GetBinContent(binnumber+1) );            
  binnumber = postfit_rareMET->GetSize()-2;   postfit_rareMET->SetBinContent(binnumber,   postfit_rareMET->GetBinContent(binnumber) +    postfit_rareMET->GetBinContent(binnumber+1) );              
  binnumber = postfit_rareHT->GetSize()-2;    postfit_rareHT->SetBinContent(binnumber,    postfit_rareHT->GetBinContent(binnumber) +     postfit_rareHT->GetBinContent(binnumber+1) );                  
 // Errors
  binnumber = p_allPhoEt->GetSize()-2;  p_allPhoEt->SetBinError(binnumber,  p_allPhoEt->GetBinError(binnumber) +   p_allPhoEt->GetBinError(binnumber+1) );   
  binnumber = p_allMET->GetSize()-2;    p_allMET->SetBinError(binnumber,    p_allMET->GetBinError(binnumber) +     p_allMET->GetBinError(binnumber+1) );   
  binnumber = p_allHT->GetSize()-2;     p_allHT->SetBinError(binnumber,     p_allHT->GetBinError(binnumber) +      p_allHT->GetBinError(binnumber+1) );    
                                                                                                                           
  binnumber = postfit_elePhoEt->GetSize()-2;  postfit_elePhoEt->SetBinError(binnumber,  postfit_elePhoEt->GetBinError(binnumber) +   postfit_elePhoEt->GetBinError(binnumber+1) ); 
  binnumber = postfit_eleMET->GetSize()-2;    postfit_eleMET->SetBinError(binnumber,    postfit_eleMET->GetBinError(binnumber) +     postfit_eleMET->GetBinError(binnumber+1) );   
  binnumber = postfit_eleHT->GetSize()-2;     postfit_eleHT->SetBinError(binnumber,     postfit_eleHT->GetBinError(binnumber) +      postfit_eleHT->GetBinError(binnumber+1) );    
                                                                                                                                                                            
  binnumber = postfit_jetPhoEt->GetSize()-2;  postfit_jetPhoEt->SetBinError(binnumber,  postfit_jetPhoEt->GetBinError(binnumber) +   postfit_jetPhoEt->GetBinError(binnumber+1) );             
  binnumber = postfit_jetMET->GetSize()-2;    postfit_jetMET->SetBinError(binnumber,    postfit_jetMET->GetBinError(binnumber) +     postfit_jetMET->GetBinError(binnumber+1) );               
  binnumber = postfit_jetHT->GetSize()-2;     postfit_jetHT->SetBinError(binnumber,     postfit_jetHT->GetBinError(binnumber) +      postfit_jetHT->GetBinError(binnumber+1) );  
                                                                                                                                                                            
  binnumber = postfit_qcdPhoEt->GetSize()-2;  postfit_qcdPhoEt->SetBinError(binnumber,  postfit_qcdPhoEt->GetBinError(binnumber) +   postfit_qcdPhoEt->GetBinError(binnumber+1) );  
  binnumber = postfit_qcdMET->GetSize()-2;    postfit_qcdMET->SetBinError(binnumber,    postfit_qcdMET->GetBinError(binnumber) +     postfit_qcdMET->GetBinError(binnumber+1) );    
  binnumber = postfit_qcdHT->GetSize()-2;     postfit_qcdHT->SetBinError(binnumber,     postfit_qcdHT->GetBinError(binnumber) +      postfit_qcdHT->GetBinError(binnumber+1) );  
                                                                                                                                                                            
  binnumber = postfit_VGPhoEt->GetSize()-2;   postfit_VGPhoEt->SetBinError(binnumber,   postfit_VGPhoEt->GetBinError(binnumber) +    postfit_VGPhoEt->GetBinError(binnumber+1) );   
  binnumber = postfit_VGMET->GetSize()-2;     postfit_VGMET->SetBinError(binnumber,     postfit_VGMET->GetBinError(binnumber) +      postfit_VGMET->GetBinError(binnumber+1) );     
  binnumber = postfit_VGHT->GetSize()-2;      postfit_VGHT->SetBinError(binnumber,      postfit_VGHT->GetBinError(binnumber) +       postfit_VGHT->GetBinError(binnumber+1) );   
                                                                                                                                                                            
  binnumber = postfit_rarePhoEt->GetSize()-2; postfit_rarePhoEt->SetBinError(binnumber, postfit_rarePhoEt->GetBinError(binnumber) +  postfit_rarePhoEt->GetBinError(binnumber+1) );            
  binnumber = postfit_rareMET->GetSize()-2;   postfit_rareMET->SetBinError(binnumber,   postfit_rareMET->GetBinError(binnumber) +    postfit_rareMET->GetBinError(binnumber+1) );              
  binnumber = postfit_rareHT->GetSize()-2;    postfit_rareHT->SetBinError(binnumber,    postfit_rareHT->GetBinError(binnumber) +     postfit_rareHT->GetBinError(binnumber+1) );                  


	std::ostringstream etplot;  etplot.str("");
	std::ostringstream ptplot;  ptplot.str("");
	std::ostringstream metplot; metplot.str("");
	std::ostringstream mtplot;  mtplot.str("");
	std::ostringstream htplot;  htplot.str("");
	
	if(channel == 1){
			etplot << "SIGNAL_egamma_2016ReMiniAOD_pt.pdf";     
			ptplot << "SIGNAL_egamma_2016ReMiniAOD_leppt.pdf";	 
			metplot << "SIGNAL_egamma_2016ReMiniAOD_met.pdf";   
			mtplot << "SIGNAL_egamma_2016ReMiniAOD_mt.pdf";     
			htplot << "SIGNAL_egamma_2016ReMiniAOD_ht.pdf";	   
	}                                                      
	else if(channel == 2){
			etplot << "SIGNAL_mg_2016ReMiniAOD_pt.pdf";
			ptplot << "SIGNAL_mg_2016ReMiniAOD_leppt.pdf";		
			metplot << "SIGNAL_mg_2016ReMiniAOD_met.pdf";
			mtplot << "SIGNAL_mg_2016ReMiniAOD_mt.pdf";
			htplot << "SIGNAL_mg_2016ReMiniAOD_ht.pdf";	
	}

  gStyle->SetOptStat(0);
	TCanvas *c_pt = new TCanvas("Photon_Pt", "Photon P_{T}",600,600);
	setCanvas(c_pt); 
	c_pt->cd();
	TPad *pt_pad1 = new TPad("pt_pad1", "pt_pad1", 0, 0.3, 1, 1.0);
	setTopPad(pt_pad1); 
	pt_pad1->Draw();  
	pt_pad1->cd();  
	gPad->SetLogy();
	p_allPhoEt->SetMaximum(50*p_allPhoEt->GetBinContent(p_allPhoEt->GetMaximumBin()));
	p_allPhoEt->SetMinimum(0.5);
	p_allPhoEt->GetXaxis()->SetRangeUser(35,800);
	p_allPhoEt->SetLineColor(1);
	p_allPhoEt->SetMarkerStyle(20);
	p_allPhoEt->Draw("P");
	postfit_VGPhoEt->SetFillStyle(1001);
	postfit_VGPhoEt->SetLineColor(kGreen-7);
	postfit_VGPhoEt->SetFillColor(kGreen-7);
	postfit_rarePhoEt->SetFillStyle(1001);
	postfit_rarePhoEt->SetLineColor(kAzure-9);
	postfit_rarePhoEt->SetFillColor(kAzure-9);
	postfit_qcdPhoEt->SetFillStyle(1001);
	postfit_qcdPhoEt->SetLineColor(kYellow-9);
	postfit_qcdPhoEt->SetFillColor(kYellow-9);
	postfit_elePhoEt->SetFillStyle(1001);
	postfit_elePhoEt->SetLineColor(kMagenta-9);
	postfit_elePhoEt->SetFillColor(kMagenta-9);
	postfit_jetPhoEt->SetFillStyle(1001);
	postfit_jetPhoEt->SetLineColor(kOrange-9);
	postfit_jetPhoEt->SetFillColor(kOrange-9);
	postfit_elePhoEt->Add(postfit_rarePhoEt); // ele 2nd
	postfit_jetPhoEt->Add(postfit_elePhoEt);  // jet 3rd
	postfit_qcdPhoEt->Add(postfit_jetPhoEt);  // qcd 4th
	postfit_VGPhoEt->Add(postfit_qcdPhoEt);   // VG  5th
	postfit_VGPhoEt->Sumw2();
	postfit_VGPhoEt->Draw("hist same");
	postfit_qcdPhoEt->Draw("hist same");
	postfit_jetPhoEt->Draw("hist same");
	postfit_elePhoEt->Draw("hist same");
	postfit_rarePhoEt->Draw("hist same");
	p_allPhoEt->Draw("E same");
  TLatex chantex;
  chantex.SetNDC();
  chantex.SetTextFont(42);
  chantex.SetTextSize(0.05);    
  chantex.DrawLatex(0.4,0.87,"M_{T} > 100 GeV, postfit_{T}^{miss} > 120 GeV");
  chantex.SetTextSize(0.07);    
	chantex.DrawLatex(0.8,0.5,"(d)");
 	gPad->RedrawAxis();
  CMS_lumi( pt_pad1, 11 );

	c_pt->cd();
	TPad *pt_pad2 = new TPad("pt_pad2", "pt_pad2", 0, 0, 1, 0.3);
	pt_pad2->SetBottomMargin(0.4);
	pt_pad2->Draw();
	pt_pad2->cd();
  TLine *flatratio = new TLine(35,1,800,1);
	TH1F *ratio=(TH1F*)p_allPhoEt->Clone("transfer factor");
	ratio->SetMinimum(0);
	ratio->SetMaximum(2);
	ratio->SetMarkerStyle(20);
	ratio->SetLineColor(kBlack);
	ratio->Divide(postfit_VGPhoEt);
	ratio->SetTitle("");
	ratio->GetYaxis()->SetTitle("#frac{Obs.}{Bkg.}");
	ratio->GetYaxis()->SetNdivisions(504);
	ratio->Draw();
	ratio->Draw("same");
	flatratio->Draw("same");
	c_pt->SaveAs(etplot.str().c_str());

// ******** MET ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_met = new TCanvas("MET", "MET",600,600);
	setCanvas(c_met); 
	c_met->cd();
	TPad *met_pad1 = new TPad("met_pad1", "met_pad1", 0, 0.3, 1, 1.0);
	setTopPad(met_pad1); 
	met_pad1->Draw();  
	met_pad1->cd();  
	gPad->SetLogy();
	p_allMET->GetYaxis()->SetRangeUser(0.05, 50*p_allMET->GetBinContent(p_allMET->GetMaximumBin()));
	p_allMET->SetMinimum(0.5);
	p_allMET->GetXaxis()->SetRangeUser(0,600);
	p_allMET->SetLineColor(1);
	p_allMET->SetMarkerStyle(20);
	p_allMET->Draw("P");
	postfit_VGMET->SetFillStyle(1001);
	postfit_VGMET->SetLineColor(kGreen-7);
	postfit_VGMET->SetFillColor(kGreen-7);
	postfit_rareMET->SetFillStyle(1001);
	postfit_rareMET->SetLineColor(kAzure-9);
	postfit_rareMET->SetFillColor(kAzure-9);
	postfit_qcdMET->SetFillStyle(1001);
	postfit_qcdMET->SetLineColor(kYellow-9);
	postfit_qcdMET->SetFillColor(kYellow-9);
	postfit_eleMET->SetFillStyle(1001);
	postfit_eleMET->SetLineColor(kMagenta-9);
	postfit_eleMET->SetFillColor(kMagenta-9);
	postfit_jetMET->SetFillStyle(1001);
	postfit_jetMET->SetLineColor(kOrange-9);
	postfit_jetMET->SetFillColor(kOrange-9);
	postfit_eleMET->Add(postfit_rareMET); // ele 2nd
	postfit_jetMET->Add(postfit_eleMET);  // jet 3rd
	postfit_qcdMET->Add(postfit_jetMET);  // qcd 4th
	postfit_VGMET->Add(postfit_qcdMET);   // VG  5th
	postfit_VGMET->Sumw2();
	postfit_VGMET->Draw("hist same");
	postfit_qcdMET->Draw("hist same");
	postfit_jetMET->Draw("hist same");
	postfit_eleMET->Draw("hist same");
	postfit_rareMET->Draw("hist same");
	p_allMET->Draw("E same");
  chantex.SetTextSize(0.05);    
  chantex.DrawLatex(0.4,0.87,"M_{T} > 100 GeV");
  chantex.SetTextSize(0.07);    
	chantex.DrawLatex(0.8,0.5,"(b)");
 	gPad->RedrawAxis();
  CMS_lumi( met_pad1, 11 );

	c_met->cd();
	TPad *met_pad2 = new TPad("met_pad2", "met_pad2", 0, 0, 1, 0.3);
	met_pad2->SetBottomMargin(0.4);
	met_pad2->Draw();
	met_pad2->cd();
  TLine *flatratio_met = new TLine(0,1,600,1);
	TH1F *ratio_met=(TH1F*)p_allMET->Clone("transfer factor");
	ratio_met->GetXaxis()->SetRangeUser(0,600);
	ratio_met->GetYaxis()->SetNdivisions(504);
	ratio_met->SetLineColor(kBlack);
	ratio_met->SetMarkerStyle(20);
	ratio_met->Divide(postfit_VGMET);
	ratio_met->SetTitle("");
	ratio_met->GetYaxis()->SetTitle("#frac{Obs.}{Bkg.}");
	ratio_met->GetYaxis()->SetRangeUser(0,2.3);
	ratio_met->Draw();
	ratio_met->Draw("same");
	flatratio_met->Draw("same");
	c_met->SaveAs(metplot.str().c_str());


// ******** HT ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_HT = new TCanvas("HT", "HT",600,600);
	setCanvas(c_HT); 
	c_HT->cd();
	TPad *HT_pad1 = new TPad("HT_pad1", "HT_pad1", 0, 0.3, 1, 1.0);
	setTopPad(HT_pad1); 
	HT_pad1->Draw();  
	HT_pad1->cd();  
	gPad->SetLogy();
	p_allHT->SetMinimum(0.5);
	p_allHT->SetMaximum(50*p_allHT->GetBinContent(p_allHT->GetMaximumBin()));
	p_allHT->SetLineColor(1);
	p_allHT->SetMarkerStyle(20);
	p_allHT->Draw("P");
	postfit_VGHT->SetFillStyle(1001);
	postfit_VGHT->SetLineColor(kGreen-7);
	postfit_VGHT->SetFillColor(kGreen-7);
	postfit_rareHT->SetFillStyle(1001);
	postfit_rareHT->SetLineColor(kAzure-9);
	postfit_rareHT->SetFillColor(kAzure-9);
	postfit_qcdHT->SetFillStyle(1001);
	postfit_qcdHT->SetLineColor(kYellow-9);
	postfit_qcdHT->SetFillColor(kYellow-9);
	postfit_eleHT->SetFillStyle(1001);
	postfit_eleHT->SetLineColor(kMagenta-9);
	postfit_eleHT->SetFillColor(kMagenta-9);
	postfit_jetHT->SetFillStyle(1001);
	postfit_jetHT->SetLineColor(kOrange-9);
	postfit_jetHT->SetFillColor(kOrange-9);
	postfit_eleHT->Add(postfit_rareHT); // ele 2nd
	postfit_jetHT->Add(postfit_eleHT);  // jet 3rd
	postfit_qcdHT->Add(postfit_jetHT);  // qcd 4th
	postfit_VGHT->Add(postfit_qcdHT);   // VG  5th
	postfit_VGHT->Draw("hist same");
	postfit_qcdHT->Draw("hist same");
	postfit_jetHT->Draw("hist same");
	postfit_eleHT->Draw("hist same");
	postfit_rareHT->Draw("hist same");
	p_allHT->Draw("E same");
  chantex.SetTextSize(0.05);    
  chantex.DrawLatex(0.4,0.87,"M_{T} > 100 GeV, postfit_{T}^{miss} > 120 GeV");
  chantex.SetTextSize(0.07);    
	chantex.DrawLatex(0.8,0.5,"(f)");
 	gPad->RedrawAxis();
  CMS_lumi( HT_pad1, 11 );

	c_HT->cd();
	TPad *HT_pad2 = new TPad("HT_pad2", "HT_pad2", 0, 0, 1, 0.3);
	HT_pad2->SetBottomMargin(0.4);
	HT_pad2->Draw();
	HT_pad2->cd();
  TLine *flatratio_HT = new TLine(0,1,2000,1);
	TH1F *ratio_HT=(TH1F*)p_allHT->Clone("transfer factor");
	ratio_HT->SetMarkerStyle(20);
	ratio_HT->SetLineColor(kBlack);
	ratio_HT->GetXaxis()->SetRangeUser(0,2000);
	ratio_HT->GetYaxis()->SetRangeUser(0,2);
	ratio_HT->GetYaxis()->SetNdivisions(504);
	ratio_HT->SetMinimum(0);
	ratio_HT->SetMaximum(2);
	ratio_HT->Divide(postfit_VGHT);
	ratio_HT->SetTitle("");
	ratio_HT->GetYaxis()->SetTitle("#frac{Obs.}{Bkg.}");
	ratio_HT->Draw();
	flatratio_HT->Draw("same");
	c_HT->SaveAs(htplot.str().c_str());
}


