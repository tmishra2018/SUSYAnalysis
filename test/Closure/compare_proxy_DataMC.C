	// run by root -l -q "compare_proxy_DataMC.C(1,2016,1)"

	#include "../../include/analysis_commoncode.h"
	#include "TArrow.h"

	//#define NTOY 1000
	#define NTOY 1
	bool useGaussFit=false;

	#define MAXET 199
	#define MAXMT 399
	#define MAXMET 399
	#define MAXHT 399

	void compare_proxy_DataMC(int ichannel, int Year, bool ISpreVFP){
	  gROOT->SetBatch(kTRUE);
	  setTDRStyle();
	  gStyle->SetTitleXOffset(2.5);
	  gSystem->Load("../../lib/libAnaClasses.so");
	  int channelType = ichannel; // eg = 1; mg =2;

		std::string whichVFP;
		if(Year==2016 and ISpreVFP == true) whichVFP = "preVFP";
		if(Year==2016 and ISpreVFP == false) whichVFP = "postVFP";
		if(Year==2017 or  Year == 2018) whichVFP = "";


	// Signal Tree //
	//*********** histo list **********************//
	TH1D *p_PhoEt = new TH1D("p_PhoEt",";p_{T}^{#gamma} (GeV);Events / bin",nBkgEtBins,bkgEtBins);
	TH1D *p_LepPt = new TH1D("p_LepPt","",nBkgPtBins,bkgPtBins);
	TH1D *p_MET = new TH1D("p_MET","; p_{T}^{miss} (GeV);Events / bin",nBkgMETBins, bkgMETBins);
	TH1D *p_Mt = new TH1D("p_Mt","; M_{T} (GeV);Events / bin",nBkgMtBins,bkgMtBins);
	TH1D *p_HT = new TH1D("p_HT","; HT (GeV);Events / bin",nBkgHTBins, bkgHTBins); 
	TH1D *p_dPhiEleMET = new TH1D("p_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *p_PhoEta = new TH1D("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *p_LepEta = new TH1D("p_LepEta","p_LepEta",60,-3,3);
	TH1D *p_eventcount = new TH1D("p_eventcount","eventcount",9,0,9);
	TH1D *p_nJet = new TH1D("p_nJet","p_nJet",10,0,10);
	TH1D *p_nBJet = new TH1D("p_nBJet","p_nBJet",10,0,10);
	//************ Proxy Tree Data **********************//

	TChain *proxytreeData = new TChain("proxyTree");
        if(channelType==1)proxytreeData->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_%d%s.root",Year,whichVFP.c_str()));
        if(channelType==2)proxytreeData->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s_Muon20.root",Year,whichVFP.c_str()));
	float phoEt(0);
	float phoEta(0);
	float phoPhi(0);
	float lepPt(0);
	float lepEta(0);
	float lepPhi(0);
	float sigMT(0);
	float sigMET(0);
	float sigMETPhi(0);
	float dPhiLepMET(0);
	int   nVertex(0);
	float dRPhoLep(0);
	float HT(0);
	float nJetFloat(0);
	int nJetInt(0);

	int   nBJet(0);

	proxytreeData->SetBranchAddress("phoEt",     &phoEt);
	proxytreeData->SetBranchAddress("phoEta",    &phoEta);
	proxytreeData->SetBranchAddress("phoPhi",    &phoPhi);
	proxytreeData->SetBranchAddress("lepPt",     &lepPt);
	proxytreeData->SetBranchAddress("lepEta",    &lepEta);
	proxytreeData->SetBranchAddress("lepPhi",    &lepPhi);
	proxytreeData->SetBranchAddress("sigMT",     &sigMT);
	proxytreeData->SetBranchAddress("sigMET",    &sigMET);
	proxytreeData->SetBranchAddress("sigMETPhi", &sigMETPhi);
	proxytreeData->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
	proxytreeData->SetBranchAddress("nVertex",   &nVertex);
	proxytreeData->SetBranchAddress("dRPhoLep",  &dRPhoLep);
	proxytreeData->SetBranchAddress("HT",        &HT);
	proxytreeData->SetBranchAddress("nBJet",     &nBJet);
	if (channelType == 1) proxytreeData->SetBranchAddress("nJet", &nJetFloat);
	else proxytreeData->SetBranchAddress("nJet", &nJetInt);

	for (unsigned ievt(0); ievt<proxytreeData->GetEntries(); ++ievt){//loop on entries
		proxytreeData->GetEntry(ievt);
		if (channelType == 1 && nJetFloat <1 ) continue; // suggestion from convenors
		if (channelType == 2 && nJetInt <1 ) continue;
		//p_PU->Fill(nVertex);
		/** cut flow *****/
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
		if(sigMET < lowMET)continue;
		if(highMET > 0 && sigMET > highMET)continue;
		if(sigMT < lowMt)continue;
		if(highMt > 0 && sigMT > highMt)continue;
		if(lepPt < lowPt)continue;
		if(highPt > 0 && lepPt > highPt)continue;
		// different MET, MT and lep pT cuts
		p_PhoEt->Fill(phoEt);
		p_PhoEta->Fill(phoEta);
		p_MET->Fill(sigMET);
		p_Mt->Fill(sigMT);
		p_HT->Fill(HT);
		p_LepPt->Fill(lepPt);
		p_LepEta->Fill(lepEta);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET));
		if (channelType == 1) p_nJet->Fill(nJetFloat);
                if (channelType == 2) p_nJet->Fill(nJetInt);

		p_nBJet->Fill(nBJet);
		// ttbar events when nBJets >= 1
	}


	// Fake Tree //
	//*********** histo list **********************//
	std::ostringstream histname;
	TH1D *pred_PhoEt = new TH1D("pred_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *pred_LepPt = new TH1D("pred_LepPt","pred_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *pred_MET = new TH1D("pred_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *pred_Mt = new TH1D("pred_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *pred_HT = new TH1D("pred_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *pred_PhoEta = new TH1D("pred_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *pred_LepEta = new TH1D("pred_LepEta","pred_LepEta",60,-3,3);
	TH1D *pred_dPhiEleMET = new TH1D("pred_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *pred_nJet = new TH1D("pred_nJet","pred_nJet",10,0,10);

	TH1D *DY_PhoEt = new TH1D("DY_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *DY_LepPt = new TH1D("DY_LepPt","DY_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *DY_MET = new TH1D("DY_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *DY_Mt = new TH1D("DY_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *DY_HT = new TH1D("DY_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *DY_PhoEta = new TH1D("DY_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *DY_LepEta = new TH1D("DY_LepEta","DY_LepEta",60,-3,3);
	TH1D *DY_dPhiEleMET = new TH1D("DY_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *DY_nJet = new TH1D("DY_nJet","DY_nJet",10,0,10);

	//************ Proxy Tree **********************//
	TChain *proxytree = new TChain("proxyTree");
	if(channelType==1)proxytree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_DYJetsToLL_%d%s.root",Year,whichVFP.c_str()));
	if(channelType==2)proxytree->Add(Form("/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_DYJetsToLL_%d%s.root",Year,whichVFP.c_str()));

	float proxycrosssection(0);
	float proxyntotalevent(0);
	float proxyphoEt(0);
	float proxyphoEta(0);
	float proxyphoPhi(0);
	float proxylepPt(0);
	float proxylepEta(0);
	float proxylepPhi(0);
	float proxysigMT(0);
	float proxysigMET(0);
	float proxysigMETPhi(0);
	float proxydPhiLepMET(0);
	int   proxynVertex(0);
	float proxydRPhoLep(0);
	float proxyHT(0);
	float proxynJet(0);
	
	proxytree->SetBranchAddress("crosssection",&proxycrosssection);
	proxytree->SetBranchAddress("ntotalevent", &proxyntotalevent);
	proxytree->SetBranchAddress("phoEt",     	 &proxyphoEt);
	proxytree->SetBranchAddress("phoEta",    	 &proxyphoEta);
	proxytree->SetBranchAddress("phoPhi",    	 &proxyphoPhi);
	proxytree->SetBranchAddress("lepPt",     	 &proxylepPt);
	proxytree->SetBranchAddress("lepEta",    	 &proxylepEta);
	proxytree->SetBranchAddress("lepPhi",    	 &proxylepPhi);
	proxytree->SetBranchAddress("sigMT",     	 &proxysigMT);
	proxytree->SetBranchAddress("sigMET",    	 &proxysigMET);
	proxytree->SetBranchAddress("sigMETPhi", 	 &proxysigMETPhi);
	proxytree->SetBranchAddress("dPhiLepMET",	 &proxydPhiLepMET);
	proxytree->SetBranchAddress("nVertex",   	 &proxynVertex);
	proxytree->SetBranchAddress("dRPhoLep",  	 &proxydRPhoLep);
	proxytree->SetBranchAddress("HT",        	 &proxyHT);
	proxytree->SetBranchAddress("nJet",      	 &proxynJet);

	 cout<<"Year "<<Year<<endl<<endl;
	for (unsigned ievt(0); ievt<proxytree->GetEntries(); ++ievt){//loop on entries
		proxytree->GetEntry(ievt);
		
		double weight = 1;
		if(channelType == 1){
                        if(Year == 2016 and ISpreVFP == 1)             weight = lumi_2016preVFP_DoubleEG*1000*proxycrosssection/proxyntotalevent;
                        else if(Year == 2016 and ISpreVFP == 0)        weight = lumi_2016postVFP_DoubleEG*1000*proxycrosssection/proxyntotalevent;
                        else if(Year == 2017)                        weight = lumi_2017_DoubleEG*1000*proxycrosssection/proxyntotalevent;
                        else if(Year == 2018)                        weight = lumi_2018_DoubleEG*1000*proxycrosssection/proxyntotalevent;}

                else if(channelType == 2){
                        if(Year == 2016 and ISpreVFP == 1)             weight = lumi_2016preVFP_MuonEG*1000*proxycrosssection/proxyntotalevent;
                        else if(Year == 2016 and ISpreVFP == 0)        weight = lumi_2016postVFP_MuonEG*1000*proxycrosssection/proxyntotalevent;
                        else if(Year == 2017)                        weight = lumi_2017_MuonEG*1000*proxycrosssection/proxyntotalevent;
                        else if(Year == 2018)                        weight = lumi_2018_MuonEG*1000*proxycrosssection/proxyntotalevent;}
		
		if(ievt==0) cout<<"DY cross-section "<<proxycrosssection<<endl;

		if(proxyphoEt > MAXET)proxyphoEt = MAXET;
		if(proxysigMET > MAXMET)proxysigMET = MAXMET;
		if(proxysigMT > MAXMT)proxysigMT = MAXMT;
		if(proxyHT > MAXHT)proxyHT = MAXHT;	
		
		/** cut flow *****/
		if ((ichannel == 1 && (proxyphoEt < 35 || proxylepPt < 25)) || (ichannel == 2 && (proxyphoEt < 35 || proxylepPt < 20))) continue;
		if(fabs(proxyphoEta) > 1.4442 || fabs(proxylepEta) > 2.5)continue;
		pred_PhoEt->Fill(proxyphoEt,weight);
		pred_PhoEta->Fill(proxyphoEta, weight);
		pred_MET->Fill(proxysigMET, weight);
		pred_Mt->Fill(proxysigMT, weight);
		pred_HT->Fill(proxyHT, weight);
		pred_LepPt->Fill(proxylepPt, weight);
		pred_LepEta->Fill(proxylepEta, weight);
		pred_dPhiEleMET->Fill(fabs(proxydPhiLepMET), weight);
		pred_nJet->Fill(proxynJet, weight);

		DY_PhoEt->Fill(proxyphoEt,weight);
		DY_PhoEta->Fill(proxyphoEta, weight);
		DY_MET->Fill(proxysigMET, weight);
		DY_Mt->Fill(proxysigMT, weight);
		DY_HT->Fill(proxyHT, weight);
		DY_LepPt->Fill(proxylepPt, weight);
		DY_LepEta->Fill(proxylepEta, weight);
		DY_dPhiEleMET->Fill(fabs(proxydPhiLepMET), weight);
		DY_nJet->Fill(proxynJet, weight);

	}

	//************ Proxy Tree **********************//
	TChain *raretree = new TChain("proxyTree");
	// TTJets, WW, WZ rarely contribute to electron fake photon background, proxy event
	if(channelType==1)raretree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_TTJets_%d%s.root",Year,whichVFP.c_str()));
	if(channelType==1)raretree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_WW_%d%s.root",Year,whichVFP.c_str()));
	if(channelType==1)raretree->Add(Form("/eos/uscms/store/user/tmishra/egMC/resTree_egsignal_WZ_%d%s.root",Year,whichVFP.c_str()));

	if(channelType==2)raretree->Add(Form("/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_TTJets_%d%s.root",Year,whichVFP.c_str()));
	if(channelType==2)raretree->Add(Form("/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_WW_%d%s.root",Year,whichVFP.c_str()));
	if(channelType==2)raretree->Add(Form("/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_WZ_%d%s.root",Year,whichVFP.c_str()));

	float rarecrosssection(0);
	float rarentotalevent(0);
	float rarephoEt(0);
	float rarephoEta(0);
	float rarephoPhi(0);
	float rarelepPt(0);
	float rarelepEta(0);
	float rarelepPhi(0);
	float raresigMT(0);
	float raresigMET(0);
	float raresigMETPhi(0);
	float raredPhiLepMET(0);
	int   rarenVertex(0);
	float raredRPhoLep(0);
	float rareHT(0);
	float rarenJet(0);
	
	raretree->SetBranchAddress("crosssection",&rarecrosssection);
	raretree->SetBranchAddress("ntotalevent", &rarentotalevent);
	raretree->SetBranchAddress("phoEt",     	 &rarephoEt);
	raretree->SetBranchAddress("phoEta",    	 &rarephoEta);
	raretree->SetBranchAddress("phoPhi",    	 &rarephoPhi);
	raretree->SetBranchAddress("lepPt",     	 &rarelepPt);
	raretree->SetBranchAddress("lepEta",    	 &rarelepEta);
	raretree->SetBranchAddress("lepPhi",    	 &rarelepPhi);
	raretree->SetBranchAddress("sigMT",     	 &raresigMT);
	raretree->SetBranchAddress("sigMET",    	 &raresigMET);
	raretree->SetBranchAddress("sigMETPhi", 	 &raresigMETPhi);
	raretree->SetBranchAddress("dPhiLepMET",	 &raredPhiLepMET);
	raretree->SetBranchAddress("nVertex",   	 &rarenVertex);
	raretree->SetBranchAddress("dRPhoLep",  	 &raredRPhoLep);
	raretree->SetBranchAddress("HT",        	 &rareHT);
	raretree->SetBranchAddress("nJet",      	 &rarenJet);

	for (unsigned ievt(0); ievt<raretree->GetEntries(); ++ievt){//loop on entries
		raretree->GetEntry(ievt);
		
		double weight = 1;
		if(channelType == 1){
                        if(Year == 2016 and ISpreVFP == 1)             weight = lumi_2016preVFP_DoubleEG*1000*rarecrosssection/rarentotalevent;
                        else if(Year == 2016 and ISpreVFP == 0)        weight = lumi_2016postVFP_DoubleEG*1000*rarecrosssection/rarentotalevent;
                        else if(Year == 2017)                        weight = lumi_2017_DoubleEG*1000*rarecrosssection/rarentotalevent;
                        else if(Year == 2018)                        weight = lumi_2018_DoubleEG*1000*rarecrosssection/rarentotalevent;}

                else if(channelType == 2){
                        if(Year == 2016 and ISpreVFP == 1)             weight = lumi_2016preVFP_MuonEG*1000*rarecrosssection/rarentotalevent;
                        else if(Year == 2016 and ISpreVFP == 0)        weight = lumi_2016postVFP_MuonEG*1000*rarecrosssection/rarentotalevent;
                        else if(Year == 2017)                        weight = lumi_2017_MuonEG*1000*rarecrosssection/rarentotalevent;
                        else if(Year == 2018)                        weight = lumi_2018_MuonEG*1000*rarecrosssection/rarentotalevent;}
		
		if(ievt==0) cout<<"tt cross-section "<<rarecrosssection<<endl<<endl;
		if(rarephoEt > MAXET)rarephoEt = MAXET;
		if(raresigMET > MAXMET)raresigMET = MAXMET;
		if(raresigMT > MAXMT)raresigMT = MAXMT;
		if(rareHT > MAXHT)rareHT = MAXHT;	
		

		/** cut flow *****/
		if ((ichannel == 1 && (rarephoEt < 35 || rarelepPt < 25)) || (ichannel == 2 && (rarephoEt < 35 || rarelepPt < 20))) continue;
		if(fabs(rarephoEta) > 1.4442 || fabs(rarelepEta) > 2.5)continue;
		
		pred_PhoEt->Fill(rarephoEt,weight);
		pred_PhoEta->Fill(rarephoEta, weight);
		pred_MET->Fill(raresigMET, weight);
		pred_Mt->Fill(raresigMT, weight);
		pred_HT->Fill(rareHT, weight);
		pred_LepPt->Fill(rarelepPt, weight);
		pred_LepEta->Fill(rarelepEta, weight);
		pred_dPhiEleMET->Fill(fabs(raredPhiLepMET), weight);
		pred_nJet->Fill(rarenJet, weight);

	}

		
	TGraphErrors *error_PhoEt = new TGraphErrors(nBkgEtBins);
	TGraphErrors *error_MET = new TGraphErrors(nBkgMETBins); 
	TGraphErrors *error_Mt = new TGraphErrors(nBkgMtBins); 
	TGraphErrors *error_HT = new TGraphErrors(nBkgHTBins);
	TGraphErrors *ratioerror_PhoEt = new TGraphErrors(nBkgEtBins);
	TGraphErrors *ratioerror_MET = new TGraphErrors(nBkgMETBins); 
	TGraphErrors *ratioerror_Mt = new TGraphErrors(nBkgMtBins); 
	TGraphErrors *ratioerror_HT = new TGraphErrors(nBkgHTBins);
// from signal events
	TH1F *ratio=(TH1F*)p_PhoEt->Clone("transfer factor");
	TH1F *ratio_met=(TH1F*)p_MET->Clone("transfer factor");
	TH1F *ratio_mt=(TH1F*)p_Mt->Clone("transfer factor");
	TH1F *ratio_HT=(TH1F*)p_HT->Clone("transfer factor");

	TCanvas *c_pt = new TCanvas("Photon_Pt", "Photon P_{T}",600,600);  
	setCanvas(c_pt); 
	c_pt->cd();                                                             
	TPad *pt_pad1 = new TPad("pt_pad1", "pt_pad1", 0, 0.35, 1, 1.0);        
	setTopPad(pt_pad1); 
	pt_pad1->SetBottomMargin(0);
	pt_pad1->Draw();                                                        
	pt_pad1->cd();                                                          
	gPad->SetLogy();
	p_PhoEt->SetMaximum(200*p_PhoEt->GetBinContent(p_PhoEt->GetMaximumBin())); 
	p_PhoEt->GetXaxis()->SetRangeUser(35,200);                              
	p_PhoEt->GetXaxis()->SetLabelSize(0);
	p_PhoEt->Draw();                                                        
	p_PhoEt->SetLineColor(kBlack);                                          
	p_PhoEt->SetMarkerStyle(20);                                            
	pred_PhoEt->SetFillStyle(1001);                                         
	pred_PhoEt->SetLineColor(kRed);                                         
	pred_PhoEt->SetFillColor(kRed);                                         
	DY_PhoEt->SetFillStyle(1001);                                           
	DY_PhoEt->SetLineColor(kYellow);
	DY_PhoEt->SetFillColor(kYellow);
// from proxy events, predicted backgrounds
	pred_PhoEt->Draw("hist same");
	DY_PhoEt->Draw("hist same");
	TLegend *leg =  new TLegend(0.5,0.55,0.9,0.8);
	leg->SetFillStyle(0);
	gStyle->SetLegendBorderSize(1);
	gStyle->SetLegendFillColor(0);
	pred_PhoEt->SetMarkerColor(kRed);
	DY_PhoEt->SetMarkerColor(kYellow);
	error_PhoEt->SetMarkerSize(0);
	error_PhoEt->SetLineWidth(0);
	leg->AddEntry((TObject*)0, "#bf{e fake #gamma proxy events}", "");
	leg->AddEntry(p_PhoEt,"Data", "pel");
	leg->AddEntry(pred_PhoEt,"t#bar{t} / WW / WZ");
	leg->AddEntry(DY_PhoEt,"Drell#scale[0.5]{ }#minusYan");

	//leg->AddEntry(error_PhoEt, "Total uncertainty");
	leg->Draw("same");
  	TLatex chantex;
  	chantex.SetNDC();
  	chantex.SetTextFont(42);
  	chantex.SetTextSize(0.07);    
  	if(channelType==1) chantex.DrawLatex(0.58,0.82," e + #gamma");
  	if(channelType==2) chantex.DrawLatex(0.58,0.82," #mu + #gamma");
 	gPad->RedrawAxis();
	if(Year==2016 and ISpreVFP == 1)       CMS_lumi( pt_pad1, 1,ichannel, 11 );
        else if(Year==2016 and ISpreVFP == 0)  CMS_lumi( pt_pad1, 2,ichannel, 11 );
        else if(Year==2017)                  CMS_lumi( pt_pad1, 3,ichannel, 11 );
        else if(Year==2018)                  CMS_lumi( pt_pad1, 4,ichannel, 11 );
		
	p_PhoEt->Draw("same");
	for(int ibin(1); ibin < pred_PhoEt->GetSize(); ibin++){
		error_PhoEt->SetPoint(ibin-1,pred_PhoEt->GetBinCenter(ibin), pred_PhoEt->GetBinContent(ibin));
		// syst + stat error
		float prederror = pred_PhoEt->GetBinError(ibin);
		error_PhoEt->SetPointError(ibin-1,(pred_PhoEt->GetBinLowEdge(ibin+1)-pred_PhoEt->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_PhoEt->SetPoint(ibin-1,pred_PhoEt->GetBinCenter(ibin), 1); 
		ratioerror_PhoEt->SetPointError(ibin-1,(pred_PhoEt->GetBinLowEdge(ibin+1)-pred_PhoEt->GetBinLowEdge(ibin))/2, prederror/pred_PhoEt->GetBinContent(ibin)); 
	}
	error_PhoEt->SetFillColor(12);
  	error_PhoEt->SetFillStyle(3345);
	error_PhoEt->Draw("E2 same");

	c_pt->cd();
	TPad *pt_pad2 = new TPad("pt_pad2", "pt_pad2", 0, 0, 1, 0.35);
	pt_pad2->SetTopMargin(0);
	pt_pad2->SetBottomMargin(0.3);
	pt_pad2->Draw();
	pt_pad2->cd();
  	TLine *flatratio = new TLine(35,1,200,1);
	ratio->GetXaxis()->SetRangeUser(35,200);
	ratio->GetYaxis()->SetRangeUser(0.,1.9);
	ratio->SetMarkerStyle(20);
	ratio->SetLineColor(kBlack);
	// from signal events - predicted events
	ratio->Divide(pred_PhoEt);
	ratio->SetTitle("");
	ratio->GetXaxis()->SetTitleOffset(0.9);
	ratio->GetXaxis()->SetTitle("p_{T}^{#gamma} (GeV)");
	ratio->GetYaxis()->SetTitle("#frac{Data}{Simulation}  ");
	ratio->GetYaxis()->SetNdivisions(504);
	ratio->Draw();
	ratioerror_PhoEt->SetFillColor(12);
	ratioerror_PhoEt->SetFillStyle(3345);
	//ratioerror_PhoEt->Draw("E2 same");
	ratio->Draw("same");
	flatratio->Draw("same");
	TLine *ratioValue_PhoEt = new TLine(35,p_PhoEt->Integral()/pred_PhoEt->Integral(),200,p_PhoEt->Integral()/pred_PhoEt->Integral());
        ratioValue_PhoEt->SetLineColor(kRed);
        ratioValue_PhoEt->Draw("same");
	
	if(channelType==1) c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Closure/compare_elefakepho-proxy_DataMC_PhotonEt_eg_%d%s.pdf",Year,whichVFP.c_str()));
	if(channelType==2) c_pt->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Closure/compare_elefakepho-proxy_DataMC_PhotonEt_mg_%d%s.pdf",Year,whichVFP.c_str()));


// similar for other variables
// ******** MET ************************//
	TCanvas *c_met = new TCanvas("MET", "MET",600,600);
	c_met->cd();
	TPad *met_pad1 = new TPad("met_pad1", "met_pad1", 0, 0.35, 1, 1.0);
	setTopPad(met_pad1); 
	met_pad1->SetBottomMargin(0);
	met_pad1->Draw();  
	met_pad1->cd();  
	gPad->SetLogy();
	cout<<"Direct simulation "<< p_MET->Integral()<<endl;
	//cout<<"DY " <<DY_MET->Integral()<<endl;
	cout<<"t#bar{t} / WW/ WZ "<<pred_MET->Integral()<<endl;
	cout<<"Ratio "<<p_MET->Integral()/pred_MET->Integral()<<endl;
	p_MET->SetMaximum(1000*p_MET->GetBinContent(p_MET->GetMaximumBin())); 
	p_MET->SetMinimum(0.05);
	p_MET->GetXaxis()->SetRangeUser(0,400);
	p_MET->GetXaxis()->SetLabelSize(0);
	p_MET->SetLineColor(1);
	p_MET->SetMarkerStyle(20);
	p_MET->Draw("P");
	DY_MET->SetFillStyle(1001);
	DY_MET->SetLineColor(kYellow-4);
	DY_MET->SetFillColor(kYellow-4);
	pred_MET->SetFillStyle(1001);
	pred_MET->SetLineColor(kRed);
	pred_MET->SetFillColor(kRed);
	for(int ibin(1); ibin < pred_MET->GetSize(); ibin++){
		float prederror = pred_MET->GetBinError(ibin);
		// Total stat and syst error
		error_MET->SetPoint(ibin-1,pred_MET->GetBinCenter(ibin), pred_MET->GetBinContent(ibin));
		error_MET->SetPointError(ibin-1,(pred_MET->GetBinLowEdge(ibin+1)-pred_MET->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_MET->SetPoint(ibin-1,pred_MET->GetBinCenter(ibin), 1); 
		ratioerror_MET->SetPointError(ibin-1,(pred_MET->GetBinLowEdge(ibin+1)-pred_MET->GetBinLowEdge(ibin))/2, prederror/pred_MET->GetBinContent(ibin)); 
	}
	pred_MET->Draw("hist same");
	DY_MET->Draw("hist same");
  	error_MET->SetFillColor(12);
  	error_MET->SetFillStyle(3345);
	error_MET->Draw("E2 same");
	leg->Draw("same");
	p_MET->Draw("E same");
  	if(channelType==1) chantex.DrawLatex(0.58,0.82," e + #gamma");
  	if(channelType==2) chantex.DrawLatex(0.58,0.82," #mu + #gamma");
	//TLine *line_met = new TLine(70,0,70,10000);
	TLine *line_met = new TLine(70,0.05,70,1000*p_MET->GetBinContent(p_MET->GetMaximumBin()));
	line_met->SetLineStyle(2);
	line_met->Draw("same");
	TLatex* latex = new TLatex();
	latex->SetTextSize(0.05);
	latex->DrawLatex(20, 60000,"control");
	latex->DrawLatex(20, 30000,"region");
 	gPad->RedrawAxis();
	if(Year==2016 and ISpreVFP == 1)       CMS_lumi( met_pad1, 1,ichannel, 11 );
        else if(Year==2016 and ISpreVFP == 0)  CMS_lumi( met_pad1, 2,ichannel, 11 );
        else if(Year==2017)                  CMS_lumi( met_pad1, 3,ichannel, 11 );
        else if(Year==2018)                  CMS_lumi( met_pad1, 4,ichannel, 11 );

	c_met->cd();
	TPad *met_pad2 = new TPad("met_pad2", "met_pad2", 0, 0, 1, 0.35);
	met_pad2->SetTopMargin(0);
	met_pad2->SetBottomMargin(0.3);
	met_pad2->Draw();
	met_pad2->cd();
  TLine *flatratio_met = new TLine(0,1,400,1);
	ratio_met->GetXaxis()->SetRangeUser(0,400);
	ratio_met->GetYaxis()->SetNdivisions(504);
	ratio_met->SetLineColor(kBlack);
	ratio_met->SetMarkerStyle(20);
	ratio_met->Divide(pred_MET);
	ratio_met->SetTitle("");
	ratio_met->GetXaxis()->SetTitleOffset(0.9);
	ratio_met->GetXaxis()->SetTitle("p_{T}^{miss} (GeV)");
	ratio_met->GetYaxis()->SetTitle("#frac{Data}{Simulation} ");
	ratio_met->GetYaxis()->SetRangeUser(0.,1.9);
	ratio_met->Draw();
	ratioerror_MET->SetFillColor(12);
	ratioerror_MET->SetFillStyle(3345);
	//ratioerror_MET->Draw("E2 same");
	ratio_met->Draw("same");
	flatratio_met->Draw("same");
	
	TLine *ratioValue_met = new TLine(0,p_MET->Integral()/pred_MET->Integral(),400,p_MET->Integral()/pred_MET->Integral());
        ratioValue_met->SetLineColor(kRed);
        ratioValue_met->Draw("same");

	TLine *line_met_ratio = new TLine(70,0.4,70,1.7);
	line_met_ratio->SetLineStyle(2);
	line_met_ratio->Draw("same");
	double ratioValue_MET = (pred_MET->Integral() != 0) ? p_MET->Integral() / pred_MET->Integral() : 0;
	TLatex *text_MET = new TLatex(250, ratioValue_MET + 0.05, Form("Ratio = %.2f", ratioValue_MET));
	text_MET->SetTextSize(0.12);
	text_MET->SetTextColor(kRed);
	text_MET->Draw();

 	gPad->RedrawAxis();
	if(channelType==1) c_met->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Closure/compare_elefakepho-proxy_DataMC_MET_eg_%d%s.pdf",Year,whichVFP.c_str()));
	if(channelType==2) c_met->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Closure/compare_elefakepho-proxy_DataMC_MET_mg_%d%s.pdf",Year,whichVFP.c_str()));

// ******** Mt ************************//
	TCanvas *c_mt = new TCanvas("Mt", "Mt",600,600);
	c_mt->cd();
	TPad *mt_pad1 = new TPad("mt_pad1", "mt_pad1", 0, 0.35, 1, 1.0);
	setTopPad(mt_pad1); 
	mt_pad1->SetBottomMargin(0);
	mt_pad1->Draw();  
	mt_pad1->cd();  
	gPad->SetLogy();
	p_Mt->SetMinimum(0.05);
	p_Mt->SetMaximum(100*p_Mt->GetBinContent(p_Mt->GetMaximumBin())); 
	p_Mt->GetXaxis()->SetRangeUser(0,400);
	p_Mt->GetXaxis()->SetLabelSize(0);
	p_Mt->SetLineColor(1);
	p_Mt->SetMarkerStyle(20);
	p_Mt->Draw("P");
	DY_Mt->SetFillStyle(1001);
	DY_Mt->SetLineColor(kYellow-4);
	DY_Mt->SetFillColor(kYellow-4);
	pred_Mt->SetFillStyle(1001);
	pred_Mt->SetLineColor(kRed);
	pred_Mt->SetFillColor(kRed);
	for(int ibin(1); ibin < pred_Mt->GetSize(); ibin++){
		error_Mt->SetPoint(ibin-1,pred_Mt->GetBinCenter(ibin), pred_Mt->GetBinContent(ibin));
		float prederror = pred_Mt->GetBinError(ibin);
		error_Mt->SetPointError(ibin-1,(pred_Mt->GetBinLowEdge(ibin+1)-pred_Mt->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_Mt->SetPoint(ibin-1,pred_Mt->GetBinCenter(ibin), 1); 
		ratioerror_Mt->SetPointError(ibin-1,(pred_Mt->GetBinLowEdge(ibin+1)-pred_Mt->GetBinLowEdge(ibin))/2, prederror/pred_Mt->GetBinContent(ibin)); 
	}
	pred_Mt->Draw("hist same");
	DY_Mt->Draw("hist same");
  error_Mt->SetFillColor(12);
  error_Mt->SetFillStyle(3345);
	error_Mt->Draw("E2 same");
	leg->Draw("same");
	p_Mt->Draw("E same");
  	if(channelType==1) chantex.DrawLatex(0.58,0.82," e + #gamma");
  	if(channelType==2) chantex.DrawLatex(0.58,0.82," #mu + #gamma");
 	gPad->RedrawAxis();
	if(Year==2016 and ISpreVFP == 1)       CMS_lumi( mt_pad1, 1,ichannel, 11 );
        else if(Year==2016 and ISpreVFP == 0)  CMS_lumi( mt_pad1, 2,ichannel, 11 );
        else if(Year==2017)                  CMS_lumi( mt_pad1, 3,ichannel, 11 );
        else if(Year==2018)                  CMS_lumi( mt_pad1, 4,ichannel, 11 );

	c_mt->cd();
	TPad *mt_pad2 = new TPad("mt_pad2", "mt_pad2", 0, 0, 1, 0.35);
	mt_pad2->SetTopMargin(0);
	mt_pad2->SetBottomMargin(0.3);
	mt_pad2->Draw();
	mt_pad2->cd();
  TLine *flatratio_mt = new TLine(0,1,400,1);
	ratio_mt->SetMarkerStyle(20);
	ratio_mt->SetLineColor(kBlack);
	ratio_mt->GetXaxis()->SetRangeUser(0,400);
	ratio_mt->GetYaxis()->SetRangeUser(0.,1.9);
	ratio_mt->GetYaxis()->SetNdivisions(504);
	ratio_mt->Divide(pred_Mt);
	ratio_mt->SetTitle("");
	ratio_mt->GetXaxis()->SetTitleOffset(0.9);
	ratio_mt->GetXaxis()->SetTitle("M_{T} (GeV)");
	ratio_mt->GetYaxis()->SetTitle("#frac{Data}{Simulation} ");
	ratio_mt->Draw();
	ratioerror_Mt->SetFillColor(12);
	ratioerror_Mt->SetFillStyle(3345);
	//ratioerror_Mt->Draw("E2 same");
	ratio_mt->Draw("same");

	TLine *ratioValue_mt = new TLine(0,p_Mt->Integral()/pred_Mt->Integral(),400,p_Mt->Integral()/pred_Mt->Integral());
        ratioValue_mt->SetLineColor(kRed);
        ratioValue_mt->Draw("same");

	double ratioValue_MT = (pred_Mt->Integral() != 0) ? p_Mt->Integral() / pred_Mt->Integral() : 0;
	TLatex *text_MT = new TLatex(250, ratioValue_MT + 0.05, Form("Ratio = %.2f", ratioValue_MT));
	text_MT->SetTextSize(0.12);
	text_MT->SetTextColor(kRed);
	text_MT->Draw();

	flatratio_mt->Draw("same");
	if(channelType==1) c_mt->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Closure/compare_elefakepho-proxy_DataMC_MT_eg_%d%s.pdf",Year,whichVFP.c_str()));
	if(channelType==2) c_mt->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Closure/compare_elefakepho-proxy_DataMC_MT_mg_%d%s.pdf",Year,whichVFP.c_str()));

// ******** HT ************************//
	TCanvas *c_HT = new TCanvas("HT", "HT",600,600);
	c_HT->cd();
	TPad *HT_pad1 = new TPad("HT_pad1", "HT_pad1", 0, 0.35, 1, 1.0);
	setTopPad(HT_pad1); 
	HT_pad1->SetBottomMargin(0);
	HT_pad1->Draw();  
	HT_pad1->cd();  
	gPad->SetLogy();
	p_HT->SetMaximum(100*p_HT->GetBinContent(p_HT->GetMaximumBin())); 
	p_HT->GetXaxis()->SetRangeUser(0,400);
	p_HT->GetXaxis()->SetLabelSize(0);
	p_HT->SetMinimum(5);
	p_HT->SetLineColor(1);
	p_HT->SetMarkerStyle(20);
	p_HT->Draw("P");
	DY_HT->SetFillStyle(1001);
	DY_HT->SetLineColor(kYellow-4);
	DY_HT->SetFillColor(kYellow-4);
	pred_HT->SetFillStyle(1001);
	pred_HT->SetLineColor(kRed);
	pred_HT->SetFillColor(kRed);
	for(int ibin(1); ibin < pred_HT->GetSize(); ibin++){
		error_HT->SetPoint(ibin-1,pred_HT->GetBinCenter(ibin), pred_HT->GetBinContent(ibin));
		float prederror = pred_HT->GetBinError(ibin);
		error_HT->SetPointError(ibin-1,(pred_HT->GetBinLowEdge(ibin+1)-pred_HT->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_HT->SetPoint(ibin-1,pred_HT->GetBinCenter(ibin), 1); 
		ratioerror_HT->SetPointError(ibin-1,(pred_HT->GetBinLowEdge(ibin+1)-pred_HT->GetBinLowEdge(ibin))/2, prederror/pred_HT->GetBinContent(ibin)); 
	}
	pred_HT->Draw("hist same");
	DY_HT->Draw("hist same");
	error_HT->SetFillColor(12);
  	error_HT->SetFillStyle(3345);
	error_HT->Draw("E2 same");
	leg->Draw("same");
	p_HT->Draw("E same");
	  if(channelType==1) chantex.DrawLatex(0.58,0.82," e + #gamma");
	  if(channelType==2) chantex.DrawLatex(0.58,0.82," #mu + #gamma");
 	gPad->RedrawAxis();
	if(Year==2016 and ISpreVFP == 1)       CMS_lumi( HT_pad1, 1, ichannel, 11 );
        else if(Year==2016 and ISpreVFP == 0)  CMS_lumi( HT_pad1, 2, ichannel, 11 );
        else if(Year==2017)                  CMS_lumi( HT_pad1, 3, ichannel, 11 );
        else if(Year==2018)                  CMS_lumi( HT_pad1, 4, ichannel, 11 );

	c_HT->cd();
	TPad *HT_pad2 = new TPad("HT_pad2", "HT_pad2", 0, 0, 1, 0.35);
	HT_pad2->SetTopMargin(0);
	HT_pad2->SetBottomMargin(0.3);
	HT_pad2->Draw();
	HT_pad2->cd();
  	TLine *flatratio_HT = new TLine(0,1,420,1);
	ratio_HT->SetMarkerStyle(20);
	ratio_HT->SetLineColor(kBlack);
	ratio_HT->GetXaxis()->SetRangeUser(0,420);
	ratio_HT->GetYaxis()->SetRangeUser(0.,1.9);
	ratio_HT->GetYaxis()->SetNdivisions(504);
	ratio_HT->Divide(pred_HT);
	ratio_HT->SetTitle("");
	ratio_HT->GetXaxis()->SetTitleOffset(0.9);
	ratio_HT->GetXaxis()->SetTitle("H_{T} (GeV)");
	ratio_HT->GetYaxis()->SetTitle("#frac{Data}{Simulation} ");
	ratio_HT->Draw();
	ratioerror_HT->SetFillColor(12);
	ratioerror_HT->SetFillStyle(3345);
	//ratioerror_HT->Draw("E2 same");
	flatratio_HT->Draw("same");
	TLine *ratioValue_ht = new TLine(0,p_HT->Integral()/pred_HT->Integral(),420,p_HT->Integral()/pred_HT->Integral());
        ratioValue_ht->SetLineColor(kRed);
        ratioValue_ht->Draw("same");
	
	double ratioValue_HT = (pred_HT->Integral() != 0) ? p_HT->Integral() / pred_HT->Integral() : 0;
	TLatex *text_HT = new TLatex(250, ratioValue_HT + 0.05, Form("Ratio = %.2f", ratioValue_HT));
	text_HT->SetTextSize(0.12);
	text_HT->SetTextColor(kRed);
	text_HT->Draw();


	if(channelType==1) c_HT->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Closure/compare_elefakepho-proxy_DataMC_HT_eg_%d%s.pdf",Year,whichVFP.c_str()));
	if(channelType==2) c_HT->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Closure/compare_elefakepho-proxy_DataMC_HT_mg_%d%s.pdf",Year,whichVFP.c_str()));

}
