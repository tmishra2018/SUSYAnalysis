float roundoff(float value, int prec)
{
  float pow_10 = pow(10.0f, prec);
  return round(value * pow_10) / pow_10;
}
void getSignificance(){
	gROOT->SetBatch();
	int METcut[12] = {90,100,110,120,130,140,150,160,170,180,190,200};
  	int MTcut[13] = {80,90,100,110,120,130,140,150,160,170,180,190,200};

//	int METcut[7] = {90,100,110,120,130,140,150};
// 	int MTcut[8] = {80,90,100,110,120,130,140,150};
	const char* sample[11] = {"ZGToLLG","WGJet130","WGToLNuG","DYJetsToLL","TTJets","WGJet40","WWG","WZG","TTGJets","WW","WZ"};
	TH1D *histo_bkg;
	TH1D *histo_sig, *histo_sigNoCut;
	TH2F* signif2D = new TH2F("signif2D","",12,90,210,13,80,210);
	TH2F* signalFrac2D = new TH2F("signalFrac2D","",12,90,210,13,80,210);
	//signif2D->SetTitle( "Significance (s/#sqrt{b}) for TChiWg, m_{#tilde{#chi}_{1}^{0}}  = m_{#tilde{#chi}_{1}^{+/-}} = 1000 GeV, m_{#tilde{G}} = 1 GeV ");
	signif2D->SetTitle( "Significance (s/#sqrt{b}) for T5WG, m_{#tilde{g}} = 1800 GeV, m_{#tilde{G}} = 800 GeV");
    	signif2D->SetXTitle("p_{T}^{miss} [GeV] cut at         ");
    	signif2D->SetYTitle("m_{T} [GeV] cut at          ");
	
	signalFrac2D->SetTitle( "Signal fraction (Signal with cut/Total signal) for TChiWg, m_{#tilde{#chi}_{1}^{0}}  = m_{#tilde{#chi}_{1}^{+/-}} = 1000 GeV, m_{#tilde{G}} = 1 GeV ");
	//signalFrac2D->SetTitle( "Signal fraction (Signal with cut/Total signal) for T5WG, m_{#tilde{g}} = 1800 GeV, m_{#tilde{G}} = 800 GeV");
    	signalFrac2D->SetXTitle("p_{T}^{miss} [GeV] cut at         ");
    	signalFrac2D->SetYTitle("m_{T} [GeV] cut at          ");
	for(int p=0; p<12 ; p++){                
		cout<<endl;                                                                                                              
        	for(int q=0; q<13 ; q++){
			float bkgSum=0, signif=0, signalfrac;
        		for(int s=0; s<11 ; s++){
				TFile *fbkg = TFile::Open(Form ("/eos/uscms/store/user/tmishra/egMC/Hists_eg_%s_2016_MET%d_MT%d.root",sample[s],METcut[p],MTcut[q]),"READ");
				histo_bkg = (TH1D *)fbkg->Get("hsigMET");
				bkgSum = bkgSum+histo_bkg->Integral();	
				fbkg->Close();
				delete fbkg;
			}
			
			
			//TFile *fsigNoCut = TFile::Open("/eos/uscms/store/user/tmishra/egMC/Hists_eg_T5Wg_mG-1800_mLSP-800_2016.root","READ");
			//TFile *fsig = TFile::Open(Form ("/eos/uscms/store/user/tmishra/egMC/Hists_eg_T5Wg_mG-1800_mLSP-800_2016_MET%d_MT%d.root",METcut[p],MTcut[q]),"READ");
			TFile *fsigNoCut = TFile::Open("/eos/uscms/store/user/tmishra/egMC/Hists_eg_SMS-TChiWg_mChi-1000_mLSP-1_2016.root","READ");
			TFile *fsig = TFile::Open(Form ("/eos/uscms/store/user/tmishra/egMC/Hists_eg_SMS-TChiWg_mChi-1000_mLSP-1_2016_MET%d_MT%d.root",METcut[p],MTcut[q]),"READ");
			histo_sigNoCut = (TH1D *)fsigNoCut->Get("hsigMET");
			histo_sig = (TH1D *)fsig->Get("hsigMET");
			histo_sig->Scale(100);
			
			signalfrac = histo_sig->Integral()/histo_sigNoCut->Integral();
			signalFrac2D-> SetBinContent(p+1,q+1,roundoff(signalfrac,4));
			//cout<<"("<<METcut[p]<<","<<MTcut[q]<< "):"<< roundoff(signalfrac,4)<<" ";
			cout<<signalfrac<<" ";


			signif = histo_sig->Integral()/sqrt(bkgSum);
			signif2D-> SetBinContent(p+1,q+1,roundoff(signif,4));
//			cout<<"("<<METcut[p]<<","<<MTcut[q]<< "):"<< roundoff(signif,4)<<" ";
			fsig->Close();
			fsigNoCut->Close();
			delete fsig,fsigNoCut;
		}
	}
	cout<<endl;
	TCanvas *c1 = new TCanvas("c1", "c1",65,24,2480,1356);
	c1->cd();
	signalFrac2D->SetStats(0);
	signalFrac2D->Draw("colz text");
//	signif2D->SetStats(0);
//	signif2D->Draw("colz text");
	//c1->SaveAs("2D_signif_T5Wg_mG-1800_mLSP-800.png");
	//c1->SaveAs("2D_signif_TChiWg_mChi-1000_mLSP-1.png");
//	c1->SaveAs("2D_signalFrac_T5Wg_mG-1800_mLSP-800.png");
	c1->SaveAs("2D_signalFrac_TChiWg_mChi-1000_mLSP-1.png");
	
}
