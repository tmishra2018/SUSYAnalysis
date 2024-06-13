#define MultiJet_cxx
#include "MultiJet.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <fstream>
#include "btag/BTagCorrector.h"

using namespace std;

int main(int argc, char* argv[])
{

  if (argc < 2) {
    cerr << "Please give 3 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" << endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];

  MultiJet ana(inputFileList, outFileName, data);
  cout << "dataset " << data << " " << endl;

  ana.EventLoop(data,inputFileList);

  return 0;
}

void MultiJet::EventLoop(const char *data,const char *inputFileList) {
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;

  TString s_data=data;

  Long64_t nbytes = 0, nb = 0;
  int decade = 0;
  bool do_AB_reweighting=0;
  int evtSurvived=0;
  TFile *f_HLR=TFile::Open("HLR_gjets_qcd.root");
  // TH1D *h_HLratio=(TH1D*)f_HLR->Get("HLratio_1D");
  TH2D *h2_HLRatio;

  TFile* pufile = TFile::Open("PileupHistograms_0121_69p2mb_pm4p6.root","READ");
  //choose central, up, or down
  TH1* puhist = (TH1*)pufile->Get("pu_weights_down");
  //----------- btags SFs-----------------
  bool applybTagSFs=1;
  int fListIndxOld=-1;
  double prob0=-100,prob1=-100;
  vector<TString> inFileName;
  TString sampleName;
  string str1;
  ifstream runListFile(inputFileList);
  TFile *currFile;
  while (std::getline(runListFile, str1)) {
    inFileName.push_back(str1);
  }runListFile.close();
  cout<<"applying b-tag SFs? "<<applybTagSFs<<endl;
  BTagCorrector btagcorr;
  //if fastsim
  if(s_data.Contains("FastSim"))
    btagcorr.SetFastSim(true);
  //--------------------------------------
  //-----------ISR rewighting for signal ------------------
  bool isrReweight = false;
  TFile *fISR;  TH2D *h2_isrWtCorr;
  vector<double> isrwt_arr={1., 0.920, 0.821, 0.715, 0.662, 0.561, 0.511};
  vector<double> isrwtUnc_arr={0., 0.04, 0.09, 0.143, 0.169, 0.219, 0.244};//if applying unc
  if(s_data.Contains("FastSim")){
    isrReweight = 1;//set here for doing ISR weighting. Reset for not applying.
    applybTagSFs = 1;
    cout<<"Applying ISR weights to these signal samples? "<<isrReweight<<endl;
    cout<<"applying b-tag SFs for signal?(ignore last message on bTag SFs) "<<applybTagSFs<<endl;
    if(s_data.Contains("T5bbbb")) fISR = new TFile("T5bbbbZg_MassScan.root");
    else if(s_data.Contains("T5tttt")) fISR = new TFile("T5ttttZg_MassScan.root");
    else if(s_data.Contains("T5qqqq")) fISR = new TFile("T5qqqqHg_MassScan.root");
    else if(s_data.Contains("T6tt")) fISR = new TFile("T6ttZg_MassScan.root");
    h2_isrWtCorr = (TH2D*)fISR->Get("Nevts_NoISRWt");
    h2_isrWtCorr->Divide((TH2D*)fISR->Get("Nevts_ISRWt"));
  }
  //-------------------------------------------------------
  cout<<"************* Applying weights for AB? "<<do_AB_reweighting<<endl;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    // ==============print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int (progress);
    if (k > decade)
      cout << 10 * k << " %" <<endl;
    decade = k;
    // cout<<"j:"<<jentry<<" fcurrent:"<<fCurrent<<endl;
    // ===============read this entry == == == == == == == == == == == 
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    //    wt=Weight*1000.0*lumiInfb;
    if(s_data.Contains("FastSim")) wt=Weight*1000.0*lumiInfb;//no PU reweighting for signal
    else  wt=Weight*1000.0*lumiInfb*(puhist->GetBinContent(puhist->GetXaxis()->FindBin(min(TrueNumInteractions,puhist->GetBinLowEdge(puhist->GetNbinsX()+1)))));

    if( !(s_data.Contains("FastSim")) ){
      if(!(CSCTightHaloFilter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && BadPFMuonFilter && NVtx > 0)) continue;  
    }
    
    if(fListIndxOld!=fCurrent){
      fListIndxOld = fCurrent;
      sampleName = inFileName[fCurrent];
      //----------- btags SFs-----------------                                                     
      if(applybTagSFs){
        currFile = TFile::Open(sampleName);
        btagcorr.SetEffs(currFile);
        btagcorr.SetCalib("btag/CSVv2_Moriond17_B_H_mod.csv");
	if(s_data.Contains("FastSim"))btagcorr.SetCalibFastSim("btag/fastsim_csvv2_ttbar_26_1_2017.csv");
      }
    }
    vector<double> prob;
    if(applybTagSFs){
      //get prob with SF up
      // btagcorr.SetBtagSFunc(1);
      // btagcorr.SetMistagSFunc(1);
      // btagcorr.SetBtagCFunc(1);
      // btagcorr.SetMistagCFunc(1);     
      prob = btagcorr.GetCorrections(Jets,Jets_hadronFlavor,Jets_HTMask);
      prob0 = prob[0]; prob1 = prob[1]+prob[2]+prob[3];
      //double corr = btagcorr.GetSimpleCorrection(Jets,Jets_hadronFlavor,Jets_HTMask,Jets_bDiscriminatorCSV);                                                 
    }
    
    //-----------ISR reweighting for signal ------------------
    double isrWt = 0,isrWtCorr = 0.;
    if(s_data.Contains("FastSim") && isrReweight){
      isrWtCorr = h2_isrWtCorr->GetBinContent(h2_isrWtCorr->GetXaxis()->FindBin(SusyMotherMass),h2_isrWtCorr->GetYaxis()->FindBin(SusyLSPMass));
      if(NJetsISR>=6) isrWt = isrwt_arr[6];//+isrwtUnc_arr[6];//if applying unc, +isrwtUnc_arr[6];
      else isrWt = isrwt_arr[NJetsISR];//+isrwtUnc_arr[6];//if applying unc, +isrwtUnc_arr[NJetsISR];
      wt = wt*isrWt*isrWtCorr;
    }
    //--------------------------------------------------------

    if( (Electrons->size() !=0) || (Muons->size() !=0) ) continue;   
    if(isoElectronTracks!=0 || isoMuonTracks!=0 || isoPionTracks!=0) continue;
    
    TLorentzVector bestPhoton=getBestPhoton();
    bool eMatchedG=check_eMatchedtoGamma();//this may not be necessary since e veto is there.
    if(bestPhoton.Pt()<0.1) continue;
    char regType='n';
    bool process=true;
    bool noFakeJet = true;
    //    bool bestPromptPho=true;
    if(s_data=="QCD"){
      if(jentry==0){cout<<"**********processing "<<s_data<<endl;}
      //      if(madMinPhotonDeltaR > 0.4 && madMinDeltaRStatus==1 && !((*Photons_nonPrompt)[bestPhotonIndxAmongPhotons]) ) continue;
      if(madMinPhotonDeltaR > 0.4 && madMinDeltaRStatus==1 && hasGenPromptPhoton) continue;
    }//QCD
    else if(s_data=="GJets"){
      if(jentry==0){cout<<"**********processing "<<s_data<<endl;}
      //      if( ((*Photons_nonPrompt)[bestPhotonIndxAmongPhotons]) ) continue;
      if(!hasGenPromptPhoton) continue;
    }//GJets
    else if(s_data=="DiPhoton"){
      if(jentry==0){cout<<"**********processing "<<s_data<<endl;}
      //      if( ((*Photons_nonPrompt)[bestPhotonIndxAmongPhotons]) ) continue;
      if(allBestPhotons.size() < 2) continue;
    }//GJets
    else if(s_data.Contains("FastSim")){
      //      //reject events with any jet pt>20, |eta|<2.5 NOT matched to a GenJet (w/in DeltaR<0.3) and chfrac < 0.1
      for(unsigned j = 0; j < Jets->size(); ++j){
      	if(Jets->at(j).Pt() <= 20 || fabs(Jets->at(j).Eta())>=2.5) continue;
      	bool genMatched = false;
      	for(unsigned g = 0; g < GenJets->size(); ++g){
      	  if(GenJets->at(g).DeltaR(Jets->at(j)) < 0.3) {
      	    genMatched = true;
      	    break;
      	  }
      	}
      	if(!genMatched && Jets_chargedHadronEnergyFraction->at(j) < 0.1){
      	  noFakeJet = false;
      	  break;
      	}
      }
    }

    bool hadJetID=true;
    int minDRindx=-100,photonMatchingJetIndx=-100,nHadJets=0;
    double minDR=99999,ST=0,myHT=0,mtPho=0;
    vector<TLorentzVector> hadJets;

    for(int i=0;i<Jets->size();i++){
      if( ((*Jets)[i].Pt() > MHT_PtCut) && (abs((*Jets)[i].Eta()) <= HT_EtaCut) ){
	double dR=bestPhoton.DeltaR((*Jets)[i]);
	if(dR<minDR){minDR=dR;minDRindx=i;}
      }
    }
    for(int i=0;i<Jets->size();i++){
      if( ((*Jets)[i].Pt() > MHT_PtCut) && (abs((*Jets)[i].Eta()) <= HT_EtaCut) ){
	if( !(minDR < 0.3 && i==minDRindx) ){
	  hadJets.push_back((*Jets)[i]);
	  if(hadJetID) hadJetID=(*Jets_ID)[i];
	}
      }
    }
    if( minDR<0.3 ) photonMatchingJetIndx=minDRindx;
    for(int i=0;i<hadJets.size();i++){
      if( (abs(hadJets[i].Eta()) < HT_EtaCut) ){ST=ST+(hadJets[i].Pt());}
      if( (abs(hadJets[i].Eta()) < Njets_EtaCut) ){nHadJets++;}
    }
    if( minDR<0.3 ){
      myHT=ST;
      ST=ST+bestPhoton.Pt();
    }
    else myHT=ST;
    sortTLorVec(&hadJets);

    if(s_data.Contains("FastSim")){
      hadJetID = true;
      if(!noFakeJet) continue;
      if(s_data.Contains("GenMET")){
	MET = GenMET;
	METPhi = GenMETPhi;
      }
    }
    double dphi1=3.8,dphi2=3.8,dphi3=3.8,dphi4=3.8,dphiG_MET=3.8;
    if(bestPhoton.Pt()>0.1) dphiG_MET=abs(DeltaPhi(METPhi,bestPhoton.Phi()));
    if(hadJets.size() > 0 ) dphi1 = abs(DeltaPhi(METPhi,(hadJets)[0].Phi()));
    if(hadJets.size() > 1 ) dphi2 = abs(DeltaPhi(METPhi,(hadJets)[1].Phi()));
    if(hadJets.size() > 2 ) dphi3 = abs(DeltaPhi(METPhi,(hadJets)[2].Phi()));
    if(hadJets.size() > 3 ) dphi4 = abs(DeltaPhi(METPhi,(hadJets)[3].Phi()));
    mtPho = sqrt(2*bestPhoton.Pt()*MET*(1-cos(DeltaPhi(METPhi,bestPhoton.Phi()))));
    if     (MET <  200 && !(dphi1 > 0.3 && dphi2 > 0.3)) regType = 'A';
    else if(MET >= 200 && !(dphi1 > 0.3 && dphi2 > 0.3)) regType = 'B';
    else if(MET <  200 &&  (dphi1 > 0.3 && dphi2 > 0.3)) regType = 'C';
    else if(MET >= 200 &&  (dphi1 > 0.3 && dphi2 > 0.3)) regType = 'D';
    else cout<<"AHHHH:Cannot assign region!!!"<<endl;

    double dphi1_genMET=3.8,dphi2_genMET=3.8;
    if(hadJets.size() > 0 ) dphi1_genMET = abs(DeltaPhi(GenMETPhi,(hadJets)[0].Phi()));
    if(hadJets.size() > 1 ) dphi2_genMET = abs(DeltaPhi(GenMETPhi,(hadJets)[1].Phi()));
    
    if(photonMatchingJetIndx>=0 && ((*Jets)[photonMatchingJetIndx].Pt())/(bestPhoton.Pt()) < 1.0) continue;
    if(photonMatchingJetIndx<0) continue;
    // int sBin4 = getBinNoV4(nHadJets),  sBin7 = getBinNoV7(nHadJets);
    // if(sBin7!=4) continue;
    if( !((ST>800 && bestPhoton.Pt()>100) || (bestPhoton.Pt()>190)) ) continue;
    if( process && !eMatchedG && bestPhoton.Pt()>=100 && (Electrons->size()==0) && (Muons->size()==0) && ST>500 && nHadJets>=2 && hadJetID && GenMET > 100 && dphi1_genMET > 0.3 && dphi2_genMET > 0.3) h_GenMETvBin_CD->Fill(GenMET,wt);

    process = process && !eMatchedG && bestPhoton.Pt()>=100 && (Electrons->size()==0) && (Muons->size()==0) && ST>500 && nHadJets>=2 && MET > 100;
    //    print(jentry);
    
    if(process){
      evtSurvived++;
      h_RunNum->Fill(RunNum);
      h_intLumi->Fill(lumiInfb,wt);
      //------------------------ Sbins----------------------------
      int searchRegion=0,sBin1=-100,m_i1=0;
      if     (nHadJets >= 2 && nHadJets <= 4 ){ searchRegion=1; sBin1 = 0;}
      else if(nHadJets == 5 || nHadJets == 6 ){ searchRegion=2; sBin1 = 7;}
      else if(nHadJets >= 7                  ){ searchRegion=3; sBin1 = 14;}
      if( searchRegion > 0 && searchRegion < 4){
	for(int i=0;i<METBinLowEdge2.size()-1;i++){
	  if(METBinLowEdge2[i]<99.99) continue;
	  m_i1++;
	  if(MET>=METBinLowEdge2[i] && MET<METBinLowEdge2[i+1]){ sBin1 = sBin1+m_i1;break; }
	  else if(MET>=METBinLowEdge2[METBinLowEdge2.size()-1]){ sBin1 = sBin1+7   ;break; }
	}
      }
      else cout<<"Event outside search region! ";
      int sBin4 = getBinNoV4(nHadJets),  sBin7 = getBinNoV7(nHadJets);
      //------------------------ Sbins----------------------------
      if(do_AB_reweighting && (regType=='A' || regType=='B')){
	double parX=MET;
	double parY=nHadJets;
	if(BTags==0) h2_HLRatio=(TH2D*)f_HLR->Get("HLRatio_0");
	else h2_HLRatio=(TH2D*)f_HLR->Get("HLRatio_1");
	//	wt=wt*(h2_HLRatio->GetBinContent(h2_HLRatio->GetXaxis()->FindBin(parX),h2_HLRatio->GetYaxis()->FindBin(parY)));
	double hlr1=h2_HLRatio->GetBinContent(h2_HLRatio->GetXaxis()->FindBin(101),h2_HLRatio->GetYaxis()->FindBin(parY));
	double doubleratio = (h2_HLRatio->GetBinContent(h2_HLRatio->GetXaxis()->FindBin(205),h2_HLRatio->GetYaxis()->FindBin(parY)))/hlr1;
	if(regType=='B') wt=wt*hlr1*doubleratio;
	else wt=wt*hlr1;
	//	MET=MET-7;
      }
      if(hadJetID){
	int r_i=-1;
	if(regType=='A') r_i=1;
	else if(regType=='B') r_i=2;
	else if(regType=='C') r_i=3;
	else if(regType=='D') r_i=4;
	r_i=r_i-1;
	
	h_ST_[r_i]->Fill(ST,wt);
	h_MET_[r_i]->Fill(MET,wt);
	h_nHadJets_[r_i]->Fill(nHadJets,wt);
	if(applybTagSFs){
	  h_BTags_[r_i]->Fill(0.0,wt*prob[0]);
	  h_BTags_[r_i]->Fill(1.0,wt*prob[1]);
	  h_BTags_[r_i]->Fill(2.0,wt*prob[2]);
	  h_BTags_[r_i]->Fill(3.0,wt*prob[3]);
	}
	else{
	  h_BTags_[r_i]->Fill(BTags,wt);
	}
	h_nVtx_[r_i]->Fill(NVtx,wt);
	h_GenMET_[r_i]->Fill(GenMET,wt);
	h_METPhi_[r_i]->Fill(METPhi,wt);
	h_myHT_[r_i]->Fill(myHT,wt);
	h_SBins_v1_[r_i]->Fill(sBin1,wt);
	h_SBins_v4_[r_i]->Fill(sBin4,wt);

	//---------------- for ISR syst -----------------
	if(s_data.Contains("FastSim")){
	  h_nHadJets_SBin_v7_[r_i]->Fill(nHadJets,wt);
	  h_nHadJets_NoISRWt_SBin_v7_[r_i]->Fill(nHadJets,wt/(isrWt*isrWtCorr));
	  h_nHadJets_ISRUncSq_SBin_v7_[r_i]->Fill(nHadJets,wt*(abs(1-isrWt)/2)*(abs(1-isrWt)/2));
	}
	//-----------------------------------------------

	if(!applybTagSFs)
	  h_SBins_v7_[r_i]->Fill(sBin7,wt);
	else{
	  int nBtagsOrg = BTags;
	  BTags = 0;
	  int binNum0 = getBinNoV7(nHadJets);
	  BTags = 1;
	  int binNum1 = getBinNoV7(nHadJets);
	  h_SBins_v7_[r_i]->Fill(binNum0,wt*prob0);
	  h_SBins_v7_[r_i]->Fill(binNum1,wt*prob1);
	  // h_SBins_v7_bTagSFup_[r_i]->Fill(binNum0,wt*probUp[0]);
	  // h_SBins_v7_bTagSFup_[r_i]->Fill(binNum1,wt*(1-probUp[0]));
	  BTags = nBtagsOrg;
	}
	h_BestPhotonPt_[r_i]->Fill( bestPhoton.Pt(),wt );	
	h_BestPhotonPtvBin_[r_i]->Fill(bestPhoton.Pt(),wt);
	h_BestPhotonEta_[r_i]->Fill(bestPhoton.Eta(),wt);
	h_BestPhotonPhi_[r_i]->Fill(bestPhoton.Phi(),wt);
	h_mTPho_[r_i]->Fill(mtPho,wt);
	h_dPhi_METBestPhoton_[r_i]->Fill(dphiG_MET,wt);
	h_dPhiPhotonJet1_[r_i]->Fill(bestPhoton.DeltaPhi(hadJets[0]),wt);
	if(photonMatchingJetIndx>=0) h_RatioJetPhoPt_[r_i]->Fill( ((*Jets)[photonMatchingJetIndx].Pt())/(bestPhoton.Pt()),wt);

	if(nHadJets<=4) h_PhoPt_nJ2to4_[r_i]->Fill(bestPhoton.Pt(),wt );
	else h_PhoPt_minNJ5_[r_i]->Fill(bestPhoton.Pt(),wt );

	h2_dPhi1dPhi2_[r_i]->Fill(dphi1,dphi2,wt);
	h2_PtPhotonvsMET_[r_i]->Fill( bestPhoton.Pt(),MET,wt);
	h2_NJST_[r_i]->Fill(nHadJets,ST,wt);
	h2_dPhiMETPho_NJ_[r_i]->Fill(dphiG_MET,nHadJets,wt);
	if(photonMatchingJetIndx>=0) h2_RatioJetPhoPtVsPhoPt_[r_i]->Fill(bestPhoton.Pt(),((*Jets)[photonMatchingJetIndx].Pt())/(bestPhoton.Pt()),wt);
	//	h2_RatioJetPhoPtVsPhoPt_[r_i]->Fill(bestPhoton.Pt(),MHT/MET,wt);
	if(photonMatchingJetIndx>=0) h2_RatioJetPhoPtVsdPhiG_[r_i]->Fill(dphiG_MET,      ((*Jets)[photonMatchingJetIndx].Pt())/(bestPhoton.Pt()),wt);
	if(photonMatchingJetIndx>=0) h2_RatioJetPhoPtVsMET_[r_i]->Fill(MET,            ((*Jets)[photonMatchingJetIndx].Pt())/(bestPhoton.Pt()),wt);
	h_mindPhi1dPhi2_[r_i]->Fill(min(dphi1,dphi2),wt);
	h_mindPhi1dPhi2_ABCD->Fill(min(dphi1,dphi2),wt);

	h_HT_[r_i]->Fill(HT,wt);
	h_MHT_[r_i]->Fill(MHT,wt);
	h_nJets_[r_i]->Fill(NJets,wt);

	h_STvBin_[r_i]->Fill(ST,wt);
	h_METvBin_[r_i]->Fill(MET,wt);
	h_jet1Pt_[r_i]->Fill(hadJets[0].Pt(),wt);
	h_jet2Pt_[r_i]->Fill(hadJets[1].Pt(),wt);

	if(dphi1 < dphi2) h_PtjetNearMET_[r_i]->Fill(hadJets[0].Pt(),wt);
	else h_PtjetNearMET_[r_i]->Fill(hadJets[1].Pt(),wt);

	if(hadJets.size() > 0 ) h_dPhi_METjet1_[r_i]->Fill( dphi1,wt );
	if(hadJets.size() > 1 ) h_dPhi_METjet2_[r_i]->Fill( dphi2,wt );
	if(hadJets.size() > 2 ) h_dPhi_METjet3_[r_i]->Fill( dphi3,wt );
	if(hadJets.size() > 3 ) h_dPhi_METjet4_[r_i]->Fill( dphi4,wt );
	if(regType=='A' || regType=='B'){
	  h_ST_AB->Fill(ST,wt);
	  h_MET_AB->Fill(MET,wt);
	  h_nHadJets_AB->Fill(nHadJets,wt);
	  if(applybTagSFs){
	    h_BTags_AB->Fill(0.0,wt*prob[0]);
	    h_BTags_AB->Fill(1.0,wt*prob[1]);
	    h_BTags_AB->Fill(2.0,wt*prob[2]);
	    h_BTags_AB->Fill(3.0,wt*prob[3]);
	  }
	  else{
	    h_BTags_AB->Fill(BTags,wt);
	  }
	  h_nVtx_AB->Fill(NVtx,wt);
	  h_GenMET_AB->Fill(GenMET,wt);
	  h_METPhi_AB->Fill(METPhi,wt);
	  h_myHT_AB->Fill(myHT,wt);
	  h_SBins_v1_AB->Fill(sBin1,wt);
	  h_SBins_v4_AB->Fill(sBin4,wt);
	  if(!applybTagSFs)
	    h_SBins_v7_AB->Fill(sBin7,wt);
	  else{
	    int nBtagsOrg = BTags;
	    BTags = 0;
	    int binNum0 = getBinNoV7(nHadJets);
	    BTags = 1;
	    int binNum1 = getBinNoV7(nHadJets);
	    h_SBins_v7_AB->Fill(binNum0,wt*prob0);
	    h_SBins_v7_AB->Fill(binNum1,wt*prob1);
	    // h_SBins_v7_bTagSFup_AB->Fill(binNum0,wt*probUp[0]);
	    // h_SBins_v7_bTagSFup_AB->Fill(binNum1,wt*(1-probUp[0]));
	    BTags = nBtagsOrg;
	  }
	  h_BestPhotonPt_AB->Fill( bestPhoton.Pt(),wt );	
	  h_BestPhotonPtvBin_AB->Fill(bestPhoton.Pt(),wt);
	  h_BestPhotonEta_AB->Fill(bestPhoton.Eta(),wt);
	  h_BestPhotonPhi_AB->Fill(bestPhoton.Phi(),wt);
	  h_mTPho_AB->Fill(mtPho,wt);
	  h_dPhi_METBestPhoton_AB->Fill(dphiG_MET,wt);
	  h_dPhiPhotonJet1_AB->Fill(bestPhoton.DeltaPhi(hadJets[0]),wt);
	  if(photonMatchingJetIndx>=0) h_RatioJetPhoPt_AB->Fill( ((*Jets)[photonMatchingJetIndx].Pt())/(bestPhoton.Pt()),wt);
	  h_mindPhi1dPhi2_AB->Fill(min(dphi1,dphi2),wt);

	  if(nHadJets<=4) h_PhoPt_nJ2to4_AB->Fill(bestPhoton.Pt(),wt );
	  else h_PhoPt_minNJ5_AB->Fill(bestPhoton.Pt(),wt );

	  h2_PtPhotonvsMET_AB->Fill( bestPhoton.Pt(),MET,wt);
	  h2_dPhi1dPhi2_AB->Fill(dphi1,dphi2,wt);
	  h2_NJST_AB->Fill(nHadJets,ST,wt);
	  h2_dPhiMETPho_NJ_AB->Fill(dphiG_MET,nHadJets,wt);
	  if(photonMatchingJetIndx>=0) h2_RatioJetPhoPtVsPhoPt_AB->Fill(bestPhoton.Pt(),((*Jets)[photonMatchingJetIndx].Pt())/(bestPhoton.Pt()),wt);
	  //h2_RatioJetPhoPtVsPhoPt_AB->Fill(bestPhoton.Pt(),MHT/MET,wt);
	  if(photonMatchingJetIndx>=0) h2_RatioJetPhoPtVsdPhiG_AB->Fill(dphiG_MET,      ((*Jets)[photonMatchingJetIndx].Pt())/(bestPhoton.Pt()),wt);
	  if(photonMatchingJetIndx>=0) h2_RatioJetPhoPtVsMET_AB->Fill(MET,              ((*Jets)[photonMatchingJetIndx].Pt())/(bestPhoton.Pt()),wt);

	  if(!applybTagSFs){
	    if(BTags==0) h2_METnHadJ_0b_AB->Fill(MET,nHadJets,wt);
	    else h2_METnHadJ_min1b_AB->Fill(MET,nHadJets,wt);
	  }
	  else{
	    h2_METnHadJ_0b_AB->Fill(MET,nHadJets,wt*prob0);
	    h2_METnHadJ_min1b_AB->Fill(MET,nHadJets,wt*prob1);
	  }

	  h_HT_AB->Fill(HT,wt);
	  h_MHT_AB->Fill(MHT,wt);
	  h_nJets_AB->Fill(NJets,wt);

	  h_STvBin_AB->Fill(ST,wt);
	  h_METvBin_AB->Fill(MET,wt);

	  h_jet1Pt_AB->Fill(hadJets[0].Pt(),wt);
	  h_jet2Pt_AB->Fill(hadJets[1].Pt(),wt);
	  if(dphi1 < dphi2) h_PtjetNearMET_AB->Fill(hadJets[0].Pt(),wt);
	  else h_PtjetNearMET_AB->Fill(hadJets[1].Pt(),wt);
  
	  if(hadJets.size() > 0 ) h_dPhi_METjet1_AB->Fill( dphi1,wt );
	  if(hadJets.size() > 1 ) h_dPhi_METjet2_AB->Fill( dphi2,wt );
	  if(hadJets.size() > 2 ) h_dPhi_METjet3_AB->Fill( dphi3,wt );
	  if(hadJets.size() > 3 ) h_dPhi_METjet4_AB->Fill( dphi4,wt );
	}
	if(regType=='C' || regType=='D'){
	  //	  print(jentry);
	  h_ST_CD->Fill(ST,wt);
	  h_MET_CD->Fill(MET,wt);
	  h_nHadJets_CD->Fill(nHadJets,wt);
	  if(applybTagSFs){
	    h_BTags_CD->Fill(0.0,wt*prob[0]);
	    h_BTags_CD->Fill(1.0,wt*prob[1]);
	    h_BTags_CD->Fill(2.0,wt*prob[2]);
	    h_BTags_CD->Fill(3.0,wt*prob[3]);
	  }
	  else{
	    h_BTags_CD->Fill(BTags,wt);
	  }
	  h_nVtx_CD->Fill(NVtx,wt);
	  h_GenMET_CD->Fill(GenMET,wt);
	  h_METPhi_CD->Fill(METPhi,wt);
	  h_myHT_CD->Fill(myHT,wt);
	  h_SBins_v1_CD->Fill(sBin1,wt);
	  h_SBins_v4_CD->Fill(sBin4,wt);
	  if(!applybTagSFs)
	    h_SBins_v7_CD->Fill(sBin7,wt);
	  else{
	    int nBtagsOrg = BTags;
	    BTags = 0;
	    int binNum0 = getBinNoV7(nHadJets);
	    BTags = 1;
	    int binNum1 = getBinNoV7(nHadJets);
	    h_SBins_v7_CD->Fill(binNum0,wt*prob0);
	    h_SBins_v7_CD->Fill(binNum1,wt*prob1);
	    // h_SBins_v7_bTagSFup_CD->Fill(binNum0,wt*probUp[0]);
	    // h_SBins_v7_bTagSFup_CD->Fill(binNum1,wt*(1-probUp[0]));
	    BTags = nBtagsOrg;
	  }
	  //---------------- for ISR syst -----------------
	  if(s_data.Contains("FastSim")){
	    h_SBins_v7_ISRwtNoBtagSF->Fill(sBin7,wt);
	    h_SBins_v7_ISRUncSq_CD->Fill(sBin7,wt*(abs(1-isrWt)/2)*(abs(1-isrWt)/2));
	    h_SBins_v7_NoISRWt_CD->Fill(sBin7,wt/(isrWt*isrWtCorr));
	    h_SBins_v7_ISRUncSqNoISRwt_CD->Fill(sBin7,wt*((abs(1-isrWt)/2)*(abs(1-isrWt)/2)) / (isrWt*isrWtCorr));
	  }
	  //-----------------------------------------------
	  h_BestPhotonPt_CD->Fill( bestPhoton.Pt(),wt );	
	  h_BestPhotonPtvBin_CD->Fill(bestPhoton.Pt(),wt);
	  h_BestPhotonEta_CD->Fill(bestPhoton.Eta(),wt);
	  h_BestPhotonPhi_CD->Fill(bestPhoton.Phi(),wt);
	  h_mTPho_CD->Fill(mtPho,wt);
	  h_dPhi_METBestPhoton_CD->Fill(dphiG_MET,wt);
	  h_dPhiPhotonJet1_CD->Fill(bestPhoton.DeltaPhi(hadJets[0]),wt);
	  if(photonMatchingJetIndx>=0) h_RatioJetPhoPt_CD->Fill( ((*Jets)[photonMatchingJetIndx].Pt())/(bestPhoton.Pt()),wt);
	  h_mindPhi1dPhi2_CD->Fill(min(dphi1,dphi2),wt);

	  if(nHadJets<=4) h_PhoPt_nJ2to4_CD->Fill(bestPhoton.Pt(),wt );
	  else h_PhoPt_minNJ5_CD->Fill(bestPhoton.Pt(),wt );

	  h2_PtPhotonvsMET_CD->Fill( bestPhoton.Pt(),MET,wt);
	  h2_dPhi1dPhi2_CD->Fill(dphi1,dphi2,wt);
	  h2_NJST_CD->Fill(nHadJets,ST,wt);
	  h2_dPhiMETPho_NJ_CD->Fill(dphiG_MET,nHadJets,wt);
	  if(photonMatchingJetIndx>=0) h2_RatioJetPhoPtVsPhoPt_CD->Fill(bestPhoton.Pt(),((*Jets)[photonMatchingJetIndx].Pt())/(bestPhoton.Pt()),wt);
	  //h2_RatioJetPhoPtVsPhoPt_CD->Fill(bestPhoton.Pt(),MHT/MET,wt);
	  if(photonMatchingJetIndx>=0) h2_RatioJetPhoPtVsdPhiG_CD->Fill(dphiG_MET,      ((*Jets)[photonMatchingJetIndx].Pt())/(bestPhoton.Pt()),wt);
	  if(photonMatchingJetIndx>=0) h2_RatioJetPhoPtVsMET_CD->Fill(MET,              ((*Jets)[photonMatchingJetIndx].Pt())/(bestPhoton.Pt()),wt);

	  if(!applybTagSFs){
	    if(BTags==0) h2_METnHadJ_0b_CD->Fill(MET,nHadJets,wt);
	    else h2_METnHadJ_min1b_CD->Fill(MET,nHadJets,wt);
	  }
	  else{
	    h2_METnHadJ_0b_CD->Fill(MET,nHadJets,wt*prob0);
	    h2_METnHadJ_min1b_CD->Fill(MET,nHadJets,wt*prob1);
	  }
	  h_HT_CD->Fill(HT,wt);
	  h_MHT_CD->Fill(MHT,wt);
	  h_nJets_CD->Fill(NJets,wt);

	  h_STvBin_CD->Fill(ST,wt);
	  h_METvBin_CD->Fill(MET,wt);

	  h_jet1Pt_CD->Fill(hadJets[0].Pt(),wt);
	  h_jet2Pt_CD->Fill(hadJets[1].Pt(),wt);

	  if(dphi1 < dphi2) h_PtjetNearMET_CD->Fill(hadJets[0].Pt(),wt);
	  else h_PtjetNearMET_CD->Fill(hadJets[1].Pt(),wt);

	  if(hadJets.size() > 0 ) h_dPhi_METjet1_CD->Fill( dphi1,wt );
	  if(hadJets.size() > 1 ) h_dPhi_METjet2_CD->Fill( dphi2,wt );
	  if(hadJets.size() > 2 ) h_dPhi_METjet3_CD->Fill( dphi3,wt );
	  if(hadJets.size() > 3 ) h_dPhi_METjet4_CD->Fill( dphi4,wt );
	}
	findObjMatchedtoG(bestPhoton);
	h_nBestPho->Fill(allBestPhotons.size(),wt);
	if(allBestPhotons.size()>1){
	  sortTLorVec(&allBestPhotons);
	  h_BestPho2Pt->Fill(allBestPhotons[1].Pt(),wt);
	  double dR=100;
	  int match=-100;
	  for(int i=0;i<GenParticles->size();i++){
	    if((*GenParticles)[i].Pt()!=0){
	      if(dR > (allBestPhotons[1].DeltaR((*GenParticles)[i])) && ( (allBestPhotons[1].Pt()/(*GenParticles)[match].Pt()> 0.5) || (allBestPhotons[1].Pt()/(*GenParticles)[match].Pt()< 1.5) ) ){
		dR=allBestPhotons[1].DeltaR((*GenParticles)[i]);
		match=i;
	      }
	    }
	  }
	  if(dR<0.2 && match>0)  h_2ndPhoGenMatch->Fill(abs((*GenParticles_PdgId)[match]),wt);
	  else h_2ndPhoGenMatch->Fill(0.0,wt);
	}
	h_madHT->Fill(madHT,wt);
	for(int i=0;i<hadJets.size() && i<4;i++){
	  h_jetPt[i] ->Fill(hadJets[i].Pt(),wt);
	  h_jetEta[i]->Fill(hadJets[i].Eta(),wt);
	  h_jetPhi[i]->Fill(hadJets[i].Phi(),wt);
	}
	for(int i=0;i<Jets->size();i++){
	  if((*Jets)[i].Pt()>30 && (*Jets_ID)[i] ){
	    h_NM->Fill((*Jets_neutralMultiplicity)[i],wt);
	    if( abs((*Jets)[i].Eta())<3){h_NMEta3->Fill((*Jets_neutralMultiplicity)[i],wt);}
	    if( abs((*Jets)[i].Eta())>=3){h_NMHF->Fill((*Jets_neutralMultiplicity)[i],wt);}
	    h_CM->Fill((*Jets_chargedMultiplicity)[i],wt);
	    h_NHF->Fill((*Jets_neutralHadronEnergyFraction)[i],wt);
	    h_NEMF->Fill((*Jets_neutralEmEnergyFraction)[i],wt);
	    h_CHF->Fill((*Jets_chargedHadronEnergyFraction)[i],wt);
	    h_CEMF->Fill((*Jets_chargedEmEnergyFraction)[i],wt);

	    h_PhoMult->Fill((*Jets_photonMultiplicity)[i],wt);
	    h_PhoFrac->Fill((*Jets_photonEnergyFraction)[i],wt);
	  }			
	}	
      }
    }//process
  }// loop over entries
  cout<<"Events Survied:"<<evtSurvived<<endl;
}


TLorentzVector MultiJet::getBestPhoton(){
  vector<TLorentzVector> goodPho;
  vector<int> goodPhoIndx;
  allBestPhotons.resize(0);
  for(int iPho=0;iPho<Photons->size();iPho++){
    if( (*Photons_fullID)[iPho] && ((*Photons_hasPixelSeed)[iPho]<0.001) ) {
      goodPho.push_back( (*Photons)[iPho] );
      goodPhoIndx.push_back(iPho);
      allBestPhotons.push_back((*Photons)[iPho]);
    }
  }
  
  int highPtIndx=-100;
  for(int i=0;i<goodPho.size();i++){
    if(i==0) highPtIndx=0;
    else if( (goodPho[highPtIndx].Pt()) < (goodPho[i].Pt()) ){highPtIndx=i;}
  }

  if(highPtIndx>=0){
    bestPhotonIndxAmongPhotons = goodPhoIndx[highPtIndx];
  }
  else bestPhotonIndxAmongPhotons = -100;
  if(highPtIndx==-100){TLorentzVector v0;return v0;}
  else return goodPho[highPtIndx];  
}

int MultiJet::getBinNoV4(int nHadJets){
  int sBin=-100,m_i=0;
  if(BTags==0){
    if(nHadJets>=2 && nHadJets<=4)     { sBin=0;}
    else if(nHadJets==5 || nHadJets==6){ sBin=8;}
    else if(nHadJets>=7)               { sBin=15;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)     { sBin=22;}
    else if(nHadJets==5 || nHadJets==6){ sBin=29;}
    else if(nHadJets>=7)               { sBin=36;}
  }
  if(sBin==0){
    for(int i=0;i<METBinLowEdgeV4_njLow.size()-1;i++){
      if(METBinLowEdgeV4_njLow[i]<99.99) continue;
      m_i++;
      if(MET >= METBinLowEdgeV4_njLow[i] && MET < METBinLowEdgeV4_njLow[i+1]){ sBin = sBin+m_i;break; }
      else if(MET >= METBinLowEdgeV4_njLow[METBinLowEdgeV4_njLow.size()-1])  { sBin = 8         ;break; }
    }
  }
  else{
    for(int i=0;i<METBinLowEdgeV4.size()-1;i++){
      if(METBinLowEdgeV4[i]<99.99) continue;
      m_i++;
      if(MET >= METBinLowEdgeV4[i] && MET < METBinLowEdgeV4[i+1]){ sBin = sBin+m_i;break; }
      else if(MET >= METBinLowEdgeV4[METBinLowEdgeV4.size()-1])  { sBin = sBin+7   ;break; }
    }
  }
  return sBin;
}

int MultiJet::getBinNoV7(int nHadJets){
  int sBin=-100,m_i=0;
  if(BTags==0){
    if(nHadJets>=2 && nHadJets<=4)     { sBin=0;}
    else if(nHadJets==5 || nHadJets==6){ sBin=6;}
    else if(nHadJets>=7)               { sBin=11;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)     { sBin=16;}
    else if(nHadJets==5 || nHadJets==6){ sBin=21;}
    else if(nHadJets>=7)               { sBin=26;}
  }
  if(sBin==0){
    for(int i=0;i<METBinLowEdgeV7_njLow.size()-1;i++){
      if(METBinLowEdgeV7_njLow[i]<99.99) continue;
      m_i++;
      if(MET >= METBinLowEdgeV7_njLow[i] && MET < METBinLowEdgeV7_njLow[i+1]){ sBin = sBin+m_i;break; }
      else if(MET >= METBinLowEdgeV7_njLow[METBinLowEdgeV7_njLow.size()-1])  { sBin = 6         ;break; }
    }
  }
  else{
    for(int i=0;i<METBinLowEdgeV7.size()-1;i++){
      if(METBinLowEdgeV7[i]<99.99) continue;
      m_i++;
      if(MET >= METBinLowEdgeV7[i] && MET < METBinLowEdgeV7[i+1]){ sBin = sBin+m_i;break; }
      else if(MET >= METBinLowEdgeV7[METBinLowEdgeV7.size()-1])  { sBin = sBin+5   ;break; }
    }
  }
  return sBin;
}

bool MultiJet::check_eMatchedtoGamma(){
  if(bestPhotonIndxAmongPhotons>=0){
    for(int i=0;i<Electrons->size();i++){
      if( (*Photons)[bestPhotonIndxAmongPhotons].DeltaR( (*Electrons)[i] ) < 0.1){
	//	cout<<(*Electrons)[i].Pt()<<" "<<(*Electrons)[i].Eta()<<" "<<(*Electrons)[i].Phi()<<" "<<(*Photons)[bestPhotonIndxAmongPhotons].Pt()<<" "<<(*Photons)[bestPhotonIndxAmongPhotons].Eta()<<" "<<(*Photons)[bestPhotonIndxAmongPhotons].Phi()<<" dR:"<<(*Photons)[bestPhotonIndxAmongPhotons].DeltaR( (*Electrons)[i])<<endl;
	return true;
      }
    }
  }
  else
    return false;
}

void  MultiJet::findObjMatchedtoG(TLorentzVector bestPhoton){
  double dR=100,mindr=1000;
  int match=-100;
  for(int i=0;i<GenParticles->size();i++){
    if((*GenParticles)[i].Pt()!=0){
      if(dR > (bestPhoton.DeltaR((*GenParticles)[i])) && ( (bestPhoton.Pt()/(*GenParticles)[i].Pt() > 0.5) && (bestPhoton.Pt()/(*GenParticles)[i].Pt()<1.5)) ){
	dR=bestPhoton.DeltaR((*GenParticles)[i]);
	match=i;
      }
    }
  }
  //cout<<"Pt: "<<bestPhoton.Pt()<<" Eta:"<<bestPhoton.Eta()<<" Phi:"<<bestPhoton.Phi()<<" PdgID:"<<(*GenParticles_PdgId)[match]<<" Pt:"<<(*GenParticles)[match].Pt()<<" Eta:"<<(*GenParticles)[match].Eta()<<" Phi:"<<(*GenParticles)[match].Phi()<<" parentId:"<<(*GenParticles_ParentId)[match]<<endl;
  if(dR<0.2){
    h_GmatchedObj->Fill(abs((*GenParticles_PdgId)[match]),wt);
    if(abs((*GenParticles_PdgId)[match])==22)  h_PdgIdPhoParent->Fill(abs((*GenParticles_ParentId)[match]),wt);
  }
  else{
    h_GmatchedObj->Fill(0.0,wt);
    h_PdgIdPhoParent->Fill(0.0,wt);
  }
}

void MultiJet::print(Long64_t jentry){
  //cout<<endl;
  TLorentzVector v1,photo;
  for(int i=0;i<GenParticles->size();i++){
    cout<<EvtNum<<" "<<jentry<<" "<<GenParticles->size()<<" "<<i<<" PdgId:"<<(*GenParticles_PdgId)[i]<<" parentId:"<<(*GenParticles_ParentId)[i]<<" parentIndx:"<<(*GenParticles_ParentIdx)[i]<<" Status:"<<(*GenParticles_Status)[i]<</*"\tPx:"<<(*GenParticles)[i].Px()<<" Py:"<<(*GenParticles)[i].Py()<<" Pz:"<<(*GenParticles)[i].Pz()<<*/"\tPt:"<<(*GenParticles)[i].Pt()<<" Eta:"<<(*GenParticles)[i].Eta()<<" Phi:"<<(*GenParticles)[i].Phi()<<" E:"<<(*GenParticles)[i].Energy()<<endl;
  }

  for(int i=0;i<GenJets->size();i++){
    if(i==0) cout<<"-------------------------------- GenJets -------------------------------------------"<<endl;
    cout<<"GenJetsPt:"<<(*GenJets)[i].Pt()<<" eta:"<<(*GenJets)[i].Eta()<<" phi:"<<(*GenJets)[i].Phi()<<endl;
  }
  
  for(int i=0;i<Photons->size();i++){
    if(i==0) cout<<"-------------------------------- Photons -------------------------------------------"<<endl;
    double dR=0;//DeltaR( bestPhoton.Eta(),bestPhoton.Phi(),(*Photons)[i].Eta(),(*Photons)[i].Phi() );
    //    cout<<jentry<<" i:"<<i<<" phoSize:"<<Photons->size()<<" Pt:"<<bestPhoton.Pt()<<" eta:"<<bestPhoton.Eta()<<" phi:"<<bestPhoton.Phi()<<" otherP:"<<(*Photons)[i].Pt()<<" eta:"<<(*Photons)[i].Eta()<<" phi:"<<(*Photons)[i].Phi()<<" dR:"<<dR<<endl;
    cout<<"Photon Pt:"<<(*Photons)[i].Pt()<<" eta:"<<(*Photons)[i].Eta()<<" phi:"<<(*Photons)[i].Phi()<<" E:"<<(*Photons)[i].Energy()<<" FullID:"<<(*Photons_fullID)[i]<<" Pixel Seed:"<<(*Photons_hasPixelSeed)[i]<<endl;
  }
  
  // for(int i=0;i<TAPElectronTracks->size();i++){
  //   cout<<"Ele Track Pt:"<<(*TAPElectronTracks)[i].Pt()<<" Eta:"<<(*TAPElectronTracks)[i].Eta()<<" Phi:"<<(*TAPElectronTracks)[i].Phi()<<endl;
  // }
  // for(int i=0;i<TAPMuonTracks->size();i++){
  //   cout<<"Mu Track Pt:"<<(*TAPMuonTracks)[i].Pt()<<" Eta:"<<(*TAPMuonTracks)[i].Eta()<<" Phi:"<<(*TAPMuonTracks)[i].Phi()<<endl;
  // }
  // for(int i=0;i<TAPPionTracks->size();i++){
  //   cout<<"Pi Track Pt:"<<(*TAPPionTracks)[i].Pt()<<" Eta:"<<(*TAPPionTracks)[i].Eta()<<" Phi:"<<(*TAPPionTracks)[i].Phi()<<endl;
  // }

  for(int i=0;i<Jets->size();i++){
    //  if( ((*Jets)[i].Pt() > MHT_PtCut) && (abs((*Jets)[i].Eta()) <= MHT_EtaCut) ){
    cout<<"JetPt:"<<(*Jets)[i].Pt()<<" JetEta:"<<(*Jets)[i].Eta()<<" JetPhi:"<<(*Jets)[i].Phi()<<" energy:"<<(*Jets)[i].Energy()<<" QMult:"<<(*Jets_chargedMultiplicity)[i]<<" CSV:"<<(*Jets_bDiscriminatorCSV)[i]<<" dPhi with photon:"<<(*Jets)[i].DeltaPhi((*Photons)[bestPhotonIndxAmongPhotons])<<endl;
  }
  
  for(int i=0;i<Muons->size();i++){
    if(i==0) cout<<"-------------------------------- Muons -------------------------------------------"<<endl;
    cout<<"MuonPt: "<<(*Muons)[i].Pt()<<" Eta: "<<(*Muons)[i].Eta()<<" Phi: "<<(*Muons)[i].Phi()<<" M: "<<(*Muons)[i].M()<<" E:"<<(*Muons)[i].Energy()<<endl;
  }
  for(int i=0;i<Electrons->size();i++){
    if(i==0) cout<<"-------------------------------- Electrons -------------------------------------------"<<endl;
    cout<<"ElectronPt: "<<(*Electrons)[i].Pt()<<" Eta: "<<(*Electrons)[i].Eta()<<" Phi: "<<(*Electrons)[i].Phi()<<" M: "<<(*Electrons)[i].M()<<" E:"<<(*Electrons)[i].Energy()<<endl;
  }
  cout<<"MET:"<<MET<<" METPhi:"<<METPhi<<endl;
  cout<<"^^^^^^^^^^^^^^^^^^ Event ends ^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<endl<<endl;
}
