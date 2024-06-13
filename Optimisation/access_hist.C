// compile using g++ `root-config --cflags` access_hist.C -o access_hist.exe `root-config --libs`
#include "TChain.h"
#include <iostream>
#include "TTree.h"
#include "TFileCollection.h"
#include "TFile.h"
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "TFile.h"
#include "TTree.h"
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
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TFileCollection.h"

using namespace std;
void access_hist(const char *inputfile, const char *outputfile){

  Float_t         phoEt, phoEta, phoPhi, lepPt, lepEta, lepPhi, sigMT, sigMET, sigMETPhi, dPhiLepMET, dRPhoLep, HT, nJet, eleR9, phoR9, MiniIso,lumiWeight;
  Int_t           nVertex, nBJet;

  TFile *f;
  f = TFile::Open(inputfile);
  TTree *t = (TTree*)f->Get("signalTree"); 

  t->SetBranchAddress("phoEt", &phoEt);
  t->SetBranchAddress("phoEta", &phoEta);
  t->SetBranchAddress("phoPhi", &phoPhi);
  t->SetBranchAddress("lepPt", &lepPt);
  t->SetBranchAddress("lepEta", &lepEta);
  t->SetBranchAddress("lepPhi", &lepPhi);
  t->SetBranchAddress("sigMT", &sigMT);
  t->SetBranchAddress("sigMET", &sigMET);
  t->SetBranchAddress("sigMETPhi", &sigMETPhi);
  t->SetBranchAddress("dPhiLepMET", &dPhiLepMET);
  t->SetBranchAddress("nVertex", &nVertex);
  t->SetBranchAddress("dRPhoLep", &dRPhoLep);
  t->SetBranchAddress("HT", &HT);
  t->SetBranchAddress("nJet", &nJet);
  t->SetBranchAddress("nBJet", &nBJet);
  t->SetBranchAddress("eleR9", &eleR9);
  t->SetBranchAddress("phoR9", &phoR9);
  t->SetBranchAddress("MiniIso", &MiniIso);
  t->SetBranchAddress("lumiWeight", &lumiWeight);

  		TH1F* hPhopt = new TH1F("hPhopt","",45,0,900);
  		TH1F* hPhoeta = new TH1F("hPhoeta","",51,-1.7,1.7);
  		TH1F* hPhophi = new TH1F("hPhophi","",40,-4.0,4.0);
  		TH1F* hLeppt = new TH1F("hLeppt","",70,0,700);
  		TH1F* hLepeta = new TH1F("hLepeta","",90,-3.0,3.0);
  		TH1F* hLepphi = new TH1F("hLepphi","",60,-4.0,4.0);
  		TH1F* hsigMT = new TH1F("hsigMT","",60,0,600);
  		TH1F* hsigMET = new TH1F("hsigMET","",60,0,600);
  		TH1F* hsigMETPhi = new TH1F("hsigMETPhi","",60,-4.0,4.0);
  		TH1F* hdPhiLepMET = new TH1F("hdPhiLepMET","",60,-4.0,4.0);
  		TH1F* hnVertex = new TH1F("hnVertex","",100,0,100);
  		TH1F* hdRPhoLep = new TH1F("hdRPhoLep","",60,0,6);
  		TH1F* hHT = new TH1F("hHT","",60,0,1200);
  		TH1F* hnJet = new TH1F("hnJet","",10,0,10);
  		TH1F* hnBJet = new TH1F("hnBJet","",10,0,10);
  		TH1F* heleR9 = new TH1F("heleR9","",32,0.4,1.2);
  		TH1F* hphoR9 = new TH1F("hphoR9","",32,0.4,1.2);
  		TH1F* hMiniIso = new TH1F("hMiniIso","",20,0,0.2);
  		TH2F* hMET_vs_MT = new TH2F("hMET_vs_MT","",21,90,300,22,80,300);
//  int METcut, MTcut, option;
  
  Int_t nentries1 = (Int_t)t->GetEntries();
  cout<< nentries1 <<endl;
  bool cut;
  int METcut[12] = {90,100,110,120,130,140,150,160,170,180,190,200};
  int MTcut[13] = {80,90,100,110,120,130,140,150,160,170,180,190,200};

  
  for(int p=0; p<12 ; p++){                        
         for(int q=0; q<13 ; q++){                    
		hPhopt->Reset();
   		hPhoeta->Reset();
   		hPhophi->Reset();
   		hLeppt->Reset();
   		hLepeta->Reset();
   		hLepphi->Reset();
   		hsigMT->Reset();
   		hsigMET->Reset();
   		hsigMETPhi->Reset();
   		hdPhiLepMET->Reset();
   		hnVertex->Reset();
   		hdRPhoLep->Reset();
   		hHT->Reset();
   		hnJet->Reset();
   		hnBJet->Reset();
   		heleR9->Reset();
   		hphoR9->Reset();
   		hMiniIso->Reset();
		hMET_vs_MT ->Reset();
 		cout<<METcut[p] << "  "<<MTcut[q] << "  "<<p<<" "<<q<<endl;
   		for(Int_t i=0; i<nentries1; i++){                                                                                                                               
                  t->GetEntry(i);
                  cut = sigMET < METcut[p] || sigMT < MTcut[q];
		  if(cut) continue;
                  hPhopt->Fill(phoEt,lumiWeight);
                  hPhoeta ->Fill(phoEta,lumiWeight);
                  hPhophi ->Fill(phoPhi,lumiWeight);
                  hLeppt ->Fill(lepPt,lumiWeight);
                  hLepeta ->Fill(lepEta,lumiWeight);
                  hLepphi ->Fill(lepPhi,lumiWeight);
                  hsigMT ->Fill(sigMT,lumiWeight);
                  hsigMET ->Fill(sigMET,lumiWeight);
                  hsigMETPhi ->Fill(sigMETPhi,lumiWeight);
                  hdPhiLepMET ->Fill(dPhiLepMET,lumiWeight);
                  hnVertex ->Fill(nVertex,lumiWeight);
                  hdRPhoLep ->Fill(dRPhoLep,lumiWeight);
                  hHT ->Fill(HT,lumiWeight);
                  hnJet ->Fill(nJet,lumiWeight);
                  hnBJet ->Fill(nBJet,lumiWeight);
                  heleR9 ->Fill(eleR9,lumiWeight);
                  hphoR9 ->Fill(phoR9,lumiWeight);
                  hMiniIso ->Fill(MiniIso,lumiWeight);
		  hMET_vs_MT ->Fill(sigMET,sigMT);
   		}

  		TFile *file;
   		file = new TFile(Form("%s_MET%d_MT%d.root",outputfile,METcut[p],MTcut[q]), "RECREATE");
		hPhopt->Write();
   		hPhoeta->Write();
   		hPhophi->Write();
   		hLeppt->Write();
   		hLepeta->Write();
   		hLepphi->Write();
   		hsigMT->Write();
   		hsigMET->Write();
   		hsigMETPhi->Write();
   		hdPhiLepMET->Write();
   		hnVertex->Write();
   		hdRPhoLep->Write();
   		hHT->Write();
   		hnJet->Write();
   		hnBJet->Write();
   		heleR9->Write();
   		hphoR9->Write();
   		hMiniIso->Write();
		hMET_vs_MT ->Write();
		file->Close();
//		delete hPhopt, hPhoeta, hPhophi, hLeppt, hLepeta, hLepphi, hsigMT, hsigMET, hsigMETPhi, hdPhiLepMET, hnVertex, hdRPhoLep, hHT, hnJet, hnBJet, heleR9, hphoR9,  hMiniIso, hMET_vs_MT; 
	}
  }
  f->Close();
}

int main(int argc, char** argv)
{
    if(argc < 3)
    std::cout << "You have to provide two arguments!!\n";
    access_hist(argv[1],argv[2]);
    return 0;
}
