#ifndef PlotLimits_H
#define PlotLimits_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "NtupleVariables.h"
#include "TH1F.h"
#include "TH2.h"
#include <TProfile.h>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TDirectory.h"

class PlotLimits : public NtupleVariables{

 public:
  PlotLimits(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data");
  ~PlotLimits();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *,const char *);
  void     BookHistogram(const char *);

  TH2D *h2_mGlmNLSP_r;
  TH2D *h2_mGlmNLSP_rUnc;
  TH2D *h2_mGlmNLSP_XsecUL,*h2_mGlmNLSP_XsecULexp;
  TH2D *h2_mGlmNLSP_XsecUL_v1;

  TH2D *h2_mGlmNLSP_r_fb;
  TH2D *h2_mGlmNLSP_XsecUL_fb;

  TH2D *h2_mGlmNLSP_r_XsecUp,*h2_mGlmNLSP_r_XsecDn;

  TH2D *h2_mGlmNLSP_median;
  TH2D *h2_mGlmNLSP_16pc;
  TH2D *h2_mGlmNLSP_84pc;
  TH2D *h2_mGlmNLSP_2p5pc;
  TH2D *h2_mGlmNLSP_97p5pc;

  TFile *oFile;
};
#endif

#ifdef PlotLimits_cxx

void PlotLimits::BookHistogram(const char *outFileName) {

  //  char hname[200], htit[200];
  //  double xlow = 0.0,  xhigh = 2000.0;
  //  int nbins = 2000;
  char name[100],title[100];
 
  oFile = new TFile(outFileName, "recreate");
  TH1::SetDefaultSumw2(1);
  h2_mGlmNLSP_r = new TH2D("mGlmNLSP_r","x: mass of gluino, y: mass of neutralino, z: r-value",100,12.5,2512.5,350,5,3505);
  h2_mGlmNLSP_XsecUL = new TH2D("mGlmNLSP_XsecUL","x: mass of gluino, y: mass of neutralino, z: 95% CL UpLimit xsec(fb)",100,12.5,2512.5,350,5,3505);
  h2_mGlmNLSP_r_fb = new TH2D("mGlmNLSP_r_fb","x: mass of gluino, y: mass of neutralino, z: r-value fine binned",3000,-0.5,2999.5,3000,-0.5,2999.5);
  h2_mGlmNLSP_XsecUL_fb = new TH2D("mGlmNLSP_XsecUL_fb","x: mass of gluino, y: mass of neutralino, z: 95% CL UpLimit xsec(fb) fine binned",3000,-0.5,2999.5,3000,-0.5,2999.5);
  ///////////////////////
  h2_mGlmNLSP_XsecUL_v1 = new TH2D("mGlmNLSP_XsecUL_v1","x: mass of gluino, y: mass of neutralino, z: 95% CL UpLimit xsec(fb)",60,25,3025,350,5,3505);
  //////////////////////
  h2_mGlmNLSP_rUnc = new TH2D("mGlmNLSP_rUnc","x: mass of gluino, y: mass of neutralino, z: unc r-value",100,12.5,2512.5,350,5,3505);
  h2_mGlmNLSP_r_XsecUp = new TH2D("mGlmNLSP_r_XsecUp","x: mass of gluino, y: mass of neutralino, z: r-value with #sigma+#sigma_{up}(fb)",100,12.5,2512.5,350,5,3505);
  h2_mGlmNLSP_r_XsecDn = new TH2D("mGlmNLSP_r_XsecDn","x: mass of gluino, y: mass of neutralino, z: r-value with #sigma-#sigma_{down}(fb)",100,12.5,2512.5,350,5,3505);

  h2_mGlmNLSP_XsecULexp = new TH2D("mGlmNLSP_XsecULexp","x: mass of gluino, y: mass of neutralino, z: Expected 95% CL UpLimit xsec(fb)",100,12.5,2512.5,350,5,3505);
  h2_mGlmNLSP_median = new TH2D("mGlmNLSP_median","x: mass of gluino, y: mass of neutralino, z: median exp limit xsec(fb)",100,12.5,2512.5,350,5,3505);
  h2_mGlmNLSP_16pc = new TH2D("mGlmNLSP_16pc","x: mass of gluino, y: mass of neutralino, z: 16pc exp limit xsec(fb)",100,12.5,2512.5,350,5,3505);
  h2_mGlmNLSP_84pc = new TH2D("mGlmNLSP_84pc","x: mass of gluino, y: mass of neutralino, z: 84pc exp limit xsec(fb)",100,12.5,2512.5,350,5,3505);
  h2_mGlmNLSP_2p5pc = new TH2D("mGlmNLSP_2p5pc","x: mass of gluino, y: mass of neutralino, z: 2.5pc exp limit xsec(fb)",100,12.5,2512.5,350,5,3505);
  h2_mGlmNLSP_97p5pc = new TH2D("mGlmNLSP_97p5pc","x: mass of gluino, y: mass of neutralino, z: 97p5pc exp limit xsec(fb)",100,12.5,2512.5,350,5,3505);
}


PlotLimits::PlotLimits(const TString &inputFileList, const char *outFileName, const char* dataset) {
  string nameData=dataset;//vvv
  TChain *tree = new TChain("limit");
  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
  }

  NtupleVariables::Init(tree,nameData);

  BookHistogram(outFileName);
  
}

Bool_t PlotLimits::FillChain(TChain *chain, const TString &inputFileList) {

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    //std::cout << "Adding tree from " << buffer.c_str() << std::endl;                                                              
    chain->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  return kTRUE;
}

Long64_t PlotLimits::LoadTree(Long64_t entry) {
  // Set the environment to read one entry                                                                                          
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }
  return centry;
}

PlotLimits::~PlotLimits() { 

  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  oFile->Close();

}

#endif

