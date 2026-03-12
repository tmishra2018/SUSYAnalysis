#include "../../include/analysis_commoncode.h"

void quantify_doubleFake(){

  SetRunConfig();
  int channelType = ichannel; // 1=eg, 2=mg
  std::string whichVFP;
  if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
  if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
  if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

  std::string baseDir = "/uscms_data/d3/tmishra/Background/";

  std::string treeType;
  switch(anatype){
        case 0: treeType = "controlTree_"; break;
        case 1: treeType = "bkgTree_";     break;
        case 2: treeType = "validTree_";   break;
        case 3: treeType = "signalTree_";  break;
        default:
        std::cout << "ERROR: invalid anatype" << std::endl;
        return;
  }

  std::string chlabel = (channelType==1 ? "egamma" : "mg");

  std::string suffix="";
  if(anatype==0){
        std::ostringstream tmp;
        tmp << "_met" << lowMET << "_" << highMET
        << "_pt"  << lowPt  << "_" << highPt;
        suffix = tmp.str();
  }

  auto makeName = [&](std::string sample){
  std::ostringstream fname;
  fname << baseDir
        << treeType
        << chlabel << "_"
        << sample  << "_"
        << RunYear << whichVFP
        << suffix
        << ".root";
  return fname.str();
  };

  std::string eleName  = makeName("eleBkg");
  std::string jetName  = makeName("jetbkg");
  std::string qcdName  = makeName("qcd");
  std::string dfName   = makeName("doubleFake");
  std::string vgName   = makeName("VGBkg");
  std::string rareName = makeName("rareBkg");

  TFile *f_ele = TFile::Open(eleName.c_str());
  TFile *f_jet = TFile::Open(jetName.c_str());
  TFile *f_qcd = TFile::Open(qcdName.c_str());
  TFile *f_vg  = TFile::Open(vgName.c_str());
  TFile *f_rare= TFile::Open(rareName.c_str());
  TFile *f_df  = TFile::Open(dfName.c_str());

  if(!f_ele || !f_jet || !f_qcd || !f_vg || !f_rare || !f_df){
    std::cout << "ERROR: failed to open one or more files." << std::endl;
    std::cout << "  " << eleName.c_str() << std::endl;
    std::cout << "  " << jetName.c_str() << std::endl;
    std::cout << "  " << qcdName.c_str() << std::endl;
    std::cout << "  " << vgName.c_str() << std::endl;
    std::cout << "  " << rareName.c_str() << std::endl;
    std::cout << "  " << dfName.c_str() << std::endl;
    return;
  }

  TH1D *h_ele = (TH1D*)f_ele->Get("p_MET");
  TH1D *h_jet = (TH1D*)f_jet->Get("p_MET");
  TH1D *h_qcd = (TH1D*)f_qcd->Get("p_MET");
  TH1D *h_vg  = (TH1D*)f_vg->Get("p_MET");
  TH1D *h_rare= (TH1D*)f_rare->Get("p_MET");
  TH1D *h_df  = (TH1D*)f_df->Get("p_MET");

  if(!h_ele || !h_jet || !h_qcd || !h_vg || !h_rare || !h_df){
    std::cout << "ERROR: p_MET histogram missing in one of the files." << std::endl;
    return;
  }

  double y_ele = h_ele->Integral(0, h_ele->GetNbinsX()+1);
  double y_jet = h_jet->Integral(0, h_jet->GetNbinsX()+1);
  double y_qcd = h_qcd->Integral(0, h_qcd->GetNbinsX()+1);
  double y_vg  = h_vg->Integral(0, h_vg->GetNbinsX()+1);
  double y_rare= h_rare->Integral(0, h_rare->GetNbinsX()+1);
  double y_df  = h_df->Integral(0, h_df->GetNbinsX()+1);

  double y_singleSum = y_ele + y_jet + y_qcd;
  double y_corrected = y_singleSum - y_df;
  double overlapFrac = (y_singleSum>0 ? y_df/y_singleSum : 0);
  double y_smRaw = y_singleSum + y_vg + y_rare;
  double y_smCorrected = y_corrected + y_vg + y_rare;
  double dfFracSMRaw = (y_smRaw > 0 ? y_df/y_smRaw : 0);
  double dfFracSMCorrected = (y_smCorrected > 0 ? y_df/y_smCorrected : 0);

  std::cout << "================ Double-fake overlap summary ================" << std::endl;
  std::cout << "channel        : " << chlabel << std::endl;
  std::cout << "anatype        : " << anatype << std::endl;
  std::cout << "single sum     : ele + jet + qcd = " << y_singleSum << std::endl;
  std::cout << "double fake    : " << y_df << std::endl;
  std::cout << "VGamma         : " << y_vg << std::endl;
  std::cout << "Rare           : " << y_rare << std::endl;
  std::cout << "corrected total: single - double = " << y_corrected << std::endl;
  std::cout << "overlap fraction (double/single): " << overlapFrac << std::endl;
  std::cout << "overlap percentage (double/single): " << overlapFrac*100.0 << " %" << std::endl;
  std::cout << "SM total (raw) : single + VGamma + rare = " << y_smRaw << std::endl;
  std::cout << "SM total (corr): (single-double) + VGamma + rare = " << y_smCorrected << std::endl;
  std::cout << "double-fake contribution out of total SM (raw): " << dfFracSMRaw*100.0 << " %" << std::endl;
  std::cout << "double-fake contribution out of total SM (corrected): " << dfFracSMCorrected*100.0 << " %" << std::endl;

  std::cout << "\nPer-bin overlap fraction in MET:" << std::endl;
  for(int ib=1; ib<=h_ele->GetNbinsX(); ib++){
    double bsingle = h_ele->GetBinContent(ib) + h_jet->GetBinContent(ib) + h_qcd->GetBinContent(ib);
    double bdf = h_df->GetBinContent(ib);
    double frac = (bsingle>0 ? bdf/bsingle : 0);
    std::cout << "  bin " << ib << " [" << h_ele->GetXaxis()->GetBinLowEdge(ib) << ", " << h_ele->GetXaxis()->GetBinUpEdge(ib)
              << "] : single=" << bsingle << ", double=" << bdf << ", frac=" << frac << ", percent=" << frac*100.0 << " %" << std::endl;
  }


  // Write LaTeX table snippet for overlap summary
  std::ostringstream texName;
  texName << treeType << chlabel << "_doubleFakeOverlapTable";
  if(anatype==0){
    texName << "_met" << lowMET << "_" << highMET << "_pt" << lowPt << "_" << highPt << "_iso" << lepIso;
  }
  texName << ".tex";

  std::ofstream texout(texName.str().c_str());
  if(texout.is_open()){
    texout << "\\begin{table}[htbp]\n";
    texout << "\\centering\n";
    texout << "\\caption{Double-fake overlap summary for " << chlabel << ", anatype=" << anatype << ".}\n";
    texout << "\\begin{tabular}{lccc}\n";
    texout << "\\hline\n";
    texout << "Region & Single sum & Double fake & Overlap (\\%) \\\\ \n";
    texout << "\\hline\n";
    texout << "Inclusive & " << y_singleSum << " & " << y_df << " & " << overlapFrac*100.0 << " \\\\ \n";
    texout << "SM total (raw) & " << y_smRaw << " & " << y_df << " & " << dfFracSMRaw*100.0 << " \\\\ \n";
    texout << "SM total (corrected) & " << y_smCorrected << " & " << y_df << " & " << dfFracSMCorrected*100.0 << " \\\\ \n";
    texout << "\\hline\n";
    texout << "\\end{tabular}\n";
    texout << "\\vspace{0.3cm}\n";
    texout << "\\begin{tabular}{lccc}\n";
    texout << "\\hline\n";
    texout << "MET bin [GeV] & Single sum & Double fake & Overlap (\\%) \\\\ \n";
    texout << "\\hline\n";
    for(int ib=1; ib<=h_ele->GetNbinsX(); ib++){
      double bsingle = h_ele->GetBinContent(ib) + h_jet->GetBinContent(ib) + h_qcd->GetBinContent(ib);
      double bdf = h_df->GetBinContent(ib);
      double frac = (bsingle>0 ? bdf/bsingle : 0);
      texout << "[" << h_ele->GetXaxis()->GetBinLowEdge(ib) << ", " << h_ele->GetXaxis()->GetBinUpEdge(ib)
             << "] & " << bsingle << " & " << bdf << " & " << frac*100.0 << " \\\\ \n";
    }
    texout << "\\hline\n";
    texout << "\\end{tabular}\n";
    texout << "\\end{table}\n";
    texout.close();
    std::cout << "Wrote LaTeX overlap table to: " << texName.str() << std::endl;
  }

  TH1D *h_elePart = (TH1D*)f_df->Get("p_elePhoPart");
  TH1D *h_jetPart = (TH1D*)f_df->Get("p_jetPhoPart");
  if(h_elePart && h_jetPart){
    std::cout << "\nDouble-fake decomposition:" << std::endl;
    std::cout << "  e->gamma part : " << h_elePart->Integral(0, h_elePart->GetNbinsX()+1) << std::endl;
    std::cout << "  jet->gamma part: " << h_jetPart->Integral(0, h_jetPart->GetNbinsX()+1) << std::endl;
  }

  // Write corrected histograms for direct use in analysis/plotting
  std::vector<std::string> hnames;
  hnames.push_back("p_PhoEt");
  hnames.push_back("p_LepPt");
  hnames.push_back("p_MET");
  hnames.push_back("p_Mt");
  hnames.push_back("p_HT");
  hnames.push_back("p_PhoEta");
  hnames.push_back("p_LepEta");
  hnames.push_back("p_dPhiEleMET");
  hnames.push_back("p_nJet");
  hnames.push_back("p_nBJet");
  hnames.push_back("p_PhoEt_TT");
  hnames.push_back("p_MET_TT");
  hnames.push_back("p_Mt_TT");
  hnames.push_back("p_HT_TT");

  std::ostringstream outName;
  outName << treeType << chlabel << "_doubleFakeCorrected";
  if(anatype==0){
    outName << "_met" << lowMET << "_" << highMET << "_pt" << lowPt << "_" << highPt << "_iso" << lepIso;
  }
  outName << ".root";

  TFile *f_out = TFile::Open(outName.str().c_str(), "RECREATE");
  for(unsigned ih(0); ih<hnames.size(); ih++){
    TH1D *he = (TH1D*)f_ele->Get(hnames[ih].c_str());
    TH1D *hj = (TH1D*)f_jet->Get(hnames[ih].c_str());
    TH1D *hq = (TH1D*)f_qcd->Get(hnames[ih].c_str());
    TH1D *hd = (TH1D*)f_df->Get(hnames[ih].c_str());
    if(!he || !hj || !hq || !hd)continue;

    TH1D *h_single = (TH1D*)he->Clone((hnames[ih] + std::string("_singleSum")).c_str());
    h_single->Add(hj);
    h_single->Add(hq);
    TH1D *h_corr = (TH1D*)h_single->Clone((hnames[ih] + std::string("_corrected")).c_str());
    h_corr->Add(hd, -1.0);

    h_single->Write();
    h_corr->Write();
  }
  f_out->Write();
  f_out->Close();
  std::cout << "\nWrote corrected histograms to: " << outName.str() << std::endl;
  std::cout << "==============================================================" << std::endl;
}
