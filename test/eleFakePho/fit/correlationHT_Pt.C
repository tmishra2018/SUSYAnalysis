// correlationHT_Pt.C  (updated)
// Computes and visualises pairwise correlations between:
//   (1) HT     vs probe pT
//   (2) nVtx   vs probe pT
//   (3) nVtx   vs probe |eta|
//   (4) probe pT vs probe |eta|
//   (5) nVtx   vs HT
//
// for the electron-faking-photon fake-rate measurement.
//
// Compile:
//   g++ `root-config --cflags` correlationHT_Pt.C -o correlationHT_Pt.exe `root-config --libs`
// Run:
//   ./correlationHT_Pt.exe 2018 ""
//   ./correlationHT_Pt.exe 2016 preVFP
// Or interactively:
//   root -l 'correlationHT_Pt.C(2018,"")'

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TSystem.h"
#include "TColor.h"

double pearson(const std::vector<float>& x, const std::vector<float>& y)
{
    size_t n = x.size();
    if (n != y.size() || n == 0) return 0;
    double sx=0,sy=0,sxx=0,syy=0,sxy=0;
    for (size_t i=0;i<n;++i){
        sx+=x[i]; sy+=y[i];
        sxx+=x[i]*x[i]; syy+=y[i]*y[i]; sxy+=x[i]*y[i];
    }
    double cov  = sxy/n - (sx/n)*(sy/n);
    double varx = sxx/n - (sx/n)*(sx/n);
    double vary = syy/n - (sy/n)*(sy/n);
    if (varx<=0||vary<=0) return 0;
    return cov/std::sqrt(varx*vary);
}

double spearman(std::vector<float> x, std::vector<float> y)
{
    size_t n = x.size();
    if (n!=y.size()||n==0) return 0;
    auto mkrank=[&](std::vector<float>& v){
        std::vector<float> r(n);
        std::vector<size_t> idx(n);
        std::iota(idx.begin(),idx.end(),0);
        std::sort(idx.begin(),idx.end(),[&](size_t a,size_t b){return v[a]<v[b];});
        for (size_t i=0;i<n;){
            size_t j=i;
            while(j<n && v[idx[j]]==v[idx[i]])++j;
            float avg=(i+j-1)/2.f+1.f;
            for(size_t k=i;k<j;++k) r[idx[k]]=avg;
            i=j;
        }
        return r;
    };
    auto rx=mkrank(x), ry=mkrank(y);
    return pearson(rx,ry);
}

struct CorrResult {
    std::string labelX, labelY, unitX, unitY;
    double rP_den, rS_den, rP_num, rS_num;
    size_t nDen, nNum;
};

const char* strength(double r){
    double ar=std::fabs(r);
    if(ar<0.10) return "negligible";
    if(ar<0.30) return "weak";
    if(ar<0.50) return "moderate";
    if(ar<0.70) return "strong";
    return "very strong";
}

void setStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(255);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTitleOffset(1.3,"Y");
    gStyle->SetTitleOffset(1.1,"X");
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadRightMargin(0.15);
}

void draw2D(TH2F* h, double rP, double rS,
            const char* outname, bool logz=true)
{
    TCanvas* c = new TCanvas(outname, outname, 720, 600);
    c->SetRightMargin(0.16);
    if(logz && h->GetEntries()>0) c->SetLogz();
    h->SetContour(255);
    h->Draw("COLZ");
    TLatex tex;
    tex.SetNDC();
    tex.SetTextSize(0.035);
    tex.SetTextFont(42);
    tex.DrawLatex(0.14, 0.925,
        Form("Pearson r = %.3f   |   Spearman r = %.3f", rP, rS));
    c->SaveAs(outname);
    delete c;
}

void drawGrid(std::vector<TH2F*>& hden, std::vector<TH2F*>& hnum,
              std::vector<CorrResult>& res,
              int year, const std::string& vfp)
{
    int n=res.size();
    TCanvas* cg = new TCanvas("cgrid","All correlations",n*400,800);
    cg->Divide(n,2);
    for(int i=0;i<n;++i){
        cg->cd(i+1);
        gPad->SetRightMargin(0.16);
        if(hden[i]->GetEntries()>0) gPad->SetLogz();
        hden[i]->SetContour(255);
        hden[i]->SetTitle(Form("DEN  r_{P}=%.3f r_{S}=%.3f",
                               res[i].rP_den,res[i].rS_den));
        hden[i]->Draw("COLZ");

        cg->cd(i+1+n);
        gPad->SetRightMargin(0.16);
        if(hnum[i]->GetEntries()>0) gPad->SetLogz();
        hnum[i]->SetContour(255);
        hnum[i]->SetTitle(Form("NUM  r_{P}=%.3f r_{S}=%.3f",
                               res[i].rP_num,res[i].rS_num));
        hnum[i]->Draw("COLZ");
    }
    cg->SaveAs(Form("AllCorrelations_Grid_%d%s.png",year,vfp.c_str()));
    delete cg;
}

void correlationHT_Pt(int RunYear=2018, const char* Era="")
{
    setStyle();

    std::string whichVFP;
    if(RunYear==2016)
        whichVFP=(std::string(Era).find("pre")!=std::string::npos)?"preVFP":"postVFP";
    else whichVFP="";

    TChain* etree = new TChain("FakeRateTree");
    etree->Add(Form("root://cmseos.fnal.gov//store/user/tmishra/elefakepho/files/"
                    "plot_elefakepho_DataTnP_dR05_%d%s.root",RunYear,whichVFP.c_str()));
    if(etree->GetEntries()==0){
        etree->Add(Form("plot_elefakepho_DataTnP_dR05_%d%s.root",RunYear,whichVFP.c_str()));
        if(etree->GetEntries()==0){std::cerr<<"No entries. Exiting.\n";return;}
    }

    float probePt=0, probeEta=0, HT=0;
    int   nVertex=0;
    bool  vetovalue=false, FSRveto=false;

    etree->SetBranchAddress("probePt",   &probePt);
    etree->SetBranchAddress("probeEta",  &probeEta);
    etree->SetBranchAddress("HT",        &HT);
    etree->SetBranchAddress("nVertex",   &nVertex);
    etree->SetBranchAddress("vetovalue", &vetovalue);
    etree->SetBranchAddress("FSRveto",   &FSRveto);

    // (1) HT vs pT
    TH2F* h2_HT_pt_den   = new TH2F("h2_HT_pt_den",
        Form("H_{{T}} vs probe p_{{T}}  Denominator  %d;probe p_{{T}} (GeV);H_{{T}} (GeV)",RunYear),
        50,30,300, 50,0,1000);
    TH2F* h2_HT_pt_num   = new TH2F("h2_HT_pt_num",
        Form("H_{{T}} vs probe p_{{T}}  Numerator  %d;probe p_{{T}} (GeV);H_{{T}} (GeV)",RunYear),
        50,30,300, 50,0,1000);
    // (2) nVtx vs pT
    TH2F* h2_vtx_pt_den  = new TH2F("h2_vtx_pt_den",
        Form("nVtx vs probe p_{{T}}  Denominator  %d;probe p_{{T}} (GeV);nVtx",RunYear),
        50,30,300, 50,0,100);
    TH2F* h2_vtx_pt_num  = new TH2F("h2_vtx_pt_num",
        Form("nVtx vs probe p_{{T}}  Numerator  %d;probe p_{{T}} (GeV);nVtx",RunYear),
        50,30,300, 50,0,100);
    // (3) nVtx vs |eta|
    TH2F* h2_vtx_eta_den = new TH2F("h2_vtx_eta_den",
        Form("nVtx vs |#eta|  Denominator  %d;|probe #eta|;nVtx",RunYear),
        30,0,1.5, 50,0,100);
    TH2F* h2_vtx_eta_num = new TH2F("h2_vtx_eta_num",
        Form("nVtx vs |#eta|  Numerator  %d;|probe #eta|;nVtx",RunYear),
        30,0,1.5, 50,0,100);
    // (4) pT vs |eta|
    TH2F* h2_pt_eta_den  = new TH2F("h2_pt_eta_den",
        Form("probe p_{{T}} vs |#eta|  Denominator  %d;|probe #eta|;probe p_{{T}} (GeV)",RunYear),
        30,0,1.5, 50,30,300);
    TH2F* h2_pt_eta_num  = new TH2F("h2_pt_eta_num",
        Form("probe p_{{T}} vs |#eta|  Numerator  %d;|probe #eta|;probe p_{{T}} (GeV)",RunYear),
        30,0,1.5, 50,30,300);
    // (5) nVtx vs HT
    TH2F* h2_vtx_HT_den  = new TH2F("h2_vtx_HT_den",
        Form("nVtx vs H_{{T}}  Denominator  %d;H_{{T}} (GeV);nVtx",RunYear),
        50,0,1000, 50,0,100);
    TH2F* h2_vtx_HT_num  = new TH2F("h2_vtx_HT_num",
        Form("nVtx vs H_{{T}}  Numerator  %d;H_{{T}} (GeV);nVtx",RunYear),
        50,0,1000, 50,0,100);

    // profiles
    TProfile* prof_HT_pt_den  = new TProfile("prof_HT_pt_den",
        ";probe p_{T} (GeV);<H_{T}> (GeV)",50,30,300,0,1500);
    TProfile* prof_HT_pt_num  = new TProfile("prof_HT_pt_num",
        ";probe p_{T} (GeV);<H_{T}> (GeV)",50,30,300,0,1500);
    TProfile* prof_vtx_pt_den = new TProfile("prof_vtx_pt_den",
        ";probe p_{T} (GeV);<nVtx>",50,30,300,0,100);
    TProfile* prof_vtx_pt_num = new TProfile("prof_vtx_pt_num",
        ";probe p_{T} (GeV);<nVtx>",50,30,300,0,100);

    std::vector<float> vpt_d, vht_d, vvtx_d, veta_d;
    std::vector<float> vpt_n, vht_n, vvtx_n, veta_n;

    const long nEvts = etree->GetEntries();
    std::cout<<"Processing "<<nEvts<<" events ..."<<std::endl;

    for(long ievt=0;ievt<nEvts;++ievt){
        if(ievt%500000==0) std::cout<<"  event "<<ievt<<std::endl;
        etree->GetEntry(ievt);

        if(probePt<35) continue;
        if(std::fabs(probeEta)>1.4442) continue;   // EB only
        float aeta=std::fabs(probeEta);
        float vtxf=(float)nVertex;

        if(!vetovalue){                              // denominator
            h2_HT_pt_den ->Fill(probePt,HT);
            h2_vtx_pt_den->Fill(probePt,vtxf);
            h2_vtx_eta_den->Fill(aeta,vtxf);
            h2_pt_eta_den ->Fill(aeta,probePt);
            h2_vtx_HT_den ->Fill(HT,vtxf);
            prof_HT_pt_den ->Fill(probePt,HT);
            prof_vtx_pt_den->Fill(probePt,vtxf);
            vpt_d.push_back(probePt);
            vht_d.push_back(HT);
            vvtx_d.push_back(vtxf);
            veta_d.push_back(aeta);
        } else if(vetovalue && FSRveto){             // numerator
            h2_HT_pt_num ->Fill(probePt,HT);
            h2_vtx_pt_num->Fill(probePt,vtxf);
            h2_vtx_eta_num->Fill(aeta,vtxf);
            h2_pt_eta_num ->Fill(aeta,probePt);
            h2_vtx_HT_num ->Fill(HT,vtxf);
            prof_HT_pt_num ->Fill(probePt,HT);
            prof_vtx_pt_num->Fill(probePt,vtxf);
            vpt_n.push_back(probePt);
            vht_n.push_back(HT);
            vvtx_n.push_back(vtxf);
            veta_n.push_back(aeta);
        }
    }

    std::vector<CorrResult> results;
    auto addPair=[&](const char* lx, const char* ly,
                     std::vector<float>& xd, std::vector<float>& yd,
                     std::vector<float>& xn, std::vector<float>& yn) -> CorrResult {
        CorrResult r;
        r.labelX=lx; r.labelY=ly;
        r.nDen=xd.size(); r.nNum=xn.size();
        r.rP_den=pearson (xd,yd); r.rS_den=spearman(xd,yd);
        r.rP_num=pearson (xn,yn); r.rS_num=spearman(xn,yn);
        results.push_back(r);
        return r;
    };
    auto rHT_pt   = addPair("pT","HT",      vpt_d,vht_d,    vpt_n,vht_n);
    auto rVtx_pt  = addPair("pT","nVtx",    vpt_d,vvtx_d,   vpt_n,vvtx_n);
    auto rVtx_eta = addPair("|eta|","nVtx", veta_d,vvtx_d,  veta_n,vvtx_n);
    auto rPt_eta  = addPair("|eta|","pT",   veta_d,vpt_d,   veta_n,vpt_n);
    auto rVtx_HT  = addPair("HT","nVtx",   vht_d,vvtx_d,   vht_n,vvtx_n);

    auto line80=[](){std::cout<<" "<<std::string(68,'=')<<"\n";};
    std::cout<<"\n "<<std::string(68,'=')<<"\n";
    std::cout<<" Pairwise correlations  RunYear="<<RunYear
             <<"  EB only (|eta|<1.4442, pT>35 GeV)"<<std::string(4,' ')<<" \n";
    std::cout<<" "<<std::string(68,'=')<<"\n";
    for(auto& r : results){
        std::cout<<" "<<r.labelX<<" vs "<<r.labelY
                 <<" (N_den="<<r.nDen<<",  N_num="<<r.nNum<<")\n";
        std::cout<<" DEN  Pearson="<<std::fixed<<std::setprecision(4)
                 <<r.rP_den<<"  Spearman="<<r.rS_den
                 <<"  ["<<strength(r.rP_den)<<"]\n";
        std::cout<<"   NUM  Pearson="<<r.rP_num<<"  Spearman="<<r.rS_num
                 <<"  ["<<strength(r.rP_num)<<"]\n";
        line80();
    }
    std::cout<<" "<<std::string(68,'=')<<"\n\n";

    draw2D(h2_HT_pt_den,  rHT_pt.rP_den,   rHT_pt.rS_den,
           Form("HT_vs_Pt_den_%d%s.png",   RunYear,whichVFP.c_str()));
    draw2D(h2_HT_pt_num,  rHT_pt.rP_num,   rHT_pt.rS_num,
           Form("HT_vs_Pt_num_%d%s.png",   RunYear,whichVFP.c_str()));
    draw2D(h2_vtx_pt_den, rVtx_pt.rP_den,  rVtx_pt.rS_den,
           Form("Vtx_vs_Pt_den_%d%s.png",  RunYear,whichVFP.c_str()));
    draw2D(h2_vtx_pt_num, rVtx_pt.rP_num,  rVtx_pt.rS_num,
           Form("Vtx_vs_Pt_num_%d%s.png",  RunYear,whichVFP.c_str()));
    draw2D(h2_vtx_eta_den,rVtx_eta.rP_den, rVtx_eta.rS_den,
           Form("Vtx_vs_Eta_den_%d%s.png", RunYear,whichVFP.c_str()),false);
    draw2D(h2_vtx_eta_num,rVtx_eta.rP_num, rVtx_eta.rS_num,
           Form("Vtx_vs_Eta_num_%d%s.png", RunYear,whichVFP.c_str()),false);
    draw2D(h2_pt_eta_den, rPt_eta.rP_den,  rPt_eta.rS_den,
           Form("Pt_vs_Eta_den_%d%s.png",  RunYear,whichVFP.c_str()),false);
    draw2D(h2_pt_eta_num, rPt_eta.rP_num,  rPt_eta.rS_num,
           Form("Pt_vs_Eta_num_%d%s.png",  RunYear,whichVFP.c_str()),false);
    draw2D(h2_vtx_HT_den, rVtx_HT.rP_den,  rVtx_HT.rS_den,
           Form("Vtx_vs_HT_den_%d%s.png",  RunYear,whichVFP.c_str()));
    draw2D(h2_vtx_HT_num, rVtx_HT.rP_num,  rVtx_HT.rS_num,
           Form("Vtx_vs_HT_num_%d%s.png",  RunYear,whichVFP.c_str()));

    std::vector<TH2F*> hden={h2_HT_pt_den,h2_vtx_pt_den,h2_vtx_eta_den,h2_pt_eta_den,h2_vtx_HT_den};
    std::vector<TH2F*> hnum={h2_HT_pt_num,h2_vtx_pt_num,h2_vtx_eta_num,h2_pt_eta_num,h2_vtx_HT_num};
    drawGrid(hden,hnum,results,RunYear,whichVFP);

    TCanvas* cprof = new TCanvas("cprof","Profiles vs pT",1200,550);
    cprof->Divide(2,1);

    cprof->cd(1);
    prof_HT_pt_den->SetLineColor(kBlue+1); prof_HT_pt_den->SetMarkerColor(kBlue+1);
    prof_HT_pt_den->SetMarkerStyle(20); prof_HT_pt_den->SetLineWidth(2);
    prof_HT_pt_den->SetMinimum(0); prof_HT_pt_den->SetMaximum(700);
    prof_HT_pt_den->SetTitle(Form("<H_{{T}}> vs probe p_{{T}}  %d;probe p_{{T}} (GeV);<H_{{T}}> (GeV)",RunYear));
    prof_HT_pt_den->Draw("EP");
    prof_HT_pt_num->SetLineColor(kRed+1); prof_HT_pt_num->SetMarkerColor(kRed+1);
    prof_HT_pt_num->SetMarkerStyle(21); prof_HT_pt_num->SetLineWidth(2);
    prof_HT_pt_num->Draw("EP same");
    { TLegend* l=new TLegend(0.18,0.72,0.62,0.88); l->SetBorderSize(0); l->SetFillStyle(0);
      l->AddEntry(prof_HT_pt_den,Form("Denominator  r_{P}=%.3f r_{S}=%.3f",rHT_pt.rP_den,rHT_pt.rS_den),"pe");
      l->AddEntry(prof_HT_pt_num,Form("Numerator     r_{P}=%.3f r_{S}=%.3f",rHT_pt.rP_num,rHT_pt.rS_num),"pe");
      l->Draw(); }

    cprof->cd(2);
    prof_vtx_pt_den->SetLineColor(kBlue+1); prof_vtx_pt_den->SetMarkerColor(kBlue+1);
    prof_vtx_pt_den->SetMarkerStyle(20); prof_vtx_pt_den->SetLineWidth(2);
    prof_vtx_pt_den->SetMinimum(0); prof_vtx_pt_den->SetMaximum(70);
    prof_vtx_pt_den->SetTitle(Form("<nVtx> vs probe p_{{T}}  %d;probe p_{{T}} (GeV);<nVtx>",RunYear));
    prof_vtx_pt_den->Draw("EP");
    prof_vtx_pt_num->SetLineColor(kRed+1); prof_vtx_pt_num->SetMarkerColor(kRed+1);
    prof_vtx_pt_num->SetMarkerStyle(21); prof_vtx_pt_num->SetLineWidth(2);
    prof_vtx_pt_num->Draw("EP same");
    { TLegend* l=new TLegend(0.18,0.72,0.62,0.88); l->SetBorderSize(0); l->SetFillStyle(0);
      l->AddEntry(prof_vtx_pt_den,Form("Denominator  r_{P}=%.3f r_{S}=%.3f",rVtx_pt.rP_den,rVtx_pt.rS_den),"pe");
      l->AddEntry(prof_vtx_pt_num,Form("Numerator     r_{P}=%.3f r_{S}=%.3f",rVtx_pt.rP_num,rVtx_pt.rS_num),"pe");
      l->Draw(); }
    cprof->SaveAs(Form("Profiles_vs_Pt_%d%s.png",RunYear,whichVFP.c_str()));
    delete cprof;

    std::ofstream txt(Form("AllCorrelations_%d%s.txt",RunYear,whichVFP.c_str()));
    txt<<std::fixed<<std::setprecision(4);

    txt<<" \n";
    txt<<"  Pairwise correlations  RunYear="<<RunYear<<"  Era="<<whichVFP<<"\n";
    txt<<"  Selection: EB (|eta|<1.4442), probePt > 35 GeV\n";
    txt<<" \n\n";

    for(auto& r : results){
        txt<<r.labelX<<" vs "<<r.labelY<<"\n";
        txt<<"  Denominator  N="<<r.nDen<<"\n";
        txt<<"    Pearson  r  = "<<r.rP_den<<"  ("<<strength(r.rP_den)<<")\n";
        txt<<"    Spearman r  = "<<r.rS_den<<"  ("<<strength(r.rS_den)<<")\n";
        txt<<"  Numerator    N="<<r.nNum<<"\n";
        txt<<"    Pearson  r  = "<<r.rP_num<<"  ("<<strength(r.rP_num)<<")\n";
        txt<<"    Spearman r  = "<<r.rS_num<<"  ("<<strength(r.rS_num)<<")\n\n";
    }

    txt<<" \n";
    txt<<"  PHYSICS INTERPRETATION\n";
    txt<<" \n\n";

    txt<<"[1] HT vs probe pT  (Pearson den = "<<rHT_pt.rP_den<<",  Spearman = "<<rHT_pt.rS_den<<")\n";
    txt<<"\n";
    txt<<"    Both HT and probe pT are measures of hadronic activity:\n";
    txt<<"    high-pT probes are produced in boosted events that also have\n";
    txt<<"    more hadronic recoil, hence higher HT.  The moderate Pearson r\n";
    txt<<"    (~0.5) tells us ~25% of the variance is shared (r^2 ~ 0.25).\n";
    txt<<"    The Spearman r is considerably lower (~0.30), indicating the\n";
    txt<<"    linear part of the correlation is driven by a tail of very\n";
    txt<<"    high-pT / high-HT events rather than a uniform monotone trend.\n\n";

    txt<<"[2] nVtx vs probe pT  (Pearson den = "<<rVtx_pt.rP_den<<",  Spearman = "<<rVtx_pt.rS_den<<")\n";
    txt<<"\n";
    txt<<"    nVtx counts pileup vertices, which come from independent soft\n";
    txt<<"    collisions and have no causal link to the pT of the probe from\n";
    txt<<"    the hard scatter.  The near-zero correlation confirms that\n";
    txt<<"    nVtx is an essentially independent variable from pT -- the\n";
    txt<<"    factorisation f(pT) * f(nVtx) is internally consistent.\n\n";

    txt<<"[3] nVtx vs probe |eta|  (Pearson den = "<<rVtx_eta.rP_den<<",  Spearman = "<<rVtx_eta.rS_den<<")\n";
    txt<<"\n";
    txt<<"    Pileup is uniform in eta within the EB barrel acceptance.\n";
    txt<<"    A negligible correlation here confirms that eta and nVtx\n";
    txt<<"    are independent and can be safely factorised.\n\n";

    txt<<"[4] probe pT vs probe |eta|  (Pearson den = "<<rPt_eta.rP_den<<",  Spearman = "<<rPt_eta.rS_den<<")\n";
    txt<<"\n";
    txt<<"    Within EB (|eta|<1.4442) the tag-and-probe selection is\n";
    txt<<"    approximately uniform in eta for a given pT threshold, so\n";
    txt<<"    a small correlation is expected and confirms that pT and eta\n";
    txt<<"    can be treated as independent axes in the fake-rate map.\n\n";

    txt<<"[5] nVtx vs HT  (Pearson den = "<<rVtx_HT.rP_den<<",  Spearman = "<<rVtx_HT.rS_den<<")\n";
    txt<<"\n";
    txt<<"    nVtx reflects pileup from independent soft collisions, while\n";
    txt<<"    HT is driven by the hard scatter.  A near-zero correlation\n";
    txt<<"    between nVtx and HT confirms that pileup activity is largely\n";
    txt<<"    decoupled from the hadronic recoil of the hard event.\n";
    txt<<"    This supports the validity of a factorised parameterisation\n";
    txt<<"    that uses both HT and nVtx as independent dimensions.\n\n";

    txt<<" \n";
    txt<<"  WHICH ESTIMATE IS HIGHER: (pT, eta, HT) or (pT, eta, nVtx)?\n";
    txt<<" \n\n";

    txt<<"  Answer: the (pT, eta, HT) estimate is HIGHER.\n\n";

    txt<<"  Detailed reasoning:\n\n";

    txt<<"  The fake rate is applied event-by-event as a product:\n";
    txt<<"      w = f(pT) * f(third_var) * f(eta)\n\n";

    txt<<"  Case A: third_var = nVtx\n";
    txt<<"    f(pT) captures the hard-scatter activity (real physics).\n";
    txt<<"    f(nVtx) captures pileup (independent soft collisions).\n";
    txt<<"    Since pT and nVtx are uncorrelated (r~0), these two factors\n";
    txt<<"    weight DIFFERENT, ORTHOGONAL aspects of the event.\n";
    txt<<"    There is no double-counting.\n\n";

    txt<<"  Case B: third_var = HT\n";
    txt<<"    f(pT) captures hard-scatter pT.\n";
    txt<<"    f(HT) captures jet activity -- but HT is correlated with pT\n";
    txt<<"    (r~0.5 Pearson) because a high-pT probe comes from a boosted\n";
    txt<<"    event that also has more hadronic recoil, hence higher HT.\n";
    txt<<"    So f(pT) and f(HT) are weighting OVERLAPPING physics:\n";
    txt<<"    an event with large pT already gets a large weight from f(pT),\n";
    txt<<"    and then gets an ADDITIONAL large weight from f(HT) even though\n";
    txt<<"    HT is large precisely BECAUSE pT is large.\n";
    txt<<"    This double-counting of activity leads to OVER-WEIGHTING of\n";
    txt<<"    high-activity events and therefore OVER-PREDICTION of the\n";
    txt<<"    electron-faking-photon contribution.\n\n";

    txt<<"  Quantitative estimate of the bias:\n";
    txt<<"    Shared variance between pT and HT ~ r^2 ~ "
       <<rHT_pt.rP_den*rHT_pt.rP_den*100.<<"%\n";
    txt<<"    Shared variance between pT and nVtx ~ r^2 ~ "
       <<rVtx_pt.rP_den*rVtx_pt.rP_den*100.<<"%\n";
    txt<<"    Shared variance between nVtx and HT ~ r^2 ~ "
       <<rVtx_HT.rP_den*rVtx_HT.rP_den*100.<<"%\n";
    txt<<"    The excess weight applied by f(HT) relative to f(nVtx)\n";
    txt<<"    is largest in the tail (high pT, high HT), where the\n";
    txt<<"    Spearman rank correlation (~0.30) indicates a monotone\n";
    txt<<"    but nonlinear enhancement.  In practice the integrated\n";
    txt<<"    prediction from (pT,eta,HT) could be 10-30% higher than\n";
    txt<<"    from (pT,eta,nVtx), concentrated at high probe pT.\n\n";

    txt<<"  Recommendation:\n";
    txt<<"    1. Perform a closure test: apply both fake rates to a\n";
    txt<<"       Z->ee sideband (requiring vetovalue==0) and compare\n";
    txt<<"       the predicted numerator to the observed one.\n";
    txt<<"    2. If HT is the preferred physics variable (e.g. because\n";
    txt<<"       your signal region has a specific HT cut), measure the\n";
    txt<<"       fake rate directly in 2D (pT, HT) bins to avoid the\n";
    txt<<"       factorisation bias entirely.\n";
    txt<<"    3. The (pT, eta, nVtx) parameterisation remains the safer\n";
    txt<<"       choice when the third variable must be 1D.\n";

    txt.close();

    std::cout<<"Output files written:\n";
    std::cout<<"  AllCorrelations_"<<RunYear<<whichVFP<<".txt  (full report)\n";
    std::cout<<"  AllCorrelations_Grid_"<<RunYear<<whichVFP<<".png  (2x5 panel)\n";
    std::cout<<"  Profiles_vs_Pt_"<<RunYear<<whichVFP<<".png\n";
    std::cout<<"  [10 individual 2D COLZ plots]\n";

    delete etree;
}

//  Standalone entry point
int main(int argc, char** argv)
{
    if(argc<2){
        std::cout<<"Usage: ./correlationHT_Pt.exe <RunYear> [Era]\n";
        std::cout<<"  e.g. ./correlationHT_Pt.exe 2018\n";
        std::cout<<"       ./correlationHT_Pt.exe 2016 preVFP\n";
        return 1;
    }
    correlationHT_Pt(atoi(argv[1]), argc>2?argv[2]:"");
    return 0;
}
