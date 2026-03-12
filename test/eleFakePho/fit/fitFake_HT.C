#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TMatrixDSym.h"
#include "TRandom3.h"
#include "../../../include/tdrstyle.C"
#include "TROOT.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooAbsReal.h"
#include "RooKeysPdf.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooNumIntConfig.h"
#include "RooFFTConvPdf.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooHistPdf.h"
#include "RooMCStudy.h"
#include "RooChi2MCSModule.h"
#include "RooFitResult.h"
#include "TFitResult.h"
#include "TVirtualFitter.h"
#include "RooMultiVarGaussian.h"
#include "../../../include/RooCMSShape.h"
#include "../../../include/RooDCBShape.h"
#include "../../../include/RooUserPoly.h"
// #include "../../../include/analysis_fakes.h"
#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_jet.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_mcData.h"
#define NTOY 100

#define MAXVTX 48

int MINVTX = 18;
bool doEB = true;
const char*processName ="Data"; bool doDrellYan = false;
void fitFake_HT(int RunYear, bool preVFP)
{
    gROOT->SetBatch(1);
    setTDRStyle();
    gStyle->SetOptStat(0);

    std::string whichVFP="";
    if(RunYear==2016 && preVFP)  whichVFP="preVFP";
    if(RunYear==2016 && !preVFP) whichVFP="postVFP";

    // ---------- File paths ----------
    std::string f_HT0_100   = Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s_HT0-100/EleFakeRate-DATA-Bw-ker-vtx-60-120.txt",RunYear,whichVFP.c_str());
    std::string f_HT100_400 = Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s_HT100-400/EleFakeRate-DATA-Bw-ker-vtx-60-120.txt",RunYear,whichVFP.c_str());
    std::string f_HTgt400   = Form("/eos/uscms/store/user/tmishra/elefakepho/DATAResult%d%s_HT-gt400/EleFakeRate-DATA-Bw-ker-vtx-60-120.txt",RunYear,whichVFP.c_str());

    std::vector<std::string> files  = {f_HT0_100,f_HT100_400,f_HTgt400};
    std::vector<int> colors         = {kBlue,kRed,kGreen+2};
    std::vector<std::string> labels = {"HT 0-100","HT 100-400","HT >400"};

    std::vector<TGraphErrors*> graphs;
    std::vector<TF1*> fits;

    // ---------- Loop over HT files ----------
    for(int ifile=0;ifile<3;ifile++){

        std::ifstream fin(files[ifile]);
        if(!fin.is_open()){
            std::cout<<"Cannot open "<<files[ifile]<<std::endl;
            continue;
        }

        std::vector<double> x, den, num, denErr, numErr;

        std::string bintype,numtype;
        float lowcut, signal1, error1, fitmean, fitrms;

        std::string line;
        int nLines=0;
        while(std::getline(fin,line)) nLines++;
        fin.clear(); fin.seekg(0);

        int nBins = nLines/2;

        for(int i=0;i<nBins;i++){
            fin>>bintype>>numtype>>lowcut>>signal1>>error1>>fitmean>>fitrms;
            den.push_back(fitmean);
            denErr.push_back(fitrms);
            x.push_back(lowcut);
        }
        for(int i=0;i<nBins;i++){
            fin>>bintype>>numtype>>lowcut>>signal1>>error1>>fitmean>>fitrms;
            num.push_back(fitmean);
            numErr.push_back(fitrms);
        }

        TGraphErrors *gr = new TGraphErrors(nBins);

        for(int i=0;i<nBins;i++){
            double fr = num[i]/den[i];
            double er = sqrt( pow(numErr[i]/den[i],2)
                            + pow(num[i]*denErr[i]/(den[i]*den[i]),2) );

            gr->SetPoint(i,x[i],fr);
            gr->SetPointError(i,0,er);
        }

        gr->SetMarkerStyle(20);
        gr->SetMarkerColor(colors[ifile]);
        gr->SetLineColor(colors[ifile]);
        gr->SetLineWidth(2);

        // ---------- Straight line fit ----------
        TF1 *fit = new TF1(Form("fit_%d",ifile),"pol1",MINVTX,46);
        fit->SetLineColor(colors[ifile]);
        fit->SetParameters(0.02,0.0);

        TFitResultPtr r = gr->Fit(fit,"S");

        double p0     = fit->GetParameter(0);
        double p1     = fit->GetParameter(1);
        double p0err  = fit->GetParError(0);
        double p1err  = fit->GetParError(1);
        double chi2ndf= fit->GetChisquare()/fit->GetNDF();

        std::cout<<"===================================="<<std::endl;
        std::cout<<"Fit results : "<<labels[ifile]<<std::endl;
        std::cout<<"Intercept p0 = "<<p0<<" +/- "<<p0err<<std::endl;
        std::cout<<"Slope     p1 = "<<p1<<" +/- "<<p1err<<std::endl;
        std::cout<<"chi2/ndf     = "<<chi2ndf<<std::endl;

        graphs.push_back(gr);
        fits.push_back(fit);
    }

    // ---------- Plot ----------
    TCanvas *c = new TCanvas("c","",700,600);

    TH1D *dummy = new TH1D("dummy",";nVtx;fake rate",MAXVTX,0,MAXVTX);
    dummy->SetMaximum(0.07);
    dummy->Draw();

    TLegend *leg = new TLegend(0.45,0.7,0.78,0.88);

    for(size_t i=0; i<graphs.size(); i++){
        graphs[i]->Draw("EP same");
        fits[i]->Draw("same");
        leg->AddEntry(graphs[i],labels[i].c_str(),"pe");
    }

    leg->Draw();

    c->SaveAs(Form("/eos/uscms/store/user/tmishra/elefakepho/Plots/elefake_vtx_HTcompare_%d%s.png",
                   RunYear,whichVFP.c_str()));
}
