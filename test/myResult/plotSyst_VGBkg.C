// plotSyst_VGBkg.C
//
// Signal-bin plots use pre-computed syserr_* histograms from the ROOT file
// to construct up = nom + syserr and down = nom - syserr per bin.
// This is identical to the syserr/norm convention used in plot_eventct.C's
// LaTeX table, so the ratio panel directly reflects the table values.
//
// For QCDscale: syserr_qcdmuR and syserr_qcdmuF are enveloped bin-by-bin
//   (max of the two) before constructing up/down, matching plot_eventct.C
//   lines 798-800.  VGamma has only syserr_qcdmuF (no muR stored).
//
// Kinematic-variable plots (MET, HT, phoEt, MT) continue to use the raw
// up/down varied histograms — these are correct for shape visualisation.
//
// Y-axis shows absolute event yields (no area normalisation).
//
// Usage:
//   root -l -b -q 'plotSyst_VGBkg.C("eg.root","mg.root","2018")'
//   root -l -b -q 'plotSyst_VGBkg.C("eg.root","mg.root","2016preVFP","myplots")'

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <iomanip>

void setStyle(){
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    // Use Helvetica (font 42) everywhere
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleFont(42, "xyz");
}

TH1D* getH(TFile* f, const std::string& name){
    TH1D* h = (TH1D*)f->Get(name.c_str());
    if(!h){
        std::cerr << "WARNING: histogram \"" << name
                  << "\" not found in " << f->GetName()
                  << " — using empty placeholder." << std::endl;
        h = new TH1D(("_ph_"+name).c_str(), "", 10, 0, 1);
    }
    h->SetDirectory(nullptr);   // detach so it survives after file close
    return h;
}

std::string fmtPct(double frac){
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(1) << frac * 100.0 << "%";
    return ss.str();
}

// hsyserr: optional pre-computed symmetric uncertainty histogram (absolute yield units).
//   When provided (signal-bin plots), a 4th histogram is drawn:
//     upper panel: nom+syserr (green solid) and nom-syserr (violet solid)
//     ratio panel: 1 ± syserr/nom band
//   When nullptr (kinematic-variable plots), only Nominal/Up/Down are drawn.
void drawOne(TH1D* hnom,
             TH1D* hup,
             TH1D* hdown,
             TH1D* hsyserr,           // pass nullptr for kin-var plots
             const std::string& xTitle,
             const std::string& nomLabel,
             const std::string& upLabel,
             const std::string& doLabel,
             const std::string& syserrLabel, // label for 4th histogram (ignored if hsyserr==nullptr)
             const std::string& chLabel,
             const std::string& srcLabel,
             const std::string& era,
             const std::string& outPDF)
{
    TH1D* hn = (TH1D*)hnom ->Clone("_hn");
    TH1D* hu = (TH1D*)hup  ->Clone("_hu");
    TH1D* hd = (TH1D*)hdown->Clone("_hd");

    // No normalisation: keep absolute yields so ratio = Var/Nom in yield space.

    // 4th histograms: nom ± syserr (only when hsyserr provided)
    TH1D* hSysUp = nullptr;
    TH1D* hSysDo = nullptr;
    if(hsyserr){
        hSysUp = (TH1D*)hnom->Clone("_hSysUp");
        hSysDo = (TH1D*)hnom->Clone("_hSysDo");
        for(int ib = 1; ib <= hnom->GetNbinsX(); ++ib){
            double n = hnom->GetBinContent(ib);
            double s = hsyserr->GetBinContent(ib);
            hSysUp->SetBinContent(ib, n + s);
            hSysDo->SetBinContent(ib, std::max(0.0, n - s));
        }
    }

    TH1D* rU = (TH1D*)hn->Clone("_rU");
    TH1D* rD = (TH1D*)hn->Clone("_rD");
    rU->Reset(); rD->Reset();

    // Ratio histograms for syserr band (1 ± syserr/nom)
    TH1D* rSysUp = nullptr;
    TH1D* rSysDo = nullptr;
    if(hsyserr){
        rSysUp = (TH1D*)hn->Clone("_rSysUp");
        rSysDo = (TH1D*)hn->Clone("_rSysDo");
        rSysUp->Reset(); rSysDo->Reset();
    }

    double rMin =  9999.0;
    double rMax = -9999.0;
    double maxDevUp  = 0.0;
    double maxDevDo  = 0.0;
    double maxDevSys = 0.0;

    for(int ib = 1; ib <= hn->GetNbinsX(); ++ib){
        double nc = hn->GetBinContent(ib);
        if(nc > 1e-10){
            double ru = hu->GetBinContent(ib) / nc;
            double rd = hd->GetBinContent(ib) / nc;
            rU->SetBinContent(ib, ru);
            rD->SetBinContent(ib, rd);
            rU->SetBinError(ib, 0);
            rD->SetBinError(ib, 0);
            rMin = std::min({rMin, ru, rd});
            rMax = std::max({rMax, ru, rd});
            maxDevUp = std::max(maxDevUp, std::fabs(ru - 1.0));
            maxDevDo = std::max(maxDevDo, std::fabs(rd - 1.0));
            if(hsyserr && rSysUp && rSysDo){
                double s = hsyserr->GetBinContent(ib);
                double rSU = 1.0 + s / nc;
                double rSD = 1.0 - s / nc;
                rSysUp->SetBinContent(ib, rSU);
                rSysDo->SetBinContent(ib, rSD);
                rSysUp->SetBinError(ib, 0);
                rSysDo->SetBinError(ib, 0);
                rMin = std::min({rMin, rSU, rSD});
                rMax = std::max({rMax, rSU, rSD});
                maxDevSys = std::max(maxDevSys, s / nc);
            }
        } else {
            rU->SetBinContent(ib, 1.0);
            rD->SetBinContent(ib, 1.0);
            if(rSysUp){ rSysUp->SetBinContent(ib, 1.0); rSysDo->SetBinContent(ib, 1.0); }
        }
    }
    if(rMin > rMax){ rMin = 0.95; rMax = 1.05; }

    // 15% headroom, symmetric about 1, minimum ±5%
    double spread   = rMax - rMin;
    double headroom = std::max(0.04, spread * 0.15);
    double halfRange = std::max({std::fabs(rMax + headroom - 1.0),
                                 std::fabs(rMin - headroom - 1.0), 0.05});
    double rAxisMin = std::round((1.0 - halfRange) * 10.0) / 10.0;
    double rAxisMax = std::round((1.0 + halfRange) * 10.0) / 10.0;

    std::string concUp  = "Up: "     + fmtPct(maxDevUp)  + " max dev.";
    std::string concDo  = "Down: "   + fmtPct(maxDevDo)  + " max dev.";
    std::string concSys = "Syserr: " + fmtPct(maxDevSys) + " max dev.";

    // Colours and styles
    hn->SetFillColorAlpha(kGray+2, 0.40);
    hn->SetLineColor(kBlack);    hn->SetLineWidth(2);  hn->SetMarkerSize(0);
    hu->SetFillStyle(0);
    hu->SetLineColor(kRed+1);    hu->SetLineStyle(2);  hu->SetLineWidth(2);  hu->SetMarkerSize(0);
    hd->SetFillStyle(0);
    hd->SetLineColor(kAzure+2);  hd->SetLineStyle(2);  hd->SetLineWidth(2);  hd->SetMarkerSize(0);
    if(hSysUp){
        hSysUp->SetFillStyle(0);
        hSysUp->SetLineColor(kGreen+2); hSysUp->SetLineStyle(1); hSysUp->SetLineWidth(2); hSysUp->SetMarkerSize(0);
        hSysDo->SetFillStyle(0);
        hSysDo->SetLineColor(kViolet+1); hSysDo->SetLineStyle(1); hSysDo->SetLineWidth(2); hSysDo->SetMarkerSize(0);
    }
    rU->SetFillStyle(0);
    rU->SetLineColor(kRed+1);    rU->SetLineWidth(2);  rU->SetMarkerSize(0);
    rD->SetFillStyle(0);
    rD->SetLineColor(kAzure+2);  rD->SetLineWidth(2);  rD->SetMarkerSize(0);
    if(rSysUp){
        rSysUp->SetFillStyle(0);
        rSysUp->SetLineColor(kGreen+2);  rSysUp->SetLineWidth(2); rSysUp->SetMarkerSize(0);
        rSysDo->SetFillStyle(0);
        rSysDo->SetLineColor(kViolet+1); rSysDo->SetLineWidth(2); rSysDo->SetMarkerSize(0);
    }

    static int gCanvIdx = 0;
    std::string cname = Form("_c%d_", gCanvIdx++);
    TCanvas* c = new TCanvas(cname.c_str(), outPDF.c_str(), 720, 800);
    c->SetFillColor(0);

    // Upper pad (top 68% of canvas)
    TPad* pTop = new TPad("_pTop","", 0.0, 0.32, 1.0, 1.0);
    pTop->SetBottomMargin(0.018);
    pTop->SetLeftMargin(0.16);
    pTop->SetRightMargin(0.05);
    pTop->SetTopMargin(0.095);
    pTop->Draw();
    pTop->cd();

    double ymax = std::max({hn->GetMaximum(), hu->GetMaximum(), hd->GetMaximum()});
    if(hSysUp) ymax = std::max(ymax, hSysUp->GetMaximum());
    hn->GetYaxis()->SetRangeUser(0.0, ymax * 1.65);

    hn->GetYaxis()->SetTitle("Events / bin");
    hn->GetYaxis()->SetTitleFont(42);
    hn->GetYaxis()->SetTitleSize(0.060);
    hn->GetYaxis()->SetTitleOffset(1.15);
    hn->GetYaxis()->SetLabelFont(42);
    hn->GetYaxis()->SetLabelSize(0.054);
    hn->GetYaxis()->SetNdivisions(506);

    hn->GetXaxis()->SetLabelSize(0.0);
    hn->GetXaxis()->SetTitleSize(0.0);
    hn->GetXaxis()->SetTickLength(0.035);

    hn->Draw("HIST");
    hu->Draw("HIST SAME");
    hd->Draw("HIST SAME");
    if(hSysUp){ hSysUp->Draw("HIST SAME"); hSysDo->Draw("HIST SAME"); }

    // Legend — shrink box to fit 4 entries when needed
    double legY1 = hsyserr ? 0.42 : 0.52;
    TLegend* leg = new TLegend(0.57, legY1, 0.93, 0.80);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.040);
    leg->AddEntry(hn, nomLabel.c_str(), "f");
    leg->AddEntry(hu, upLabel.c_str(),  "l");
    leg->AddEntry(hd, doLabel.c_str(),  "l");
    if(hSysUp){
        leg->AddEntry(hSysUp, (syserrLabel + " [Nom +]").c_str(), "l");
        leg->AddEntry(hSysDo, (syserrLabel + " [Nom #minus]").c_str(), "l");
    }
    leg->Draw();

    // Labels in upper pad
    TLatex lat;
    lat.SetNDC();
    // "CMS" bold
    lat.SetTextFont(62); lat.SetTextSize(0.058); lat.SetTextColor(kBlack);
    lat.DrawLatex(0.17, 0.920, "CMS");
    // "Simulation" italic
    lat.SetTextFont(52); lat.SetTextSize(0.044);
    lat.DrawLatex(0.29, 0.920, "Simulation");
    // Channel + era — top right
    lat.SetTextFont(42); lat.SetTextSize(0.042);
    //lat.DrawLatex(0.48, 0.920, Form("  V#gamma  %s |  (Run 2) ", chLabel.c_str()));
    lat.DrawLatex(0.48, 0.920, Form(" electroweak  %s |  (Run 2) ", chLabel.c_str()));

    // Uncertainty source label below CMS
    lat.SetTextFont(42); lat.SetTextSize(0.046);
    lat.DrawLatex(0.17, 0.845, Form("Source: %s", srcLabel.c_str()));

    // Lower pad: ratio Var/Nom (bottom 32% of canvas)
    c->cd();
    TPad* pBot = new TPad("_pBot","", 0.0, 0.0, 1.0, 0.32);
    pBot->SetTopMargin(0.012);
    pBot->SetBottomMargin(0.40);    // large bottom margin for x-axis title
    pBot->SetLeftMargin(0.16);
    pBot->SetRightMargin(0.05);
    pBot->Draw();
    pBot->cd();

    // Y axis — ratio pad
    rU->GetYaxis()->SetRangeUser(rAxisMin, rAxisMax);
    rU->GetYaxis()->SetTitle("Var / Nom");
    rU->GetYaxis()->SetTitleFont(42);
    rU->GetYaxis()->SetTitleSize(0.148);
    rU->GetYaxis()->SetTitleOffset(0.40);
    rU->GetYaxis()->SetLabelFont(42);
    rU->GetYaxis()->SetLabelSize(0.130);
    rU->GetYaxis()->SetNdivisions(504);
    rU->GetYaxis()->SetTickLength(0.022);

    rU->GetXaxis()->SetTitle(xTitle.c_str());
    rU->GetXaxis()->SetTitleFont(42);
    rU->GetXaxis()->SetTitleSize(0.155);
    rU->GetXaxis()->SetTitleOffset(1.00);
    rU->GetXaxis()->SetLabelFont(42);
    rU->GetXaxis()->SetLabelSize(0.11);
    rU->GetXaxis()->SetLabelOffset(0.012);
    rU->GetXaxis()->SetTickLength(0.060);

    rU->Draw("HIST");
    rD->Draw("HIST SAME");
    if(rSysUp){ rSysUp->Draw("HIST SAME"); rSysDo->Draw("HIST SAME"); }

    // Reference and guide lines (only draw if within visible Y range)
    double xlo = rU->GetXaxis()->GetXmin();
    double xhi = rU->GetXaxis()->GetXmax();

    auto drawLine = [&](double y, int col, int style, int width){
        if(y >= rAxisMin - 1e-6 && y <= rAxisMax + 1e-6){
            TLine* l = new TLine(xlo, y, xhi, y);
            l->SetLineColor(col); l->SetLineStyle(style); l->SetLineWidth(width);
            l->Draw();
        }
    };
    drawLine(1.00, kBlack,  2, 2);
    drawLine(1.05, kGray+1, 3, 1);
    drawLine(0.95, kGray+1, 3, 1);
    drawLine(1.10, kGray,   3, 1);
    drawLine(0.90, kGray,   3, 1);

    // Conclusion text — 3 lines when syserr present, 2 otherwise.
    // Positions shift upward to stay within bottom margin space.
    TLatex rLat;
    rLat.SetNDC();
    rLat.SetTextFont(42);
    rLat.SetTextSize(0.115);

    if(hsyserr){
        rLat.SetTextColor(kRed+1);
        rLat.DrawLatex(0.17, 0.185, concUp.c_str());
        rLat.SetTextColor(kAzure+2);
        rLat.DrawLatex(0.17, 0.100, concDo.c_str());
        rLat.SetTextColor(kGreen+2);
        rLat.DrawLatex(0.17, 0.015, concSys.c_str());
    } else {
        rLat.SetTextColor(kRed+1);
        rLat.DrawLatex(0.17, 0.128, concUp.c_str());
        rLat.SetTextColor(kAzure+2);
        rLat.DrawLatex(0.17, 0.028, concDo.c_str());
    }

    // Save
    c->SaveAs(outPDF.c_str());
    std::cout << "  saved: " << outPDF << std::endl;

    // Clean up clones
    delete hn; delete hu; delete hd;
    if(hSysUp){ delete hSysUp; delete hSysDo; }
    delete rU; delete rD;
    if(rSysUp){ delete rSysUp; delete rSysDo; }
    delete c;
}
void drawOneISR(TH1D* hnom,          // nominal histogram
                TH1D* hisrup,        // ISR-up varied histogram (isrup_*)
                TH1D* hsyserr,       // absolute syserr histogram (syserr_isr or nullptr)
                const std::string& xTitle,
                const std::string& nomLabel,
                const std::string& isrLabel,
                const std::string& chLabel,
                const std::string& srcLabel,
                const std::string& era,
                const std::string& outPDF)
{
    TH1D* hn  = (TH1D*)hnom  ->Clone("_isrHn");
    TH1D* hu  = (TH1D*)hisrup->Clone("_isrHu");
 
    // No normalisation: keep absolute yields so ratio = Var/Nom in yield space.
 
    TH1D* rU = (TH1D*)hn->Clone("_isrRu");
    rU->Reset();
 
    TH1D* rBandUp = (TH1D*)hn->Clone("_isrBu");
    TH1D* rBandDo = (TH1D*)hn->Clone("_isrBd");
    rBandUp->Reset(); rBandDo->Reset();
 
    double rMin = 9999.0, rMax = -9999.0;
    double maxDevUp = 0.0;
 
    for(int ib = 1; ib <= hn->GetNbinsX(); ++ib){
        double nc = hn->GetBinContent(ib);
        if(nc > 1e-10){
            double ru = hu->GetBinContent(ib) / nc;
            rU->SetBinContent(ib, ru);
            rU->SetBinError(ib, 0);
            rMin = std::min(rMin, ru);
            rMax = std::max(rMax, ru);
            maxDevUp = std::max(maxDevUp, std::fabs(ru - 1.0));
            if(hsyserr){
                double relErr = hsyserr->GetBinContent(ib) / nc;
                rBandUp->SetBinContent(ib, 1.0 + relErr);
                rBandDo->SetBinContent(ib, 1.0 - relErr);
            } else {
                rBandUp->SetBinContent(ib, 1.0);
                rBandDo->SetBinContent(ib, 1.0);
            }
        } else {
            rU->SetBinContent(ib, 1.0);
            rBandUp->SetBinContent(ib, 1.0);
            rBandDo->SetBinContent(ib, 1.0);
        }
    }
    if(rMin > rMax){ rMin = 0.95; rMax = 1.05; }
 
    // Also account for the syserr band extent
    double bandMax = rBandUp->GetMaximum();
    double bandMin = rBandDo->GetMinimum();
    rMax = std::max(rMax, bandMax);
    rMin = std::min(rMin, bandMin);
 
    double spread   = rMax - rMin;
    double headroom = std::max(0.04, spread * 0.15);
    double halfRange = std::max(std::fabs(rMax + headroom - 1.0),
                                std::fabs(rMin - headroom - 1.0));
    halfRange = std::max(halfRange, 0.05);
    double rAxisMin = std::round((1.0 - halfRange) * 10.0) / 10.0;
    double rAxisMax = std::round((1.0 + halfRange) * 10.0) / 10.0;

    std::string concUp = fmtPct(maxDevUp) + " max dev.";
 
    hn->SetFillColorAlpha(kGray+2, 0.40);
    hn->SetLineColor(kBlack);  hn->SetLineWidth(2); hn->SetMarkerSize(0);
    hu->SetFillStyle(0);
    hu->SetLineColor(kOrange+7); hu->SetLineStyle(2); hu->SetLineWidth(2); hu->SetMarkerSize(0);
    rU->SetFillStyle(0);
    rU->SetLineColor(kOrange+7); rU->SetLineWidth(2); rU->SetMarkerSize(0);
    // Yellow band for syserr
    rBandUp->SetFillColorAlpha(kYellow-9, 0.70);
    rBandUp->SetLineColor(0); rBandUp->SetMarkerSize(0);
    rBandDo->SetFillColorAlpha(kYellow-9, 0.70);
    rBandDo->SetLineColor(0); rBandDo->SetMarkerSize(0);
 
    static int gISRIdx = 0;
    std::string cname = Form("_cisr%d_", gISRIdx++);
    TCanvas* c = new TCanvas(cname.c_str(), outPDF.c_str(), 720, 800);
    c->SetFillColor(0);
 
    TPad* pTop = new TPad("_isrTop","", 0.0, 0.32, 1.0, 1.0);
    pTop->SetBottomMargin(0.018);
    pTop->SetLeftMargin(0.16);
    pTop->SetRightMargin(0.05);
    pTop->SetTopMargin(0.095);
    pTop->Draw(); pTop->cd();
 
    double ymax = std::max(hn->GetMaximum(), hu->GetMaximum());
    hn->GetYaxis()->SetRangeUser(0.0, ymax * 1.65);
    hn->GetYaxis()->SetTitle("Events / bin");
    hn->GetYaxis()->SetTitleFont(42); hn->GetYaxis()->SetTitleSize(0.060);
    hn->GetYaxis()->SetTitleOffset(1.15);
    hn->GetYaxis()->SetLabelFont(42); hn->GetYaxis()->SetLabelSize(0.054);
    hn->GetYaxis()->SetNdivisions(506);
    hn->GetXaxis()->SetLabelSize(0.0);
    hn->GetXaxis()->SetTitleSize(0.0);
    hn->GetXaxis()->SetTickLength(0.035);
    hn->Draw("HIST");
    hu->Draw("HIST SAME");
 
    TLegend* leg = new TLegend(0.50, 0.62, 0.93, 0.86);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->SetTextFont(42);  leg->SetTextSize(0.048);
    leg->AddEntry(hn, nomLabel.c_str(), "f");
    leg->AddEntry(hu, isrLabel.c_str(), "l");
    if(hsyserr){
        TH1D* dummy = (TH1D*)hn->Clone("_isrDummy");
        dummy->SetFillColorAlpha(kYellow-9, 0.70); dummy->SetLineColor(0);
        leg->AddEntry(dummy, "Syst. unc. band", "f");
    }
    leg->Draw();
 
    TLatex lat; lat.SetNDC();
    lat.SetTextFont(62); lat.SetTextSize(0.058); lat.SetTextColor(kBlack);
    lat.DrawLatex(0.17, 0.920, "CMS");
    lat.SetTextFont(52); lat.SetTextSize(0.044);
    lat.DrawLatex(0.29, 0.920, "Simulation");
    lat.SetTextFont(42); lat.SetTextSize(0.042);
    //lat.DrawLatex(0.48, 0.920, Form(" V#gamma   %s |  (Run 2) ", chLabel.c_str()));
    lat.DrawLatex(0.48, 0.920, Form(" electroweak  %s |  (Run 2) ", chLabel.c_str()));
    lat.SetTextSize(0.046);
    lat.DrawLatex(0.17, 0.845, Form("Source: %s", srcLabel.c_str()));
 
    c->cd();
    TPad* pBot = new TPad("_isrBot","", 0.0, 0.0, 1.0, 0.32);
    pBot->SetTopMargin(0.012);
    pBot->SetBottomMargin(0.40);
    pBot->SetLeftMargin(0.16);
    pBot->SetRightMargin(0.05);
    pBot->Draw(); pBot->cd();
 
    rU->GetYaxis()->SetRangeUser(rAxisMin, rAxisMax);
    rU->GetYaxis()->SetTitle("no ISR / Nom");
    rU->GetYaxis()->SetTitleFont(42); rU->GetYaxis()->SetTitleSize(0.148);
    rU->GetYaxis()->SetTitleOffset(0.40);
    rU->GetYaxis()->SetLabelFont(42); rU->GetYaxis()->SetLabelSize(0.130);
    rU->GetYaxis()->SetNdivisions(504);
    rU->GetYaxis()->SetTickLength(0.022);
    rU->GetXaxis()->SetTitle(xTitle.c_str());
    rU->GetXaxis()->SetTitleFont(42); rU->GetXaxis()->SetTitleSize(0.155);
    rU->GetXaxis()->SetTitleOffset(1.00);
    rU->GetXaxis()->SetLabelFont(42); rU->GetXaxis()->SetLabelSize(0.138);
    rU->GetXaxis()->SetLabelOffset(0.012);
    rU->GetXaxis()->SetTickLength(0.060);
 
    // Draw syserr band first (behind ratio line)
    rU->Draw("HIST");   // sets the axes frame
    if(hsyserr){
        // Build a proper filled band between rBandDo and rBandUp
        // using TH1 error bars as band: set bin content = mid, error = half-width
        TH1D* hBand = (TH1D*)rBandUp->Clone("_isrBand");
        for(int ib=1; ib<=hBand->GetNbinsX(); ++ib){
            double up = rBandUp->GetBinContent(ib);
            double dn = rBandDo->GetBinContent(ib);
            double mid = 0.5*(up+dn);
            double err = 0.5*(up-dn);
            hBand->SetBinContent(ib, mid);
            hBand->SetBinError(ib, err);
        }
        hBand->SetFillColorAlpha(kYellow-9, 0.70);
        hBand->SetLineColor(0);
        hBand->SetMarkerSize(0);
        hBand->Draw("E2 SAME");   // E2 = filled error band
    }
    rU->Draw("HIST SAME");   // ratio line on top of band
 
    // Reference lines
    double xlo = rU->GetXaxis()->GetXmin();
    double xhi = rU->GetXaxis()->GetXmax();
    auto drawLine = [&](double y, int col, int style, int width){
        if(y >= rAxisMin-1e-6 && y <= rAxisMax+1e-6){
            TLine* l = new TLine(xlo,y,xhi,y);
            l->SetLineColor(col); l->SetLineStyle(style); l->SetLineWidth(width);
            l->Draw();
        }
    };
    drawLine(1.00, kBlack,  2, 2);
    drawLine(1.05, kGray+1, 3, 1);
    drawLine(0.95, kGray+1, 3, 1);
    drawLine(1.10, kGray,   3, 1);
    drawLine(0.90, kGray,   3, 1);
 
    // Conclusion text
    TLatex rLat; rLat.SetNDC();
    rLat.SetTextFont(42); rLat.SetTextSize(0.115);
    rLat.SetTextColor(kOrange+7);
    rLat.DrawLatex(0.17, 0.128, concUp.c_str());

    c->SaveAs(outPDF.c_str());
    std::cout << "  saved: " << outPDF << std::endl;
 
    delete hn; delete hu; delete rU;
    delete rBandUp; delete rBandDo;
    delete c;
}
// ---------------------------------------------------------------------------
// Helper: extract background type token from a ROOT filename.
// Looks for "VGBkg" → returns "VGamma", or "rareBkg" → returns "rare".
// Falls back to "VGamma" if neither is found.
// ---------------------------------------------------------------------------
std::string extractBkgType(const std::string& fname){
    if(fname.find("rareBkg") != std::string::npos) return "rare";
    if(fname.find("VGBkg")   != std::string::npos) return "VGamma";
    std::cerr << "WARNING: Cannot determine bkg type from filename \""
              << fname << "\". Defaulting to \"VGamma\"." << std::endl;
    return "VGamma";
}

// Main entry point
void plotSyst_VGBkg(const char* egFile  = "signalTree_egamma_VGBkg_2018.root",
                    const char* mgFile  = "signalTree_mg_VGBkg_2018.root",
                    const char* era     = "2018",
                    const char* outdir  = nullptr)
{
    setStyle();

    // --- Determine background type from the input filename ---
    std::string bkgType = extractBkgType(std::string(mgFile));
    std::cout << "Detected background type: " << bkgType << std::endl;

    // --- Build output directory name if not provided ---
    std::string outdirStr;
    if(outdir && std::string(outdir).size() > 0){
        outdirStr = outdir;
    } else {
        outdirStr = std::string("syst_plots_") + bkgType + "_" + era;
    }
    gSystem->mkdir(outdirStr.c_str(), true);

    TFile* fEG = TFile::Open(egFile, "READ");
    TFile* fMG = TFile::Open(mgFile, "READ");
    if(!fEG || fEG->IsZombie()){ std::cerr << "Cannot open " << egFile << std::endl; return; }
    if(!fMG || fMG->IsZombie()){ std::cerr << "Cannot open " << mgFile << std::endl; return; }
    std::cout << "Opened:\n  " << egFile << "\n  " << mgFile << std::endl;

    struct Chan {
        TFile*      f;
        std::string prefix;  // histogram name prefix in file ("eg" or "mg")
        std::string label;   // TLatex string for plot title
        std::string tag;     // used in output filename
    };
    std::vector<Chan> channels = {
        {fEG, "eg", "e#gamma",   "eg"},
        {fMG, "mg", "#mu#gamma", "mg"}
    };

    struct Source {
        std::string tag;
        std::string label;         // shown on canvas
        std::string nomLabel;      // legend: nominal
        std::string upLabel;       // legend: up variation
        std::string doLabel;       // legend: down variation
        std::string syserrLabel;   // legend: assigned uncertainty (4th histogram pair)
        // Histogram names in the ROOT file for each observable
        std::string met_up,    met_do;
	std::string mt_up,    mt_do;
        std::string ht_up,     ht_do;
        std::string phoEt_up,  phoEt_do;
        // Suffix for signal-bin histograms: {ch}_{bkgType}_{suffix}
        std::string sigBin_up_suffix, sigBin_do_suffix;
    };

    std::vector<Source> sources = {
        {
            "PU",
            "PU reweighting (#pm4.6% #sigma_{mb})",
            "Nominal",
            "PU Up (+4.6%)",
            "PU Down (#minus4.6%)",
            "Assigned unc.",
            "puup_MET",      "pudo_MET",
            "puup_Mt",      "pudo_Mt",
            "puup_HT",       "pudo_HT",
            "puup_PhoEt",    "pudo_PhoEt",
            "puUp",          "puDown"
        },
        {
            "JEC",
            "Jet Energy Corrections (JEC)",
            "Nominal",
            "JEC Up (+1#sigma)",
            "JEC Down (#minus1#sigma)",
            "Assigned unc.",
            "jesup_MET",     "jesdo_MET",
            "jesup_Mt",     "jesdo_Mt",
            "jesup_HT",      "jesdo_HT",
            "p_PhoEt",       "p_PhoEt",
            "jesUp",         "jesDown"
        },
        {
            "JER",
            "Jet Energy Resolution (JER)",
            "Nominal",
            "JER Up (+1#sigma)",
            "JER Down (#minus1#sigma)",
            "Assigned unc.",
            "jerup_MET",     "jerdo_MET",
            "jerup_Mt",     "jerdo_Mt",
            "p_HT",          "p_HT",
            "p_PhoEt",       "p_PhoEt",
            "jerUp",         "jerDown"
        },
        {
            "PSW_ISR",
            "PS weight: ISR (#mu_{R} #times 2 / 0.5)",
            "Nominal",
            "ISR Up (#mu_{R} #times 2)",
            "ISR Down (#mu_{R} #times 0.5)",
            "Assigned unc.",
            "pswIsrUp_MET",  "pswIsrDo_MET",
            "pswIsrUp_Mt",  "pswIsrDo_Mt",
            "pswIsrUp_HT",   "pswIsrDo_HT",
            "pswIsrUp_PhoEt","pswIsrDo_PhoEt",
            "pswIsrUp",      "pswIsrDown"
        },
        {
            "PSW_FSR",
            "PS weight: FSR (#mu_{R} #times 2 / 0.5)",
            "Nominal",
            "FSR Up (#mu_{R} #times 2)",
            "FSR Down (#mu_{R} #times 0.5)",
            "Assigned unc.",
            "pswFsrUp_MET",  "pswFsrDo_MET",
            "pswFsrUp_Mt",  "pswFsrDo_Mt",
            "pswFsrUp_HT",   "pswFsrDo_HT",
            "pswFsrUp_PhoEt","pswFsrDo_PhoEt",
            "pswFsrUp",      "pswFsrDown"
        },
        {
            "QCDscale_muF",
            "QCD scale (#mu_{F} variation)",
            "Nominal",
            "QCD #mu_{F} Up (#times 2)",
            "QCD #mu_{F} Down (#times 0.5)",
            "Assigned unc.",
            "muFup_MET",   "muFdo_MET",
            "muFup_Mt",    "muFdo_Mt",
            "muFup_HT",    "muFdo_HT",
            "muFup_PhoEt", "muFdo_PhoEt",
            "muFup",       "muFdown"
        },
        {
            "QCDscale_muR",
            "QCD scale (#mu_{R} variation)",
            "Nominal",
            "QCD #mu_{R} Up (#times 2)",
            "QCD #mu_{R} Down (#times 0.5)",
            "Assigned unc.",
            "muRup_MET",   "muRdo_MET",
            "muRup_Mt",    "muRdo_Mt",
            "muRup_HT",    "muRdo_HT",
            "muRup_PhoEt", "muRdo_PhoEt",
            "muRup",       "muRdown"
        },
        {
            "PDF",
            "PDF uncertainty (Hessian RMS)",
            "Nominal",
            "PDF Up (+1#sigma)",
            "PDF Down (#minus1#sigma)",
            "Assigned unc.",
            "pdfup_MET",     "pdfdo_MET",
            "pdfup_Mt",      "pdfdo_Mt",
            "pdfup_HT",      "pdfdo_HT",
            "pdfup_PhoEt",   "pdfdo_PhoEt",
            "pdfUp",         "pdfDown"
        }
    };

    // nomHistName, full x-axis label with units, filename tag
    struct Var {
        std::string nom;      // histogram name in file
        std::string xTitle;   // x-axis label (with LaTeX + units)
        std::string tag;      // used in output filename
    };
    std::vector<Var> kinVars = {
        {"p_MET",   "p_{T}^{miss} (GeV)",         "MET"  },
	{"p_Mt",    "M_{T} (GeV)",           "MT"   },
        {"p_HT",    "H_{T} (GeV)",                 "HT"   },
        {"p_PhoEt", "p_{T}^{#gamma} (GeV)",        "phoEt"}
    };

    int nPlots = 0;

    // Main loop:
    //   5 sources × 2 channels × (3 kin vars + 1 sigbin) = 35 plots
    for(const auto& src : sources){
        std::cout << "\n--- Source: " << src.tag << " ---" << std::endl;
        for(const auto& ch : channels){
            for(const auto& var : kinVars){
                std::string upName, doName;
                if     (var.tag == "MET")  { upName=src.met_up;   doName=src.met_do;   }
		else if(var.tag == "MT")   { upName=src.mt_up;    doName=src.mt_do;    }
                else if(var.tag == "HT")   { upName=src.ht_up;    doName=src.ht_do;    }
                else                       { upName=src.phoEt_up; doName=src.phoEt_do; }

                TH1D* hnom = getH(ch.f, var.nom);
                TH1D* hup  = getH(ch.f, upName);
                TH1D* hdo  = getH(ch.f, doName);

                std::string outPDF = outdirStr + "/" + bkgType + "_" +
                    ch.tag + "_" + var.tag + "_" + src.tag + "_" + era + ".png";

                drawOne(hnom, hup, hdo,
                        nullptr,          // no syserr for kinematic-variable plots
                        var.xTitle,
                        src.nomLabel, src.upLabel, src.doLabel,
                        "",               // syserrLabel unused when nullptr
                        ch.label, src.label, era, outPDF);
                ++nPlots;
            }

            // --- Signal-bin plot ---
            // Draws 4 histograms:
            //   Nominal, Up (raw varied), Down (raw varied), ± syserr band
            // The raw up/down show the actual variation shape; the syserr band
            // shows the pre-computed symmetric uncertainty used in the LaTeX table.

            std::string nomName = ch.prefix + "_" + bkgType + "_norm";
            std::string upName  = ch.prefix + "_" + bkgType + "_" + src.sigBin_up_suffix;
            std::string doName  = ch.prefix + "_" + bkgType + "_" + src.sigBin_do_suffix;

            TH1D* hnom_sigbin = getH(ch.f, nomName);
            TH1D* hup_sigbin  = getH(ch.f, upName);
            TH1D* hdo_sigbin  = getH(ch.f, doName);

            if(!hnom_sigbin){ ++nPlots; continue; }

            // Load syserr histogram(s); envelope muR+muF for QCDscale
            TH1D* hSysEnv = nullptr;
            bool isQCDscale = (src.tag == "QCDscale_muF" || src.tag == "QCDscale_muR");

            if(isQCDscale){
                std::string muFname = ch.prefix + "_" + bkgType + "_syserr_qcdmuF";
                std::string muRname = ch.prefix + "_" + bkgType + "_syserr_qcdmuR";
                TH1D* hmuF = (TH1D*)ch.f->Get(muFname.c_str());
                TH1D* hmuR = (TH1D*)ch.f->Get(muRname.c_str());
                if(hmuF){
                    hSysEnv = (TH1D*)hmuF->Clone("_sysEnvQCD");
                    if(hmuR){
                        for(int ib = 1; ib <= hSysEnv->GetNbinsX(); ++ib)
                            hSysEnv->SetBinContent(ib,
                                std::max(hmuF->GetBinContent(ib),
                                         hmuR->GetBinContent(ib)));
                    }
                }
            } else if(src.tag == "PDF"){
                TH1D* h = (TH1D*)ch.f->Get((ch.prefix+"_"+bkgType+"_syserr_pdf").c_str());
                if(h) hSysEnv = (TH1D*)h->Clone("_sysEnv");
            } else if(src.tag == "PU"){
                TH1D* h = (TH1D*)ch.f->Get((ch.prefix+"_"+bkgType+"_syserr_pu").c_str());
                if(h) hSysEnv = (TH1D*)h->Clone("_sysEnv");
            } else if(src.tag == "JEC"){
                TH1D* h = (TH1D*)ch.f->Get((ch.prefix+"_"+bkgType+"_syserr_jes").c_str());
                if(h) hSysEnv = (TH1D*)h->Clone("_sysEnv");
            } else if(src.tag == "JER"){
                TH1D* h = (TH1D*)ch.f->Get((ch.prefix+"_"+bkgType+"_syserr_jer").c_str());
                if(h) hSysEnv = (TH1D*)h->Clone("_sysEnv");
            } else if(src.tag == "PSW_ISR"){
                TH1D* h = (TH1D*)ch.f->Get((ch.prefix+"_"+bkgType+"_syserr_isr").c_str());
                if(h) hSysEnv = (TH1D*)h->Clone("_sysEnv");
            } else if(src.tag == "PSW_FSR"){
                TH1D* h = (TH1D*)ch.f->Get((ch.prefix+"_"+bkgType+"_syserr_fsr").c_str());
                if(h) hSysEnv = (TH1D*)h->Clone("_sysEnv");
            }

            // syserrLabel shown in the legend for the 4th histogram pair
            std::string outPDF = outdirStr + "/" + bkgType + "_" +
                ch.tag + "_SignalBin_" + src.tag + "_" + era + ".png";

            drawOne(hnom_sigbin, hup_sigbin, hdo_sigbin,
                    hSysEnv,                // 4th histogram: pre-computed syserr
                    "Signal Bin Index",
                    src.nomLabel, src.upLabel, src.doLabel,
                    src.syserrLabel,        // per-source physically meaningful label
                    ch.label, src.label, era, outPDF);

            if(hSysEnv) delete hSysEnv;
            ++nPlots;
        }
    }
    std::cout << "\n--- Source: ISR_ALT (single-sided) ---" << std::endl;

    // ISR kinematic variable definitions: {nomName, isrupName, xTitle, tag}
    struct ISRVar {
        std::string nom;
        std::string isrup;
        std::string xTitle;
        std::string tag;
    };
    std::vector<ISRVar> isrVars = {
        {"p_MET",   "isrup_MET",   "p_{T}^{miss} (GeV)",  "MET"  },
        {"p_Mt",    "isrup_Mt",    "M_{T} (GeV)",          "MT"   },
        {"p_HT",    "isrup_HT",    "H_{T} (GeV)",          "HT"   },
        {"p_PhoEt", "isrup_PhoEt", "p_{T}^{#gamma} (GeV)", "phoEt"}
    };

    for(const auto& ch : channels){
        for(const auto& iv : isrVars){
            TH1D* hnom   = getH(ch.f, iv.nom);
            TH1D* hisrup = getH(ch.f, iv.isrup);
            TH1D* hsyserr = nullptr;

            std::string outPDF = outdirStr + "/" + bkgType + "_" +
                ch.tag + "_" + iv.tag + "_ISR_" + era + ".png";

            drawOneISR(hnom, hisrup, hsyserr,
                       iv.xTitle,
                       "ISR weighted (Nominal)",
                       "No ISR reweighting",
                       ch.label,
                       "ISR reweighting",
                       era, outPDF);
            ++nPlots;
        }

        std::string nomName    = ch.prefix + "_" + bkgType + "_norm";
        std::string isrupName  = ch.prefix + "_" + bkgType + "_isrAlter";
        std::string syserrName = ch.prefix + "_" + bkgType + "_syserr_isr";

        TH1D* hnom    = getH(ch.f, nomName);
        TH1D* hisrup  = getH(ch.f, isrupName);
        TH1D* hsyserr = getH(ch.f, syserrName);

        std::string outPDF = outdirStr + "/" + bkgType + "_" + ch.tag + "_SignalBin_ISR_" + era + ".png";

        drawOneISR(hnom, hisrup, hsyserr,
                   "Signal Bin Index",
                   "ISR weighted (Nominal)",
                   "No ISR reweighting",
                   ch.label,
                   "ISR reweighting",
                   era, outPDF);
        ++nPlots;
    }

    fEG->Close();
    fMG->Close();

    std::cout << "\n========================================" << std::endl;
    std::cout << "  Done:  " << nPlots << " plots  (expected 49)" << std::endl;
    std::cout << "  Saved to: " << outdirStr << "/" << std::endl;
    std::cout << "  Pattern: " << bkgType << "_{ch}_{var}_{src}_" << era << ".png" << std::endl;
    std::cout << "  ch  : eg, mg" << std::endl;
    std::cout << "  var : MET, MT, HT, phoEt, SignalBin" << std::endl;
    std::cout << "  src : PU, JEC, JER, PSW_ISR, PSW_FSR, QCDscale, PDF" << std::endl;
    std::cout << "========================================" << std::endl;
}
