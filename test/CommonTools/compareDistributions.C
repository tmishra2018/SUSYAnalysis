#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

void compareDistributions(int RunYear, std::string whichVFP) {
    TFile *file_t5_tri = TFile::Open(Form("/uscms/home/tmishra/nobackup/signal_trees/signalTree_T5WG_%d%s.root", RunYear, whichVFP.c_str()));
    TFile *file_t5 = TFile::Open("/uscms/homes/m/mengleis/work/SUSY2016/SUSYAnalysis/test/Result/signalTree_T5WG.root");
    
    TH1D *h_MET_t5_tri = (TH1D*)file_t5_tri->Get("p_t5wg_MET_signal_1700_1000_mg");
    TH1D *h_HT_t5_tri = (TH1D*)file_t5_tri->Get("p_t5wg_HT_signal_1700_1000_mg");
    TH1D *h_PhoEt_t5_tri = (TH1D*)file_t5_tri->Get("p_t5wg_PhoEt_signal_1700_1000_mg");
    
    TH1D *h_MET_t5 = (TH1D*)file_t5->Get("p_t5wg_MET_signal_1700_1000_mg");
    TH1D *h_HT_t5 = (TH1D*)file_t5->Get("p_t5wg_HT_signal_1700_1000_mg");
    TH1D *h_PhoEt_t5 = (TH1D*)file_t5->Get("p_t5wg_PhoEt_signal_1700_1000_mg");
    
    if (!h_MET_t5_tri || !h_HT_t5_tri || !h_PhoEt_t5_tri || !h_MET_t5 || !h_HT_t5 || !h_PhoEt_t5) {
        printf("Error: One or more histograms not found!\n");
        return;
    }
    
    h_MET_t5_tri->SetLineColor(kRed);
    h_MET_t5->SetLineColor(kBlue);
    
    h_HT_t5_tri->SetLineColor(kRed);
    h_HT_t5->SetLineColor(kBlue);
    
    h_PhoEt_t5_tri->SetLineColor(kRed);
    h_PhoEt_t5->SetLineColor(kBlue);
    
    h_MET_t5_tri->GetYaxis()->SetRangeUser(0, 1.2 * std::max(h_MET_t5_tri->GetMaximum(), h_MET_t5->GetMaximum()));
    h_HT_t5_tri->GetYaxis()->SetRangeUser(0, 1.2 * std::max(h_HT_t5_tri->GetMaximum(), h_HT_t5->GetMaximum()));
    h_PhoEt_t5_tri->GetYaxis()->SetRangeUser(0, 1.2 * std::max(h_PhoEt_t5_tri->GetMaximum(), h_PhoEt_t5->GetMaximum()));
    
    TCanvas *c = new TCanvas(Form("c_%d_%s", RunYear, whichVFP.c_str()), Form("Comparison of Distributions %d %s", RunYear, whichVFP.c_str()), 1200, 400);
    c->Divide(3, 1);
    
    c->cd(1);
    h_MET_t5_tri->Draw("HIST");
    h_MET_t5->Draw("HIST SAME");
    TLegend *leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg1->AddEntry(h_MET_t5_tri, "T5 Tri", "l");
    leg1->AddEntry(h_MET_t5, "T5", "l");
    leg1->Draw();
    
    c->cd(2);
    h_HT_t5_tri->Draw("HIST");
    h_HT_t5->Draw("HIST SAME");
    TLegend *leg2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg2->AddEntry(h_HT_t5_tri, "T5 Tri", "l");
    leg2->AddEntry(h_HT_t5, "T5", "l");
    leg2->Draw();
    
    c->cd(3);
    h_PhoEt_t5_tri->Draw("HIST");
    h_PhoEt_t5->Draw("HIST SAME");
    TLegend *leg3 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg3->AddEntry(h_PhoEt_t5_tri, "T5 Tri", "l");
    leg3->AddEntry(h_PhoEt_t5, "T5", "l");
    leg3->Draw();
    
    c->SaveAs(Form("comparison_plot_%d_%s.png", RunYear, whichVFP.c_str()));
    file_t5_tri->Close();
    file_t5->Close();
}

