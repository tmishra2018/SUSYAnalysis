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
#include "TChain.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "../../include/tdrstyle.C"
int RunYear = 2016;
bool preVFP = true;
bool getPUweight = false;

float getPUESF(int nvertex){
	
	float pileupweights[100];
    if(RunYear == 2016 and preVFP == true){
        pileupweights[0] = 0;
        pileupweights[1] = 0;
        pileupweights[2] = 1.10726;
        pileupweights[3] = 0.851737;
        pileupweights[4] = 0.616496;
        pileupweights[5] = 0.711674;
        pileupweights[6] = 0.591054;
        pileupweights[7] = 0.834372;
        pileupweights[8] = 0.770663;
        pileupweights[9] = 0.828467;
        pileupweights[10] = 0.908953;
        pileupweights[11] = 0.963917;
        pileupweights[12] = 0.970074;
        pileupweights[13] = 1.04855;
        pileupweights[14] = 1.0199;
        pileupweights[15] = 1.03302;
        pileupweights[16] = 1.06986;
        pileupweights[17] = 1.00919;
        pileupweights[18] = 1.03648;
        pileupweights[19] = 1.05633;
        pileupweights[20] = 1.07702;
        pileupweights[21] = 0.987433;
        pileupweights[22] = 0.987196;
        pileupweights[23] = 0.973624;
        pileupweights[24] = 0.947854;
        pileupweights[25] = 0.952608;
        pileupweights[26] = 0.942894;
        pileupweights[27] = 0.944227;
        pileupweights[28] = 0.919677;
        pileupweights[29] = 0.866373;
        pileupweights[30] = 0.952522;
        pileupweights[31] = 0.998731;
        pileupweights[32] = 1.06274;
        pileupweights[33] = 1.06342;
        pileupweights[34] = 1.04084;
        pileupweights[35] = 1.23387;
        pileupweights[36] = 1.52291;
        pileupweights[37] = 1.15981;
        pileupweights[38] = 0.929168;
        pileupweights[39] = 2.17193;
        pileupweights[40] = 2.10969;
        pileupweights[41] = 1.68928;
        pileupweights[42] = 2.00158;
        pileupweights[43] = 5.33755;
        pileupweights[44] = 1.34149;
        pileupweights[45] = 2.83912;
        pileupweights[46] = 7.15459;
        pileupweights[47] = 2.21452;
        pileupweights[48] = 0;
        pileupweights[49] = 0;
        pileupweights[50] = 0;
        pileupweights[51] = 4.42903;
        pileupweights[52] = 1.36278;
        pileupweights[53] = 0;
        pileupweights[54] = 0;
        pileupweights[55] = 0;
        pileupweights[56] = 0;
        pileupweights[57] = 0;
        pileupweights[58] = 0;
        pileupweights[59] = 0;
        pileupweights[60] = 0;
        pileupweights[61] = 0;
        pileupweights[62] = 0;
        pileupweights[63] = 0;
        pileupweights[64] = 0;
        pileupweights[65] = 0;
        pileupweights[66] = 0;
        pileupweights[67] = 0;
        pileupweights[68] = 0;
        pileupweights[69] = 0;
        pileupweights[70] = 0;
        pileupweights[71] = 0;
        pileupweights[72] = 0;
        pileupweights[73] = 0;
        pileupweights[74] = 0;
        pileupweights[75] = 0;
        pileupweights[76] = 0;
        pileupweights[77] = 0;
        pileupweights[78] = 0;
        pileupweights[79] = 0;
        pileupweights[80] = 0;
        pileupweights[81] = 0;
        pileupweights[82] = 0;
        pileupweights[83] = 0;
        pileupweights[84] = 0;
        pileupweights[85] = 0;
        pileupweights[86] = 0;
        pileupweights[87] = 0;
        pileupweights[88] = 0;
        pileupweights[89] = 0;
        pileupweights[90] = 0;
        pileupweights[91] = 0;
        pileupweights[92] = 0;
        pileupweights[93] = 0;
        pileupweights[94] = 0;
        pileupweights[95] = 0;
        pileupweights[96] = 0;
        pileupweights[97] = 0;
        pileupweights[98] = 0;
        pileupweights[99] = 0;}

   if(RunYear == 2016 and preVFP == false){
        pileupweights[0] = 0;
        pileupweights[1] = 0;
        pileupweights[2] = 0.294853;
        pileupweights[3] = 0.147426;
        pileupweights[4] = 0.234148;
        pileupweights[5] = 0.239158;
        pileupweights[6] = 0.23329;
        pileupweights[7] = 0.279377;
        pileupweights[8] = 0.373123;
        pileupweights[9] = 0.447634;
        pileupweights[10] = 0.507057;
        pileupweights[11] = 0.506309;
        pileupweights[12] = 0.588541;
        pileupweights[13] = 0.649759;
        pileupweights[14] = 0.698702;
        pileupweights[15] = 0.726097;
        pileupweights[16] = 0.797024;
        pileupweights[17] = 0.883971;
        pileupweights[18] = 0.913238;
        pileupweights[19] = 0.917722;
        pileupweights[20] = 1.02381;
        pileupweights[21] = 1.1397;
        pileupweights[22] = 1.15752;
        pileupweights[23] = 1.25488;
        pileupweights[24] = 1.32716;
        pileupweights[25] = 1.42432;
        pileupweights[26] = 1.52131;
        pileupweights[27] = 1.47958;
        pileupweights[28] = 1.88607;
        pileupweights[29] = 1.81664;
        pileupweights[30] = 1.78508;
        pileupweights[31] = 1.88761;
        pileupweights[32] = 1.98657;
        pileupweights[33] = 2.3814;
        pileupweights[34] = 2.35535;
        pileupweights[35] = 2.11436;
        pileupweights[36] = 2.51853;
        pileupweights[37] = 3.60922;
        pileupweights[38] = 2.47676;
        pileupweights[39] = 2.86318;
        pileupweights[40] = 5.06778;
        pileupweights[41] = 5.23364;
        pileupweights[42] = 6.83076;
        pileupweights[43] = 4.27536;
        pileupweights[44] = 14.4478;
        pileupweights[45] = 11.4993;
        pileupweights[46] = 8.25588;
        pileupweights[47] = 2.35882;
        pileupweights[48] = 3.98051;
        pileupweights[49] = 0;
        pileupweights[50] = 2.8011;
        pileupweights[51] = 0;
        pileupweights[52] = 0;
        pileupweights[53] = 0;
        pileupweights[54] = 0;
        pileupweights[55] = 0;
        pileupweights[56] = 0;
        pileupweights[57] = 0;
        pileupweights[58] = 0;
        pileupweights[59] = 0;
        pileupweights[60] = 0;
        pileupweights[61] = 0;
        pileupweights[62] = 0;
        pileupweights[63] = 0;
        pileupweights[64] = 0;
        pileupweights[65] = 0;
        pileupweights[66] = 0;
        pileupweights[67] = 0;
        pileupweights[68] = 0;
        pileupweights[69] = 0;
        pileupweights[70] = 0;
        pileupweights[71] = 0;
        pileupweights[72] = 0;
        pileupweights[73] = 0;
        pileupweights[74] = 0;
        pileupweights[75] = 0;
        pileupweights[76] = 0;
        pileupweights[77] = 0;
        pileupweights[78] = 0;
        pileupweights[79] = 0;
        pileupweights[80] = 0;
        pileupweights[81] = 0;
        pileupweights[82] = 0;
        pileupweights[83] = 0;
        pileupweights[84] = 0;
        pileupweights[85] = 0;
        pileupweights[86] = 0;
        pileupweights[87] = 0;
        pileupweights[88] = 0;
        pileupweights[89] = 0;
        pileupweights[90] = 0;
        pileupweights[91] = 0;
        pileupweights[92] = 0;
        pileupweights[93] = 0;
        pileupweights[94] = 0;
        pileupweights[95] = 0;
        pileupweights[96] = 0;
        pileupweights[97] = 0;
        pileupweights[98] = 0;
        pileupweights[99] = 0;}


   if(RunYear == 2017){
	           pileupweights[0] = 0;
        pileupweights[1] = 0;
        pileupweights[2] = 0.589099;
        pileupweights[3] = 0.883649;
        pileupweights[4] = 0.621467;
        pileupweights[5] = 0.53937;
        pileupweights[6] = 0.366083;
        pileupweights[7] = 0.380591;
        pileupweights[8] = 0.428996;
        pileupweights[9] = 0.410105;
        pileupweights[10] = 0.407469;
        pileupweights[11] = 0.447313;
        pileupweights[12] = 0.46196;
        pileupweights[13] = 0.463202;
        pileupweights[14] = 0.488909;
        pileupweights[15] = 0.497674;
        pileupweights[16] = 0.512438;
        pileupweights[17] = 0.540785;
        pileupweights[18] = 0.576647;
        pileupweights[19] = 0.599647;
        pileupweights[20] = 0.610656;
        pileupweights[21] = 0.632355;
        pileupweights[22] = 0.687481;
        pileupweights[23] = 0.729983;
        pileupweights[24] = 0.743171;
        pileupweights[25] = 0.760318;
        pileupweights[26] = 0.793052;
        pileupweights[27] = 0.832838;
        pileupweights[28] = 0.882234;
        pileupweights[29] = 0.920284;
        pileupweights[30] = 0.957036;
        pileupweights[31] = 0.99634;
        pileupweights[32] = 1.10119;
        pileupweights[33] = 1.12411;
        pileupweights[34] = 1.17196;
        pileupweights[35] = 1.25766;
        pileupweights[36] = 1.25956;
        pileupweights[37] = 1.43347;
        pileupweights[38] = 1.45999;
        pileupweights[39] = 1.67431;
        pileupweights[40] = 1.68499;
        pileupweights[41] = 1.73574;
        pileupweights[42] = 1.8537;
        pileupweights[43] = 1.79025;
        pileupweights[44] = 2.09781;
        pileupweights[45] = 2.35716;
        pileupweights[46] = 2.35598;
        pileupweights[47] = 2.34743;
        pileupweights[48] = 2.49629;
        pileupweights[49] = 3.01763;
        pileupweights[50] = 2.79909;
        pileupweights[51] = 3.61932;
        pileupweights[52] = 3.87594;
        pileupweights[53] = 4.38064;
        pileupweights[54] = 4.8741;
        pileupweights[55] = 4.08162;
        pileupweights[56] = 5.10182;
        pileupweights[57] = 6.05931;
        pileupweights[58] = 7.28221;
        pileupweights[59] = 10.446;
        pileupweights[60] = 6.66524;
        pileupweights[61] = 8.75233;
        pileupweights[62] = 6.81672;
        pileupweights[63] = 6.35912;
        pileupweights[64] = 9.5939;
        pileupweights[65] = 16.0067;
        pileupweights[66] = 11.9166;
        pileupweights[67] = 30.9277;
        pileupweights[68] = 13.381;
        pileupweights[69] = 15.1483;
        pileupweights[70] = 0;
        pileupweights[71] = 0;
        pileupweights[72] = 30.0441;
        pileupweights[73] = 0;
        pileupweights[74] = 0;
        pileupweights[75] = 10.6038;
        pileupweights[76] = 0;
        pileupweights[77] = 0;
        pileupweights[78] = 0;
        pileupweights[79] = 11.8661;
        pileupweights[80] = 7.57413;
        pileupweights[81] = 0;
        pileupweights[82] = 0;
        pileupweights[83] = 0;
        pileupweights[84] = 0;
        pileupweights[85] = 0;
        pileupweights[86] = 0;
        pileupweights[87] = 0;
        pileupweights[88] = 0;
        pileupweights[89] = 2.01977;
        pileupweights[90] = 0;
        pileupweights[91] = 0;
        pileupweights[92] = 0;
        pileupweights[93] = 0;
        pileupweights[94] = 0;
        pileupweights[95] = 0;
        pileupweights[96] = 0;
        pileupweights[97] = 0;
        pileupweights[98] = 0;
        pileupweights[99] = 0;}

if(RunYear == 2018){
        pileupweights[0] = 0;
        pileupweights[1] = 0;
        pileupweights[2] = 0.676206;
        pileupweights[3] = 0.855201;
        pileupweights[4] = 0.784171;
        pileupweights[5] = 0.556875;
        pileupweights[6] = 0.557444;
        pileupweights[7] = 0.617633;
        pileupweights[8] = 0.649433;
        pileupweights[9] = 0.652437;
        pileupweights[10] = 0.634621;
        pileupweights[11] = 0.677451;
        pileupweights[12] = 0.77526;
        pileupweights[13] = 0.707019;
        pileupweights[14] = 0.747273;
        pileupweights[15] = 0.780067;
        pileupweights[16] = 0.77217;
        pileupweights[17] = 0.775911;
        pileupweights[18] = 0.800261;
        pileupweights[19] = 0.819243;
        pileupweights[20] = 0.832591;
        pileupweights[21] = 0.820051;
        pileupweights[22] = 0.86297;
        pileupweights[23] = 0.895992;
        pileupweights[24] = 0.898081;
        pileupweights[25] = 0.878921;
        pileupweights[26] = 0.891497;
        pileupweights[27] = 0.90168;
        pileupweights[28] = 0.925917;
        pileupweights[29] = 0.91614;
        pileupweights[30] = 0.929198;
        pileupweights[31] = 0.983414;
        pileupweights[32] = 0.995418;
        pileupweights[33] = 1.02254;
        pileupweights[34] = 1.02474;
        pileupweights[35] = 1.13139;
        pileupweights[36] = 1.20221;
        pileupweights[37] = 1.15919;
        pileupweights[38] = 1.2064;
        pileupweights[39] = 1.32074;
        pileupweights[40] = 1.4429;
        pileupweights[41] = 1.54274;
        pileupweights[42] = 1.5031;
        pileupweights[43] = 1.80117;
        pileupweights[44] = 1.9134;
        pileupweights[45] = 2.4657;
        pileupweights[46] = 2.43828;
        pileupweights[47] = 2.47887;
        pileupweights[48] = 2.7438;
        pileupweights[49] = 3.30075;
        pileupweights[50] = 3.30364;
        pileupweights[51] = 3.37535;
        pileupweights[52] = 4.88012;
        pileupweights[53] = 4.13148;
        pileupweights[54] = 5.12375;
        pileupweights[55] = 6.47096;
        pileupweights[56] = 6.91454;
        pileupweights[57] = 8.83045;
        pileupweights[58] = 6.45975;
        pileupweights[59] = 11.8989;
        pileupweights[60] = 6.03161;
        pileupweights[61] = 12.1478;
        pileupweights[62] = 9.23618;
        pileupweights[63] = 5.92177;
        pileupweights[64] = 12.6092;
        pileupweights[65] = 16.0499;
        pileupweights[66] = 7.66698;
        pileupweights[67] = 14.26;
        pileupweights[68] = 20.1668;
        pileupweights[69] = 20.1668;
        pileupweights[70] = 0;
        pileupweights[71] = 7.39849;
        pileupweights[72] = 0;
        pileupweights[73] = 6.20518;
        pileupweights[74] = 10.5011;
        pileupweights[75] = 8.94978;
        pileupweights[76] = 0;
        pileupweights[77] = 0;
        pileupweights[78] = 0;
        pileupweights[79] = 0;
        pileupweights[80] = 0;
        pileupweights[81] = 0;
        pileupweights[82] = 0;
        pileupweights[83] = 0;
        pileupweights[84] = 0;
        pileupweights[85] = 0;
        pileupweights[86] = 2.50594;
        pileupweights[87] = 0;
        pileupweights[88] = 0;
        pileupweights[89] = 0;
        pileupweights[90] = 1.43197;
        pileupweights[91] = 0;
        pileupweights[92] = 0;
        pileupweights[93] = 0;
        pileupweights[94] = 0;
        pileupweights[95] = 0;
        pileupweights[96] = 0;
        pileupweights[97] = 0;
        pileupweights[98] = 0;
        pileupweights[99] = 0;}



	if(nvertex > 99)return 0;
	else return pileupweights[nvertex]; 

}
void plot_PUReweight(){
	setTDRStyle();
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);
        gROOT->SetBatch(kTRUE);	
	TH1D *p_PU_data = new TH1D("p_PU_data",";N_{vtx};",100,0,100);
	TH1D *p_PU_MC = new TH1D("p_PU_MC","",100,0,100);
	TH1D *p_PU_MC_raw = new TH1D("p_PU_MC_raw","",100,0,100);

  	std::string whichVFP;
  	if(RunYear==2016 and preVFP == 1) whichVFP = "preVFP";
  	if(RunYear==2016 and preVFP == 0) whichVFP = "postVFP";
  	if(RunYear==2017 or  RunYear == 2018) whichVFP = "";

	TChain *sigtree = new TChain("signalTree");
	sigtree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s_Muon20.root",RunYear,whichVFP.c_str()));
	sigtree->Draw("nVertex >> p_PU_data");
 
	TChain *mctree;
  	mctree = new TChain("mgTree","mgTree");
	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_DYJetsToLL_%d%s.root",RunYear,whichVFP.c_str()));

	int   nVertex(0);
	float PUweight(1);
	float puWeight;
  	mctree->SetBranchAddress("nVertex",   &nVertex);
	mctree->SetBranchAddress("PUweight",  &PUweight);
	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);
		puWeight = getPUESF(nVertex);
		p_PU_MC->Fill(nVertex,puWeight);
		p_PU_MC_raw->Fill(nVertex);
	}

	
	p_PU_data->Sumw2();
	p_PU_data->Scale(1.0/p_PU_data->Integral(1,100));
	p_PU_MC->Scale(1.0/p_PU_MC->Integral(1,100));
	p_PU_MC_raw->Scale(1.0/p_PU_MC_raw->Integral(1,100));
	
	gStyle->SetOptStat(0);
	TCanvas *c_PU = new TCanvas("PU", "PU",600,600);
	c_PU->cd();
	TPad *PU_pad1 = new TPad("PU_pad1", "PU_pad1", 0, 0.3, 1, 0.98);
	PU_pad1->SetBottomMargin(0);
	PU_pad1->Draw();  
	PU_pad1->cd(); 
	p_PU_data->SetMaximum(1.4*p_PU_data->GetBinContent(p_PU_data->GetMaximumBin()));
	p_PU_data->SetLineColor(1);
	p_PU_data->SetMarkerStyle(20);
	p_PU_data->GetXaxis()->SetTitle("N_{vtx}");
	p_PU_data->Draw("P ");
	p_PU_MC_raw->Draw("hist same ");
	p_PU_MC->SetLineColor(kMagenta);
	p_PU_MC_raw->SetLineColor(kBlue);
	p_PU_MC->Draw("hist same");
	TLegend *leg_PU =  new TLegend(0.45,0.7,0.8,0.9);
	leg_PU->SetFillStyle(0);
	leg_PU->AddEntry(p_PU_data,"Data");
	leg_PU->AddEntry(p_PU_MC_raw,"Simulation (no correction)");
	leg_PU->AddEntry(p_PU_MC,"Simulation (corrected)");
	leg_PU->Draw("same");
	//if(RunYear==2016)lumitex->Draw();
	//if(RunYear==2017)lumitex17->Draw();
	//if(RunYear==2018)lumitex18->Draw();

	c_PU->cd();
	TPad *PU_pad2 = new TPad("PU_pad2", "PU_pad2", 0, 0.0, 1, 0.3);
	PU_pad2->SetTopMargin(0);
	PU_pad2->SetBottomMargin(0.4);
	PU_pad2->Draw();
	PU_pad2->cd();
  	TLine *flatratio_PU = new TLine(0,1,100,1);
	TH1F *ratio_PU=(TH1F*)p_PU_data->Clone("transfer factor");
	ratio_PU->SetMarkerStyle(20);
	ratio_PU->SetLineColor(kBlack);
	ratio_PU->GetXaxis()->SetRangeUser(0,100);
	ratio_PU->GetYaxis()->SetRangeUser(0,4);
	ratio_PU->SetMinimum(0);
	ratio_PU->SetMaximum(2.7);
	if(getPUweight) ratio_PU->Divide(p_PU_MC_raw); //to get the pileupweights
	if(!getPUweight) ratio_PU->Divide(p_PU_MC);
	
	cout<< RunYear <<endl;
	for(int bin=0;bin<100;bin++){
		if(getPUweight) cout<<"        pileupweights["<<bin<<"] = "<<ratio_PU->GetBinContent(bin)<<";\n";
		else 		cout<<"Pileup corrected now"<<"\n"<<"        pileupweights["<<bin<<"] = "<<ratio_PU->GetBinContent(bin)<<";\n";
	}
	cout<<endl;
	ratio_PU->SetTitle("");
	ratio_PU->GetYaxis()->SetTitle("Data/MC");
	ratio_PU->Draw();
	flatratio_PU->Draw("same");
        if(RunYear==2016 and preVFP == 1)       CMS_lumi(c_PU, 1, 11);
        else if(RunYear==2016 and preVFP == 0)  CMS_lumi(c_PU, 2, 11);
        else if(RunYear==2017)                  CMS_lumi(c_PU, 3, 11);
        else if(RunYear==2018)                  CMS_lumi(c_PU, 4, 11);
	if(!getPUweight)	c_PU->SaveAs(Form("/eos/uscms/store/user/tmishra/Plots_myAN/PLOT_PUreweight_%d%s.pdf",RunYear,whichVFP.c_str()));
}
