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
int RunYear = 2018;
bool preVFP = false;
bool getPUweight = false; // for making PU plots
//bool getPUweight = true;

float getPUESF(int nvertex){
	
	float pileupweights[100];
    if(RunYear == 2016 and preVFP == true){
	        pileupweights[1] = 0;
        pileupweights[2] = 1.06566;
        pileupweights[3] = 0.852528;
        pileupweights[4] = 0.634321;
        pileupweights[5] = 0.698599;
        pileupweights[6] = 0.600106;
        pileupweights[7] = 0.84494;
        pileupweights[8] = 0.776627;
        pileupweights[9] = 0.841392;
        pileupweights[10] = 0.916784;
        pileupweights[11] = 0.966177;
        pileupweights[12] = 0.977609;
        pileupweights[13] = 1.05143;
        pileupweights[14] = 1.02428;
        pileupweights[15] = 1.03655;
        pileupweights[16] = 1.0732;
        pileupweights[17] = 1.0132;
        pileupweights[18] = 1.03824;
        pileupweights[19] = 1.05213;
        pileupweights[20] = 1.07704;
        pileupweights[21] = 0.987459;
        pileupweights[22] = 0.978869;
        pileupweights[23] = 0.973576;
        pileupweights[24] = 0.940438;
        pileupweights[25] = 0.940926;
        pileupweights[26] = 0.936716;
        pileupweights[27] = 0.936811;
        pileupweights[28] = 0.917883;
        pileupweights[29] = 0.853291;
        pileupweights[30] = 0.951073;
        pileupweights[31] = 0.980436;
        pileupweights[32] = 1.05052;
        pileupweights[33] = 1.06566;
        pileupweights[34] = 1.02818;
        pileupweights[35] = 1.24327;
        pileupweights[36] = 1.51324;
        pileupweights[37] = 1.14124;
        pileupweights[38] = 0.924379;
        pileupweights[39] = 2.16462;
        pileupweights[40] = 2.07667;
        pileupweights[41] = 1.7021;
        pileupweights[42] = 1.95371;
        pileupweights[43] = 4.97308;
        pileupweights[44] = 1.26547;
        pileupweights[45] = 2.66415;
        pileupweights[46] = 6.92679;
        pileupweights[47] = 2.13132;
        pileupweights[48] = 0;
        pileupweights[49] = 0;
        pileupweights[50] = 0;
        pileupweights[51] = 4.44025;
        pileupweights[52] = 1.42088;
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
        pileupweights[2] = 0.308456;
        pileupweights[3] = 0.154228;
        pileupweights[4] = 0.226806;
        pileupweights[5] = 0.243337;
        pileupweights[6] = 0.22541;
        pileupweights[7] = 0.279324;
        pileupweights[8] = 0.370708;
        pileupweights[9] = 0.451854;
        pileupweights[10] = 0.51652;
        pileupweights[11] = 0.509773;
        pileupweights[12] = 0.593442;
        pileupweights[13] = 0.655468;
        pileupweights[14] = 0.703156;
        pileupweights[15] = 0.728414;
        pileupweights[16] = 0.800497;
        pileupweights[17] = 0.890329;
        pileupweights[18] = 0.913117;
        pileupweights[19] = 0.919977;
        pileupweights[20] = 1.02277;
        pileupweights[21] = 1.13619;
        pileupweights[22] = 1.15348;
        pileupweights[23] = 1.24963;
        pileupweights[24] = 1.32115;
        pileupweights[25] = 1.42108;
        pileupweights[26] = 1.5223;
        pileupweights[27] = 1.48465;
        pileupweights[28] = 1.87111;
        pileupweights[29] = 1.81841;
        pileupweights[30] = 1.76388;
        pileupweights[31] = 1.87812;
        pileupweights[32] = 1.97198;
        pileupweights[33] = 2.36066;
        pileupweights[34] = 2.34426;
        pileupweights[35] = 2.07524;
        pileupweights[36] = 2.54155;
        pileupweights[37] = 3.50726;
        pileupweights[38] = 2.4368;
        pileupweights[39] = 2.8248;
        pileupweights[40] = 4.89673;
        pileupweights[41] = 5.24375;
        pileupweights[42] = 6.78603;
        pileupweights[43] = 4.27982;
        pileupweights[44] = 14.0347;
        pileupweights[45] = 11.1044;
        pileupweights[46] = 8.01985;
        pileupweights[47] = 2.41624;
        pileupweights[48] = 3.8557;
        pileupweights[49] = 0;
        pileupweights[50] = 2.7761;
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
        pileupweights[2] = 0.6179;
        pileupweights[3] = 0.882714;
        pileupweights[4] = 0.63148;
        pileupweights[5] = 0.565739;
        pileupweights[6] = 0.37736;
        pileupweights[7] = 0.38734;
        pileupweights[8] = 0.432745;
        pileupweights[9] = 0.404615;
        pileupweights[10] = 0.415658;
        pileupweights[11] = 0.450588;
        pileupweights[12] = 0.466161;
        pileupweights[13] = 0.472166;
        pileupweights[14] = 0.494203;
        pileupweights[15] = 0.504551;
        pileupweights[16] = 0.518071;
        pileupweights[17] = 0.543537;
        pileupweights[18] = 0.581046;
        pileupweights[19] = 0.604227;
        pileupweights[20] = 0.616439;
        pileupweights[21] = 0.636305;
        pileupweights[22] = 0.69373;
        pileupweights[23] = 0.734144;
        pileupweights[24] = 0.749157;
        pileupweights[25] = 0.766277;
        pileupweights[26] = 0.795866;
        pileupweights[27] = 0.835479;
        pileupweights[28] = 0.885063;
        pileupweights[29] = 0.921053;
        pileupweights[30] = 0.959103;
        pileupweights[31] = 0.994048;
        pileupweights[32] = 1.09906;
        pileupweights[33] = 1.12441;
        pileupweights[34] = 1.17567;
        pileupweights[35] = 1.25765;
        pileupweights[36] = 1.25387;
        pileupweights[37] = 1.43409;
        pileupweights[38] = 1.45452;
        pileupweights[39] = 1.66797;
        pileupweights[40] = 1.67394;
        pileupweights[41] = 1.72591;
        pileupweights[42] = 1.83773;
        pileupweights[43] = 1.78216;
        pileupweights[44] = 2.08957;
        pileupweights[45] = 2.33969;
        pileupweights[46] = 2.32932;
        pileupweights[47] = 2.32311;
        pileupweights[48] = 2.46955;
        pileupweights[49] = 3.00237;
        pileupweights[50] = 2.78357;
        pileupweights[51] = 3.57233;
        pileupweights[52] = 3.7879;
        pileupweights[53] = 4.27928;
        pileupweights[54] = 4.76666;
        pileupweights[55] = 3.98919;
        pileupweights[56] = 5.02647;
        pileupweights[57] = 6.03776;
        pileupweights[58] = 7.39771;
        pileupweights[59] = 10.1953;
        pileupweights[60] = 6.5568;
        pileupweights[61] = 8.79183;
        pileupweights[62] = 6.59093;
        pileupweights[63] = 6.28934;
        pileupweights[64] = 9.34794;
        pileupweights[65] = 15.2533;
        pileupweights[66] = 11.9696;
        pileupweights[67] = 31.5129;
        pileupweights[68] = 13.7041;
        pileupweights[69] = 14.7413;
        pileupweights[70] = 0;
        pileupweights[71] = 0;
        pileupweights[72] = 29.3944;
        pileupweights[73] = 0;
        pileupweights[74] = 0;
        pileupweights[75] = 10.0629;
        pileupweights[76] = 0;
        pileupweights[77] = 0;
        pileupweights[78] = 0;
        pileupweights[79] = 11.387;
        pileupweights[80] = 7.67961;
        pileupweights[81] = 0;
        pileupweights[82] = 0;
        pileupweights[83] = 0;
        pileupweights[84] = 0;
        pileupweights[85] = 0;
        pileupweights[86] = 0;
        pileupweights[87] = 0;
        pileupweights[88] = 0;
        pileupweights[89] = 1.8537;
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
        pileupweights[2] = 0.67751;
        pileupweights[3] = 0.856851;
        pileupweights[4] = 0.828384;
        pileupweights[5] = 0.564081;
        pileupweights[6] = 0.565351;
        pileupweights[7] = 0.621452;
        pileupweights[8] = 0.654518;
        pileupweights[9] = 0.660033;
        pileupweights[10] = 0.63146;
        pileupweights[11] = 0.675482;
        pileupweights[12] = 0.775103;
        pileupweights[13] = 0.702737;
        pileupweights[14] = 0.745734;
        pileupweights[15] = 0.772471;
        pileupweights[16] = 0.765172;
        pileupweights[17] = 0.769519;
        pileupweights[18] = 0.798855;
        pileupweights[19] = 0.819572;
        pileupweights[20] = 0.834141;
        pileupweights[21] = 0.819299;
        pileupweights[22] = 0.861469;
        pileupweights[23] = 0.895489;
        pileupweights[24] = 0.9012;
        pileupweights[25] = 0.881325;
        pileupweights[26] = 0.896319;
        pileupweights[27] = 0.903338;
        pileupweights[28] = 0.928523;
        pileupweights[29] = 0.918563;
        pileupweights[30] = 0.930576;
        pileupweights[31] = 0.981967;
        pileupweights[32] = 0.995202;
        pileupweights[33] = 1.01947;
        pileupweights[34] = 1.0227;
        pileupweights[35] = 1.12535;
        pileupweights[36] = 1.20962;
        pileupweights[37] = 1.15884;
        pileupweights[38] = 1.20439;
        pileupweights[39] = 1.31233;
        pileupweights[40] = 1.45323;
        pileupweights[41] = 1.53299;
        pileupweights[42] = 1.50747;
        pileupweights[43] = 1.80464;
        pileupweights[44] = 1.92413;
        pileupweights[45] = 2.46721;
        pileupweights[46] = 2.44886;
        pileupweights[47] = 2.49157;
        pileupweights[48] = 2.73156;
        pileupweights[49] = 3.30053;
        pileupweights[50] = 3.3373;
        pileupweights[51] = 3.39273;
        pileupweights[52] = 4.82477;
        pileupweights[53] = 4.18994;
        pileupweights[54] = 4.98418;
        pileupweights[55] = 6.48345;
        pileupweights[56] = 6.9478;
        pileupweights[57] = 8.77391;
        pileupweights[58] = 6.52801;
        pileupweights[59] = 11.8023;
        pileupweights[60] = 5.96716;
        pileupweights[61] = 12.1952;
        pileupweights[62] = 9.56485;
        pileupweights[63] = 5.9033;
        pileupweights[64] = 12.3147;
        pileupweights[65] = 16.0211;
        pileupweights[66] = 7.65188;
        pileupweights[67] = 13.4506;
        pileupweights[68] = 20.4449;
        pileupweights[69] = 20.4449;
        pileupweights[70] = 0;
        pileupweights[71] = 7.17364;
        pileupweights[72] = 0;
        pileupweights[73] = 6.63561;
        pileupweights[74] = 10.0431;
        pileupweights[75] = 9.32573;
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
        pileupweights[86] = 2.27165;
        pileupweights[87] = 0;
        pileupweights[88] = 0;
        pileupweights[89] = 0;
        pileupweights[90] = 1.55429;
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
	sigtree->Add(Form("/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_%d%s.root",RunYear,whichVFP.c_str()));
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
	if(RunYear==2016 and preVFP == 1)       CMS_lumi( PU_pad1,1, 1, 11 );
        else if(RunYear==2016 and preVFP == 0)  CMS_lumi( PU_pad1,2, 1, 11 );
        else if(RunYear==2017)                  CMS_lumi( PU_pad1,3, 1, 11 );
        else if(RunYear==2018)                  CMS_lumi( PU_pad1,4, 1, 11 );

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
//        if(RunYear==2016 and preVFP == 1)       CMS_lumi(c_PU, 1, 1, 11);
//        else if(RunYear==2016 and preVFP == 0)  CMS_lumi(c_PU, 2, 1, 11);
//        else if(RunYear==2017)                  CMS_lumi(c_PU, 3, 1, 11);
//        else if(RunYear==2018)                  CMS_lumi(c_PU, 4, 1,11);
	if(!getPUweight)	c_PU->SaveAs(Form("/eos/uscms/store/user/tmishra/Plots_myAN/PLOT_PUreweight_%d%s.pdf",RunYear,whichVFP.c_str()));
}
