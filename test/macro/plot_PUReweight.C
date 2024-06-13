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

float getPUESF(int nvertex){
	
	float pileupweights[100];
    if(RunYear == 2016 and preVFP == true){
	pileupweights[0] = 0;
        pileupweights[1] = 0;
        pileupweights[2] = 0.91083;
        pileupweights[3] = 0.850108;
        pileupweights[4] = 0.663444;
        pileupweights[5] = 0.788321;
        pileupweights[6] = 0.62194;
        pileupweights[7] = 0.834791;
        pileupweights[8] = 0.801725;
        pileupweights[9] = 0.842834;
        pileupweights[10] = 0.945161;
        pileupweights[11] = 1.00213;
        pileupweights[12] = 0.968471;
        pileupweights[13] = 1.06945;
        pileupweights[14] = 1.02802;
        pileupweights[15] = 1.04851;
        pileupweights[16] = 1.06629;
        pileupweights[17] = 1.00153;
        pileupweights[18] = 1.02247;
        pileupweights[19] = 1.06958;
        pileupweights[20] = 1.06102;
        pileupweights[21] = 0.986573;
        pileupweights[22] = 0.982752;
        pileupweights[23] = 0.979845;
        pileupweights[24] = 0.922098;
        pileupweights[25] = 0.94373;
        pileupweights[26] = 0.935909;
        pileupweights[27] = 0.909424;
        pileupweights[28] = 0.91083;
        pileupweights[29] = 0.842551;
        pileupweights[30] = 0.919567;
        pileupweights[31] = 0.976259;
        pileupweights[32] = 0.976832;
        pileupweights[33] = 1.07258;
        pileupweights[34] = 1.0708;
        pileupweights[35] = 1.23001;
        pileupweights[36] = 1.48125;
        pileupweights[37] = 1.21444;
        pileupweights[38] = 0.937462;
        pileupweights[39] = 2.24204;
        pileupweights[40] = 1.89951;
        pileupweights[41] = 1.43371;
        pileupweights[42] = 1.92286;
        pileupweights[43] = 4.52042;
        pileupweights[44] = 1.26504;
        pileupweights[45] = 3.74452;
        pileupweights[46] = 6.0722;
        pileupweights[47] = 1.68672;
        pileupweights[48] = 0;
        pileupweights[49] = 0;
        pileupweights[50] = 0;
        pileupweights[51] = 3.64332;
        pileupweights[52] = 1.61925;
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
        pileupweights[2] = 0.355255;
        pileupweights[3] = 0.177628;
        pileupweights[4] = 0.189469;
        pileupweights[5] = 0.2505;
        pileupweights[6] = 0.248223;
        pileupweights[7] = 0.288226;
        pileupweights[8] = 0.367382;
        pileupweights[9] = 0.439095;
        pileupweights[10] = 0.522513;
        pileupweights[11] = 0.495065;
        pileupweights[12] = 0.601674;
        pileupweights[13] = 0.663256;
        pileupweights[14] = 0.712658;
        pileupweights[15] = 0.7244;
        pileupweights[16] = 0.798803;
        pileupweights[17] = 0.889869;
        pileupweights[18] = 0.905278;
        pileupweights[19] = 0.923172;
        pileupweights[20] = 1.05572;
        pileupweights[21] = 1.12308;
        pileupweights[22] = 1.16928;
        pileupweights[23] = 1.25401;
        pileupweights[24] = 1.33411;
        pileupweights[25] = 1.39892;
        pileupweights[26] = 1.53307;
        pileupweights[27] = 1.46587;
        pileupweights[28] = 1.84712;
        pileupweights[29] = 1.7436;
        pileupweights[30] = 1.84671;
        pileupweights[31] = 2.03409;
        pileupweights[32] = 1.92409;
        pileupweights[33] = 2.39517;
        pileupweights[34] = 2.38564;
        pileupweights[35] = 2.26544;
        pileupweights[36] = 2.48256;
        pileupweights[37] = 3.57476;
        pileupweights[38] = 2.41164;
        pileupweights[39] = 2.41771;
        pileupweights[40] = 5.0497;
        pileupweights[41] = 4.5295;
        pileupweights[42] = 9.32545;
        pileupweights[43] = 3.64137;
        pileupweights[44] = 13.855;
        pileupweights[45] = 8.52613;
        pileupweights[46] = 7.63799;
        pileupweights[47] = 1.71707;
        pileupweights[48] = 3.55255;
        pileupweights[49] = 0;
        pileupweights[50] = 2.84204;
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
        pileupweights[2] = 0.752414;
        pileupweights[3] = 1.20386;
        pileupweights[4] = 0.769134;
        pileupweights[5] = 0.742381;
        pileupweights[6] = 0.423233;
        pileupweights[7] = 0.416335;
        pileupweights[8] = 0.467837;
        pileupweights[9] = 0.396949;
        pileupweights[10] = 0.438457;
        pileupweights[11] = 0.472729;
        pileupweights[12] = 0.480629;
        pileupweights[13] = 0.480112;
        pileupweights[14] = 0.505426;
        pileupweights[15] = 0.506403;
        pileupweights[16] = 0.511837;
        pileupweights[17] = 0.547386;
        pileupweights[18] = 0.58658;
        pileupweights[19] = 0.608867;
        pileupweights[20] = 0.622004;
        pileupweights[21] = 0.641974;
        pileupweights[22] = 0.688321;
        pileupweights[23] = 0.739286;
        pileupweights[24] = 0.737365;
        pileupweights[25] = 0.742806;
        pileupweights[26] = 0.797297;
        pileupweights[27] = 0.830932;
        pileupweights[28] = 0.879505;
        pileupweights[29] = 0.925101;
        pileupweights[30] = 0.940651;
        pileupweights[31] = 0.982011;
        pileupweights[32] = 1.1115;
        pileupweights[33] = 1.14709;
        pileupweights[34] = 1.16921;
        pileupweights[35] = 1.27611;
        pileupweights[36] = 1.25675;
        pileupweights[37] = 1.43619;
        pileupweights[38] = 1.44172;
        pileupweights[39] = 1.63863;
        pileupweights[40] = 1.68784;
        pileupweights[41] = 1.76483;
        pileupweights[42] = 1.85339;
        pileupweights[43] = 1.78295;
        pileupweights[44] = 2.03227;
        pileupweights[45] = 2.3974;
        pileupweights[46] = 2.38789;
        pileupweights[47] = 2.33702;
        pileupweights[48] = 2.4359;
        pileupweights[49] = 3.03524;
        pileupweights[50] = 2.68272;
        pileupweights[51] = 3.57738;
        pileupweights[52] = 3.66485;
        pileupweights[53] = 4.84843;
        pileupweights[54] = 4.64413;
        pileupweights[55] = 3.85821;
        pileupweights[56] = 5.26689;
        pileupweights[57] = 5.85593;
        pileupweights[58] = 7.74407;
        pileupweights[59] = 8.74681;
        pileupweights[60] = 5.99423;
        pileupweights[61] = 8.57751;
        pileupweights[62] = 5.91309;
        pileupweights[63] = 6.30021;
        pileupweights[64] = 10.1576;
        pileupweights[65] = 13.0619;
        pileupweights[66] = 13.1672;
        pileupweights[67] = 26.184;
        pileupweights[68] = 11.6624;
        pileupweights[69] = 12.5402;
        pileupweights[70] = 0;
        pileupweights[71] = 0;
        pileupweights[72] = 26.184;
        pileupweights[73] = 0;
        pileupweights[74] = 0;
        pileupweights[75] = 9.78138;
        pileupweights[76] = 0;
        pileupweights[77] = 0;
        pileupweights[78] = 0;
        pileupweights[79] = 0;
        pileupweights[80] = 6.9222;
        pileupweights[81] = 0;
        pileupweights[82] = 0;
        pileupweights[83] = 0;
        pileupweights[84] = 0;
        pileupweights[85] = 0;
        pileupweights[86] = 0;
        pileupweights[87] = 0;
        pileupweights[88] = 0;
        pileupweights[89] = 1.80579;
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
        pileupweights[2] = 0.765237;
        pileupweights[3] = 0.832196;
        pileupweights[4] = 0.772595;
        pileupweights[5] = 0.603468;
        pileupweights[6] = 0.585689;
        pileupweights[7] = 0.689118;
        pileupweights[8] = 0.64109;
        pileupweights[9] = 0.648843;
        pileupweights[10] = 0.649582;
        pileupweights[11] = 0.69035;
        pileupweights[12] = 0.797717;
        pileupweights[13] = 0.726607;
        pileupweights[14] = 0.766009;
        pileupweights[15] = 0.791764;
        pileupweights[16] = 0.781238;
        pileupweights[17] = 0.786906;
        pileupweights[18] = 0.817605;
        pileupweights[19] = 0.830885;
        pileupweights[20] = 0.838369;
        pileupweights[21] = 0.826052;
        pileupweights[22] = 0.878224;
        pileupweights[23] = 0.911633;
        pileupweights[24] = 0.901887;
        pileupweights[25] = 0.887659;
        pileupweights[26] = 0.885957;
        pileupweights[27] = 0.895242;
        pileupweights[28] = 0.94038;
        pileupweights[29] = 0.908796;
        pileupweights[30] = 0.927954;
        pileupweights[31] = 0.97809;
        pileupweights[32] = 0.995497;
        pileupweights[33] = 1.01097;
        pileupweights[34] = 1.00282;
        pileupweights[35] = 1.13299;
        pileupweights[36] = 1.21651;
        pileupweights[37] = 1.15839;
        pileupweights[38] = 1.16259;
        pileupweights[39] = 1.29246;
        pileupweights[40] = 1.41637;
        pileupweights[41] = 1.53915;
        pileupweights[42] = 1.43431;
        pileupweights[43] = 1.72249;
        pileupweights[44] = 1.87819;
        pileupweights[45] = 2.36869;
        pileupweights[46] = 2.51663;
        pileupweights[47] = 2.47825;
        pileupweights[48] = 2.84686;
        pileupweights[49] = 3.36;
        pileupweights[50] = 3.23184;
        pileupweights[51] = 3.35623;
        pileupweights[52] = 4.6348;
        pileupweights[53] = 4.2618;
        pileupweights[54] = 4.64696;
        pileupweights[55] = 6.93496;
        pileupweights[56] = 7.22193;
        pileupweights[57] = 7.65973;
        pileupweights[58] = 6.18076;
        pileupweights[59] = 10.5357;
        pileupweights[60] = 5.76798;
        pileupweights[61] = 10.9907;
        pileupweights[62] = 10.1514;
        pileupweights[63] = 5.90326;
        pileupweights[64] = 10.4742;
        pileupweights[65] = 13.3438;
        pileupweights[66] = 6.45669;
        pileupweights[67] = 12.4829;
        pileupweights[68] = 18.9396;
        pileupweights[69] = 17.2178;
        pileupweights[70] = 0;
        pileupweights[71] = 5.73928;
        pileupweights[72] = 0;
        pileupweights[73] = 5.16535;
        pileupweights[74] = 9.32633;
        pileupweights[75] = 7.74803;
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
        pileupweights[86] = 1.86527;
        pileupweights[87] = 0;
        pileupweights[88] = 0;
        pileupweights[89] = 0;
        pileupweights[90] = 1.43482;
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
	//ratio_PU->Divide(p_PU_MC_raw); //to get the pileupweights
	ratio_PU->Divide(p_PU_MC);
	
	cout<< RunYear <<endl;
	for(int bin=0;bin<100;bin++){
		cout<<"        pileupweights["<<bin<<"] = "<<ratio_PU->GetBinContent(bin)<<";\n";
		//cout<<ratio_PU->GetBinContent(bin)<<"; ";
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
	c_PU->SaveAs(Form("/eos/uscms/store/user/tmishra/Plots_myAN/PLOT_PUreweight_%d%s.pdf",RunYear,whichVFP.c_str()));
}
