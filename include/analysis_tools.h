// https://github.com/scodella/LatinoAnalysis/blob/UL_worker/NanoGardener/python/framework/samples/samplesCrossSections2016.py
// https://github.com/TreeMaker/TreeMaker/blob/Run2_UL/WeightProducer/python/MCSampleValues.py
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
// https://twiki.cern.ch/twiki/bin/view/CMS/HowToGenXSecAnalyzer#Automated_scripts_to_compute_the



enum MCType{
  NOMC = 0,
  WGJetInclusive = 1,
  WGJet40 = 2,
  WGJet130 = 3,
	ZGInclusive = 4,
	DYLL50 = 5,
	TTG = 6,
	WWG = 7,
	WZG = 8,
	DYLL10 = 9, // not used
	ZG130 = 10, // not used
	TT = 11,
	WW = 12,
	WZ = 13,
  W  = 14,
	QCDEM30 = 15, // not used
	QCDEM40 = 16, // not used 
	QCDMU = 17, // not used 
	GJet = 18,  

	T5WgmG1800mLSP800 = 19,
	TChiWgmChi1000mLSP1 = 20,
	T5WgmG1800mLSP1600 = 21,
	T5WgmG1500mLSP1 = 22,
  generalMC = 23,
  ZZ = 24,
  nMCType = 25
};

// XSec for T5Wg_mG-1800_mLSP-800  = .001027
// XSec for TChiWg_mChi-1000_mLSP-1 = .001019
// multiplied a k-factor of 1.34 on WGToLNuG, WGJet samples

double MC_XS[25] = {1, 551.008 /*WGToLNuG*/, 26.5454/*WGJet40*/, 1.083658/*WGJet130*/, 51.53 /*ZGInclusive*/, 6424.0/*DY*/, 3.774/*TTG*/, 0.3369/*WWG*/, 0.07876/*WZG*/, 18610, 0.143, 734.6/*TTBar*/, 75.95/*WW*/, 27.55/*WZ*/, 61526.7/*W*/, 108000000/*QCDEM30*/, 113100.0 /*QCDEM40*/, 1/*MU*/, 874.2/*GJet*/, 0.001027 /*T5Wg_mG-1800_mLSP-800*/, 0.001019 /*TChiWg_mChi-1000_mLSP-1*/, 0.001057/*T5Wg_mG-1800_mLSP-1600*/,   0.005763 /*T5Wg_mG-1500_mLSP-1*/, 1, 12.23/*ZZ*/};


/* 117400.0 QCDEM40 for 2017*/

//0.143 // LO ZG130
//0.1404 // NLO ZG130


// Define luminosity values (in fb^-1)
float lumi_2016preVFP_DoubleEG = 19.499256;
float lumi_2016postVFP_DoubleEG = 16.810813;
float lumi_2017_DoubleEG = 41.4613;
float lumi_2018_DoubleEG = 59.816229;
float lumi_678_DoubleEG = 137.5876;

float lumi_2016preVFP_MuonEG = 19.495077;
float lumi_2016postVFP_MuonEG = 16.810813;
float lumi_2017_MuonEG = 27.100211; // considering only DEF eras, because of availabilty of trigger
float lumi_2018_MuonEG = 59.810698;
float lumi_678_MuonEG = 123.2168;


float getEvtWeight(int ichannel, int year, bool preVFP, double XSec, double nEvents_MC) {
    // Initialize luminosity to handle invalid cases
    float luminosity = 0.0;

    // Determine luminosity based on channel and year
    if (ichannel == 1) { // DoubleEG channel
        if (year == 2016 && preVFP) luminosity = lumi_2016preVFP_DoubleEG*1000;
        else if (year == 2016 && !preVFP) luminosity = lumi_2016postVFP_DoubleEG*1000;
        else if (year == 2017) luminosity = lumi_2017_DoubleEG*1000;
        else if (year == 2018) luminosity = lumi_2018_DoubleEG*1000;
    } else if (ichannel == 2) { // MuonEG channel
        if (year == 2016 && preVFP) luminosity = lumi_2016preVFP_MuonEG*1000;
        else if (year == 2016 && !preVFP) luminosity = lumi_2016postVFP_MuonEG*1000;
        else if (year == 2017) luminosity = lumi_2017_MuonEG*1000;
        else if (year == 2018) luminosity = lumi_2018_MuonEG*1000;
    }

    // Check for invalid inputs
    if (luminosity == 0.0) {
        std::cerr << "Error: Invalid ichannel (" << ichannel << ") or year (" << year << ")." << std::endl;
        return 0.0; // Return weight as 0 to indicate an error
    }

    // Calculate event weight
    double evtWeight = XSec * luminosity / nEvents_MC;
    return evtWeight;
}




float DeltaPhi(float phi1, float phi2){
  float deltaPhi = phi1 - phi2;
  if(fabs(deltaPhi) > TMath::Pi()){
    if(deltaPhi > 0)deltaPhi = -1.0*(TMath::TwoPi() - fabs(deltaPhi));
    else deltaPhi = TMath::TwoPi() - fabs(deltaPhi);
  }
  return deltaPhi;
}

float DeltaR(float eta1,float phi1,float eta2,float phi2)
{
	float deltaPhi = TMath::Abs(phi1-phi2);
	float deltaEta = eta1-eta2;
	if(deltaPhi > TMath::Pi())
	deltaPhi = TMath::TwoPi() - deltaPhi;
		return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}
// prompt photon if momID
// 1-6 for quarks, 11, 13, 15 leptons, 21, 22, 24 for gluons, photons, W+
// check for photon (PID 22) and electron(PID 11)
bool isHad(int PID, int momID){
   bool isFakePho;
   if(fabs(PID) == 22 || fabs(PID) == 11){
       switch(momID){
         case 1: isFakePho = false; break;
         case 2: isFakePho = false; break;
         case 3: isFakePho = false; break;
         case 4: isFakePho = false; break;
         case 5: isFakePho = false; break;
         case 6: isFakePho = false; break;
         case 999: isFakePho = false; break;
         case 11: isFakePho = false; break;
         case 13: isFakePho = false; break;
         case 15: isFakePho = false; break;
         case 22: isFakePho = false; break;
         case 24: isFakePho = false; break;
         case 21: isFakePho = false; break;
         default: isFakePho = true; break;
       }
  }
  else isFakePho = true;

  return isFakePho;
}

template<class recoParticle>
bool Veto(recoParticle reco, std::vector<recoPhoton>::iterator itpho, float dRcut){
  bool passVeto(true);
  typename recoParticle::iterator ire;
  for(ire = reco.begin(); ire != reco.end(); ire++)
    if(ire->getEt() > 2.0 && DeltaR(itpho->getEta(), itpho->getPhi(), ire->getEta(), ire->getPhi()) < dRcut)passVeto = false;
  return passVeto;
}


Double_t mybw(Double_t* x, Double_t* par)
{
  Double_t arg1 = 14.0/22.0; // 2 over pi
  Double_t arg2 = par[1]*par[1]*par[2]*par[2]; //Gamma=par[1]  M=par[2]
  Double_t arg3 = ((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
  Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[1]*par[1])/(par[2]*par[2]));
  Double_t arg5 = exp(par[3]+par[4]*x[0]);
  return par[0]*arg1*arg2/(arg3 + arg4)+ arg5;
}


Double_t mybwpol(Double_t* x, Double_t* par)
{
  Double_t arg1 = 14.0/22.0; // 2 over pi
  Double_t arg2 = par[1]*par[1]*par[2]*par[2]; //Gamma=par[1]  M=par[2]
  Double_t arg3 = ((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
  Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[1]*par[1])/(par[2]*par[2]));
  Double_t arg5 = par[3]*x[0]*x[0] + par[4]*x[0]+par[5]+par[6];
  return par[0]*arg1*arg2/(arg3 + arg4)+ arg5;
}

void DrawHisto(TObjArray *obj_list,TCanvas *canvas, int iPad, bool setLog){
  std::ostringstream name;
  int ihist(0);
  TObjArrayIter itr(obj_list);
  TH1F* hist(0);

  canvas->cd(iPad);
  if(setLog)gPad->SetLogy();
  TLegend *leg =  new TLegend(0.5,0.6,0.96,0.75);
  while( (hist = static_cast<TH1F*>(itr.Next())) ){
     ihist += 1;
     hist->SetLineColor(ihist+1);
     hist->SetLineWidth(2);
     hist->SetMarkerStyle(8);
     hist->SetMarkerColor(ihist+1);
     name.str("");
     name << hist->GetName();

     if(name.str().find("loose")!= std::string::npos){hist->SetLineColor(1);hist->SetMarkerColor(1);}
     else if(name.str().find("medium")!= std::string::npos){hist->SetLineColor(8);hist->SetMarkerColor(8);}
     else if(name.str().find("tight")!= std::string::npos){hist->SetLineColor(4);hist->SetMarkerColor(4);}
     leg->AddEntry(hist,name.str().c_str());
     if(ihist==1)hist->Draw();
     else hist->Draw("same");
  }
  leg->Draw("same");
}




float PrefireWeight(int mcType, rawData raw){
	if(mcType==0) return 1.0;
	else return raw.L1ECALPrefire;
}

int StartHEM = 319077;
double passHEMVeto(int mcType, rawData raw){
    if(mcType==0 && raw.run < StartHEM) return true;
    for(int iJet(0); iJet < raw.nJet; iJet++){
        if(recoJet(raw, iJet).getPt() > 30 && recoJet(raw, iJet).getEta() > -3.2 && recoJet(raw, iJet).getEta() < -1.2 && recoJet(raw, iJet).getPhi() > -1.77 && recoJet(raw, iJet).getPhi() < -0.67 && DeltaPhi(recoJet(raw, iJet).getPhi(),raw.pfMETPhi)<0.5) 
	return false;
    }

    for(int iEle(0); iEle < raw.nEle; iEle++){
	if(recoEle(raw, iEle).getCalibEt() > 30 && recoEle(raw, iEle).getEta() > -3.0 && recoEle(raw, iEle).getEta() < -1.4 && recoEle(raw, iEle).getPhi() > -1.57 && recoEle(raw, iEle).getPhi() < -0.87) 
	return false;
    }
        return true;   
}

double passPixelIssue17(rawData raw){
    for(int iEle(0); iEle < raw.nEle; iEle++){
	if(1.0 < recoEle(raw, iEle).getEta() && recoEle(raw, iEle).getEta() < 1.5 && recoEle(raw, iEle).getPhi() > 2.7) 
	return false;
    }
    for(int iPho(0); iPho < raw.nPho; iPho++){
	if(1.0 < recoPhoton(raw, iPho).getEta() && recoPhoton(raw, iPho).getEta() < 1.5 && recoPhoton(raw, iPho).getPhi() > 2.7) 
	return false;
    }
    for(int iMu(0); iMu < raw.nMu; iMu++){
	if(1.0 < recoMuon(raw, iMu).getEta() && recoMuon(raw, iMu).getEta() < 1.5 && recoMuon(raw, iMu).getPhi() > 2.7) 
	return false;
    }
	return true;
}

double passPixelIssue18(rawData raw){
    for(int iEle(0); iEle < raw.nEle; iEle++){
	if(0.3 < recoEle(raw, iEle).getEta() && recoEle(raw, iEle).getEta() < 1.2 && 0.4 < recoEle(raw, iEle).getPhi() && recoEle(raw, iEle).getPhi() < 0.8) 
	return false;
    }
    for(int iPho(0); iPho < raw.nPho; iPho++){
	if(0.3 < recoPhoton(raw, iPho).getEta() && recoPhoton(raw, iPho).getEta() < 1.2 && 0.4 < recoPhoton(raw, iPho).getPhi() && recoPhoton(raw, iPho).getPhi() < 0.8) 
	return false;
    }
    for(int iMu(0); iMu < raw.nMu; iMu++){
	if(0.3 < recoMuon(raw, iMu).getEta() && recoMuon(raw, iMu).getEta() < 1.2 && 0.4 < recoMuon(raw, iMu).getPhi() && recoMuon(raw, iMu).getPhi() < 0.8) 
	return false;
    }
	return true;
}



unsigned findIndex(float* array, float kinevar,unsigned len){
  unsigned Index(0);
  for(unsigned i(0); i< len; i++){
    if(i != len-1){
      if(kinevar >= *(array+i) && kinevar < *(array+i+1))Index = i;}
    else{
      if(kinevar >= *(array+i))Index = i;
    } 
  }
  return Index;
}


float getPUESF16preVFP(int nvertex){
	
	float pileupweights[100];
	        pileupweights[0] = 0;
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
        pileupweights[99] = 0;

	if(nvertex > 99)return 0;
	else return pileupweights[nvertex]; 
}

float getPUESF16(int nvertex){
	
	float pileupweights[100];
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
        pileupweights[99] = 0;

	if(nvertex > 99)return 0;
	else return pileupweights[nvertex]; 
}
float getPUESF17(int nvertex){
	
	float pileupweights[100];	
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
        pileupweights[99] = 0;

	if(nvertex > 99)return 0;
	else return pileupweights[nvertex]; 
}

float getPUESF18(int nvertex){
	
	float pileupweights[100];
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
        pileupweights[99] = 0;

	
	if(nvertex > 99)return 0;
	else return pileupweights[nvertex]; 
}

//Double_t bkgEtBins[]={35,40,50,60,70,80,90,100,110,120,130,140,150,160,170,185,200,215,230,250,275,290, 305,325,345,370,400,500,800};
Double_t bkgEtBins[]={35,50,100,150,200,250,300,500,800};
int nBkgEtBins= sizeof(bkgEtBins)/sizeof(bkgEtBins[0]) -1;

Double_t bkgPtBins[]={25,50,75,100,125,150,200,400,800};
int nBkgPtBins= sizeof(bkgPtBins)/sizeof(bkgPtBins[0])-1;
Double_t bkgMETBins[]={0,40,60,80,100,120,140,160,180,210,240,280,320,400,600,1000};
int nBkgMETBins= sizeof(bkgMETBins)/sizeof(bkgMETBins[0]) -1;
Double_t bkgMtBins[]={0,20,40,60,80,100,120,140,160,180,200,300,400,500,1000};
int nBkgMtBins= sizeof(bkgMtBins)/sizeof(bkgMtBins[0]) -1;
Double_t bkgHTBins[]={0,40,60,80,100,120,140,160,180,200,225,250,275,300,340,380,420,500,600,1000};
int nBkgHTBins = sizeof(bkgHTBins)/sizeof(bkgHTBins[0]) -1;

//Double_t sigEtBins[]={35,50,100,150,200,300,500,800};
//int nSigEtBins= sizeof(sigEtBins)/sizeof(sigEtBins[0]) -1;
//Double_t sigPtBins[]={25,50,100,150,200,300,500,800};
//int nSigPtBins= sizeof(sigPtBins)/sizeof(sigPtBins[0])-1;
//Double_t sigMETBins[]={120,200,300,400,550,1000};
//int nSigMETBins= sizeof(sigMETBins)/sizeof(sigMETBins[0]) -1;
//Double_t sigMtBins[]={100,200,300,400,600,1000};
//int nSigMtBins= sizeof(sigMtBins)/sizeof(sigMtBins[0]) -1;
//Double_t sigHTBins[]={0,400,800,1500,2000};
//int nSigHTBins= sizeof(sigHTBins)/sizeof(sigHTBins[0]) -1;
//
Double_t sigEtBins[]={35,50,100,150,200,300,500,800};
int nSigEtBins= sizeof(sigEtBins)/sizeof(sigEtBins[0]) -1;
Double_t sigPtBins[]={25,50,100,150,200,300,500,800};
int nSigPtBins= sizeof(sigPtBins)/sizeof(sigPtBins[0])-1;
Double_t sigMETBins[]={0,20,40,60,80,100,120,150,200,250,300,350,400,600};
int nSigMETBins= sizeof(sigMETBins)/sizeof(sigMETBins[0]) -1;
Double_t sigMtBins[]={100,200,300,400,600,1000};
int nSigMtBins= sizeof(sigMtBins)/sizeof(sigMtBins[0]) -1;
Double_t sigHTBins[]={0,100,200,300,400,800,1500,2000};
int nSigHTBins= sizeof(sigHTBins)/sizeof(sigHTBins[0]) -1;
