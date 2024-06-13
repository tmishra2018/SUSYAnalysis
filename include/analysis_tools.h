// https://github.com/scodella/LatinoAnalysis/blob/UL_worker/NanoGardener/python/framework/samples/samplesCrossSections2016.py
// https://github.com/TreeMaker/TreeMaker/blob/Run2_UL/WeightProducer/python/MCSampleValues.py
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
// https://twiki.cern.ch/twiki/bin/view/CMS/HowToGenXSecAnalyzer#Automated_scripts_to_compute_the


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
  nMCType = 24
};

// XSec for T5Wg_mG-1800_mLSP-800  = .001027
// XSec for TChiWg_mChi-1000_mLSP-1 = .001019


// multiplied a k-factor of 1.34 on WGToLNuG, WGJet samples

// menglei's crosssection
//double MC_XS[20] = {1, 489.0/*WGJetInclusive*/, 17.018/*WGJet40*/, 0.87971/*WGJet130*/, 55.48/*ZGInclusive*/, 5670/*DY*/, 3.697/*TTG*/, 0.21/*WWG*/,0.04/*WZG*/, 18610, 0.143, 750/*TTBar*/,63.21/*WW*/,22.82/*WZ*/, 61526/*W*/, 108000000/*QCDEM30*/, 54120000/*QCDEM40*/, 1/*MU*/, 16792/*GJet*/,1};



// Using GenXSecAnalyzer using full datasets
// WGToLNuG is from CMS-XSDB as GenXSecAnalyser didn't work

//double MC_XS[25] = {1, 412.7/*WGToLNuG*/, 19.78 /*WGJet40*/, 0.8094/*WGJet130*/, 51.37 /* 51.43 ZGInclusive*/, 6417/*DY*/, 3.782/*TTG*/, 0.3369/*WWG*/, 0.07876/*WZG*/, 18610, 0.143, 754.6/*TTBar*/, 75.9/*WW*/, 27.57/*WZ*/, 53870/*W*/, 108000000/*QCDEM30*/, 113100.0 /*QCDEM40*/, 1/*MU*/, 874.2/*GJet*/, 0.001027 /*T5Wg_mG-1800_mLSP-800*/, 0.001019 /*TChiWg_mChi-1000_mLSP-1*/, 0.001057/*T5Wg_mG-1800_mLSP-1600*/,   0.005763 /*T5Wg_mG-1500_mLSP-1*/, 1};

// from XSDB; using k-factor = 1.34 for Wgamma samples

double MC_XS[25] = {1, 551.008/*WGToLNuG*/, 26.5454/*WGJet40*/, 1.083658/*WGJet130*/, 51.53 /*ZGInclusive*/, 6424.0/*DY*/, 3.774/*TTG*/, 0.3369/*WWG*/, 0.07876/*WZG*/, 18610, 0.143, 734.6/*TTBar*/, 75.95/*WW*/, 27.59/*WZ*/, 53870/*W*/, 108000000/*QCDEM30*/, 113100.0 /*QCDEM40*/, 1/*MU*/, 874.2/*GJet*/, 0.001027 /*T5Wg_mG-1800_mLSP-800*/, 0.001019 /*TChiWg_mChi-1000_mLSP-1*/, 0.001057/*T5Wg_mG-1800_mLSP-1600*/,   0.005763 /*T5Wg_mG-1500_mLSP-1*/, 1};


/* 117400.0 QCDEM40 for 2017*/

//0.143 // LO ZG130
//0.1404 // NLO ZG130


float getEvtWeight(int year, double XSec, double nEvents_MC){
    
    float luminosity;
    if(year==2016) luminosity = 35870;
    if(year==2017) luminosity = 41540;
    //if(year==2017) luminosity = 27130;
    if(year==2018) luminosity = 59960;
    double evtWeight = 1.;
	
    evtWeight = XSec * luminosity / nEvents_MC;
    //cout << "Using event weight " << evtWeight << endl;
    //cout << "XS = " << evtWeight/luminosity*nEvents_MC << endl;
    //cout << "lumi = " << luminosity << endl;
    //cout << "nEvents_MC = " << nEvents_MC << endl;
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
        pileupweights[99] = 0;

	if(nvertex > 99)return 0;
	else return pileupweights[nvertex]; 
}

float getPUESF16(int nvertex){
	
	float pileupweights[100];
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
        pileupweights[99] = 0;
	if(nvertex > 99)return 0;
	else return pileupweights[nvertex]; 
}
float getPUESF17(int nvertex){
	
	float pileupweights[100];	
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
        pileupweights[99] = 0;
	if(nvertex > 99)return 0;
	else return pileupweights[nvertex]; 
}

float getPUESF18(int nvertex){
	
	float pileupweights[100];
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
        pileupweights[99] = 0;
	
	if(nvertex > 99)return 0;
	else return pileupweights[nvertex]; 
}

//Double_t bkgEtBins[]={35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400,500,800};
Double_t bkgEtBins[]={35,40,50,60,70,80,90,100,110,120,130,140,150,160,170,185,200,215,230,250,275,290, 305,325,345,370,400,500,800};
int nBkgEtBins= sizeof(bkgEtBins)/sizeof(bkgEtBins[0]) -1;
//Double_t bkgPtBins[]={25,30,35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400,500,800};
Double_t bkgPtBins[]={25,50,75,100,125,150,200,400,800};
int nBkgPtBins= sizeof(bkgPtBins)/sizeof(bkgPtBins[0])-1;
Double_t bkgMETBins[]={0,40,60,80,100,120,140,160,180,210,240,280,320,400,600,1000};
int nBkgMETBins= sizeof(bkgMETBins)/sizeof(bkgMETBins[0]) -1;
Double_t bkgMtBins[]={0,20,40,60,80,100,120,140,160,180,200,300,400,500,1000};
int nBkgMtBins= sizeof(bkgMtBins)/sizeof(bkgMtBins[0]) -1;
Double_t bkgHTBins[]={0,40,60,80,100,120,140,160,180,200,225,250,275,300,340,380,420,500,600,1000};
int nBkgHTBins= sizeof(bkgHTBins)/sizeof(bkgHTBins[0]) -1;


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
