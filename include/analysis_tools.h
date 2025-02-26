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

double MC_XS[25] = {1, 551.008/*WGToLNuG*/, 26.5454/*WGJet40*/, 1.083658/*WGJet130*/, 51.53 /*ZGInclusive*/, 6424.0/*DY*/, 3.774/*TTG*/, 0.3369/*WWG*/, 0.07876/*WZG*/, 18610, 0.143, 734.6/*TTBar*/, 75.95/*WW*/, 27.55/*WZ*/, 61526.7/*W*/, 108000000/*QCDEM30*/, 113100.0 /*QCDEM40*/, 1/*MU*/, 874.2/*GJet*/, 0.001027 /*T5Wg_mG-1800_mLSP-800*/, 0.001019 /*TChiWg_mChi-1000_mLSP-1*/, 0.001057/*T5Wg_mG-1800_mLSP-1600*/,   0.005763 /*T5Wg_mG-1500_mLSP-1*/, 1};


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
        pileupweights[99] = 0;

	if(nvertex > 99)return 0;
	else return pileupweights[nvertex]; 
}

float getPUESF16(int nvertex){
	
	float pileupweights[100];
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
        pileupweights[99] = 0;

	if(nvertex > 99)return 0;
	else return pileupweights[nvertex]; 
}
float getPUESF17(int nvertex){
	
	float pileupweights[100];	
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
        pileupweights[99] = 0;

	if(nvertex > 99)return 0;
	else return pileupweights[nvertex]; 
}

float getPUESF18(int nvertex){
	
	float pileupweights[100];
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
