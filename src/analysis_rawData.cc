#include "../include/analysis_rawData.h"

bool rawData::passHLT(){
  bool passDataHLT(false);
  
  switch(runtype_){
    case MC: passDataHLT = true; break;

    case DoubleEG2016: if(((HLTPho >> 14) &1) !=0)passDataHLT = true; break; // HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90
    case DoubleEG2017: if(((HLTPho >> 14) &1) !=0)passDataHLT = true; break; // HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90
    case DoubleEG2018: if(((HLTPho >> 14) &1) !=0)passDataHLT = true; break; // HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90
    case MCDoubleEG2016: if(((HLTPho >> 14) &1) !=0)passDataHLT = true; break;
    case MCDoubleEG2017: if(((HLTPho >> 14) &1) !=0)passDataHLT = true; break;
    case MCDoubleEG2018: if(((HLTPho >> 14) &1) !=0)passDataHLT = true; break;
    
    case MuonEG2016: if(((HLTEleMuX >> 8) &1) !=0 || ((HLTEleMuX >> 51) &1) !=0)passDataHLT = true; break;  // HLT_Mu17_Photon30_CaloIdL_L1ISO_v || HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL
    case MuonEG2017: if(((HLTEleMuX >> 8) &1) !=0 || ((HLTEleMuX >> 57) &1) !=0)passDataHLT = true; break;  // HLT_Mu17_Photon30_IsoCaloId_v || HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL
    case MuonEG2018: if(((HLTEleMuX >> 8) &1) !=0 || ((HLTEleMuX >> 57) &1) !=0)passDataHLT = true; break;  // HLT_Mu17_Photon30_IsoCaloId_v || HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL
    case MCMuonEG2016: if(((HLTEleMuX >> 8) &1) !=0 || ((HLTEleMuX >> 51) &1) !=0)passDataHLT = true; break;
    case MCMuonEG2017: if(((HLTEleMuX >> 8) &1) !=0 || ((HLTEleMuX >> 57) &1) !=0)passDataHLT = true; break;
    case MCMuonEG2018: if(((HLTEleMuX >> 8) &1) !=0 || ((HLTEleMuX >> 57) &1) !=0)passDataHLT = true; break;

    // https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary#2016
    // https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary#2017
    // https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary#2018
    case SingleElectron2016:if(((HLTEleMuX >> 4) &1) !=0)passDataHLT = true; break;  // HLT_Ele27_WPTight_Gsf_v
    case SingleElectron2017:if(((HLTEleMuX >> 3) &1) !=0)passDataHLT = true; break;  // HLT_Ele35_WPTight_Gsf_v
    case SingleElectron2018:if(((HLTEleMuX >> 55) &1) !=0)passDataHLT = true; break; // HLT_Ele32_WPTight_Gsf_v
    
    case MCSingleElectron2016:if(((HLTEleMuX >> 4) &1) !=0)passDataHLT = true; break;  // HLT_Ele27_WPTight_Gsf_v
    case MCSingleElectron2017:if(((HLTEleMuX >> 3) &1) !=0)passDataHLT = true; break;  // HLT_Ele35_WPTight_Gsf_v
    case MCSingleElectron2018:if(((HLTEleMuX >> 55) &1) !=0)passDataHLT = true; break; // HLT_Ele32_WPTight_Gsf_v
    
  // https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2016
  // https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2017
  // https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2018
    case SingleMuon2016: if(((HLTEleMuX >> 19) &1)!=0 || ((HLTEleMuX >> 20) &1)!=0)passDataHLT = true; break; // HLT_IsoMu24_v || HLT_IsoTkMu24_v
    case SingleMuon2017: if(((HLTEleMuX >> 56) &1)!=0)passDataHLT = true; break;	// HLT_IsoMu27_v
    case SingleMuon2018: if(((HLTEleMuX >> 19) &1)!=0)passDataHLT = true; break;	// HLT_IsoMu24_v
    
    case MCSingleMuon2016: if(((HLTEleMuX >> 19) &1)!=0 || ((HLTEleMuX >> 20) &1)!=0)passDataHLT = true; break;  // HLT_IsoMu24_v || HLT_IsoTkMu24_v
    case MCSingleMuon2017: if(((HLTEleMuX >> 56) &1)!=0)passDataHLT = true; break;	// HLT_IsoMu27_v
    case MCSingleMuon2018: if(((HLTEleMuX >> 19) &1)!=0)passDataHLT = true; break;	// HLT_IsoMu24_v
    
    case DoubleMuon2016: if(((HLTEleMuX >>14 ) &1) !=0 || ((HLTEleMuX >>15 ) &1)!=0 || ((HLTEleMuX >>16 ) &1)!=0 )passDataHLT = true; break;
    case DoubleMuon2017: if(((HLTEleMuX >>14 ) &1) !=0 || ((HLTEleMuX >>15 ) &1)!=0 || ((HLTEleMuX >>16 ) &1)!=0 )passDataHLT = true; break;
    case DoubleMuon2018: if(((HLTEleMuX >>14 ) &1) !=0 || ((HLTEleMuX >>15 ) &1)!=0 || ((HLTEleMuX >>16 ) &1)!=0 )passDataHLT = true; break;
    case MCDoubleMuon2016: if(((HLTEleMuX >>14 ) &1) !=0 || ((HLTEleMuX >>15 ) &1)!=0 || ((HLTEleMuX >>16 ) &1)!=0 )passDataHLT = true; break;
    case MCDoubleMuon2017: if(((HLTEleMuX >>14 ) &1) !=0 || ((HLTEleMuX >>15 ) &1)!=0 || ((HLTEleMuX >>16 ) &1)!=0 )passDataHLT = true; break;
    case MCDoubleMuon2018: if(((HLTEleMuX >>14 ) &1) !=0 || ((HLTEleMuX >>15 ) &1)!=0 || ((HLTEleMuX >>16 ) &1)!=0 )passDataHLT = true; break;

    case MET2016: if(((HLTJet >>8 ) &1) !=0 || ((HLTJet >>24 ) &1)!=0 || ((HLTJet >>25 ) &1)!=0)passDataHLT = true; break;
    case MET2017: if(((HLTJet >>8 ) &1) !=0 || ((HLTJet >>24 ) &1)!=0 || ((HLTJet >>25 ) &1)!=0)passDataHLT = true; break;
    case MET2018: if(((HLTJet >>8 ) &1) !=0 || ((HLTJet >>24 ) &1)!=0 || ((HLTJet >>25 ) &1)!=0)passDataHLT = true; break;
    case MCMET2016: if(((HLTJet >>8 ) &1) !=0 || ((HLTJet >>24 ) &1)!=0 || ((HLTJet >>25 ) &1)!=0)passDataHLT = true; break;
    case MCMET2017: if(((HLTJet >>8 ) &1) !=0 || ((HLTJet >>24 ) &1)!=0 || ((HLTJet >>25 ) &1)!=0)passDataHLT = true; break;
    case MCMET2018: if(((HLTJet >>8 ) &1) !=0 || ((HLTJet >>24 ) &1)!=0 || ((HLTJet >>25 ) &1)!=0)passDataHLT = true; break;
    
    default: break;
  }
  
  return passDataHLT;
    
}

bool rawData::passMETFilter(int RunYear, int filter) {
    bool passfilter = true;
    for (int im = 1; im <= 9; im++) {  
        // Skip "Flag_ecalBadCalibFilter" (index 8) for year 2016
        if (im == 8 && RunYear == 2016) continue;

        if (((filter >> im) & 1) != 0) { 
            passfilter = false; 
            return passfilter;  
        }
    }
    return passfilter; 
}

void rawData::failFilterStepHistogram(int RunYear, int filter, TH1D* h_failFilters) {
    for (int im = 1; im <= 9; im++) {
        if (im == 8 && RunYear == 2016) continue;  // Skip filter #8 for 2016 data
        
        if ((filter >> im) & 1) {
            h_failFilters->Fill(im);
        }
    }
}

int rawData::failFilterStep(int RunYear, int filter) {
    int failStep = -1;
    for (int im = 1; im <= 9; im++) {  
        // Skip "Flag_ecalBadCalibFilter" (index 8) for year 2016
        if (im == 8 && RunYear == 2016) continue;

        if (((filter >> im) & 1) != 0) { 
            failStep = im;  
            return failStep;  // Return the first failed filter index
        }
    }
    return failStep;  // If no filter fails, return -1
}
