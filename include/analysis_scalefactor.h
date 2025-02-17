#ifndef ANALYSIS_SCALEFACTOR 
#define ANALYSIS_SCALEFACTOR

#ifndef ROOT_TFile
#include "TFile.h"
#endif
#ifndef ROOT_TH2F
#include "TH2F.h"
#endif
#ifndef ROOT_TAxis
#include "TAxis.h"
#endif
#include<iostream>
#include<fstream>
#include <filesystem>
namespace fs = std::filesystem;

class esfScaleFactor{
  public:

    esfScaleFactor(){
		int RunYear=0;
		bool preVFP=0;

		std::string current_path = std::filesystem::current_path().string();
        	std::string configFilePath;

        	if (current_path.find("/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/test/myResult") != std::string::npos) {
            		configFilePath = "/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/test/myResult/SigConfig.txt";}
        	else if (current_path.find("/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/test/Background") != std::string::npos) {
            		configFilePath = "/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/test/Background/BkgPredConfig.txt";} 
		else {	std::cerr << "Error: Unknown execution path!" << std::endl;
            		return;}
        	std::ifstream configfile(configFilePath);
        	std::string conftype;
        	double confvalue;

        	if(configfile.is_open()){
        		for(int i(0); i<11; i++){
                        	configfile >> conftype >> confvalue;
                        	if(conftype.find("RunYear")!=std::string::npos)RunYear = confvalue;
                        	if(conftype.find("preVFP")!=std::string::npos)preVFP = confvalue;
          		}
        	}
        	configfile.close(); 

		TFile *electronIDFile=0;
		TFile *photonIDFile=0;
		TFile *muonIDFile=0;
		TFile *photonTRGFile=0;
		TFile *electronTRGFile=0;
		TFile *muonTRGFile=0;
		
		if((RunYear==2016 and preVFP==1) or RunYear==678){
			electronIDFile = TFile::Open("/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/SFs/2016preVFP/egammaEffi.txt_Ele_Medium_EGM2D_UL16_preVFP.root");
			photonIDFile   = TFile::Open("/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/SFs/2016preVFP/egammaEffi.txt_EGM2D_Pho_Loose_UL16_preVFP.root");
			muonIDFile     = TFile::Open("/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/SFs/2016preVFP/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root");
			photonTRGFile  = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/diphoton_pholeg_2016preVFP.root");
			electronTRGFile = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/diphoton_eleg_2016preVFP.root");
			muonTRGFile = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/muonphoton_trigger_2016preVFP.root");
		}
		else if(RunYear==2016 and preVFP==0){
			// pt range 500 GeV and eta range 2.5
			electronIDFile = TFile::Open("/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/SFs/2016/egammaEffi.txt_Ele_Medium_EGM2D_UL16_postVFP.root");
			// pt range 500 GeV and eta range 2.5
			photonIDFile   = TFile::Open("/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/SFs/2016/egammaEffi.txt_EGM2D_Pho_Loose_UL16_postVFP.root");
			// pt range 120 GeV and eta range 2.4
			muonIDFile     = TFile::Open("/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/SFs/2016/Efficiencies_muon_generalTracks_Z_Run2016_UL_ID.root");
			photonTRGFile  = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/diphoton_pholeg_2016postVFP.root");
			electronTRGFile = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/diphoton_eleg_2016postVFP.root");
			muonTRGFile = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/muonphoton_trigger_2016postVFP.root");
		}
			
		else if(RunYear==2017){
			electronIDFile = TFile::Open("/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/SFs/2017/egammaEffi.txt_Ele_Medium_EGM2D_UL17.root");
			photonIDFile   = TFile::Open("/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/SFs/2017/egammaEffi.txt_EGM2D_PHO_Loose_UL17.root");
			muonIDFile     = TFile::Open("/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/SFs/2017/Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root");
			photonTRGFile  = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/diphoton_pholeg_2017.root");
			electronTRGFile = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/diphoton_eleg_2017.root");
			muonTRGFile = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/muonphoton_trigger_2016postVFP.root");
		}
		else if(RunYear==2018){
			electronIDFile = TFile::Open("/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/SFs/2018/egammaEffi.txt_Ele_Medium_EGM2D_UL18.root");
			photonIDFile   = TFile::Open("/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/SFs/2018/egammaEffi.txt_EGM2D_Pho_Loose_UL18.root");
			muonIDFile     = TFile::Open("/uscms/homes/t/tmishra/work/CMSSW_14_0_7/src/SUSYAnalysis/SFs/2018/Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root");
			photonTRGFile  = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/diphoton_pholeg_2018.root");
			electronTRGFile = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/diphoton_eleg_2018.root");
			muonTRGFile = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/muonphoton_trigger_2016postVFP.root");
		}
	
		TFile *muonIsoFile    = TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/TnP_NUM_MiniIsoTight_DENOM_MediumID_VAR_map_pt_eta.root");
		TFile *muonIPFile     = TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/TnP_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root");
		TFile *eleR9File          = TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/R9_SF.root");
		
		TFile *fastelectronIDFile = TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/sf_el_mediumCB.root");
		TFile *fastelectronIsoFile= TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/sf_el_mini01.root");
		TFile *fastmuonIDFile     = TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/sf_mu_mediumID.root");
		TFile *fastmuonIsoFile    = TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/sf_mu_mediumID_mini02.root");
		TFile *fastmuonIPFile     = TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/sf_mu_mediumID_tightIP2D.root");

		if(electronIDFile->IsZombie() || fastelectronIDFile->IsZombie() || fastelectronIsoFile->IsZombie()){
			std::cout << "no electron ESF file" << std::endl;
		}
		else{
			electronIDESF = (TH2F*)electronIDFile->Get("EGamma_SF2D"); // Run2016_CutBasedVetoNoIso94XV2
//			electronISOESF= (TH2F*)electronIDFile->Get("MVAVLooseElectronToMini");  // Run2016_Mini
			electronIDFast= (TH2F*)fastelectronIDFile->Get("histo2D");
			electronISOFast=(TH2F*)fastelectronIsoFile->Get("histo2D");
		}
		
		if(eleR9File ->IsZombie()){
			std::cout << "no electron R9 file" << std::endl;
		}
		else{
			electronR9ESF = (TH2F*)eleR9File->Get("EGamma_SF2D");
		}

		if(photonIDFile->IsZombie()){
			std::cout << "no photon ESF file" << std::endl;
		}
		else{
			photonIDESF   = (TH2F*)photonIDFile->Get("EGamma_SF2D");
		}

		if(muonIDFile->IsZombie() || fastmuonIDFile->IsZombie() || fastmuonIsoFile->IsZombie() || fastmuonIPFile->IsZombie()){
      std::cout << "no muon ESF file" << std::endl;
    }
		else{
			muonIDESF    = (TH2F*)muonIDFile->Get("NUM_MediumPromptID_DEN_TrackerMuons_abseta_pt");
			muonIDFast   = (TH2F*)fastmuonIDFile->Get("histo2D");
			muonISOFast  = (TH2F*)fastmuonIsoFile->Get("histo2D");
			muonIPFast   = (TH2F*)fastmuonIPFile->Get("histo2D");
		}

		if(muonIsoFile->IsZombie()){
      std::cout << "no muon Iso ESF file" << std::endl;
    }
		else{
			muonISOESF    = (TH2F*)muonIsoFile->Get("SF");
			//muonISOESF    = (TH2F*)muonIsoFile->Get("cr");
		}

		if(muonIPFile->IsZombie()){
      std::cout << "no muon ESF file" << std::endl;
    }
		else{
			muonIPESF    = (TH2F*)muonIPFile->Get("SF");
		}


		if(photonTRGFile->IsZombie()){
			std::cout << "no photon TRG ESF file" << std::endl;
		}
		else{
			photonTRGESF   = (TH2F*)photonTRGFile->Get("ESF_Lead");
			photonTRGEFF   = (TH2F*)photonTRGFile->Get("Leadeff_data");
		}

		if(electronTRGFile->IsZombie()){
			std::cout << "no electron TRG ESF file" << std::endl;
		}
		else{
			electronTRGESF   = (TH2F*)electronTRGFile->Get("p_TriggerScale");
			electronTRGEFF   = (TH2F*)electronTRGFile->Get("Traileff_data");
		}

		if(muonTRGFile->IsZombie()){
			std::cout << "no muon TRG ESF file" << std::endl;
		}
		else{
			muonTRGESF   = (TH2F*)muonTRGFile->Get("p_mgESF");
			muonTRGEFF   = (TH2F*)muonTRGFile->Get("p_mgEff");
		}
	}

	~esfScaleFactor(){
	};

	float getElectronESF(float pt, float eta);
	float getPhotonESF(float pt, float eta);
	float getMuonESF(float pt, float eta);
	float getElectronESFError(float pt, float eta);
	float getPhotonESFError(float pt, float eta);
	float getMuonESFError(float pt, float eta);
	float getElectronTRGESF(float pt, float eta);
	float getegPhotonTRGESF(float pt, float eta);
	float getMuonEGTRGESF(float pt1, float pt2);
	float getegPhotonTRGESFError(float pt, float eta);
	float getElectronTRGESFError(float pt, float eta);
	float getMuonEGTRGESFError(float pt1, float pt2);

	float getFastElectronESF(float pt, float eta);
	float getFastPhotonESF(float pt, float eta);
	float getFastMuonESF(float pt, float eta);
	float getFastElectronESFError(float pt, float eta);
	float getFastPhotonESFError(float pt, float eta);
	float getFastMuonESFError(float pt, float eta);
	float getFastElectronTRGESF(float pt, float eta);
	float getFastegPhotonTRGESF(float pt, float eta);
	float getFastMuonEGTRGESF(float pt1, float pt2);
	float getFastegPhotonTRGESFError(float pt, float eta);
	float getFastElectronTRGESFError(float pt, float eta);
	float getFastMuonEGTRGESFError(float pt1, float pt2);
	float getR9ESF(float pt, float eta);
	float getR9ESFError(float pt, float eta);

	private:
		TH2F *electronIDESF;
	//	TH2F *electronISOESF;
		TH2F *photonIDESF;
		TH2F *muonIDESF;
		TH2F *muonISOESF;
		TH2F *muonIPESF;
		TH2F *electronIDFast;
		TH2F *electronISOFast;
		TH2F *muonIDFast;
		TH2F *muonISOFast;
		TH2F *muonIPFast;
		TH2F *electronTRGESF;
		TH2F *photonTRGESF;
		TH2F *muonTRGESF;
		TH2F *electronTRGEFF;
		TH2F *photonTRGEFF;
		TH2F *muonTRGEFF;
		TH2F *electronR9ESF;
};

#endif
