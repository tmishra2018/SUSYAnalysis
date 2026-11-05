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


        	// config chosen by the chain directory we run from: only the path suffix is matched, and the
        	// analysis base directory (anaBase) is taken from the same path, so the area can be copied or moved
        	std::string anaBase;
        	const std::pair<std::string,std::string> confs[] = {
        		{"/test/myResult", "SigConfig.txt"}, {"/test/Background", "BkgPredConfig.txt"}, {"/test/ISRReweight", "ISRConfig.txt"}};
        	for(const auto &c : confs){
        		size_t pos = current_path.find(c.first);
        		if(pos != std::string::npos){
        			anaBase = current_path.substr(0, pos);
        			configFilePath = anaBase + c.first + "/" + c.second;
        			break;
        		}
        	}
		if(configFilePath.empty()){ std::cerr << "Error: Unknown execution path!" << std::endl;
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
			electronIDFile = TFile::Open((anaBase + "/SFs/2016preVFP/egammaEffi.txt_Ele_Medium_EGM2D_UL16_preVFP.root").c_str());
			photonIDFile   = TFile::Open((anaBase + "/SFs/2016preVFP/egammaEffi.txt_EGM2D_Pho_Loose_UL16_preVFP.root").c_str());
			muonIDFile     = TFile::Open((anaBase + "/SFs/2016preVFP/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root").c_str());
			photonTRGFile  = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/diphoton_pholeg_2016preVFP.root");
			electronTRGFile = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/diphoton_eleg_2016preVFP.root");
			muonTRGFile = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/muonphoton_trigger_2016preVFP.root");
		}
		else if(RunYear==2016 and preVFP==0){
			// pt range 500 GeV and eta range 2.5
			electronIDFile = TFile::Open((anaBase + "/SFs/2016/egammaEffi.txt_Ele_Medium_EGM2D_UL16_postVFP.root").c_str());
			// pt range 500 GeV and eta range 2.5
			photonIDFile   = TFile::Open((anaBase + "/SFs/2016/egammaEffi.txt_EGM2D_Pho_Loose_UL16_postVFP.root").c_str());
			// pt range 120 GeV and eta range 2.4
			muonIDFile     = TFile::Open((anaBase + "/SFs/2016/Efficiencies_muon_generalTracks_Z_Run2016_UL_ID.root").c_str());
			photonTRGFile  = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/diphoton_pholeg_2016postVFP.root");
			electronTRGFile = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/diphoton_eleg_2016postVFP.root");
			muonTRGFile = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/muonphoton_trigger_2016postVFP.root");
		}
			
		else if(RunYear==2017){
			electronIDFile = TFile::Open((anaBase + "/SFs/2017/egammaEffi.txt_Ele_Medium_EGM2D_UL17.root").c_str());
			photonIDFile   = TFile::Open((anaBase + "/SFs/2017/egammaEffi.txt_EGM2D_PHO_Loose_UL17.root").c_str());
			muonIDFile     = TFile::Open((anaBase + "/SFs/2017/Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root").c_str());
			photonTRGFile  = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/diphoton_pholeg_2017.root");
			electronTRGFile = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/diphoton_eleg_2017.root");
			muonTRGFile = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/muonphoton_trigger_2017.root");
		}
		else if(RunYear==2018){
			electronIDFile = TFile::Open((anaBase + "/SFs/2018/egammaEffi.txt_Ele_Medium_EGM2D_UL18.root").c_str());
			photonIDFile   = TFile::Open((anaBase + "/SFs/2018/egammaEffi.txt_EGM2D_Pho_Loose_UL18.root").c_str());
			muonIDFile     = TFile::Open((anaBase + "/SFs/2018/Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root").c_str());
			photonTRGFile  = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/diphoton_pholeg_2018.root");
			electronTRGFile = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/diphoton_eleg_2018.root");
			muonTRGFile = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/files/muonphoton_trigger_2018.root");
		}
	
		// UL electron R9 and miniIso SFs (egamma-tnp, test/ScaleFactor/tnp)
		std::string sfTag = ((RunYear==2016 and preVFP==1) or RunYear==678) ? "2016preVFP" : (RunYear==2016 ? "2016postVFP" : std::to_string(RunYear));
		std::string tnpDir = anaBase + "/test/ScaleFactor/tnp/results/" + sfTag + "/";
		TFile *eleR9File          = TFile::Open((tnpDir + "SF_R9_" + sfTag + "_hybrid.root").c_str());
		TFile *eleMiniIsoFile     = TFile::Open((tnpDir + "SF_miniIso_" + sfTag + "_hybrid.root").c_str());
		// UL photon R9 > 0.5 SF (Z->ee probes reconstructed as photons, tnp/run_phoR9.sh)
		TFile *phoR9File          = TFile::Open((tnpDir + "SF_phoR9_" + sfTag + "_hybrid.root").c_str());
		
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

		if(!eleMiniIsoFile || eleMiniIsoFile->IsZombie()){
			std::cout << "no electron miniIso file" << std::endl;
		}
		else{
			electronMiniIsoESF = (TH2F*)eleMiniIsoFile->Get("EGamma_SF2D");
		}

		if(!phoR9File || phoR9File->IsZombie()){
			std::cout << "no photon R9 file" << std::endl;
		}
		else{
			photonR9ESF = (TH2F*)phoR9File->Get("EGamma_SF2D");
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
	float getMiniIsoESF(float pt, float eta);
	float getMiniIsoESFError(float pt, float eta);
	float getPhotonR9ESF(float pt, float eta);
	float getPhotonR9ESFError(float pt, float eta);

	private:
		TH2F *electronIDESF;
	//	TH2F *electronISOESF;
		TH2F *photonIDESF;
		TH2F *muonIDESF;
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
		TH2F *electronMiniIsoESF;
		TH2F *photonR9ESF;
		float tnpSFBin(TH2F *h, float pt, float eta, bool err);
};

#endif
