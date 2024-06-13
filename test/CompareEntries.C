void CompareEntries()
{
	char egfile[100] = "/uscms_data/d3/mengleis/FullStatusOct/resTree_egsignal_DoubleEG_ReMiniAOD_FullEcal_newEta.root";
	char mgfile[100] = "/uscms_data/d3/mengleis/FullStatusOct/resTree_mgsignal_MuonEG_FullEcal.root";

	char ULegfile[100] = "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2016.root";
	char ULmgfile[100] = "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016.root";
	//char ULegfile[100] = "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2017.root";
	//char ULmgfile[100] = "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2017_April.root";

	TFile egFile(egfile);
	TFile mgFile(mgfile);
	TFile ULegFile(ULegfile);
	TFile ULmgFile(ULmgfile);

	TTree* EGsignalTree;
	TTree* EGjetTree;
	TTree* EGhadronTree;
	TTree* EGproxyTree;
	TTree* EGfakeLepTree;
	EGsignalTree = (TTree*)egFile.Get("signalTree");
	EGjetTree = (TTree*)egFile.Get("jetTree");
	EGhadronTree = (TTree*)egFile.Get("hadronTree");
	EGproxyTree = (TTree*)egFile.Get("proxyTree");
	EGfakeLepTree = (TTree*)egFile.Get("fakeLepTree");

	float SUS17_egSignal, SUS17_egJet, SUS17_egHadron, SUS17_egProxy, SUS17_egFakeLep;	
	SUS17_egSignal = EGsignalTree->GetEntries();
	SUS17_egJet = EGjetTree->GetEntries();
	SUS17_egHadron = EGhadronTree->GetEntries();
	SUS17_egProxy = EGproxyTree->GetEntries();
	SUS17_egFakeLep = EGfakeLepTree->GetEntries();

	TTree* MGsignalTree;
	TTree* MGjetTree;
	TTree* MGhadronTree;
	TTree* MGproxyTree;
	TTree* MGfakeLepTree;
	MGsignalTree = (TTree*)mgFile.Get("signalTree");
	MGjetTree = (TTree*)mgFile.Get("jetTree");
	MGhadronTree = (TTree*)mgFile.Get("hadronTree");
	MGproxyTree = (TTree*)mgFile.Get("proxyTree");
	MGfakeLepTree = (TTree*)mgFile.Get("fakeLepTree");

	float SUS17_mgSignal, SUS17_mgJet, SUS17_mgHadron, SUS17_mgProxy, SUS17_mgFakeLep;	
	SUS17_mgSignal = MGsignalTree->GetEntries();
	SUS17_mgJet = MGjetTree->GetEntries();
	SUS17_mgHadron = MGhadronTree->GetEntries();
	SUS17_mgProxy = MGproxyTree->GetEntries();
	SUS17_mgFakeLep = MGfakeLepTree->GetEntries();

/*******************************************************************************************************************/	
	TTree* ULEGsignalTree;
	TTree* ULEGjetTree;
	TTree* ULEGhadronTree;
	TTree* ULEGproxyTree;
	TTree* ULEGfakeLepTree;
	ULEGsignalTree = (TTree*)ULegFile.Get("signalTree");
	ULEGjetTree = (TTree*)ULegFile.Get("jetTree");
	ULEGhadronTree = (TTree*)ULegFile.Get("hadronTree");
	ULEGproxyTree = (TTree*)ULegFile.Get("proxyTree");
	ULEGfakeLepTree = (TTree*)ULegFile.Get("fakeLepTree");

	float UL_egSignal, UL_egJet, UL_egHadron, UL_egProxy, UL_egFakeLep;	
	UL_egSignal = ULEGsignalTree->GetEntries();
	UL_egJet = ULEGjetTree->GetEntries();
	UL_egHadron = ULEGhadronTree->GetEntries();
	UL_egProxy = ULEGproxyTree->GetEntries();
	UL_egFakeLep = ULEGfakeLepTree->GetEntries();
	
	TTree* ULMGsignalTree;
	TTree* ULMGjetTree;
	TTree* ULMGhadronTree;
	TTree* ULMGproxyTree;
	TTree* ULMGfakeLepTree;
	ULMGsignalTree = (TTree*)ULmgFile.Get("signalTree");
	ULMGjetTree = (TTree*)ULmgFile.Get("jetTree");
	ULMGhadronTree = (TTree*)ULmgFile.Get("hadronTree");
	ULMGproxyTree = (TTree*)ULmgFile.Get("proxyTree");
	ULMGfakeLepTree = (TTree*)ULmgFile.Get("fakeLepTree");

	float UL_mgSignal, UL_mgJet, UL_mgHadron, UL_mgProxy, UL_mgFakeLep;	
	UL_mgSignal = ULMGsignalTree->GetEntries();
	UL_mgJet = ULMGjetTree->GetEntries();
	UL_mgHadron = ULMGhadronTree->GetEntries();
	UL_mgProxy = ULMGproxyTree->GetEntries();
	UL_mgFakeLep = ULMGfakeLepTree->GetEntries();

	cout<<"eg channel :: SUS17 versus UL"<<endl<<endl;
	cout<<"signal tree :  	"<<EGsignalTree->GetEntries()<<"\t"<<ULEGsignalTree->GetEntries()<<"\t"<<1.*(EGsignalTree->GetEntries()-ULEGsignalTree->GetEntries())/EGsignalTree->GetEntries()<<endl;
	cout<<"jet tree :  	"<<EGjetTree->GetEntries()<<"\t"<<ULEGjetTree->GetEntries()<<"\t"<<1.*(EGjetTree->GetEntries()-ULEGjetTree->GetEntries())/EGjetTree->GetEntries()<<endl;
	cout<<"proxy tree :  	"<<EGproxyTree->GetEntries()<<"\t"<<ULEGproxyTree->GetEntries()<<"\t"<<1.*(EGproxyTree->GetEntries()-ULEGproxyTree->GetEntries())/EGproxyTree->GetEntries()<<endl;
	cout<<"fakeLep tree :  "<<EGfakeLepTree->GetEntries()<<"\t"<<ULEGfakeLepTree->GetEntries()<<"\t"<<1.*(EGfakeLepTree->GetEntries()-ULEGfakeLepTree->GetEntries())/EGfakeLepTree->GetEntries()<<endl;
	cout<<"hadron tree :  	"<<EGhadronTree->GetEntries()<<"\t"<<ULEGhadronTree->GetEntries()<<"\t"<<1.*(EGhadronTree->GetEntries()-ULEGhadronTree->GetEntries())/EGhadronTree->GetEntries()<<endl;

	cout<<"\n\n"<<"mg channel :: SUS17 versus UL"<<endl<<endl;
	cout<<"signal tree :  	"<<MGsignalTree->GetEntries()<<"\t"<<ULMGsignalTree->GetEntries()<<"\t"<<1.*(MGsignalTree->GetEntries()-ULMGsignalTree->GetEntries())/MGsignalTree->GetEntries()<<endl;
	cout<<"jet tree :  	"<<MGjetTree->GetEntries()<<"\t"<<ULMGjetTree->GetEntries()<<"\t"<<1.*(MGjetTree->GetEntries()-ULMGjetTree->GetEntries())/MGjetTree->GetEntries()<<endl;
	cout<<"proxy tree :  	"<<MGproxyTree->GetEntries()<<"\t"<<ULMGproxyTree->GetEntries()<<"\t"<<1.*(MGproxyTree->GetEntries()-ULMGproxyTree->GetEntries())/MGproxyTree->GetEntries()<<endl;
	cout<<"fakeLep tree :  "<<MGfakeLepTree->GetEntries()<<"\t"<<ULMGfakeLepTree->GetEntries()<<"\t"<<1.*(MGfakeLepTree->GetEntries()-ULMGfakeLepTree->GetEntries())/MGfakeLepTree->GetEntries()<<endl;
	cout<<"hadron tree :  	"<<MGhadronTree->GetEntries()<<"\t"<<ULMGhadronTree->GetEntries()<<"\t"<<1.*(MGhadronTree->GetEntries()-ULMGhadronTree->GetEntries())/MGhadronTree->GetEntries()<<endl;
}
