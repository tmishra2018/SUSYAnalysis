void modify_XSec(){

       TString process[11] = {"DYJetsToLL", "TTJets", "TTGJets", "WWG", "WZG", "WW", "WZ", "WGToLNuG", "WGJet40", "WGJet130", "ZGToLLG"};
       //float  ModifiedXS[11] = {5670, 750, 3.697, 0.2147, 0.04345, 75.94, 27.564, 489.0, 17.018, 0.87971,  117.8}; // menglei
       float  ModifiedXS[11] = {6424.0, 734.6, 3.774, 0.3369, 0.07876, 75.95, 27.59, 551.008, 26.5454, 1.083658,  51.53}; // XSDB
       //float  ModifiedXS[11] = {5670.0, 750., 3.697, 0.21, 0.04, 75.9, 27.57, 489., 17.018, 0.87971, 51.43}; // changing WWG, WZG, TTJets, TTGJets
       
       //float  ModifiedXS[11] = {6417.0, 754.6, 3.782, 0.3369, 0.07876, 75.9, 27.57, 412.7, 19.78, 0.8094, 51.43}; // GenXSecAnalyser
       for (int i=0; i<=10; i++)
       {
         TFile oldfile("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_"+process[i]+"_2017.root");
         TTree *oldtree_egTree = (TTree*)oldfile.Get("egTree");
         TTree *oldtree_mgTree = (TTree*)oldfile.Get("mgTree");
 
         float crosssection = 0;
         oldtree_egTree->SetBranchAddress("crosssection",&crosssection);
         oldtree_egTree->SetBranchStatus("*",1); //NOTE: Only active branches are copied.
 
         float crosssection_ = 0;
         oldtree_mgTree->SetBranchAddress("crosssection",&crosssection_);
         oldtree_mgTree->SetBranchStatus("*",1); //NOTE: Only active branches are copied.
         
	 TFile newfile("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_" + process[i] + "_2017_modified.root","recreate");
         TTree *newtree_egTree = oldtree_egTree->CloneTree(0);
         TTree *newtree_mgTree = oldtree_mgTree->CloneTree(0);
 
         for (Int_t j=0;j<oldtree_egTree->GetEntries(); j++)
         {
             oldtree_egTree->GetEntry(j);
             crosssection = ModifiedXS[i];
             newtree_egTree->Fill();
         }
 
         for (Int_t j=0;j<oldtree_mgTree->GetEntries(); j++)
         {
             oldtree_mgTree->GetEntry(j);
             crosssection_ = ModifiedXS[i];
             newtree_mgTree->Fill();
         }
	 // newtree_egTree->Write();
	 // newtree_mgTree->Write();
         newtree_egTree->AutoSave();
         newtree_mgTree->AutoSave();
 
         //delete oldfile;
         //delete newfile;
       }
}
