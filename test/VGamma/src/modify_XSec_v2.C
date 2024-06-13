void modify_XSec_v2(){

       TString process[5] = {"DYJetsToLL", "TTJets",  "WW", "WZ", "WGToLNuG"};
       float  ModifiedXS[5] = {5670, 750,   75.94, 27.564, 489.0 }; // menglei
       //float  ModifiedXS[11] = {6404.0, 750.5, 3.757, 0.2147, 0.04345, 75.95, 27.59, 412.7*1.25, 19.75*1.25, 0.8099*1.25,  51.1}; // XSDB
       //float  ModifiedXS[11] = {6077.22, 831.76, 3.697, 0.2147, 0.04345, 75.94, 27.564, 489.0, 17.018, 0.87971, 55.49};
       for (int i=0; i<=4; i++)
       {
         TFile oldfile("/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_"+process[i]+"_2017.root");
         TTree *oldtree_egTree = (TTree*)oldfile.Get("signalTree");
         TTree *oldtree_mgTree = (TTree*)oldfile.Get("proxyTree");
 
         float crosssection = 0;
         oldtree_egTree->SetBranchAddress("crosssection",&crosssection);
         oldtree_egTree->SetBranchStatus("*",1); //NOTE: Only active branches are copied.
 
         float crosssection_ = 0;
         oldtree_mgTree->SetBranchAddress("crosssection",&crosssection_);
         oldtree_mgTree->SetBranchStatus("*",1); //NOTE: Only active branches are copied.
         
         TFile newfile("/eos/uscms/store/user/tmishra/mgMC/resTree_mgsignal_"+process[i]+"_2017_modified.root","recreate");
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
