void plotNBJet() {
	    gROOT->SetBatch(true);
    TFile *file = TFile::Open("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_WZ_2018.root");

    // Check if the file was opened successfully
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Get the TTree "egTree"
    TTree *tree = (TTree*)file->Get("egTree");

    // Check if the tree was found
    if (!tree) {
        std::cerr << "Error finding TTree 'egTree'!" << std::endl;
        file->Close();
        return;
    }

    // Declare a variable to hold the value of the nBJet branch
    int nBJet = 0;

    // Set up the branch to read the nBJet variable
    tree->SetBranchAddress("nBJet", &nBJet);

    // Create a histogram to store the distribution of nBJet
    TH1F *hist = new TH1F("hist_nBJet", "Distribution of nBJet", 10, 0, 10);  // 10 bins from 0 to 10

    // Loop over the entries in the tree and fill the histogram
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        hist->Fill(nBJet);
    }

    // Create a canvas to draw the histogram
    TCanvas *canvas = new TCanvas("canvas", "nBJet Distribution", 800, 600);

    // Draw the histogram
    hist->Draw();

    // Save the canvas to a file
    canvas->SaveAs("nBJet_distribution.png");

    // Clean up
 //   file->Close();
 //   delete canvas;
 //   delete hist;
}
