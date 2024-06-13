void getMedianExp_fromTXT(){
  ifstream infile("MedianExp_T5bbbbZg_SBinV7.txt");
  string line;
  vector<float> mGl,mNLSP,med;
  while (getline(infile, line)){
    std::istringstream iss(line);
    float t1;
    iss>>t1;
    mGl.push_back(t1);
    iss>>t1;
    mNLSP.push_back(t1);
    iss>>t1;
    med.push_back(t1);
  }
  infile.close();
  cout<<mGl.size()<<" "<<mNLSP.size()<<" "<<med.size()<<endl;
  cout<<mGl[0]<<" "<<mNLSP[0]<<" "<<med[0]<<endl;

  TH2D *h2_mGlmNLSP_medProfL = new TH2D("mGlmNLSP_medProfL","x: mass of gluino, y: mass of neutralino, z: median exp ProfileLikelyToy400",100,12.5,2512.5,250,5,2505);
  TFile *fout = new TFile("T5bbbbZg_MedExp_ProfileLikelyToy400_SbinV7.root","recreate");
  for(int i=0;i<mGl.size();i++){
    h2_mGlmNLSP_medProfL->Fill(mGl[i],mNLSP[i],med[i]);
  }
  h2_mGlmNLSP_medProfL->Write();
}
