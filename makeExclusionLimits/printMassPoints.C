void printMassPoints(TString fName){
  TFile *f1=new TFile(fName);
  TH2D *h2=(TH2D*)f1->FindObjectAny("MGlMNLSP");
  for(int nx=1;nx<=h2->GetNbinsX();nx++){
    for(int ny=1;ny<=h2->GetNbinsY();ny++){
      if(h2->GetBinContent(nx,ny) < 0.1) continue;
      cout<<h2->GetXaxis()->GetBinCenter(nx)<<" "<<h2->GetYaxis()->GetBinCenter(ny)<<endl;
    }
  }

}
