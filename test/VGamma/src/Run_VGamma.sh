g++ `root-config --cflags` ../../../lib/libAnaClasses.so analysis_VGamma.C -o analysis_VGamma.exe `root-config --libs`

./analysis_VGamma.exe 2016 TTJets
./analysis_VGamma.exe 2016 TTGJets
./analysis_VGamma.exe 2016 WWG
./analysis_VGamma.exe 2016 WZG
./analysis_VGamma.exe 2016 WW
./analysis_VGamma.exe 2016 WZ

./analysis_VGamma.exe 2016 WGToLNuG
./analysis_VGamma.exe 2016 WGJet40
./analysis_VGamma.exe 2016 WGJet130
./analysis_VGamma.exe 2016 DYJetsToLL
./analysis_VGamma.exe 2016 ZGToLLG

