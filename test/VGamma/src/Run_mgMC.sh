g++ `root-config --cflags` ../../../lib/libAnaClasses.so analysis_mgMC.C -o analysis_mgMC.exe `root-config --libs`

./analysis_mgMC.exe 2016 DYJetsToLL
./analysis_mgMC.exe 2016 TTJets
./analysis_mgMC.exe 2016 WW
./analysis_mgMC.exe 2016 WZ
./analysis_mgMC.exe 2016 WJetsToLNu
./analysis_mgMC.exe 2016 WGToLNuG

#./analysis_mgMC.exe 2017  WJetsToLNu
#./analysis_mgMC.exe 2018  WJetsToLNu
#./analysis_mgMC.exe 2016 TTGJets
#./analysis_mgMC.exe 2016 WGJet130
#./analysis_mgMC.exe 2016 WGJet40
#./analysis_mgMC.exe 2016 WWG
#./analysis_mgMC.exe 2016 WZG
#./analysis_mgMC.exe 2016 ZGToLLG
