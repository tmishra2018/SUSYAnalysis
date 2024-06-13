g++ `root-config --cflags` ../../../lib/libAnaClasses.so analysis_egMC.C -o analysis_egMC.exe `root-config --libs`

./analysis_egMC.exe 2016 DYJetsToLL
./analysis_egMC.exe 2016 TTJets
./analysis_egMC.exe 2016 WW
./analysis_egMC.exe 2016 WZ

./analysis_egMC.exe 2016 WGToLNuG
./analysis_egMC.exe 2016 WJetsToLNu
#./analysis_egMC.exe 2017 WJetsToLNu
#./analysis_egMC.exe 2018 WJetsToLNu
#./analysis_egMC.exe 2016 WGJet40
#./analysis_egMC.exe 2016 WGJet130
#./analysis_egMC.exe 2016 ZGToLLG
#
#./analysis_egMC.exe 2016 WWG
#./analysis_egMC.exe 2016 WZG
#./analysis_egMC.exe 2016 TTGJets

#./analysis_egMC.exe 2016 WGToLNuG
#./analysis_egMC.exe 2016 DYJetsToLL
#./analysis_egMC.exe 2017 TTJets
#./analysis_egMC.exe 2017 WW
#./analysis_egMC.exe 2017 WZ

#./analysis_egMC.exe 2018 DYJetsToLL
#./analysis_egMC.exe 2017 DYJetsToLL
#./analysis_egMC.exe 2018 WW
#./analysis_egMC.exe 2018 WZ
