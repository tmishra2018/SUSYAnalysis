#!/bin/bash
# analysis_egMC(int RunYear, bool preVFP, const char *Sample)

g++ `root-config --cflags` ../../lib/libAnaClasses.so analysis_mgMC.C -o analysis_mgMC.exe `root-config --libs`

./analysis_mgMC.exe 2016 1 DYJetsToLL
./analysis_mgMC.exe 2016 1 TTJets
./analysis_mgMC.exe 2016 1 WW
./analysis_mgMC.exe 2016 1 WZ
./analysis_mgMC.exe 2016 1 WJetsToLNu
./analysis_mgMC.exe 2016 1 WGToLNuG

./analysis_mgMC.exe 2016 0 DYJetsToLL
./analysis_mgMC.exe 2016 0 TTJets
./analysis_mgMC.exe 2016 0 WW
./analysis_mgMC.exe 2016 0 WZ
./analysis_mgMC.exe 2016 0 WJetsToLNu
./analysis_mgMC.exe 2016 0 WGToLNuG

./analysis_mgMC.exe 2017 1 DYJetsToLL
./analysis_mgMC.exe 2017 1 TTJets
./analysis_mgMC.exe 2017 1 WW
./analysis_mgMC.exe 2017 1 WZ
./analysis_mgMC.exe 2017 1 WJetsToLNu
./analysis_mgMC.exe 2017 1 WGToLNuG

./analysis_mgMC.exe 2018 1 DYJetsToLL
./analysis_mgMC.exe 2018 1 TTJets
./analysis_mgMC.exe 2018 1 WW
./analysis_mgMC.exe 2018 1 WZ
./analysis_mgMC.exe 2018 1 WJetsToLNu
./analysis_mgMC.exe 2018 1 WGToLNuG

g++ `root-config --cflags` ../../lib/libAnaClasses.so analysis_egMC.C -o analysis_egMC.exe `root-config --libs`

./analysis_egMC.exe 2016 1 DYJetsToLL
./analysis_egMC.exe 2016 1 TTJets
./analysis_egMC.exe 2016 1 WW
./analysis_egMC.exe 2016 1 WZ
./analysis_egMC.exe 2016 1 WJetsToLNu
./analysis_egMC.exe 2016 1 WGToLNuG

./analysis_egMC.exe 2016 0 DYJetsToLL
./analysis_egMC.exe 2016 0 TTJets
./analysis_egMC.exe 2016 0 WW
./analysis_egMC.exe 2016 0 WZ
./analysis_egMC.exe 2016 0 WJetsToLNu
./analysis_egMC.exe 2016 0 WGToLNuG

./analysis_egMC.exe 2017 1 DYJetsToLL
./analysis_egMC.exe 2017 1 TTJets
./analysis_egMC.exe 2017 1 WW
./analysis_egMC.exe 2017 1 WZ
./analysis_egMC.exe 2017 1 WJetsToLNu
./analysis_egMC.exe 2017 1 WGToLNuG

./analysis_egMC.exe 2018 1 DYJetsToLL
./analysis_egMC.exe 2018 1 TTJets
./analysis_egMC.exe 2018 1 WW
./analysis_egMC.exe 2018 1 WZ
./analysis_egMC.exe 2018 1 WJetsToLNu
./analysis_egMC.exe 2018 1 WGToLNuG
