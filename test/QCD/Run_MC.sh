#!/bin/bash
# analysis_egMC(bool isQCD, bool isGJet, int RunYear, bool preVFP)

g++ `root-config --cflags` ../../lib/libAnaClasses.so analysis_mgMC.C -o analysis_mgMC.exe `root-config --libs`
g++ `root-config --cflags` ../../lib/libAnaClasses.so analysis_egMC.C -o analysis_egMC.exe `root-config --libs`
./analysis_egMC.exe 1 0 2016 1 &
./analysis_egMC.exe 1 0 2016 0 &
./analysis_egMC.exe 1 0 2017 0 &
./analysis_egMC.exe 1 0 2018 0 &

./analysis_mgMC.exe 1 0 2016 1 &
./analysis_mgMC.exe 1 0 2016 0 &
./analysis_mgMC.exe 1 0 2017 0 &
./analysis_mgMC.exe 1 0 2018 0 &
