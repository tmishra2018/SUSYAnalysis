g++ `root-config --cflags` ../../../lib/libAnaClasses.so analysis_mgGJet.C -o analysis_mgGJet.exe `root-config --libs`
./analysis_mgGJet.exe 2016 1
./analysis_mgGJet.exe 2016 0
./analysis_mgGJet.exe 2017 1
./analysis_mgGJet.exe 2018 1
