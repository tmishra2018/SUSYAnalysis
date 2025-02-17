g++ `root-config --cflags` ../../../lib/libAnaClasses.so analysis_HadronMC.C -o analysis_HadronMC.exe `root-config --libs`
./analysis_HadronMC.exe 2016 1
./analysis_HadronMC.exe 2016 0
./analysis_HadronMC.exe 2017 1
./analysis_HadronMC.exe 2018 1
