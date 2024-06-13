g++ `root-config --cflags` ../../lib/libAnaClasses.so analysis_ISRMC.C -o analysis_ISRMC_postVFP.exe `root-config --libs`
./analysis_ISRMC_postVFP.exe 2016 ZGToLLG
./analysis_ISRMC_postVFP.exe 2016 TTGJets
./analysis_ISRMC_postVFP.exe 2016 TTJets
./analysis_ISRMC_postVFP.exe 2016 WWG
./analysis_ISRMC_postVFP.exe 2016 WZG
