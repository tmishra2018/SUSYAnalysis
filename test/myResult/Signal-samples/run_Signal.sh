g++ `root-config --cflags` ../../../lib/libAnaClasses.so analysis_SUSY.C -o analysis_SUSY.exe `root-config --libs`
./analysis_SUSY.exe 2016 1 TChiWG
./analysis_SUSY.exe 2016 0 TChiWG
./analysis_SUSY.exe 2017 1 TChiWG
./analysis_SUSY.exe 2018 1 TChiWG

./analysis_SUSY.exe 2016 1 T5Wg
./analysis_SUSY.exe 2016 0 T5Wg
./analysis_SUSY.exe 2017 1 T5Wg
./analysis_SUSY.exe 2018 1 T5Wg
#g++ `root-config --cflags` ../../../lib/libAnaClasses.so analysis_TChiWG.C -o analysis_TChiWG.exe `root-config --libs`
#./analysis_TChiWG.exe 2016 1
#./analysis_TChiWG.exe 2016 0
#./analysis_TChiWG.exe 2017 1
#./analysis_TChiWG.exe 2018 1
