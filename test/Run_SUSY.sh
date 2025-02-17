g++ `root-config --cflags` ../lib/libAnaClasses.so analysis_SUSY.C -o analysis_SUSY.exe `root-config --libs`
./analysis_SUSY.exe 2016 1
./analysis_SUSY.exe 2016 0
./analysis_SUSY.exe 2017 1
./analysis_SUSY.exe 2018 1
