# for 2016 Data and MC
# analysis_mgTrigger(bool useData, int RunYear, const char *Era)

g++ `root-config --cflags` analysis_mgTrigger.C -o analysis_mgTrigger.exe `root-config --libs`

./analysis_mgTrigger.exe 0 2016 preVFP
./analysis_mgTrigger.exe 0 2016 postVFP

./analysis_mgTrigger.exe 1 2016 B-ver1
./analysis_mgTrigger.exe 1 2016 B-ver2
./analysis_mgTrigger.exe 1 2016 C
./analysis_mgTrigger.exe 1 2016 D
./analysis_mgTrigger.exe 1 2016 E
./analysis_mgTrigger.exe 1 2016 F_preVFP
./analysis_mgTrigger.exe 1 2016 F_postVFP
./analysis_mgTrigger.exe 1 2016 G
./analysis_mgTrigger.exe 1 2016 H
