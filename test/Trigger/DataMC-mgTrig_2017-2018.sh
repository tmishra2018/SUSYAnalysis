# for 2017 and 2018 both Data MC with one script

g++ `root-config --cflags` analysis_mgTrigger2017_2018.C -o analysis_mgTrigger2017_2018.exe `root-config --libs`

./analysis_mgTrigger2017_2018.exe 0 2017 ""
./analysis_mgTrigger2017_2018.exe 0 2018 ""

./analysis_mgTrigger2017_2018.exe 1 2018 A
./analysis_mgTrigger2017_2018.exe 1 2018 B 
./analysis_mgTrigger2017_2018.exe 1 2018 C
./analysis_mgTrigger2017_2018.exe 1 2018 D

./analysis_mgTrigger2017_2018.exe 1 2017 D
./analysis_mgTrigger2017_2018.exe 1 2017 E
./analysis_mgTrigger2017_2018.exe 1 2017 F

