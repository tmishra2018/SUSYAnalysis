#g++ `root-config --cflags` analysis_mgTriggerMC.C -o analysis_mgTriggerMC.exe `root-config --libs`
#./analysis_mgTriggerMC.exe 2016 1
#./analysis_mgTriggerMC.exe 2016 0
#./analysis_mgTriggerMC.exe 2017 1
#./analysis_mgTriggerMC.exe 2018 1

g++ `root-config --cflags` analysis_mgTrigger_MET.C -o analysis_mgTrigger_MET.exe `root-config --libs`
#./analysis_mgTrigger_MET.exe 2016 0
./analysis_mgTrigger_MET.exe 2017 1
#./analysis_mgTrigger_MET.exe 2016 1
#./analysis_mgTrigger_MET.exe 2018 1

