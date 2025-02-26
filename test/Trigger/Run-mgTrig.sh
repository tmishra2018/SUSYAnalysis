g++ `root-config --cflags` analysis_mgTriggerMC.C -o analysis_mgTriggerMC.exe `root-config --libs`
./analysis_mgTriggerMC.exe 2016 1
./analysis_mgTriggerMC.exe 2016 0
./analysis_mgTriggerMC.exe 2017 1
./analysis_mgTriggerMC.exe 2018 1


# for 2016 Data
# analysis_mgTrigger (int RunYear, const char *Era)

g++ `root-config --cflags` analysis_mgTrigger.C -o analysis_mgTrigger.exe `root-config --libs`

./analysis_mgTrigger.exe  2016 B-ver1
./analysis_mgTrigger.exe  2016 B-ver2
./analysis_mgTrigger.exe  2016 C
./analysis_mgTrigger.exe  2016 D
./analysis_mgTrigger.exe  2016 E
./analysis_mgTrigger.exe  2016 F_preVFP
./analysis_mgTrigger.exe  2016 F_postVFP
./analysis_mgTrigger.exe  2016 G
./analysis_mgTrigger.exe  2016 H

# for 2017 Data
./analysis_mgTrigger.exe  2017 D
./analysis_mgTrigger.exe  2017 E
./analysis_mgTrigger.exe  2017 F

# for 2018 Data
./analysis_mgTrigger.exe  2018 A
./analysis_mgTrigger.exe  2018 B
./analysis_mgTrigger.exe  2018 C
./analysis_mgTrigger.exe  2018 D
