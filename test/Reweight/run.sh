#!/bin/bash
g++ `root-config --cflags` ../../lib/libAnaClasses.so analysis_ISRMC.C  -o analysis_ISRMC.exe `root-config --libs`
g++ `root-config --cflags` ../../lib/libAnaClasses.so analysis_ISR.C -o analysis_ISR.exe `root-config --libs`
g++ `root-config --cflags` ../../lib/libAnaClasses.so analysis_mgHadron.C -o analysis_mgHadron.exe `root-config --libs` 


#########################################

./analysis_ISR.exe 2017 D
./analysis_ISR.exe 2017 E
./analysis_ISR.exe 2017 F

./analysis_ISRMC.exe 2017 TTGJets
./analysis_ISRMC.exe 2017 WWG
./analysis_ISRMC.exe 2017 WZG
./analysis_ISRMC.exe 2017 TTJets
./analysis_ISRMC.exe 2017 ZGToLLG

./analysis_mgHadron.exe 2017 D
./analysis_mgHadron.exe 2017 E
./analysis_mgHadron.exe 2017 F

#########################################

./analysis_ISR.exe 2018 A
./analysis_ISR.exe 2018 B
./analysis_ISR.exe 2018 C
./analysis_ISR.exe 2018 D 

./analysis_ISRMC.exe 2018 TTGJets
./analysis_ISRMC.exe 2018 WWG
./analysis_ISRMC.exe 2018 WZG
./analysis_ISRMC.exe 2018 TTJets
./analysis_ISRMC.exe 2018 ZGToLLG

./analysis_mgHadron.exe 2018 A
./analysis_mgHadron.exe 2018 B
./analysis_mgHadron.exe 2018 C
./analysis_mgHadron.exe 2018 D

#########################################

#./analysis_ISRMC.exe 2016 TTGJets
#./analysis_ISRMC.exe 2016 WWG
#./analysis_ISRMC.exe 2016 WZG
#./analysis_ISRMC.exe 2016 TTJets
#./analysis_ISRMC.exe 2016 ZGToLLG
#
#./analysis_ISRMC_postVFP.exe 2016 TTGJets
#./analysis_ISRMC_postVFP.exe 2016 WWG
#./analysis_ISRMC_postVFP.exe 2016 WZG
#./analysis_ISRMC_postVFP.exe 2016 TTJets
#./analysis_ISRMC_postVFP.exe 2016 ZGToLLG
#
#
#
#./analysis_ISR.exe 2016 B-ver2
#./analysis_ISR.exe 2016 C
#./analysis_ISR.exe 2016 D
#./analysis_ISR.exe 2016 E
#./analysis_ISR.exe 2016 F_preVFP
#./analysis_ISR.exe 2016 F_postVFP
#./analysis_ISR.exe 2016 G
#./analysis_ISR.exe 2016 H
#
#./analysis_mgHadron.exe 2016 B-ver2
#./analysis_mgHadron.exe 2016 C
#./analysis_mgHadron.exe 2016 D
#./analysis_mgHadron.exe 2016 E
#./analysis_mgHadron.exe 2016 F_preVFP
#./analysis_mgHadron.exe 2016 F_postVFP
#./analysis_mgHadron.exe 2016 G
#./analysis_mgHadron.exe 2016 H
