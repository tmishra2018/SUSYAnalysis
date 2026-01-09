#!/bin/bash
g++ `root-config --cflags` ../lib/libAnaClasses.so analysis_mg.C -o analysis_mg.exe `root-config --libs`
./analysis_mg.exe 2017 D
./analysis_mg.exe 2017 E
./analysis_mg.exe 2017 F
./analysis_mg.exe 2017 B
./analysis_mg.exe 2017 C

./analysis_mg.exe 2016 B-ver1
./analysis_mg.exe 2016 B-ver2
./analysis_mg.exe 2016 C
./analysis_mg.exe 2016 D
./analysis_mg.exe 2016 E
./analysis_mg.exe 2016 F_preVFP
./analysis_mg.exe 2016 F_postVFP
./analysis_mg.exe 2016 G
./analysis_mg.exe 2016 H

./analysis_mg.exe 2018 A
./analysis_mg.exe 2018 B
./analysis_mg.exe 2018 C
./analysis_mg.exe 2018 D
