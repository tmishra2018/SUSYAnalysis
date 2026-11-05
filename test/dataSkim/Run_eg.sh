#!/bin/bash
echo 'start'
g++ `root-config --cflags` ../../lib/libAnaClasses.so analysis_eg.C -o analysis_eg.exe `root-config --libs`
./analysis_eg.exe 2018 A
./analysis_eg.exe 2018 B
./analysis_eg.exe 2018 C
./analysis_eg.exe 2018 D

./analysis_eg.exe 2017 B
./analysis_eg.exe 2017 C
./analysis_eg.exe 2017 D
./analysis_eg.exe 2017 E
./analysis_eg.exe 2017 F

./analysis_eg.exe 2016 B-ver1
./analysis_eg.exe 2016 B-ver2
./analysis_eg.exe 2016 C
./analysis_eg.exe 2016 D
./analysis_eg.exe 2016 E
./analysis_eg.exe 2016 F_preVFP
./analysis_eg.exe 2016 F_postVFP
./analysis_eg.exe 2016 G
./analysis_eg.exe 2016 H
echo 'finish'
