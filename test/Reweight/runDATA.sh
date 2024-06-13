#!/bin/bash
g++ `root-config --cflags` ../../lib/libAnaClasses.so analysis_ISR.C -o analysis_ISR.exe `root-config --libs`
#
./analysis_ISR.exe 2016 B-ver2
./analysis_ISR.exe 2016 C
./analysis_ISR.exe 2016 D
./analysis_ISR.exe 2016 E
./analysis_ISR.exe 2016 F_preVFP
./analysis_ISR.exe 2016 F_postVFP
./analysis_ISR.exe 2016 G
./analysis_ISR.exe 2016 H
