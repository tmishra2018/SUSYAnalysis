#!/bin/bash
g++ `root-config --cflags` Data_analysis_elefakepho.C -o Data_analysis_elefakepho.exe `root-config --libs`

./Data_analysis_elefakepho.exe 2018 A
./Data_analysis_elefakepho.exe 2018 B
./Data_analysis_elefakepho.exe 2018 C
./Data_analysis_elefakepho.exe 2018 D

#./Data_analysis_elefakepho.exe 2017 B
#./Data_analysis_elefakepho.exe 2017 C
#./Data_analysis_elefakepho.exe 2017 D
#./Data_analysis_elefakepho.exe 2017 E
#./Data_analysis_elefakepho.exe 2017 F
