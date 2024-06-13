#!/bin/bash
g++ `root-config --cflags` Data_analysis_elefakepho.C -o Data_analysis_elefakepho.exe `root-config --libs`
./Data_analysis_elefakepho.exe 2016 B
./Data_analysis_elefakepho.exe 2016 C
./Data_analysis_elefakepho.exe 2016 D
./Data_analysis_elefakepho.exe 2016 E
./Data_analysis_elefakepho.exe 2016 F
./Data_analysis_elefakepho.exe 2016 G
./Data_analysis_elefakepho.exe 2016 H
