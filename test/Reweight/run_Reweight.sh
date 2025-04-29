#!/bin/bash

g++ $(root-config --cflags) ../../lib/libAnaClasses.so analysis_ISR.C -o analysis_ISR.exe $(root-config --libs)
g++ $(root-config --cflags) ../../lib/libAnaClasses.so analysis_ISRMC.C -o analysis_ISRMC.exe $(root-config --libs)
g++ $(root-config --cflags) ../../lib/libAnaClasses.so analysis_mgHadron.C -o analysis_mgHadron.exe $(root-config --libs)

years=("2017" "2018" "2016")
isr_runs_2017=("D" "E" "F")
isr_runs_2018=("A" "B" "C" "D")
isr_runs_2016=("B-ver1" "B-ver2" "C" "D" "E" "F_preVFP" "F_postVFP" "G" "H")
mc_samples=("TTGJets" "WWG" "WZG" "TTJets" "ZGToLLG")

for year in 2017 2018; do
    if [ $year -eq 2017 ]; then
        runs=("${isr_runs_2017[@]}")
    else
        runs=("${isr_runs_2018[@]}")
    fi
    for run in "${runs[@]}"; do
        ./analysis_ISR.exe $year $run
    done
done

for year in 2017 2018; do
    for sample in "${mc_samples[@]}"; do
        ./analysis_ISRMC.exe $year 0 $sample
    done
done

for year in 2017 2018; do
    if [ $year -eq 2017 ]; then
        runs=("${isr_runs_2017[@]}")
    else
        runs=("${isr_runs_2018[@]}")
    fi
    for run in "${runs[@]}"; do
        ./analysis_mgHadron.exe $year $run
    done
done

for sample in "${mc_samples[@]}"; do
    for version in 0 1; do
        ./analysis_ISRMC.exe 2016 $version $sample
    done
done

for run in "${isr_runs_2016[@]}"; do
    ./analysis_ISR.exe 2016 $run
    ./analysis_mgHadron.exe 2016 $run
done

rm analysis_ISRMC.exe analysis_ISR.exe analysis_mgHadron.exe
