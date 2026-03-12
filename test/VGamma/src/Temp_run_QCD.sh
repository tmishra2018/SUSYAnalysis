#!/bin/bash

g++ `root-config --cflags` ../../../lib/libAnaClasses.so analysis_VGamma_Run2.C \
    -o analysis_VGamma_Run2.exe `root-config --libs`

samples=(
    WJetsToLNu
    GJets_DoubleEM
#    QCD_HT300to500
#    QCD_HT500to700
#    QCD_HT700to1000
#    QCD_HT1000to1500
#    QCD_HT1500to2000
#    QCD_HT2000toInf
#    QCD_HT200to300
#    QCD_HT50to100
#    QCD_HT100to200
)

for sample in "${samples[@]}"; do
#    ./analysis_VGamma_Run2.exe 2018 0 ${sample} &
     ./analysis_VGamma_Run2.exe 2016 1 ${sample} &
#     ./analysis_VGamma_Run2.exe 2016 0 ${sample} &
#    ./analysis_VGamma_Run2.exe 2017 0 ${sample} &
done
wait
echo "All QCD jobs finished."
