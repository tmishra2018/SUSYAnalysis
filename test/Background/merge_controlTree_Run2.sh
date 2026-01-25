#!/bin/bash

OUTDIR=/eos/uscms/store/user/tmishra/Background

years=("2016preVFP" "2016postVFP" "2017" "2018")
channels=("egamma" "mg")
samples=("VGBkg" "eleBkg" "jetbkg" "qcd" "rareBkg" "signal")

for ch in "${channels[@]}"; do
  for samp in "${samples[@]}"; do

    inputs=""
    for y in "${years[@]}"; do
      inputs+=" ${OUTDIR}/controlTree_${ch}_${samp}_met40_70_pt0_1000_${y}.root"
    done

    outfile=${OUTDIR}/controlTree_${ch}_${samp}_met40_70_pt0_1000_678.root

    echo "Merging → ${outfile}"
    echo "  Inputs:${inputs}"

    hadd -f ${outfile} ${inputs}

  done
done

echo "Run-2 controlTree merging completed."
