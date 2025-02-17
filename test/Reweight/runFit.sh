#!/bin/bash

ptbins=( 35 40 50 60 70 80 100 150 200 250 300 800 )

declare -A year_preVFP
year_preVFP[2016,0]=false
year_preVFP[2016,1]=true
year_preVFP[2017,0]=false
year_preVFP[2018,0]=false

for detType in 1 2; do
    if [[ $detType -eq 1 ]]; then
        whichRegion="EB"
    elif [[ $detType -eq 2 ]]; then
        whichRegion="EE"
    fi
        # Loop for 2016
    for preVFP in 0 1; do
        year=2016
        if [[ $preVFP -eq 1 ]]; then
            whichVFP="preVFP"
        else
            whichVFP="postVFP"
        fi
        echo "Removing file: /eos/uscms/store/user/tmishra/jetfakepho/files/JetFakeRate-ISR-${whichRegion}_${year}${whichVFP}.txt"
        rm -f "/eos/uscms/store/user/tmishra/jetfakepho/files/JetFakeRate-ISR-${whichRegion}_${year}${whichVFP}.txt"
        bool_preVFP=${year_preVFP[$year,$preVFP]}

        for (( i = 0; i < ${#ptbins[@]} - 1; i++ )); do
            j=$((i + 1))  # Get the next ptbin
            root -b -q "FitJetFake.C+(${ptbins[$i]},${ptbins[$j]},${detType},${year},${bool_preVFP})"
        done
    done

    # Loop for 2017 and 2018
    for year in 2017 2018; do
        whichVFP=""
        echo "Removing file: /eos/uscms/store/user/tmishra/jetfakepho/files/JetFakeRate-ISR-${whichRegion}_${year}${whichVFP}.txt"
        rm -f "/eos/uscms/store/user/tmishra/jetfakepho/files/JetFakeRate-ISR-${whichRegion}_${year}${whichVFP}.txt"

        for (( i = 0; i < ${#ptbins[@]} - 1; i++ )); do
            j=$((i + 1))  # Get the next ptbin
            root -b -q "FitJetFake.C+(${ptbins[$i]},${ptbins[$j]},${detType},${year},0)"
        done
    done
done
