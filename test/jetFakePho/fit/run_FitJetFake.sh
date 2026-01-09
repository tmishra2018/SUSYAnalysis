#!/bin/bash
ptbins=(30 32 34 36 38 40 45 50 55 60 65 70 80 90 120 150 180 220 300)
for RunYear in 2016 2017 2018; do
    if [ "$RunYear" -eq 2016 ]; then
        ISpreVFP_list=(0 1)
    else
        ISpreVFP_list=(0)
    fi
    for ISpreVFP in "${ISpreVFP_list[@]}"; do
        for eventType in 1 2; do
            for (( i = 0; i < ${#ptbins[@]} - 1; i++ )); do
                j=$((i + 1))
                ptmin=${ptbins[$i]}
                ptmax=${ptbins[$j]}
                echo "Running for Year=$RunYear, ISpreVFP=$ISpreVFP, eventType=$eventType, pt=$ptmin-$ptmax"
                root -b -q "FitJetFake.C+($RunYear, $ISpreVFP, $eventType, $ptmin, $ptmax, 1, 0, 1.4442)"
            done
        done
    done
done
