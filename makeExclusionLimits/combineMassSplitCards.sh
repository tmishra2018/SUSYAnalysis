#!/bin/sh
selType1="_65To105ZH_"
selType2="_105To135ZH_"
#rm dataCards/*.txt
for i in TChiWH_300_1_MC2018 TChiWH_600_100_MC2018 TChiWH_600_200_MC2018 TChiWH_600_400_MC2018 TChiWH_800_100_MC2018 TChiWH_800_200_MC2018 TChiWH_800_400_MC2018 TChiWH_800_600_MC2018 TChiWW_800_100_MC2018 TChiWZ_1000_100_MC2018 TChiWZ_1000_1_MC2018 TChiWZ_300_1_MC2018 TChiWZ_600_100_MC2018 TChiWZ_600_1_MC2018 TChiWZ_800_100_MC2018 TChiWZ_800_1_MC2018 TChiWZ_800_400_MC2018
do
    # combineCards.py dataCards/combined_${i}${selType1}*allB.txt dataCards/combined_${i}${selType2}*allB.txt > dataCards/combined_${i}${selType1}${selType2}_allB.txt
    # nlspM=$(echo $i | grep -o -E '[0-9]+' |sed -n '1p')
    # lspM=$(echo $i | grep -o -E '[0-9]+' |sed -n '2p')
    # mass=$(echo "$nlspM+$lspM/10000" |bc -l)
    # echo $mass
    # combine -M AsymptoticLimits dataCards/combined_${i}${selType1}${selType2}_allB.txt -n ${selType1}${selType2}${i}_allB -m $mass -t -1 
    
    echo "=================================================="
done
for i in TChiWZ TChiWH
do
    ls -1 higgsCombine${selType1}${selType2}${i}_*allB*.root > runList_${i}_allB${selType1}${selType2}.txt
    ./plotlimit runList_${i}_allB${selType1}${selType2}.txt allB_SplitMass_${i}.root a
done
