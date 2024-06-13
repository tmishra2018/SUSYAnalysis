#!/bin/sh
echo "Exp limits" > a.txt
for i in TChiWH_800_100_MC2018 TChiWH_800_200_MC2018 TChiWH_800_400_MC2018 TChiWH_800_600_MC2018 TChiWH_600_100_MC2018 TChiWH_600_200_MC2018 TChiWH_600_400_MC2018 TChiWH_300_1_MC2018
do
    for h in METvBin_1Wt1Wm1Ht1Hm #METvBin_0-1Wm1Hm
    do
	root -l 'makeDatacard_Simple.C("'$i'","")' -q
	combineCards.py dataCards/$i*$h*.txt > dataCards/combined_$i$h.txt
	rm dataCards/$i*.txt
        echo $i |tee -a a.txt
	combine -M AsymptoticLimits dataCards/combined_${i}${h}.txt -n ${i}_${h}.txt -t -1 |tee -a a.txt
    done 
    echo "=================================================="    
    # echo $i |tee -a a.txt
    # combineCards.py dataCards/combined_${i}*METvBin_0-1Wm1Hm*.txt dataCards/combined_${i}*METvBin_1Wt1Wm1Ht1Hm*.txt > dataCards/combined_${i}_FB_WH_FR_WH.txt
    # combine -M AsymptoticLimits dataCards/combined_${i}_FB_WH_FR_WH.txt -t -1 |tee -a a.txt
done
grep 'TChi\|50.0' a.txt > opLimitCalc.txt
rm a.txt