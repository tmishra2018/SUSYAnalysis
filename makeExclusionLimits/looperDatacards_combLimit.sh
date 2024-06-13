#!/bin/sh
selType="_NomMass_"
rm dataCards/*.txt
for i in TChiWH_300_1_MC2018 TChiWH_600_100_MC2018 TChiWH_600_200_MC2018 TChiWH_600_400_MC2018 TChiWH_800_100_MC2018 TChiWH_800_200_MC2018 TChiWH_800_400_MC2018 TChiWH_800_600_MC2018 TChiWH_1000_100_MC2018 TChiWW_800_100_MC2018 TChiWZ_1000_100_MC2018 TChiWZ_1000_1_MC2018 TChiWZ_300_1_MC2018 TChiWZ_600_100_MC2018 TChiWZ_600_1_MC2018 TChiWZ_800_100_MC2018 TChiWZ_800_1_MC2018 TChiWZ_800_400_MC2018
do
    for h in METvBin_2T2M METvBin_1T2M METvBin_1T1M METvBin_FBWZ METvBin_FBWH METvBin_FBW METvBin_FBZ METvBin_FBH
    do
    	root -l 'makeDatacard_Simple.C("'$i'","'$h'")' -q
    	combineCards.py dataCards/$i*$h*.txt > dataCards/combined_${i}${selType}$h.txt
    	rm dataCards/$i*$h*.txt
	nlspM=$(echo $i | grep -o -E '[0-9]+' |sed -n '1p')
	lspM=$(echo $i | grep -o -E '[0-9]+' |sed -n '2p')
	mass=$(echo "$nlspM+$lspM/10000" |bc -l)
	echo $mass
    	combine -M AsymptoticLimits dataCards/combined_${i}${selType}${h}.txt -n $selType${i}_${h} -m $mass -t -1 
    done 
    echo "=================================================="
    # b-veto
    echo $i "bVeto"
    combineCards.py dataCards/combined_${i}${selType}*METvBin_2T2M.txt dataCards/combined_${i}${selType}*METvBin_1T2M.txt dataCards/combined_${i}${selType}*METvBin_1T1M.txt > dataCards/combined_${i}${selType}_bVeto.txt
    combine -M AsymptoticLimits dataCards/combined_${i}${selType}_bVeto.txt -n $selType${i}_bVeto -m $mass -t -1
    # >=1 b-tagged
    echo $i "bTag"
#    combineCards.py dataCards/combined_${i}${selType}*METvBin_FBWH.txt dataCards/combined_${i}${selType}*METvBin_FBW.txt dataCards/combined_${i}${selType}*METvBin_FBH.txt dataCards/combined_${i}${selType}*METvBin_FBWZ.txt dataCards/combined_${i}${selType}*METvBin_FBZ.txt > dataCards/combined_${i}${selType}_bTag.txt
    combineCards.py dataCards/combined_${i}${selType}*METvBin_FBWH.txt dataCards/combined_${i}${selType}*METvBin_FBW.txt dataCards/combined_${i}${selType}*METvBin_FBH.txt > dataCards/combined_${i}${selType}_bTag.txt
    combine -M AsymptoticLimits dataCards/combined_${i}${selType}_bTag.txt -n $selType${i}_bTag -m $mass -t -1
    # comb 0b & >=1b
    echo $i "bVeto & bTag"
    combineCards.py dataCards/combined_${i}${selType}_bVeto.txt dataCards/combined_${i}${selType}_bTag.txt > dataCards/combined_${i}${selType}_allB.txt
    combine -M AsymptoticLimits dataCards/combined_${i}${selType}_allB.txt -n $selType${i}_allB -m $mass -t -1
done

for i in TChiWZ TChiWH
do
    ls -1 higgsCombine${selType}${i}_*_MC2018_bVeto*.root > runList_${i}_bVeto${selType}.txt
    ./plotlimit runList_${i}_bVeto${selType}.txt bVeto${selType}${i}.root a

    ls -1 higgsCombine${selType}${i}_*_MC2018_bTag*.root > runList_${i}_bTag${selType}.txt
    ./plotlimit runList_${i}_bTag${selType}.txt bTag${selType}${i}.root a

    ls -1 higgsCombine${selType}${i}_*_MC2018_allB*.root > runList_${i}_allB${selType}.txt
    ./plotlimit runList_${i}_allB${selType}.txt allB${selType}${i}.root a
done

