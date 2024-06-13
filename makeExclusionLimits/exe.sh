#!/bin/sh
for i in T5bbbbZg T5qqqqHg T5ttttZg T6ttZg
#for i in T5bbbbZg
do
    echo "********************** "$i" ***********************"
#    ls ../Limits_Signif_$i\_SbinV7_v11/higgsCombine$i\_*Asymptotic*.root > runList_$i\_SbinV7_Obs_v11.txt
#    ls ../Limits_Signif_$i\_SbinV7_v11/higgsCombineExpSignif_$i*.root > runList_signif_$i\_SbinV7_Exp_v11.txt
#    ls ../Limits_Signif_$i\_SbinV7_v11/higgsCombineObsSignif_$i*.root > runList_signif_$i\_SbinV7_Obs_v11.txt
    ./plotlimit runList_$i\_SbinV7_Obs_v11.txt $i\_LimitPlots_Obs_V7_v11.root $i
    ./plotlimit runList_signif_$i\_SbinV7_Obs_v11.txt $i\_Signif_Obs_V7_v11.root $i
    ./plotlimit runList_signif_$i\_SbinV7_Exp_v11.txt $i\_Signif_Exp_V7_v11.root $i
done