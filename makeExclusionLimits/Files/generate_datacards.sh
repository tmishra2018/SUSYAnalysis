
# i=TChiWg_0_800_FastSim_v17
# for i in T5bbbbZg_1800_150_FastSim_v17 T5bbbbZg_1800_1750_FastSim_v17 TChiWg_0_800_FastSim_v17
for k in 400 800 1200
do
    for i in TChiWg_0_${k}_FastSim_v17 #TChiWg_0_800_FastSim_v17_njet2 TChiWg_0_1200_FastSim_v17_njet	     2
# for k in 150 1000 1750
# do
#     for i in T5bbbbZg_1800_${k}_FastSim_v17 #TChiWg_0_800_FastSim_v17_njet2 TChiWg_0_1200_FastSim_v17
	     # _njet2
    do
	for j in MET_nj1 MET_nj2to4 MET_nj5to7 MET_nj8toinf #METvBin2 # METvBin1 METvBin2
    # for j in  BestPhotonPt #PhovBin PhovBin1 PhoBin2
	do
	    echo "	rm -rf dataCards/${i}_${j}"
	    rm -rf dataCards/${i}_${j}
	    echo "	mkdir dataCards/${i}_${j}"
	    mkdir dataCards/${i}_${j}
	    echo "root -l -q makeDatacard_Simple_bybhumi.C+'(0,${k},"${i}.root","${j}")'"
	    root -l -q makeDatacard_Simple_bybhumi.C+'(0,'${k}',"'${i}'.root","'${j}'")'
	    # .q
	done
	echo "================================"
    done
done
    # .L makeDatacard_Simple_bybhumi.C
    # 'makeDatacard_Simple_bybhumi(0,800,"TChiWg_0_800_FastSim_v17.root","${j}")'
    

# done
# echo "================================="
# for i in T5bbbbZg_1800_150_FastSim_v17 T5bbbbZg_1800_1750_FastSim_v17 TChiWg_0_800_FastSim_v17
# do
#     for j in METvBin METvBin1 METvBin2
# 	     do
# 		 echo "${i}_${j}"
# 		 mkdir ${i}_${j}
#     done
# done
# for 
