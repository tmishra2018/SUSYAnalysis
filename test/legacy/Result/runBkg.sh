#!/bin/bash
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#anatype=0
##ch=1
#lmt=0
#hmt=-1
#lmet=40
#hmet=70
#iso=4
#lpt=0
#hpt=1000
#declare -a array=(0 30 50 70 200)
#
#for ch in {1,2}
#do
##	lpt=${array[$i]}
##	hpt=${array[$i+1]}
#
#	rm BkgPredConfig.txt
#	echo 'ichannel' $ch  >> BkgPredConfig.txt
#	echo 'anatype'  $anatype >>  BkgPredConfig.txt
#	echo 'lowMt'    $lmt >> BkgPredConfig.txt
#	echo 'highMt'   $hmt >> BkgPredConfig.txt
#	echo 'lowMET'   $lmet >>BkgPredConfig.txt
#	echo 'highMET'  $hmet >>BkgPredConfig.txt
#	echo 'lowPt'    $lpt >> BkgPredConfig.txt
#	echo 'highPt'   $hpt >> BkgPredConfig.txt
#	echo 'lepIso'   $iso    >> BkgPredConfig.txt

for METbin1 in {200,220,240,260,300}
do
	for METbin2 in {320,340,360,380,400,420,440,460,500}
	do

		rm binConfig.txt
		echo 'METbin1' $METbin1 >>  binConfig.txt
		echo 'METbin2' $METbin2 >>  binConfig.txt

		ch=1
		anatype=3
		lmt=100
		hmt=-1
		lmet=120
		hmet=-1
		iso=4
		lpt=0
		hpt=-1
		rm BkgPredConfig.txt
		echo 'ichannel' $ch  >> BkgPredConfig.txt
		echo 'anatype'  $anatype >>  BkgPredConfig.txt
		echo 'lowMt'    $lmt >> BkgPredConfig.txt
		echo 'highMt'   $hmt >> BkgPredConfig.txt
		echo 'lowMET'   $lmet >>BkgPredConfig.txt
		echo 'highMET'  $hmet >>BkgPredConfig.txt
		echo 'lowPt'    $lpt >> BkgPredConfig.txt
		echo 'highPt'   $hpt >> BkgPredConfig.txt
		echo 'lepIso'   $iso    >> BkgPredConfig.txt
		root -q analysis_VGBkg.C++
		root -q analysis_eleBkg.C++
		root -q analysis_jetBkg.C++
		root -q analysis_qcdBkg.C++
		root -q analysis_rareBkg.C++

		ch=2
		anatype=3
		lmt=100
		hmt=-1
		lmet=120
		hmet=-1
		iso=4
		lpt=0
		hpt=-1
		rm BkgPredConfig.txt
		echo 'ichannel' $ch  >> BkgPredConfig.txt
		echo 'anatype'  $anatype >>  BkgPredConfig.txt
		echo 'lowMt'    $lmt >> BkgPredConfig.txt
		echo 'highMt'   $hmt >> BkgPredConfig.txt
		echo 'lowMET'   $lmet >>BkgPredConfig.txt
		echo 'highMET'  $hmet >>BkgPredConfig.txt
		echo 'lowPt'    $lpt >> BkgPredConfig.txt
		echo 'highPt'   $hpt >> BkgPredConfig.txt
		echo 'lepIso'   $iso    >> BkgPredConfig.txt
		root -q analysis_VGBkg.C++
		root -q analysis_eleBkg.C++
		root -q analysis_jetBkg.C++
		root -q analysis_qcdBkg.C++
		root -q analysis_rareBkg.C++
		root -b -q "plot_eventct.C+(36)"
		root -q analysis_TChiWG.C++
		mv SignalSystematic_egmg.root SignalSystematic_${METbin1}_${METbin2}.root
		mv signalTree_TChiWG.root signalTree_TChiWG_${METbin1}_${METbin2}.root
		mv signalTree_T5WG.root signalTree_T5WG_${METbin1}_${METbin2}.root
	  python writeTChiWGcard.py --MET1 $METbin1 --MET2 $METbin2
	  python writeT5WGcard.py --MET1 $METbin1 --MET2 $METbin2
	done
done 
#done
#


#rm BkgPredConfig.txt
#cp BkgPred_bkgConfig_eg.txt  BkgPredConfig.txt
root -q analysis_VGBkg.C++ 
root -q analysis_eleBkg.C++
root -q analysis_jetBkg.C++
root -q analysis_qcdBkg.C++
root -q analysis_rareBkg.C++
root -q analysis_sig.C++


#rm BkgPredConfig.txt
#cp BkgPred_signalConfig_eg.txt BkgPredConfig.txt  
#root -q analysis_VGBkg.C++
#root -q analysis_eleBkg.C++
#root -q analysis_jetBkg.C++
#root -q analysis_qcdBkg.C++
#root -q analysis_rareBkg.C++
#root -q analysis_sig.C++


#rm BkgPredConfig.txt
#cp BkgPred_validConfig_eg.txt BkgPredConfig.txt
#root -q analysis_VGBkg.C++
#root -q analysis_eleBkg.C++
#root -q analysis_jetBkg.C++
#root -q analysis_qcdBkg.C++
#root -q analysis_rareBkg.C++
#root -q analysis_sig.C++


#rm BkgPredConfig.txt
#cp BkgPred_bkgConfig.txt  BkgPredConfig.txt
#root -q analysis_VGBkg.C++
##root -q analysis_eleBkg.C++
##root -q analysis_jetBkg.C++
#root -q analysis_qcdBkg.C++
##root -q analysis_rareBkg.C++
##root -q analysis_sig.C++
##
##
#rm BkgPredConfig.txt
#cp BkgPred_signalConfig.txt BkgPredConfig.txt   
#root -q analysis_VGBkg.C++
##root -q analysis_eleBkg.C++
##root -q analysis_jetBkg.C++
#root -q analysis_qcdBkg.C++
##root -q analysis_rareBkg.C++
##root -q analysis_sig.C++
##
##
#rm BkgPredConfig.txt
#cp BkgPred_validConfig.txt BkgPredConfig.txt
#root -q analysis_VGBkg.C++
##root -q analysis_eleBkg.C++
##root -q analysis_jetBkg.C++
#root -q analysis_qcdBkg.C++
##root -q analysis_rareBkg.C++
##root -q analysis_sig.C++
##
##
root -q pred_VGBkg.C++
root -q pred_eleBkg.C++
root -q pred_jetBkg.C++
root -q pred_qcdBkg.C++
root -q pred_rareBkg.C++
root -q pred_sig.C++
