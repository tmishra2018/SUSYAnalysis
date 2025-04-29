#!/bin/bash
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

METbin1=200
METbin2=400

for HTbin1 in {100,200,300,400}
do
	for HTbin2 in {400,500,600}
	do
	for PHOETbin in {100,120,150,200}
		do

		mkdir t5wg_${HTbin1}_${HTbin2}_${PHOETbin}
		rm binConfig.txt
		echo 'METbin1 200' >>  binConfig.txt
		echo 'METbin2 400' >>  binConfig.txt
		echo 'HTbin1' $HTbin1  >>  binConfig.txt
		echo 'HTbin2' $HTbin2  >>  binConfig.txt
		echo 'PHOETbin' $PHOETbin >>  binConfig.txt

		ch=1
		anatype=3
		lmt=100
		hmt=-1
		lmet=120
		hmet=-1
		iso=4
		lpt=0
		hpt=-1
		rm SigConfig.txt
		echo 'ichannel' $ch  >> SigConfig.txt
		echo 'anatype'  $anatype >>  SigConfig.txt
		echo 'lowMt'    $lmt >> SigConfig.txt
		echo 'highMt'   $hmt >> SigConfig.txt
		echo 'lowMET'   $lmet >>SigConfig.txt
		echo 'highMET'  $hmet >>SigConfig.txt
		echo 'lowPt'    $lpt >> SigConfig.txt
		echo 'highPt'   $hpt >> SigConfig.txt
		echo 'lepIso'   $iso    >> SigConfig.txt
		root -q pred_VGBkg.C++
		root -q pred_eleBkg.C++
		root -q pred_jetBkg.C++
		root -q pred_qcdBkg.C++
		root -q pred_rareBkg.C++

		ch=2
		anatype=3
		lmt=100
		hmt=-1
		lmet=120
		hmet=-1
		iso=4
		lpt=0
		hpt=-1
		rm SigConfig.txt
		echo 'ichannel' $ch  >> SigConfig.txt
		echo 'anatype'  $anatype >>  SigConfig.txt
		echo 'lowMt'    $lmt >> SigConfig.txt
		echo 'highMt'   $hmt >> SigConfig.txt
		echo 'lowMET'   $lmet >>SigConfig.txt
		echo 'highMET'  $hmet >>SigConfig.txt
		echo 'lowPt'    $lpt >> SigConfig.txt
		echo 'highPt'   $hpt >> SigConfig.txt
		echo 'lepIso'   $iso    >> SigConfig.txt
		root -q pred_VGBkg.C++
		root -q pred_eleBkg.C++
		root -q pred_jetBkg.C++
		root -q pred_qcdBkg.C++
		root -q pred_rareBkg.C++

		root -b -q plot_eventct.C++
		root -q analysis_TChiWG.C++
		mv SignalSystematic_egmg.root SignalSystematic_${METbin1}_${METbin2}.root
		mv signalTree_T5WG.root signalTree_T5WG_${METbin1}_${METbin2}.root
	  python writeNewBin.py --MET1 $METbin1 --MET2 $METbin2
		mv t5wg/* t5wg_${HTbin1}_${HTbin2}_${PHOETbin}/.
	done
done 
done
#

