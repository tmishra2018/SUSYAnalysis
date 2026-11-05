#!/bin/bash
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

METbin1=200
METbin2=400
HTbin1=100
HTbin2=400

	for PHOETbin in {100,150,200,250,300}
		do

		mkdir t5wg_${HTbin1}_${HTbin2}_${PHOETbin}
		rm binConfig.txt
		if [ "$HTbin1" -lt "$HTbin2" ]; then
		  echo 'NBIN'  18 >> binConfig.txt
		  echo 'METbin1 200' >>  binConfig.txt
		  echo 'METbin2 400' >>  binConfig.txt
		  echo 'HTbin1' $HTbin1  >>  binConfig.txt
		  echo 'HTbin2' $HTbin2  >>  binConfig.txt
		  echo 'PHOETbin' $PHOETbin >>  binConfig.txt
		else
			HTbin1=$HTbin2
		  echo 'NBIN'  12 >> binConfig.txt
		  echo 'METbin1 200' >>  binConfig.txt
		  echo 'METbin2 400' >>  binConfig.txt
		  echo 'HTbin1' $HTbin1  >>  binConfig.txt
		  echo 'HTbin2' $HTbin2  >>  binConfig.txt
		  echo 'PHOETbin' $PHOETbin >>  binConfig.txt
		fi


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
		root -q pred_VGBkg.C++   > VG_eg.log
		root -q pred_eleBkg.C++  > ele_eg.log
		root -q pred_jetBkg.C++  > jet_eg.log
		root -q pred_qcdBkg.C++  > qcd_eg.log
		root -q pred_rareBkg.C++ > rare_eg.log
		root -q pred_sig.C++

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
		root -q pred_VGBkg.C++  > VG_mg.log
		root -q pred_eleBkg.C++ > ele_mg.log
		root -q pred_jetBkg.C++ > jet_mg.log
		root -q pred_qcdBkg.C++ > qcd_mg.log
		root -q pred_rareBkg.C++ > rare_mg.log
		root -q pred_sig.C++

		if [ "$HTbin1" -lt "$HTbin2" ]; then
			root -b -q "plot_eventct.C+(18)"
		else
			root -b -q "plot_eventct.C+(12)"
		fi
		root -q analysis_TChiWG.C++
		if [ "$HTbin1" -lt "$HTbin2" ]; then
		  python writeT5WGcard.py
		else
	    python writeNewT5WG.py
		fi
		mv t5wg/* t5wg_${HTbin1}_${HTbin2}_${PHOETbin}/.
		mv *.log t5wg_${HTbin1}_${HTbin2}_${PHOETbin}/.
	done

