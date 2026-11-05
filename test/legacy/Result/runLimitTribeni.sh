#!/bin/bash
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

NBIN=18
METbin1=200
METbin2=400
METbin3=500

HTbin1=100
HTbin2=400
HTbin3=500

PHOETbin=200
PHOETBin2=200

rm binConfig.txt		  
echo 'NBIN'  $NBIN >> binConfig.txt
echo 'METbin1' $METbin1 >>  binConfig.txt
echo 'METbin2' $METbin2 >>  binConfig.txt
echo 'METbin3' $METbin3 >>  binConfig.txt
echo 'HTbin1' $HTbin1  >>  binConfig.txt
echo 'HTbin2' $HTbin2  >>  binConfig.txt
echo 'HTbin3' $HTbin3  >>  binConfig.txt
echo 'PHOETbin' $PHOETbin >>  binConfig.txt
echo 'PHOETBin2' $PHOETBin2 >>  binConfig.txt

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
root -l -q pred_VGBkg.C++   > VG_eg.log
root -l -q pred_eleBkg.C++  > ele_eg.log
root -l -q pred_jetBkg.C++  > jet_eg.log
root -l -q pred_qcdBkg.C++  > qcd_eg.log
root -l -q pred_rareBkg.C++ > rare_eg.log
root -l -q pred_sig.C++ > data_eg.log

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
root -l -q pred_VGBkg.C++  > VG_mg.log
root -l -q pred_eleBkg.C++ > ele_mg.log
root -l -q pred_jetBkg.C++ > jet_mg.log
root -l -q pred_qcdBkg.C++ > qcd_mg.log
root -l -q pred_rareBkg.C++ > rare_mg.log
root -l -q pred_sig.C++ > data_mg.log

#root -l -q "plot_eventct.C+($NBIN)" >> stat_error_MT150.txt
root -l -q "plot_eventct.C+($NBIN)" 
#cp SignalSystematic.root SignalSystematic.root
python createDataCard.py ${NBIN}

root -l -q analysis_T6WG.C++
mkdir /tmp/tribeni
python writeT5WGcard.py ${NBIN}  >> sig_MT${lmt}_MET${lmet}_${NBIN}bins_${PHOETbin}photonPt.txt
mv /tmp/tribeni/* /eos/uscms/store/user/tmishra/CombinedLimit/Sensitive_myFile/cardsOption4_MT150
		#if [ "$HTbin1" -lt "$HTbin2" ]; then
		#  python writeT5WGcard.py
		#else
	    	#  python writeNewT5WG.py
		#fi
		
		#mv t5wg/* t5wg_${HTbin1}_${HTbin2}_${PHOETbin}/.
		#mv *.log t5wg_${HTbin1}_${HTbin2}_${PHOETbin}/.

