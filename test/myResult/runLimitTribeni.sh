#!/bin/bash
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
RunYear=2016
preVFP=1

NBIN=18
METbin1=200
METbin2=400

HTbin1=100
HTbin2=400

PHOETbin=200

rm ../Result/binConfig.txt		  
echo 'NBIN'  $NBIN >> ../Result/binConfig.txt
echo 'METbin1' $METbin1 >>  ../Result/binConfig.txt
echo 'METbin2' $METbin2 >>  ../Result/binConfig.txt
echo 'HTbin1' $HTbin1  >>  ../Result/binConfig.txt
echo 'HTbin2' $HTbin2  >>  ../Result/binConfig.txt
echo 'PHOETbin' $PHOETbin >>  ../Result/binConfig.txt

ch=1
anatype=3
lmt=100
hmt=-1
lmet=120
hmet=-1
iso=4
lpt=0
hpt=-1
rm ../Result/SigConfig.txt
echo 'ichannel' $ch  >> ../Result/SigConfig.txt
echo 'anatype'  $anatype >>  ../Result/SigConfig.txt
echo 'lowMt'    $lmt >> ../Result/SigConfig.txt
echo 'highMt'   $hmt >> ../Result/SigConfig.txt
echo 'lowMET'   $lmet >> ../Result/SigConfig.txt
echo 'highMET'  $hmet >> ../Result/SigConfig.txt
echo 'lowPt'    $lpt >> ../Result/SigConfig.txt
echo 'highPt'   $hpt >> ../Result/SigConfig.txt
echo 'lepIso'   $iso    >> ../Result/SigConfig.txt
echo 'RunYear'   $RunYear    >> ../Result/SigConfig.txt
echo 'preVFP'   $preVFP  >>  ../Result/SigConfig.txt

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
rm ../Result/SigConfig.txt
echo 'ichannel' $ch  >> ../Result/SigConfig.txt
echo 'anatype'  $anatype >>  ../Result/SigConfig.txt
echo 'lowMt'    $lmt >> ../Result/SigConfig.txt
echo 'highMt'   $hmt >> ../Result/SigConfig.txt
echo 'lowMET'   $lmet >> ../Result/SigConfig.txt
echo 'highMET'  $hmet >> ../Result/SigConfig.txt
echo 'lowPt'    $lpt >> ../Result/SigConfig.txt
echo 'highPt'   $hpt >> ../Result/SigConfig.txt
echo 'lepIso'   $iso    >> ../Result/SigConfig.txt
echo 'RunYear'   $RunYear    >> ../Result/SigConfig.txt
echo 'preVFP'   $preVFP  >>  ../Result/SigConfig.txt

root -l -q pred_VGBkg.C++  > VG_mg.log
root -l -q pred_eleBkg.C++ > ele_mg.log
root -l -q pred_jetBkg.C++ > jet_mg.log
root -l -q pred_qcdBkg.C++ > qcd_mg.log
root -l -q pred_rareBkg.C++ > rare_mg.log
root -l -q pred_sig.C++ > data_mg.log


if [ ${RunYear} == 2016 -a ${preVFP} == 1 ]
then	
	rm eventcount_${RunYear}preVFP.txt
	root -l -q "plot_eventct.C+($NBIN)" >> eventcount_${RunYear}preVFP.txt
elif [ ${RunYear} == 2016 -a ${preVFP} == 0 ]	
then
	rm eventcount_${RunYear}postVFP.txt
	root -l -q "plot_eventct.C+($NBIN)" >> eventcount_${RunYear}postVFP.txt
else 
	rm eventcount_${RunYear}.txt
	root -l -q "plot_eventct.C+($NBIN)" >> eventcount_${RunYear}.txt
fi

#python createDataCard.py ${NBIN} ${RunYear}

#root -l -q analysis_T6WG.C++
#mkdir /tmp/tribeni
#python writeT5WGcard.py ${NBIN} ${RunYear}
#mv /tmp/tribeni/* /eos/uscms/store/user/tmishra/CombinedLimit/Sensitive_myFile/cardsOption4_MT150

#python writeT5WGcard.py ${NBIN}  >> sig_MT${lmt}_MET${lmet}_${NBIN}bins_${PHOETbin}photonPt.txt

# the limit plot other codes at /uscms/homes/t/tmishra/work/CMSSW_10_2_22/src/HiggsAnalysis/CombinedLimit/tools

		#if [ "$HTbin1" -lt "$HTbin2" ]; then
		#  python writeT5WGcard.py
		#else
	    	#  python writeNewT5WG.py
		#fi
		
		#mv t5wg/* t5wg_${HTbin1}_${HTbin2}_${PHOETbin}/.
		#mv *.log t5wg_${HTbin1}_${HTbin2}_${PHOETbin}/.
