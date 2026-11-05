#!/bin/bash

export jobNumber=$1
export name_list="cardsSplitted/x${1}"

cd ${_CONDOR_SCRATCH_DIR}
ls -ltrh

echo "source /cvmfs/cms.cern.ch/cmsset_default.sh"
source /cvmfs/cms.cern.ch/cmsset_default.sh

echo "scramv1 project CMSSW CMSSW_14_1_0_pre4"
scramv1 project CMSSW CMSSW_14_1_0_pre4

echo "cd CMSSW_14_1_0_pre4/src/"
cd CMSSW_14_1_0_pre4/src/

echo "eval `scramv1 runtime -sh`"
eval `scramv1 runtime -sh`
tar -zxvf ../../higgsCombine.tar.gz
scram b clean
scram b -j4

mv ../../cardsSplitted HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit

echo ' ' $name_list 

while IFS= read -r line
do 
	rm filelist
	echo "$line" > filelist
	gluino=`awk -F '[^0-9]*' '{print $3}' filelist`
	chargino=`awk -F '[^0-9]*' '{print $4}' filelist`
	echo "${gluino}"
	echo "${chargino}"
	xrdcp -f root://cmseos.fnal.gov//store/user/tmishra/CombinedLimit/Sensitive_myFile/cardsANbins/counting_t5Wg_${gluino}_${chargino}.txt .
	combine -M Significance counting_t5Wg_${gluino}_${chargino}.txt -t -1 --expectSignal=1 -n t5wg | tee signif_${gluino}_${chargino}.txt
	xrdcp -f signif_${gluino}_${chargino}.txt root://cmseos.fnal.gov//store/user/tmishra/CombinedLimit/Sensitive_myFile/signifANbin
	rm counting_t5Wg_${gluino}_${chargino}.txt
done < "$name_list"

cd ${_CONDOR_SCRATCH_DIR} 
rm -rf CMSSW_14_1_0_pre4
echo "End.. Bye now.."
