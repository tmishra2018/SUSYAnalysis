#!/bin/bash
export ARG1=$1
export ARG2=$2
export ARG3=$3
export ARG4=$4
export ARG5=$5
export ARG6=$6
export ARG7=$7

export year=2016
export preVFP=0
export isData=0

cd ${_CONDOR_SCRATCH_DIR}
echo "source /cvmfs/cms.cern.ch/cmsset_default.sh"
source /cvmfs/cms.cern.ch/cmsset_default.sh

echo "scramv1 project CMSSW CMSSW_14_0_7"
scramv1 project CMSSW CMSSW_14_0_7

echo "cd CMSSW_14_0_7/src/"
cd CMSSW_14_0_7/src/

echo "eval `scramv1 runtime -sh`"
eval `scramv1 runtime -sh`
tar -zxvf ../../files.tar.gz
cd SUSYAnalysis/include && make
cd ../test/eleFakePho/fit
bash make.sh
voms-proxy-init --voms cms --valid 168:00 -out ~/.globus/gridproxy.cert

./FitKer.exe ${ARG1} ${ARG2} ${ARG3} ${ARG4} ${ARG5} ${ARG6} ${ARG7}

if [ ${year} == 2016 -a ${isData} == 1 -a ${preVFP} == 1 ]
then
        xrdcp -f *.png         root://cmseos.fnal.gov//store/user/tmishra/elefakepho/DataFitting2016preVFP
        xrdcp -f Ele*.txt      root://cmseos.fnal.gov//store/user/tmishra/elefakepho/DataResult16preVFP
elif [ ${year} == 2016 -a ${isData} == 1 -a ${preVFP} == 0 ]
then
        xrdcp -f *.png         root://cmseos.fnal.gov//store/user/tmishra/elefakepho/DataFitting2016postVFP
        xrdcp -f Ele*.txt      root://cmseos.fnal.gov//store/user/tmishra/elefakepho/DataResult16postVFP
elif [ ${year} == 2017 -a ${isData} == 1 ]
then
	xrdcp -f *.png        root://cmseos.fnal.gov//store/user/tmishra/elefakepho/DataFitting2017
	xrdcp -f Ele*.txt     root://cmseos.fnal.gov//store/user/tmishra/elefakepho/DataResult17
elif [ ${year} == 2018 -a ${isData} == 1 ]
then
	xrdcp -f *.png        root://cmseos.fnal.gov//store/user/tmishra/elefakepho/DataFitting2018
	xrdcp -f Ele*.txt     root://cmseos.fnal.gov//store/user/tmishra/elefakepho/DataResult18
fi

if [ ${year} == 2016 -a ${isData} == 0 -a ${preVFP} == 1 ]
then
	xrdcp -f *.png        root://cmseos.fnal.gov//store/user/tmishra/elefakepho/DYFitting16preVFP
	xrdcp -f Ele*.txt     root://cmseos.fnal.gov//store/user/tmishra/elefakepho/DYResult16preVFP
elif [ ${year} == 2016 -a ${isData} == 0 -a ${preVFP} == 0 ]
then
	xrdcp -f *.png        root://cmseos.fnal.gov//store/user/tmishra/elefakepho/DYFitting16postVFP
	xrdcp -f Ele*.txt     root://cmseos.fnal.gov//store/user/tmishra/elefakepho/DYResult16postVFP
elif [ ${year} == 2017 -a ${isData} == 0 ]
then
	xrdcp -f *.png        root://cmseos.fnal.gov//store/user/tmishra/elefakepho/DYFitting17
	xrdcp -f Ele*.txt     root://cmseos.fnal.gov//store/user/tmishra/elefakepho/DYResult17
elif [ ${year} == 2018 -a ${isData} == 0 ]
then
	xrdcp -f *.png        root://cmseos.fnal.gov//store/user/tmishra/elefakepho/DYFitting18
	xrdcp -f Ele*.txt     root://cmseos.fnal.gov//store/user/tmishra/elefakepho/DYResult18
fi

rm *.txt  *.png
cd ${_CONDOR_SCRATCH_DIR} && rm -rf CMSSW_14_0_7
