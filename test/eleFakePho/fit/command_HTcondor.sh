#!/bin/bash
# Arguments: year preVFP isData doEB, then the FitKer arguments
export year=$1
export preVFP=$2
export isData=$3
export doEB=$4
shift 4

export ARG1=$1
export ARG2=$2
export ARG3=$3
export ARG4=$4
export ARG5=$5
export ARG6=$6
export ARG7=$7
export ARG8=$8
export ARG9=$9

CROSS_TAG=""
if [ -n "${ARG8}" ] && [ -n "${ARG9}" ]
then
        CROSS_LOW=$(printf '%.0f' "${ARG8}")
        CROSS_HIGH=$(printf '%.0f' "${ARG9}")

        if [ "${ARG2}" == "0" ]
        then
                if [ "${CROSS_HIGH}" -ge 10000 ]
                then
                        CROSS_TAG="_nVtx-gt${CROSS_LOW}"
                else
                        CROSS_TAG="_nVtx${CROSS_LOW}-${CROSS_HIGH}"
                fi
        elif [ "${ARG2}" == "2" ]
        then
                if [ "${CROSS_HIGH}" -ge 10000 ]
                then
                        CROSS_TAG="_pT-gt${CROSS_LOW}"
                else
                        CROSS_TAG="_pT${CROSS_LOW}-${CROSS_HIGH}"
                fi
        fi
fi
echo "CROSS_TAG = ${CROSS_TAG}"

cd ${_CONDOR_SCRATCH_DIR}
echo "source /cvmfs/cms.cern.ch/cmsset_default.sh"
source /cvmfs/cms.cern.ch/cmsset_default.sh

echo "scramv1 project CMSSW CMSSW_14_1_0_pre4"
scramv1 project CMSSW CMSSW_14_1_0_pre4

echo "cd CMSSW_14_1_0_pre4/src/"
cd CMSSW_14_1_0_pre4/src/

echo "eval `scramv1 runtime -sh`"
eval `scramv1 runtime -sh`
tar -zxvf ../../files.tar.gz
cd SUSYAnalysis/include && make
cd ../test/eleFakePho/fit
bash make.sh
voms-proxy-init --voms cms --valid 168:00 -out ~/.globus/gridproxy.cert

if [ -n "${ARG8}" ] && [ -n "${ARG9}" ]
then
        echo "./FitKer.exe ${year} ${preVFP} ${isData} ${doEB} ${ARG1} ${ARG2} ${ARG3} ${ARG4} ${ARG5} ${ARG6} ${ARG7} ${ARG8} ${ARG9}"
        ./FitKer.exe ${year} ${preVFP} ${isData} ${doEB} ${ARG1} ${ARG2} ${ARG3} ${ARG4} ${ARG5} ${ARG6} ${ARG7} ${ARG8} ${ARG9}
else
        echo "./FitKer.exe ${year} ${preVFP} ${isData} ${doEB} ${ARG1} ${ARG2} ${ARG3} ${ARG4} ${ARG5} ${ARG6} ${ARG7}"
        ./FitKer.exe ${year} ${preVFP} ${isData} ${doEB} ${ARG1} ${ARG2} ${ARG3} ${ARG4} ${ARG5} ${ARG6} ${ARG7}
fi

if [ ${year} == 2016 -a ${isData} == 1 -a ${preVFP} == 1 ]
then
        DEST_FIT=DATAFitting2016preVFP${CROSS_TAG}
        DEST_RES=DATAResult2016preVFP${CROSS_TAG}
elif [ ${year} == 2016 -a ${isData} == 1 -a ${preVFP} == 0 ]
then
        DEST_FIT=DATAFitting2016postVFP${CROSS_TAG}
        DEST_RES=DATAResult2016postVFP${CROSS_TAG}
elif [ ${year} == 2017 -a ${isData} == 1 ]
then
        DEST_FIT=DATAFitting2017${CROSS_TAG}
        DEST_RES=DATAResult2017${CROSS_TAG}
elif [ ${year} == 2018 -a ${isData} == 1 ]
then
        DEST_FIT=DATAFitting2018${CROSS_TAG}
        DEST_RES=DATAResult2018${CROSS_TAG}
fi

if [ ${year} == 2016 -a ${isData} == 0 -a ${preVFP} == 1 ]
then
        DEST_FIT=DYFitting2016preVFP${CROSS_TAG}
        DEST_RES=DYResult2016preVFP${CROSS_TAG}
elif [ ${year} == 2016 -a ${isData} == 0 -a ${preVFP} == 0 ]
then
        DEST_FIT=DYFitting2016postVFP${CROSS_TAG}
        DEST_RES=DYResult2016postVFP${CROSS_TAG}
elif [ ${year} == 2017 -a ${isData} == 0 ]
then
        DEST_FIT=DYFitting2017${CROSS_TAG}
        DEST_RES=DYResult2017${CROSS_TAG}
elif [ ${year} == 2018 -a ${isData} == 0 ]
then
        DEST_FIT=DYFitting2018${CROSS_TAG}
        DEST_RES=DYResult2018${CROSS_TAG}
fi

echo "Destination folders: ${DEST_FIT}  ${DEST_RES}"
xrdfs cmseos.fnal.gov mkdir -p /store/user/tmishra/elefakepho/${DEST_FIT} 2>/dev/null
xrdfs cmseos.fnal.gov mkdir -p /store/user/tmishra/elefakepho/${DEST_RES} 2>/dev/null

xrdcp -f *.png    root://cmseos.fnal.gov//store/user/tmishra/elefakepho/${DEST_FIT}
xrdcp -f Ele*.txt root://cmseos.fnal.gov//store/user/tmishra/elefakepho/${DEST_RES}

rm *.txt  *.png
cd ${_CONDOR_SCRATCH_DIR} && rm -rf CMSSW_14_1_0_pre4
