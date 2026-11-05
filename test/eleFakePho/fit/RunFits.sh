#!/bin/bash
# Submit the e->gamma TnP fits (FitKer.C) to condor.
# usage: bash RunFits.sh [eras] [samples] [doEB]
#   eras    : space-separated list from 2016preVFP 2016postVFP 2017 2018
#             (default "2017 2018 2016postVFP")
#   samples : "DY", "Data" or "DY Data" (default "DY Data")
#   doEB    : 1 = barrel, 0 = endcap (default 1)
# e.g. bash RunFits.sh "2017 2018 2016postVFP" "DY Data" 1
set -e

ERAS=${1:-"2017 2018 2016postVFP"}
SAMPLES=${2:-"DY Data"}
DOEB=${3:-1}

# Builds FitKer.exe, writes RunnerScript.txt and the fit configs, makes files.tar.gz
bash makeAll.sh

for era in ${ERAS}; do
  case ${era} in
    2016preVFP)  year=2016; preVFP=1 ;;
    2016postVFP) year=2016; preVFP=0 ;;
    2017|2018)   year=${era}; preVFP=0 ;;
    *) echo "Unknown era ${era}"; exit 1 ;;
  esac
  for sample in ${SAMPLES}; do
    case ${sample} in
      DY)   isData=0 ;;
      Data) isData=1 ;;
      *) echo "Unknown sample ${sample}"; exit 1 ;;
    esac
    echo "Submitting ${era} ${sample} doEB=${DOEB}"
    python3 command_HTcondor.py ${year} ${preVFP} ${isData} ${DOEB}
  done
done
