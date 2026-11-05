#!/bin/bash
# Fit the e->gamma fake-rate pT/nVtx dependence (fitFakeFunc_pT_nVtx.C) from the merged makeTxt.py files.
# usage: bash runFakeFunc.sh [era]     era: 2016preVFP 2016postVFP 2017 2018 (default: all four)
# Log: runFakeFunc_<era>_<date>.log
set -e

ERAS=${1:-"2016preVFP 2016postVFP 2017 2018"}

for era in ${ERAS}; do
  case ${era} in
    2016preVFP)  args="2016,1" ;;
    2016postVFP) args="2016,0" ;;
    2017|2018)   args="${era},1" ;;
    *) echo "Unknown era ${era}"; exit 1 ;;
  esac
  log=runFakeFunc_${era}_$(date +%F).log
  echo "Running fitFakeFunc_pT_nVtx.C++(${args}) -> ${log}"
  root -l -b -q "fitFakeFunc_pT_nVtx.C++(${args})" > ${log} 2>&1
  if grep -q "Cannot open" ${log}; then
    echo "WARNING: missing input files, see ${log}:"; grep "Cannot open" ${log}
  fi
done


#root -l -q "fitFakeFunc.C++(2016,1)"
#root -l -q "fitFakeFunc.C++(2016,0)"
#root -l -q "fitFakeFunc.C++(2017,1)"
#root -l -q "fitFakeFunc.C++(2018,1)"
#root -l -q "FakeRates_Data_vs_DY.C++(2016,1)"
#root -l -q "FakeRates_Data_vs_DY.C++(2016,0)"
#root -l -q "FakeRates_Data_vs_DY.C++(2017,1)"
#root -l -q "FakeRates_Data_vs_DY.C++(2018,1)"

#root -l -q "fitFake_HT.C++(2016,1)"
#root -l -q "fitFake_HT.C++(2016,0)"
#root -l -q "fitFake_HT.C++(2017,1)"
#root -l -q "fitFake_HT.C++(2018,1)"

#root -l -q "fitFakeFunc_HT.C++(2016,1)"
#root -l -q "fitFakeFunc_HT.C++(2016,0)"
#root -l -q "fitFakeFunc_HT.C++(2017,1)"
#root -l -q "fitFakeFunc_HT.C++(2018,1)"
