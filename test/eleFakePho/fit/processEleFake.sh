#!/bin/bash
# After the RunFits.sh condor jobs of one era are done:
#   1. checkFits.py  : every pt/eta/vtx bin has a fresh output on EOS
#   2. makeTxt.py    : merge the per-bin files (main folder + 6 pT/nVtx slices)
#   3. runFakeFunc.sh: fit the fake-rate pT/nVtx dependence
# fitFakeFunc_pT_nVtx.C reads only DATA (doDrellYan = false), so DATA must be complete;
# DY is merged too when it is complete, otherwise it is skipped with a warning.
# usage: bash processEleFake.sh <era>    era: 2016preVFP 2016postVFP 2017 2018
set -e

era=$1
case ${era} in
  2016preVFP)  tag=2016_1 ;;
  2016postVFP) tag=2016_0 ;;
  2017|2018)   tag=${era}_0 ;;
  *) echo "usage: bash processEleFake.sh <2016preVFP|2016postVFP|2017|2018>"; exit 1 ;;
esac

RANGES="pT30-50 pT50-80 pT-gt80 nVtx8-20 nVtx20-30 nVtx-gt30"

merge() {
  local type=$1 isData=$2
  # Batch submit time: skip per-bin files left over from older runs
  local since=$(date -d @$(stat -c %Y condorSubFile/condor_sub_${tag}_${isData}_1_1.txt) '+%Y-%m-%d %H:%M')
  python3 makeTxt.py --Year ${era} --Type ${type} --Since "${since}"
  for r in ${RANGES}; do
    python3 makeTxt.py --Year ${era} --Type ${type} --Range ${r} --Since "${since}"
  done
}

if ! python3 checkFits.py --Year ${era} --Type DATA; then
  echo "${era} DATA is not complete, resubmit the missing jobs first"
  exit 1
fi
merge DATA 1

if python3 checkFits.py --Year ${era} --Type DY > /dev/null; then
  merge DY 0
else
  echo "WARNING: ${era} DY is not complete, not merged (python3 checkFits.py --Year ${era} --Type DY)"
fi

bash runFakeFunc.sh ${era}
