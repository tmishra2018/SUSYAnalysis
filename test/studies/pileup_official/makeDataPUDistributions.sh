#!/bin/bash

# ============================================================
# Common settings (PU POG Run-2 recommendation)
# ============================================================

# 2026-10-18: use the UL golden JSONs (same files as test/golden/golden_<year>.json applied to data)
# instead of the PromptReco JSONs; inputs from cvmfs (the /afs certification area is not available at LPC)
PU_BASE=/cvmfs/cms-griddata.cern.ch/cat/metadata/DC
OUTDIR=PUfiles

MINBIAS_NOM=69200
MINBIAS_UP=72383
MINBIAS_DN=66017

MAXPU=200
NBINS=200

mkdir -p ${OUTDIR}

# pileupCalc.py in CMSSW_14_0_7 uses np.int (removed in numpy >= 1.24): run a patched copy
PILEUPCALC=${OUTDIR}/pileupCalc_fix.py
sed 's/astype(np\.int)/astype(int)/g' $(which pileupCalc.py) > ${PILEUPCALC}

# ------------------------------------------------------------
# Helper function
# ------------------------------------------------------------
run_pu() {
  local JSON=$1
  local PUFILE=$2
  local TAG=$3

  for VAR in NOM UP DN; do
    case ${VAR} in
      NOM) XSEC=${MINBIAS_NOM}; SUFFIX="" ;;
      UP)  XSEC=${MINBIAS_UP};  SUFFIX="_PUUp" ;;
      DN)  XSEC=${MINBIAS_DN};  SUFFIX="_PUDown" ;;
    esac

    python3 ${PILEUPCALC} -i ${JSON} --inputLumiJSON ${PUFILE} \
      --calcMode true --minBiasXsec ${XSEC} \
      --maxPileupBin ${MAXPU} --numPileupBins ${NBINS} \
      ${OUTDIR}/dataPU_${TAG}${SUFFIX}.root
  done
}

# ============================================================
# 2016 (single golden JSON, split by VFP)
# ============================================================

JSON_2016=${PU_BASE}/Collisions16/latest/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt
PU_2016=${PU_BASE}/Collisions16/latest/13TeV/PileUp/UltraLegacy/pileup_latest.txt

# Split as in the UL datasets (checked in the 2016F_preVFP/postVFP data):
# preVFP = runs <= 278768 and 278770, postVFP = run 278769 and runs >= 278771
python3 - ${JSON_2016} ${OUTDIR} <<'EOF'
import json, sys
g = json.load(open(sys.argv[1]))
pre = lambda r: r <= 278768 or r == 278770
json.dump({r: v for r, v in g.items() if pre(int(r))},     open(sys.argv[2] + "/golden_2016preVFP.json", "w"))
json.dump({r: v for r, v in g.items() if not pre(int(r))}, open(sys.argv[2] + "/golden_2016postVFP.json", "w"))
EOF

run_pu ${OUTDIR}/golden_2016preVFP.json  ${PU_2016} 2016preVFP
run_pu ${OUTDIR}/golden_2016postVFP.json ${PU_2016} 2016postVFP

# ============================================================
# 2017
# ============================================================

run_pu \
  ${PU_BASE}/Collisions17/latest/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt \
  ${PU_BASE}/Collisions17/latest/13TeV/PileUp/UltraLegacy/pileup_latest.txt \
  2017

# ============================================================
# 2018
# ============================================================

run_pu \
  ${PU_BASE}/Collisions18/latest/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt \
  ${PU_BASE}/Collisions18/latest/13TeV/PileUp/UltraLegacy/pileup_latest.txt \
  2018

echo "All Run-2 UL data PU profiles (Nominal / Up / Down) produced."
