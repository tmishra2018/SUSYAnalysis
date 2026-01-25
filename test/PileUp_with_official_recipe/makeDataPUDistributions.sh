#!/bin/bash

# ============================================================
# Common settings (PU POG Run-2 recommendation)
# ============================================================

PU_BASE=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification
OUTDIR=PUfiles

MINBIAS_NOM=69200
MINBIAS_UP=72383
MINBIAS_DN=66017

MAXPU=200
NBINS=200

mkdir -p ${OUTDIR}

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

    pileupCalc.py -i ${JSON} --inputLumiJSON ${PUFILE} \
      --calcMode true --minBiasXsec ${XSEC} \
      --maxPileupBin ${MAXPU} --numPileupBins ${NBINS} \
      ${OUTDIR}/dataPU_${TAG}${SUFFIX}.root
  done
}

# ============================================================
# 2016 (single JSON, split by VFP)
# ============================================================

JSON_2016=${PU_BASE}/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt

run_pu ${JSON_2016} \
  ${PU_BASE}/Collisions16/13TeV/PileUp/UltraLegacy/pileup_latest_preVFP.txt \
  2016preVFP

run_pu ${JSON_2016} \
  ${PU_BASE}/Collisions16/13TeV/PileUp/UltraLegacy/pileup_latest_postVFP.txt \
  2016postVFP

# ============================================================
# 2017
# ============================================================

run_pu \
  ${PU_BASE}/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt \
  ${PU_BASE}/Collisions17/13TeV/PileUp/UltraLegacy/pileup_latest.txt \
  2017

# ============================================================
# 2018
# ============================================================

run_pu \
  ${PU_BASE}/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt \
  ${PU_BASE}/Collisions18/13TeV/PileUp/UltraLegacy/pileup_latest.txt \
  2018

echo "All Run-2 UL data PU profiles (Nominal / Up / Down) produced."
