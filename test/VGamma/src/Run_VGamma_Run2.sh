#!/usr/bin/env bash
set -euo pipefail

MAX_PARALLEL="${1:-}"
if [[ -z "${MAX_PARALLEL}" ]]; then
  if command -v nproc >/dev/null 2>&1; then
    MAX_PARALLEL="$(nproc)"
  else
    MAX_PARALLEL=4
  fi
fi

if ! [[ "${MAX_PARALLEL}" =~ ^[0-9]+$ ]] || [[ "${MAX_PARALLEL}" -lt 1 ]]; then
  echo "ERROR: MAX_PARALLEL must be a positive integer." >&2
  exit 1
fi

g++ `root-config --cflags` ../../../lib/libAnaClasses.so analysis_VGamma_Run2.C -o analysis_VGamma_Run2.exe `root-config --libs`

SAMPLES=(
  TTJets
  TTGJets
  WWG
  WZG
  WW
  WZ
  ZZ
  DYJetsToLL
  ZGToLLG
  WGToLNuG
  WGJet40
  WGJet130
)

RUN_PERIODS=(
#  "2016 0"
#  "2016 1"
#  "2017 0"
  "2018 0"
)

mkdir -p logs_run2_parallel

running=0
for period in "${RUN_PERIODS[@]}"; do
  read -r year prevfp <<< "${period}"
  for sample in "${SAMPLES[@]}"; do
    log="logs_run2_parallel/${sample}_${year}_${prevfp}.log"

    echo "[LAUNCH] year=${year} preVFP=${prevfp} sample=${sample}"
    ./analysis_VGamma_Run2.exe "${year}" "${prevfp}" "${sample}" > "${log}" 2>&1 &

    ((running+=1))
    if (( running >= MAX_PARALLEL )); then
      wait -n
      ((running-=1))
    fi
  done
done

# Drain remaining jobs
wait

echo "All Run2 jobs finished. Logs are in: logs_run2_parallel/"
