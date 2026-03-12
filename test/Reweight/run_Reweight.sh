#!/usr/bin/env bash
set -euo pipefail

# Parallel ISR reweight workflow runner
# - Compiles analysis executables
# - Runs data, MC, and hadron steps in parallel batches
# - Preserves phase ordering while parallelizing jobs within each phase

MAX_JOBS="${MAX_JOBS:-8}"

echo "start"
echo "MAX_JOBS=${MAX_JOBS}"

# Avoid external CMSSW/CRAB python contamination during job execution.
unset PYTHONHOME || true
unset PYTHONPATH || true

declare -a years=("2017" "2018" "2016")
declare -a isr_runs_2017=("D" "E" "F")
declare -a isr_runs_2018=("A" "B" "C" "D")
declare -a isr_runs_2016=("B-ver1" "B-ver2" "C" "D" "E" "F_preVFP" "F_postVFP" "G" "H")
declare -a mc_samples=("TTGJets" "WWG" "WZG" "TTJets" "ZGToLLG")


if ! command -v root-config >/dev/null 2>&1; then
  echo "[ERROR] root-config not found. Please initialize a ROOT/CMSSW environment first."
  exit 1
fi

for src in analysis_ISR.C analysis_ISRMC.C analysis_mgHadron.C; do
  if [[ ! -f "${src}" ]]; then
    echo "[ERROR] Missing source file: ${src}"
    exit 1
  fi
done

active_jobs=0
failed_jobs=0

declare -a pids=()
declare -a labels=()

wait_for_one() {
  local i pid rc
  for i in "${!pids[@]}"; do
    pid="${pids[$i]}"
    if ! kill -0 "${pid}" 2>/dev/null; then
      if wait "${pid}"; then
        rc=0
        echo "[OK]   ${labels[$i]}"
      else
        rc=$?
        echo "[FAIL] ${labels[$i]} (exit ${rc})"
        failed_jobs=$((failed_jobs + 1))
      fi

      unset 'pids[i]'
      unset 'labels[i]'
      pids=("${pids[@]}")
      labels=("${labels[@]}")
      active_jobs=$((active_jobs - 1))
      return 0
    fi
  done

  sleep 1
}

wait_for_all() {
  while (( active_jobs > 0 )); do
    wait_for_one
  done
}

enqueue_job() {
  local label="$1"
  shift

  while (( active_jobs >= MAX_JOBS )); do
    wait_for_one
  done

  echo "[RUN] ${label}: $*"
  (
    unset PYTHONHOME PYTHONPATH
    "$@"
  ) &

  pids+=("$!")
  labels+=("${label}")
  active_jobs=$((active_jobs + 1))
}

compile_one() {
  local out="$1"
  local src="$2"
  echo "[BUILD] ${out} from ${src}"
  g++ $(root-config --cflags) ../../lib/libAnaClasses.so "${src}" -o "${out}" $(root-config --libs)
}

# Build executables
compile_one analysis_ISR.exe analysis_ISR.C
compile_one analysis_ISRMC.exe analysis_ISRMC.C
compile_one analysis_mgHadron.exe analysis_mgHadron.C

# Phase 1: data ISR for 2017/2018
for year in 2017 2018; do
  if [[ "${year}" == "2017" ]]; then
    runs=("${isr_runs_2017[@]}")
  else
    runs=("${isr_runs_2018[@]}")
  fi

  for run in "${runs[@]}"; do
    enqueue_job "analysis_ISR ${year} ${run}" ./analysis_ISR.exe "${year}" "${run}"
  done
done
wait_for_all

# Phase 2: MC ISR for 2017/2018
for year in 2017 2018; do
  for sample in "${mc_samples[@]}"; do
    enqueue_job "analysis_ISRMC ${year} 0 ${sample}" ./analysis_ISRMC.exe "${year}" "0" "${sample}"
  done
done
wait_for_all

# Phase 3: mg hadron for 2017/2018
for year in 2017 2018; do
  if [[ "${year}" == "2017" ]]; then
    runs=("${isr_runs_2017[@]}")
  else
    runs=("${isr_runs_2018[@]}")
  fi

  for run in "${runs[@]}"; do
    enqueue_job "analysis_mgHadron ${year} ${run}" ./analysis_mgHadron.exe "${year}" "${run}"
  done
done
wait_for_all

# Phase 4: 2016 MC versions
for sample in "${mc_samples[@]}"; do
  for version in 0 1; do
    enqueue_job "analysis_ISRMC 2016 ${version} ${sample}" ./analysis_ISRMC.exe "2016" "${version}" "${sample}"
  done
done
wait_for_all

# Phase 5: 2016 data ISR + hadron
for run in "${isr_runs_2016[@]}"; do
  enqueue_job "analysis_ISR 2016 ${run}" ./analysis_ISR.exe "2016" "${run}"
  enqueue_job "analysis_mgHadron 2016 ${run}" ./analysis_mgHadron.exe "2016" "${run}"
done
wait_for_all

# Cleanup
rm -f analysis_ISRMC.exe analysis_ISR.exe analysis_mgHadron.exe

if (( failed_jobs > 0 )); then
  echo "done with ${failed_jobs} failed job(s)"
  exit 1
fi

echo "done"
