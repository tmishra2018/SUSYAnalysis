#!/usr/bin/env bash
# run_Reweight.sh
# Compiles and runs data/MC analysis steps for both mg and egamma channels
# across all years, in parallel batches.
#
# Usage:
#   bash run_Reweight.sh [MAX_JOBS=N]
#
# Output files written per channel:
#   analysis_ISR     -> resTree_ISR_data_{mg|egamma}_{year}{era}.root
#   analysis_ISRMC   -> resTree_ISR_{sample}_{mg|egamma}_{year}{vfp}.root
#   analysis_Hadron  -> plot_hadron_ISR_{mg|egamma}_{year}{era}.root

set -euo pipefail

MAX_JOBS="${MAX_JOBS:-8}"

echo "start"
echo "MAX_JOBS=${MAX_JOBS}"

unset PYTHONHOME 2>/dev/null || true
unset PYTHONPATH 2>/dev/null || true

# ---- Era definitions -------------------------------------------------------
declare -a years=("2017" "2018" "2016")
declare -a isr_runs_2017=("D" "E" "F")
declare -a isr_runs_2018=("A" "B" "C" "D")
declare -a isr_runs_2016=("B-ver1" "B-ver2" "C" "D" "E" "F_preVFP" "F_postVFP" "G" "H")
# 2017 B,C: egamma only (MuonEG lumi/trigger covers 2017 D-F, DoubleEG lumi B-F)
declare -a isr_runs_2017_egOnly=("B" "C")
declare -a mc_samples=("TTGJets" "WWG" "WZG" "TTJets" "ZGToLLG")
declare -a channels=("mg" "egamma")

ISRDIR=/eos/uscms/store/user/tmishra/ISRweighting

# ---- Sanity checks ---------------------------------------------------------
if ! command -v root-config >/dev/null 2>&1; then
  echo "[ERROR] root-config not found. Please initialize a ROOT/CMSSW environment first."
  exit 1
fi

for src in analysis_ISR.C analysis_ISRMC.C analysis_Hadron.C; do
  if [[ ! -f "${src}" ]]; then
    echo "[ERROR] Missing source file: ${src}"
    exit 1
  fi
done

# ---- Parallel job infrastructure -------------------------------------------
active_jobs=0
failed_jobs=0
declare -a pids=()
declare -a labels=()

wait_for_one(){
  local i pid rc
  for i in "${!pids[@]}"; do
    pid="${pids[$i]}"
    if ! kill -0 "${pid}" 2>/dev/null; then
      if wait "${pid}"; then
        rc=0; echo "[OK]   ${labels[$i]}"
      else
        rc=$?; echo "[FAIL] ${labels[$i]} (exit ${rc})"
        failed_jobs=$((failed_jobs + 1))
      fi
      unset 'pids[i]' 'labels[i]'
      pids=("${pids[@]+"${pids[@]}"}")
      labels=("${labels[@]+"${labels[@]}"}")
      active_jobs=$((active_jobs - 1))
      return 0
    fi
  done
  sleep 1
}

wait_for_all(){
  while (( active_jobs > 0 )); do wait_for_one; done
}

enqueue_job(){
  local label="$1"; shift
  while (( active_jobs >= MAX_JOBS )); do wait_for_one; done
  echo "[RUN] ${label}: $*"
  ( unset PYTHONHOME PYTHONPATH; "$@" ) &
  pids+=("$!"); labels+=("${label}")
  active_jobs=$((active_jobs + 1))
}

compile_one(){
  local out="$1" src="$2"
  echo "[BUILD] ${out} from ${src}"
  g++ $(root-config --cflags) \
      ../../lib/libAnaClasses.so "${src}" \
      -o "${out}" \
      $(root-config --libs)
}

# ---- Build executables -----------------------------------------------------
compile_one analysis_ISR.exe    analysis_ISR.C
compile_one analysis_ISRMC.exe  analysis_ISRMC.C
compile_one analysis_Hadron.exe analysis_Hadron.C

# ============================================================================
# Phase 1: Data ISR for 2017/2018
#   analysis_ISR runs both channels internally, one output file per channel.
# ============================================================================
echo ""
echo "[PHASE 1] Data ISR (2017/2018)"
for year in 2017 2018; do
  declare -n runs="isr_runs_${year}"
  for run in "${runs[@]}"; do
    enqueue_job "analysis_ISR ${year} ${run}" \
      ./analysis_ISR.exe "${year}" "${run}"
  done
done
for run in "${isr_runs_2017_egOnly[@]}"; do
  enqueue_job "analysis_ISR 2017 ${run} egamma" \
    ./analysis_ISR.exe "2017" "${run}" "egamma"
done
wait_for_all

# ============================================================================
# Phase 2: MC ISR for 2017/2018 — one job per sample × channel
#   mg     -> MCMuonEG{year}   RunType, MuonEG   lumi
#   egamma -> MCEGamma2018 / MCDoubleEG{year} RunType, matching lumi
# ============================================================================
echo ""
echo "[PHASE 2] MC ISR (2017/2018)"
for year in 2017 2018; do
  for sample in "${mc_samples[@]}"; do
    for ch in "${channels[@]}"; do
      enqueue_job "analysis_ISRMC ${year} 0 ${sample} ${ch}" \
        ./analysis_ISRMC.exe "${year}" "0" "${sample}" "${ch}"
    done
  done
done
wait_for_all

# ============================================================================
# Phase 3: Hadron proxy for 2017/2018 — one job per era × channel
# ============================================================================
echo ""
echo "[PHASE 3] Hadron proxy (2017/2018)"
for year in 2017 2018; do
  declare -n runs="isr_runs_${year}"
  for run in "${runs[@]}"; do
    for ch in "${channels[@]}"; do
      enqueue_job "analysis_Hadron ${year} ${run} ${ch}" \
        ./analysis_Hadron.exe "${year}" "${run}" "${ch}"
    done
  done
done
for run in "${isr_runs_2017_egOnly[@]}"; do
  enqueue_job "analysis_Hadron 2017 ${run} egamma" \
    ./analysis_Hadron.exe "2017" "${run}" "egamma"
done
wait_for_all

# ============================================================================
# Phase 4: 2016 MC (preVFP=1 and postVFP=0) — one job per sample × version × channel
# ============================================================================
echo ""
echo "[PHASE 4] MC ISR (2016 preVFP/postVFP)"
for sample in "${mc_samples[@]}"; do
  for version in 0 1; do
    for ch in "${channels[@]}"; do
      enqueue_job "analysis_ISRMC 2016 ${version} ${sample} ${ch}" \
        ./analysis_ISRMC.exe "2016" "${version}" "${sample}" "${ch}"
    done
  done
done
wait_for_all

# ============================================================================
# Phase 5: 2016 data ISR + hadron — one ISR job per era (both channels
# internal), one Hadron job per era × channel
# ============================================================================
echo ""
echo "[PHASE 5] Data ISR + Hadron (2016)"
for run in "${isr_runs_2016[@]}"; do
  enqueue_job "analysis_ISR 2016 ${run}" \
    ./analysis_ISR.exe "2016" "${run}"
  for ch in "${channels[@]}"; do
    enqueue_job "analysis_Hadron 2016 ${run} ${ch}" \
      ./analysis_Hadron.exe "2016" "${run}" "${ch}"
  done
done
wait_for_all

# ---- Cleanup ---------------------------------------------------------------
rm -f analysis_ISR.exe analysis_ISRMC.exe analysis_Hadron.exe

if (( failed_jobs > 0 )); then
  echo "done with ${failed_jobs} failed job(s); skipping merge"
  exit 1
fi

# ============================================================================
# Phase 6: merge into the inputs read by plot_ISRweight.C (step 3)
#   ${ISRDIR}/{egamma|mgamma}/resTree_ISR_{data|<sample>}_{year}{VFP}.root
# ============================================================================
echo ""
echo "[PHASE 6] Merge per-era data and per-channel MC for step 3"
merge_one(){
  local out="$1"; shift
  echo "[HADD] ${out} <- $*"
  hadd -f "${out}" "$@"
}

for ch in "${channels[@]}"; do
  chdir=$([[ "${ch}" == "mg" ]] && echo "mgamma" || echo "egamma")
  mkdir -p "${ISRDIR}/${chdir}"

  pre=(); post=()
  for run in "${isr_runs_2016[@]}"; do
    f="${ISRDIR}/resTree_ISR_data_${ch}_2016${run}.root"
    case "${run}" in F_postVFP|G|H) post+=("${f}");; *) pre+=("${f}");; esac
  done
  merge_one "${ISRDIR}/${chdir}/resTree_ISR_data_2016preVFP.root"  "${pre[@]}"
  merge_one "${ISRDIR}/${chdir}/resTree_ISR_data_2016postVFP.root" "${post[@]}"

  runs17=("${isr_runs_2017[@]}")
  [[ "${ch}" == "egamma" ]] && runs17=("${isr_runs_2017_egOnly[@]}" "${isr_runs_2017[@]}")
  merge_one "${ISRDIR}/${chdir}/resTree_ISR_data_2017.root" \
    $(for run in "${runs17[@]}"; do echo "${ISRDIR}/resTree_ISR_data_${ch}_2017${run}.root"; done)
  merge_one "${ISRDIR}/${chdir}/resTree_ISR_data_2018.root" \
    $(for run in "${isr_runs_2018[@]}"; do echo "${ISRDIR}/resTree_ISR_data_${ch}_2018${run}.root"; done)

  for sample in "${mc_samples[@]}"; do
    for era in 2016preVFP 2016postVFP 2017 2018; do
      merge_one "${ISRDIR}/${chdir}/resTree_ISR_${sample}_${era}.root" \
        "${ISRDIR}/resTree_ISR_${sample}_${ch}_${era}.root"
    done
  done
done

echo "done"
