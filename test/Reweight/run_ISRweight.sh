#!/usr/bin/env bash
set -euo pipefail

# Build and run plot_ISRweight.exe for all eras:
#   1) derive ISR weights (applyISR=0)
#   2) apply derived ISR weights (applyISR=1)

if ! command -v root-config >/dev/null 2>&1; then
  echo "[ERROR] root-config not found. Please setup ROOT/CMSSW first."
  exit 1
fi

unset PYTHONHOME || true
unset PYTHONPATH || true

echo "[BUILD] plot_ISRweight.exe"
g++ $(root-config --cflags) ../../lib/libAnaClasses.so plot_ISRweight.C -o plot_ISRweight.exe $(root-config --libs)

for RunYear in 2016 2017 2018; do
  for preVFP in 0 1; do
    if [[ "${RunYear}" != "2016" && "${preVFP}" -eq 0 ]]; then
      continue
    fi

    if [[ "${RunYear}" -eq 2016 ]]; then
      VFP_string=$([[ "${preVFP}" -eq 1 ]] && echo "preVFP" || echo "postVFP")
    else
      VFP_string=""
    fi

    cat > ../Background/BkgPredConfig.txt <<CFG
RunYear ${RunYear}
preVFP ${preVFP}
CFG

    echo "[CFG] Wrote BkgPredConfig.txt for year=${RunYear} preVFP=${preVFP} ${VFP_string}"

    echo "[DERIVE] ./plot_ISRweight.exe ${RunYear} ${preVFP} 0"
    ./plot_ISRweight.exe "${RunYear}" "${preVFP}" 0

    echo "[APPLY ] ./plot_ISRweight.exe ${RunYear} ${preVFP} 1"
    ./plot_ISRweight.exe "${RunYear}" "${preVFP}" 1
  done
done

echo "[DONE] ISR derive+apply cycle finished for all eras."
