#!/usr/bin/env bash
set -euo pipefail

if ! command -v root-config >/dev/null 2>&1; then
  echo "[ERROR] root-config not found. Please setup ROOT/CMSSW first."
  exit 1
fi

unset PYTHONHOME || true
unset PYTHONPATH || true

echo "[BUILD] plot_ISRweight.exe"
g++ $(root-config --cflags) ../../lib/libAnaClasses.so plot_ISRweight.C -o plot_ISRweight.exe $(root-config --libs)

while IFS=':' read -r channel_id channel_name; do
  for RunYear in 2016 2017 2018; do
    for preVFP in 0 1; do
      if [[ "${RunYear}" != "2016" && "${preVFP}" -eq 1 ]]; then
        continue
      fi

      if [[ "${RunYear}" -eq 2016 ]]; then
        VFP_string=$([[ "${preVFP}" -eq 1 ]] && echo "preVFP" || echo "postVFP")
      else
        VFP_string=""
      fi

      # own config (read by esfScaleFactor in analysis_scalefactor.h) so the
      # shared ../Background/BkgPredConfig.txt cuts are not overwritten
      cat > ISRConfig.txt <<CFG
RunYear ${RunYear}
preVFP ${preVFP}
CFG

      echo "[CFG] Wrote ISRConfig.txt for channel=${channel_name} year=${RunYear} preVFP=${preVFP} ${VFP_string}"

      echo "[DERIVE] ./plot_ISRweight.exe ${RunYear} ${preVFP} 0 ${channel_id} (${channel_name})"
      ./plot_ISRweight.exe "${RunYear}" "${preVFP}" 0 "${channel_id}"

      echo "[APPLY ] ./plot_ISRweight.exe ${RunYear} ${preVFP} 1 ${channel_id} (${channel_name})"
      ./plot_ISRweight.exe "${RunYear}" "${preVFP}" 1 "${channel_id}"
    done
  done
done <<'CHANNELS'
1:eg
2:mg
CHANNELS

echo "[DONE] ISR derive+apply cycle finished for all eras and channels."
