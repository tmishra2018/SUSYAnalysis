#!/usr/bin/env bash
set -euo pipefail

# Jet->photon fake rate in the Z(ll)+gamma ISR control region, per
# channel x detector x era. Inputs come from analysis_Hadron.C (step 1) and
# the GJet MC template from ../jetFakePho/src/analysis_mgGJet.C.
# Output (read by plot_ISRweight.C):
#   /eos/uscms/store/user/tmishra/jetfakepho/files/JetFakeRate-ISR-{EB|EE}_{mg|egamma}_{year}{VFP}.txt

outdir=/eos/uscms/store/user/tmishra/jetfakepho/files

# photon pT bin edges as in the 2016 analysis (mengleis .../test/Reweight/runFit.sh);
# photon pT is capped at 799 GeV in plot_ISRweight.C, so the last bin is [300,800)
ptbins=( 35 40 50 60 70 80 100 150 200 250 300 800 )

# Detector types: 1=EB, 2=EE
detTypes=( 1 2 )

channels=( mg egamma )

# "year preVFP" pairs
eras=( "2016 1" "2016 0" "2017 0" "2018 0" )

for channel in "${channels[@]}"; do
  echo ""
  echo "===== Running FitJetFake for channel: ${channel} ====="

  for era in "${eras[@]}"; do
    read -r year preVFP <<< "${era}"
    if [[ "${year}" -eq 2016 ]]; then
      vfp=$([[ "${preVFP}" -eq 1 ]] && echo "preVFP" || echo "postVFP")
    else
      vfp=""
    fi

    for detType in "${detTypes[@]}"; do
      det_label=$([[ "${detType}" -eq 1 ]] && echo "EB" || echo "EE")
      outfile="${outdir}/JetFakeRate-ISR-${det_label}_${channel}_${year}${vfp}.txt"
      echo "  ${year}${vfp} ${det_label}: removing ${outfile}"
      rm -f "${outfile}"

      for (( i = 0; i < ${#ptbins[@]} - 1; i++ )); do
        j=$((i + 1))
        ptlo="${ptbins[$i]}"
        pthi="${ptbins[$j]}"

        echo "    pT bin: [${ptlo}, ${pthi})  channel=${channel}  det=${det_label}  era=${year}${vfp}"

        root -b -q "FitJetFake.C+(${ptlo}, ${pthi}, ${detType}, ${year}, ${preVFP}, \"${channel}\")"
      done
    done  # detType
  done    # era
done      # channel

echo ""
echo "[DONE] FitJetFake completed for all bins, detectors, eras, and channels."
