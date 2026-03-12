#!/bin/bash
set -e

dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Mode = "PTSCAN"  → scan photon pT bins, MET fixed to 40–70
# Mode = "METSCAN" → scan MET sub-bins, photon pT fixed to 0–1000

Mode="PTSCAN"
#Mode="METSCAN"

ToDeriveScale=true

if [ "$ToDeriveScale" = true ]; then
  anatype=0
  lmt=0
  hmt=-1
  iso=4

  # ---------------- Configuration by mode ----------------
  if [ "$Mode" = "PTSCAN" ]; then
    echo ">>> Running in PHOTON pT scan mode"
    ranges=("0 70" "70 1000")
    met_bins=("40 70")

  elif [ "$Mode" = "METSCAN" ]; then
    echo ">>> Running in MET sub-bin scan mode"
    ranges=("0 1000")
    met_bins=("40 50" "50 60" "60 70")

  fi

# ---------------- Main loops ----------------
for range in "${ranges[@]}"; do
  lpt=$(echo "$range" | cut -d' ' -f1)
  hpt=$(echo "$range" | cut -d' ' -f2)

  for met in "${met_bins[@]}"; do
    lmet=$(echo "$met" | cut -d' ' -f1)
    hmet=$(echo "$met" | cut -d' ' -f2)

    echo "=============================================="
    echo " pT = ${lpt}-${hpt} GeV | MET = ${lmet}-${hmet} GeV"
    echo "=============================================="

    for RunYear in 2016 2017 2018; do
      for preVFP in 0 1; do

        if [ "$RunYear" != "2016" ] && [ "$preVFP" -eq 1 ]; then
          continue
        fi

        if [ "$RunYear" -eq 2016 ]; then
          if [ "$preVFP" -eq 1 ]; then
            VFP_string="preVFP"
          else
            VFP_string="postVFP"
          fi
        else
          VFP_string=""
        fi

        for ch in 1 2; do
          echo "---- RunYear=$RunYear, VFP=$VFP_string, ch=$ch ----"

          cat > BkgPredConfig.txt <<EOF
ichannel $ch
anatype $anatype
lowMt $lmt
highMt $hmt
lowMET $lmet
highMET $hmet
lowPt $lpt
highPt $hpt
lepIso $iso
RunYear $RunYear
preVFP $preVFP
EOF
	  root -l -q analysis_VGBkg_ISRwt.C++ &
          root -l -q analysis_eleBkg_HT.C++ &
          root -l -q analysis_jetBkg_HT.C++ &
          root -l -q analysis_qcdBkg.C++ &
          root -l -q analysis_rareBkg.C++ &
          root -l -q analysis_sig.C++ &
	  wait
          echo "All 6 macros finished."

        done
      done
    done
  done
done

fi
