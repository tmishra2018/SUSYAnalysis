#!/bin/bash
set -e  # Exit immediately if any command fails
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mkdir -p logs

ToDeriveScale=true
makeCRPlots=false
makeVRPlots=false
makeSRPlots=false

# Define all regions to run
declare -a regions=()

if [ "$ToDeriveScale" = true ]; then
  regions+=("SCALE")
fi
if [ "$makeCRPlots" = true ]; then
  regions+=("CR")
fi
if [ "$makeVRPlots" = true ]; then
  regions+=("VR")
fi
if [ "$makeSRPlots" = true ]; then
  regions+=("SR")
fi

iso=4

for region in "${regions[@]}"; do

  # Set parameters based on region
  case "$region" in
    SCALE)
      anatype=0
      lmet=40; hmet=70
      lmt=0;  hmt=-1
      ranges=("0 1000")
      ;;
    CR)
      anatype=1
      lmet=0; hmet=70
      lmt=0;  hmt=-1
      ranges=("0 1000")
      ;;
    VR)
      anatype=2
      lmet=0; hmet=-1
      lmt=0;  hmt=100
      ranges=("0 1000")
      ;;
    SR)
      anatype=3
      lmet=120; hmet=-1
      lmt=100;  hmt=-1
      ranges=("0 1000")
      ;;
  esac

  echo ""
  echo "############################################"
  echo "#  Starting region: $region"
  echo "############################################"

  for range in "${ranges[@]}"; do
    lpt=$(echo "$range" | cut -d' ' -f1)
    hpt=$(echo "$range" | cut -d' ' -f2)
    echo "==== Running for lpt=$lpt, hpt=$hpt ===="

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

        # Set log file based on region
        case "$region" in
          CR)   LOGFILE="logs/BKG_${RunYear}${VFP_string}.log" ;;
          VR)   LOGFILE="logs/VALID_${RunYear}${VFP_string}.log" ;;
          *)    LOGFILE="" ;;
        esac

        if [ -n "$LOGFILE" ]; then
          : > "$LOGFILE"
        fi

        for ch in 1 2; do
          echo "---- Processing Region=$region, RunYear=$RunYear, VFP=$VFP_string, ichannel=$ch ----"

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

          # root -l -q analysis_eleBkg.C++ &
          # root -l -q analysis_jetBkg.C++ &
          
          root -l -q analysis_VGBkg.C++  &
	  root -l -q analysis_eleBkg_pT_nVtx.C++ &
          root -l -q analysis_jetBkg_pT_nVtx.C++ &
          root -l -q analysis_qcdBkg.C++ &
          root -l -q analysis_rareBkg.C++ &
          root -l -q analysis_sig.C++    &
          wait
          echo "All 6 macros finished."

          # Plot + log (CR / VR only)
          if [ -n "$LOGFILE" ]; then
            {
              echo "=============================================="
              echo "Region=$region  RunYear=$RunYear  VFP=$VFP_string  ichannel=$ch"
              echo "=============================================="
              root -l -q "plot_bkg.C($ch, $RunYear, $preVFP)"
              echo
            } >> "$LOGFILE" 2>&1
          fi

        done  # ch
      done    # preVFP
    done      # RunYear
  done        # range

  echo ""
  echo "Finished region: $region"

done  # region

echo ""
echo "Script completed successfully."
