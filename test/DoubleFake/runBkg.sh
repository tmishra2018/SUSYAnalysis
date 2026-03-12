#!/bin/bash
set -e  # Exit immediately if any command fails
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ToDeriveScale=false
makeCRPlots=false
makeVRPlots=true
makeSRPlots=false
# ---------------- Region configuration ----------------
if [ "$ToDeriveScale" = true ]; then
  anatype=0
  lmet=40
  hmet=70
  lmt=0
  hmt=-1
  ranges=("0 1000")

elif [ "$makeCRPlots" = true ]; then
  anatype=1
  lmet=0
  hmet=70
  lmt=0
  hmt=-1
  ranges=("0 1000")

elif [ "$makeVRPlots" = true ]; then
  anatype=2   # Validation region
  lmet=0
  hmet=-1
  lmt=0
  hmt=100
  ranges=("0 1000")

elif [ "$makeSRPlots" = true ]; then
  anatype=3   # Signal region
  lmet=120
  hmet=-1
  lmt=100
  hmt=-1
  ranges=("0 1000")
fi

iso=4

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

      for ch in 1 2; do
        echo "---- Processing RunYear=$RunYear, VFP=$VFP_string, ichannel=$ch ----"

        # Create config file
        cat > ../Background/BkgPredConfig.txt <<EOF
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
#	root -l -q closure_doubleFake.C++
	root -l -q closure_doubleFake_v2.C++
#       root -l -q analysis_doubleFakeBkg.C++
#	root -l -q quantify_doubleFake.C++
      done
    done
  done
done

echo "✅ Script completed successfully."
