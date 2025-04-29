#!/bin/bash

dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

temp_file="temp.txt"
rm -f $temp_file

ToDeriveScale=false
makeCRPlots=false
makeVRPlots=true
makeSRPlots=false

if [ "$ToDeriveScale" = true ]; then
  anatype=0
  lmet=40
  hmet=70
  lmt=0
  hmt=-1
elif [ "$makeCRPlots" = true ]; then
  anatype=1
  lmet=0
  hmet=70
  lmt=0
  hmt=-1
elif [ "$makeVRPlots" = true ]; then
  anatype=2 # validation region 0 < MT < 100
  lmet=0
  hmet=-1
  lmt=0
  hmt=100
elif [ "$makeSRPlots" = true ]; then
  anatype=3 # Signal region
  lmet=120
  hmet=-1
  lmt=100
  hmt=-1
else
  echo "Neither ToDeriveScale nor makeCRPlots nor makeVRPlots is set to true."
  exit 1
fi

iso=4
lpt=0
hpt=1000


for RunYear in 2016 2017 2018; do
  for preVFP in {0..1}; do
    if [ "$RunYear" != "2016" ] && [ "$preVFP" -eq 1 ]; then
      continue
    fi

    for ch in 1 2; do
      echo "Processing RunYear=$RunYear, preVFP=$preVFP, channel=$ch"

      rm -f BkgPredConfig.txt
      {
        echo "ichannel $ch"
        echo "anatype $anatype"
        echo "lowMt $lmt"
        echo "highMt $hmt"
        echo "lowMET $lmet"
        echo "highMET $hmet"
        echo "lowPt $lpt"
        echo "highPt $hpt"
        echo "lepIso $iso"
        echo "RunYear $RunYear"
        echo "preVFP $preVFP"
      } >> BkgPredConfig.txt

      root -l -q analysis_VGBkg.C++
      root -l -q analysis_eleBkg.C++
      root -l -q analysis_jetBkg.C++
      root -l -q analysis_qcdBkg.C++
      root -l -q analysis_rareBkg.C++
      root -l -q analysis_sig.C++

      if [ "$makeCRPlots" = true ] || [ "$makeVRPlots" = true ]; then
        echo "Plotting background for ch=$ch, RunYear=$RunYear, preVFP=$preVFP"
        root -l -q "plot_bkg.C($ch, $RunYear, $preVFP)" >> $temp_file
      fi
    done
  done
done

grep egamma temp.txt
grep mgamma temp.txt


if [ "$makeCRPlots" = true ]; then
    scp -r /eos/uscms/store/user/tmishra/Background/plots/* trmishra@lxplus.cern.ch:/eos/home-t/trmishra/www/Plots/SUSYAnalysis/ControlRegion
elif [ "$makeVRPlots" = true ]; then
    scp -r /eos/uscms/store/user/tmishra/Background/plots/* trmishra@lxplus.cern.ch:/eos/home-t/trmishra/www/Plots/SUSYAnalysis/ValidationRegion
fi
