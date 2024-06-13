#!/bin/bash
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

ch=2
anatype=0
lmt=0
hmt=-1
lmet=40
hmet=70
iso=4


rm ../Background/BkgPredConfig.txt
echo 'ichannel' $ch  >> ../Background/BkgPredConfig.txt
echo 'anatype'  $anatype >>  ../Background/BkgPredConfig.txt
echo 'lowMt'    $lmt >> ../Background/BkgPredConfig.txt
echo 'highMt'   $hmt >> ../Background/BkgPredConfig.txt
echo 'lowMET'   $lmet >>../Background/BkgPredConfig.txt
echo 'highMET'  $hmet >>../Background/BkgPredConfig.txt
echo 'lowPt'    0 >> ../Background/BkgPredConfig.txt
echo 'highPt'   1000 >> ../Background/BkgPredConfig.txt
echo 'lepIso'   $iso    >> ../Background/BkgPredConfig.txt
root -l -q analysis_VGBkg.C++
root -l -q analysis_eleBkg.C++
root -l -q analysis_jetBkg.C++
root -l -q analysis_qcdBkg.C++
root -l -q analysis_rareBkg.C++
root -l -q analysis_sig.C++
