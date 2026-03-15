#!/bin/bash
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

ch=2
anatype=0
lmt=0
hmt=-1
lmet=40
hmet=70
iso=4

rm BkgPredConfig.txt
echo 'ichannel' $ch  >> BkgPredConfig.txt
echo 'anatype'  $anatype >>  BkgPredConfig.txt
echo 'lowMt'    $lmt >> BkgPredConfig.txt
echo 'highMt'   $hmt >> BkgPredConfig.txt
echo 'lowMET'   $lmet >>BkgPredConfig.txt
echo 'highMET'  $hmet >>BkgPredConfig.txt
echo 'lowPt'    0 >> BkgPredConfig.txt
echo 'highPt'   50 >> BkgPredConfig.txt
echo 'lepIso'   $iso    >> BkgPredConfig.txt
root -q analysis_VGBkg.C++
root -q analysis_eleBkg.C++
root -q analysis_jetBkg.C++
root -q analysis_qcdBkg.C++
root -q analysis_doubleFakeBkg.C++
root -q analysis_rareBkg.C++
root -q analysis_sig.C++


rm BkgPredConfig.txt
echo 'ichannel' $ch  >> BkgPredConfig.txt
echo 'anatype'  $anatype >>  BkgPredConfig.txt
echo 'lowMt'    $lmt >> BkgPredConfig.txt
echo 'highMt'   $hmt >> BkgPredConfig.txt
echo 'lowMET'   $lmet >>BkgPredConfig.txt
echo 'highMET'  $hmet >>BkgPredConfig.txt
echo 'lowPt'    50 >> BkgPredConfig.txt
echo 'highPt'   70 >> BkgPredConfig.txt
echo 'lepIso'   $iso    >> BkgPredConfig.txt
root -q analysis_VGBkg.C++
root -q analysis_eleBkg.C++
root -q analysis_jetBkg.C++
root -q analysis_qcdBkg.C++
root -q analysis_doubleFakeBkg.C++
root -q analysis_rareBkg.C++
root -q analysis_sig.C++

rm BkgPredConfig.txt
echo 'ichannel' $ch  >> BkgPredConfig.txt
echo 'anatype'  $anatype >>  BkgPredConfig.txt
echo 'lowMt'    $lmt >> BkgPredConfig.txt
echo 'highMt'   $hmt >> BkgPredConfig.txt
echo 'lowMET'   $lmet >>BkgPredConfig.txt
echo 'highMET'  $hmet >>BkgPredConfig.txt
echo 'lowPt'    70 >> BkgPredConfig.txt
echo 'highPt'   100 >> BkgPredConfig.txt
echo 'lepIso'   $iso    >> BkgPredConfig.txt
root -q analysis_VGBkg.C++
root -q analysis_eleBkg.C++
root -q analysis_jetBkg.C++
root -q analysis_qcdBkg.C++
root -q analysis_doubleFakeBkg.C++
root -q analysis_rareBkg.C++
root -q analysis_sig.C++

rm BkgPredConfig.txt
echo 'ichannel' $ch  >> BkgPredConfig.txt
echo 'anatype'  $anatype >>  BkgPredConfig.txt
echo 'lowMt'    $lmt >> BkgPredConfig.txt
echo 'highMt'   $hmt >> BkgPredConfig.txt
echo 'lowMET'   $lmet >>BkgPredConfig.txt
echo 'highMET'  $hmet >>BkgPredConfig.txt
echo 'lowPt'    100 >> BkgPredConfig.txt
echo 'highPt'   1000 >> BkgPredConfig.txt
echo 'lepIso'   $iso    >> BkgPredConfig.txt
root -q analysis_VGBkg.C++
root -q analysis_eleBkg.C++
root -q analysis_jetBkg.C++
root -q analysis_qcdBkg.C++
root -q analysis_doubleFakeBkg.C++
root -q analysis_rareBkg.C++
root -q analysis_sig.C++

rm BkgPredConfig.txt
echo 'ichannel' $ch  >> BkgPredConfig.txt
echo 'anatype'  $anatype >>  BkgPredConfig.txt
echo 'lowMt'    $lmt >> BkgPredConfig.txt
echo 'highMt'   $hmt >> BkgPredConfig.txt
echo 'lowMET'   $lmet >>BkgPredConfig.txt
echo 'highMET'  $hmet >>BkgPredConfig.txt
echo 'lowPt'    0 >> BkgPredConfig.txt
echo 'highPt'   1000 >> BkgPredConfig.txt
echo 'lepIso'   $iso    >> BkgPredConfig.txt
root -q analysis_VGBkg.C++
root -q analysis_eleBkg.C++
root -q analysis_jetBkg.C++
root -q analysis_qcdBkg.C++
root -q analysis_doubleFakeBkg.C++
root -q analysis_rareBkg.C++
root -q analysis_sig.C++
