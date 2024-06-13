#!/bin/bash
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

anatype=1 # control region 0 < MET < 70
lmt=0
hmt=-1
lmet=40
hmet=70
iso=4

##################################################################

RunYear=2017
preVFP=0


ch=1

rm BkgPredConfig.txt
echo 'ichannel' $ch  >> BkgPredConfig.txt
echo 'anatype'  $anatype >>  BkgPredConfig.txt
echo 'lowMt'    $lmt >> BkgPredConfig.txt
echo 'highMt'   $hmt >> BkgPredConfig.txt
echo 'lowMET'   $lmet >> BkgPredConfig.txt
echo 'highMET'  $hmet >> BkgPredConfig.txt
echo 'lowPt'    0 >> BkgPredConfig.txt
echo 'highPt'   1000 >> BkgPredConfig.txt
echo 'lepIso'   $iso    >> BkgPredConfig.txt
echo 'RunYear'   $RunYear    >> BkgPredConfig.txt
echo 'preVFP'   $preVFP    >> BkgPredConfig.txt
root -l -q analysis_VGBkg.C++
root -l -q analysis_eleBkg.C++
root -l -q analysis_jetBkg.C++
root -l -q analysis_qcdBkg.C++
root -l -q analysis_rareBkg.C++
root -l -q analysis_sig.C++

ch=2

rm BkgPredConfig.txt
echo 'ichannel' $ch  >> BkgPredConfig.txt
echo 'anatype'  $anatype >>  BkgPredConfig.txt
echo 'lowMt'    $lmt >> BkgPredConfig.txt
echo 'highMt'   $hmt >> BkgPredConfig.txt
echo 'lowMET'   $lmet >> BkgPredConfig.txt
echo 'highMET'  $hmet >> BkgPredConfig.txt
echo 'lowPt'    0 >> BkgPredConfig.txt
echo 'highPt'   1000 >> BkgPredConfig.txt
echo 'lepIso'   $iso    >> BkgPredConfig.txt
echo 'RunYear'   $RunYear    >> BkgPredConfig.txt
echo 'preVFP'   $preVFP    >> BkgPredConfig.txt
root -l -q analysis_VGBkg.C++
root -l -q analysis_eleBkg.C++
root -l -q analysis_jetBkg.C++
root -l -q analysis_qcdBkg.C++
root -l -q analysis_rareBkg.C++
root -l -q analysis_sig.C++

##################################################################

rm temp.txt

root -l -q "plot_bkg.C(1,2017,1)" >> temp.txt
root -l -q "plot_bkg.C(2,2017,1)" >> temp.txt
