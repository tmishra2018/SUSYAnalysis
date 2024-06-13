#!/bin/bash
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

anatype=2 # validation region 0 < MT < 100
lmt=0
hmt=100
lmet=0
hmet=-1
iso=4


RunYear=2016
preVFP=1


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



############################################################


#RunYear=2016
#preVFP=0
#
#
#ch=1
#
#rm BkgPredConfig.txt
#echo 'ichannel' $ch  >> BkgPredConfig.txt
#echo 'anatype'  $anatype >>  BkgPredConfig.txt
#echo 'lowMt'    $lmt >> BkgPredConfig.txt
#echo 'highMt'   $hmt >> BkgPredConfig.txt
#echo 'lowMET'   $lmet >> BkgPredConfig.txt
#echo 'highMET'  $hmet >> BkgPredConfig.txt
#echo 'lowPt'    0 >> BkgPredConfig.txt
#echo 'highPt'   1000 >> BkgPredConfig.txt
#echo 'lepIso'   $iso    >> BkgPredConfig.txt
#echo 'RunYear'   $RunYear    >> BkgPredConfig.txt
#echo 'preVFP'   $preVFP    >> BkgPredConfig.txt
#root -l -q analysis_VGBkg.C++
#root -l -q analysis_eleBkg.C++
#root -l -q analysis_jetBkg.C++
#root -l -q analysis_qcdBkg.C++
#root -l -q analysis_rareBkg.C++
#root -l -q analysis_sig.C++
#
#ch=2
#
#rm BkgPredConfig.txt
#echo 'ichannel' $ch  >> BkgPredConfig.txt
#echo 'anatype'  $anatype >>  BkgPredConfig.txt
#echo 'lowMt'    $lmt >> BkgPredConfig.txt
#echo 'highMt'   $hmt >> BkgPredConfig.txt
#echo 'lowMET'   $lmet >> BkgPredConfig.txt
#echo 'highMET'  $hmet >> BkgPredConfig.txt
#echo 'lowPt'    0 >> BkgPredConfig.txt
#echo 'highPt'   1000 >> BkgPredConfig.txt
#echo 'lepIso'   $iso    >> BkgPredConfig.txt
#echo 'RunYear'   $RunYear    >> BkgPredConfig.txt
#echo 'preVFP'   $preVFP    >> BkgPredConfig.txt
#root -l -q analysis_VGBkg.C++
#root -l -q analysis_eleBkg.C++
#root -l -q analysis_jetBkg.C++
#root -l -q analysis_qcdBkg.C++
#root -l -q analysis_rareBkg.C++
#root -l -q analysis_sig.C++



############################################################


#RunYear=2017
#preVFP=0
#
#
#ch=1
#
#rm BkgPredConfig.txt
#echo 'ichannel' $ch  >> BkgPredConfig.txt
#echo 'anatype'  $anatype >>  BkgPredConfig.txt
#echo 'lowMt'    $lmt >> BkgPredConfig.txt
#echo 'highMt'   $hmt >> BkgPredConfig.txt
#echo 'lowMET'   $lmet >> BkgPredConfig.txt
#echo 'highMET'  $hmet >> BkgPredConfig.txt
#echo 'lowPt'    0 >> BkgPredConfig.txt
#echo 'highPt'   1000 >> BkgPredConfig.txt
#echo 'lepIso'   $iso    >> BkgPredConfig.txt
#echo 'RunYear'   $RunYear    >> BkgPredConfig.txt
#echo 'preVFP'   $preVFP    >> BkgPredConfig.txt
#root -l -q analysis_VGBkg.C++
#root -l -q analysis_eleBkg.C++
#root -l -q analysis_jetBkg.C++
#root -l -q analysis_qcdBkg.C++
#root -l -q analysis_rareBkg.C++
#root -l -q analysis_sig.C++
#
#ch=2
#
#rm BkgPredConfig.txt
#echo 'ichannel' $ch  >> BkgPredConfig.txt
#echo 'anatype'  $anatype >>  BkgPredConfig.txt
#echo 'lowMt'    $lmt >> BkgPredConfig.txt
#echo 'highMt'   $hmt >> BkgPredConfig.txt
#echo 'lowMET'   $lmet >> BkgPredConfig.txt
#echo 'highMET'  $hmet >> BkgPredConfig.txt
#echo 'lowPt'    0 >> BkgPredConfig.txt
#echo 'highPt'   1000 >> BkgPredConfig.txt
#echo 'lepIso'   $iso    >> BkgPredConfig.txt
#echo 'RunYear'   $RunYear    >> BkgPredConfig.txt
#echo 'preVFP'   $preVFP    >> BkgPredConfig.txt
#root -l -q analysis_VGBkg.C++
#root -l -q analysis_eleBkg.C++
#root -l -q analysis_jetBkg.C++
#root -l -q analysis_qcdBkg.C++
#root -l -q analysis_rareBkg.C++
#root -l -q analysis_sig.C++



############################################################

#RunYear=2018
#preVFP=0
#
#
#ch=1
#
#rm BkgPredConfig.txt
#echo 'ichannel' $ch  >> BkgPredConfig.txt
#echo 'anatype'  $anatype >>  BkgPredConfig.txt
#echo 'lowMt'    $lmt >> BkgPredConfig.txt
#echo 'highMt'   $hmt >> BkgPredConfig.txt
#echo 'lowMET'   $lmet >> BkgPredConfig.txt
#echo 'highMET'  $hmet >> BkgPredConfig.txt
#echo 'lowPt'    0 >> BkgPredConfig.txt
#echo 'highPt'   1000 >> BkgPredConfig.txt
#echo 'lepIso'   $iso    >> BkgPredConfig.txt
#echo 'RunYear'   $RunYear    >> BkgPredConfig.txt
#echo 'preVFP'   $preVFP    >> BkgPredConfig.txt
#root -l -q analysis_VGBkg.C++
#root -l -q analysis_eleBkg.C++
#root -l -q analysis_jetBkg.C++
#root -l -q analysis_qcdBkg.C++
#root -l -q analysis_rareBkg.C++
#root -l -q analysis_sig.C++
#
#ch=2
#
#rm BkgPredConfig.txt
#echo 'ichannel' $ch  >> BkgPredConfig.txt
#echo 'anatype'  $anatype >>  BkgPredConfig.txt
#echo 'lowMt'    $lmt >> BkgPredConfig.txt
#echo 'highMt'   $hmt >> BkgPredConfig.txt
#echo 'lowMET'   $lmet >> BkgPredConfig.txt
#echo 'highMET'  $hmet >> BkgPredConfig.txt
#echo 'lowPt'    0 >> BkgPredConfig.txt
#echo 'highPt'   1000 >> BkgPredConfig.txt
#echo 'lepIso'   $iso    >> BkgPredConfig.txt
#echo 'RunYear'   $RunYear    >> BkgPredConfig.txt
#echo 'preVFP'   $preVFP    >> BkgPredConfig.txt
#root -l -q analysis_VGBkg.C++
#root -l -q analysis_eleBkg.C++
#root -l -q analysis_jetBkg.C++
#root -l -q analysis_qcdBkg.C++
#root -l -q analysis_rareBkg.C++
#root -l -q analysis_sig.C++

############################################################
