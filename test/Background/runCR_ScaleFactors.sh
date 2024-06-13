#!/bin/bash
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

anatype=0 # control region 40 < MET < 70
lmt=0
hmt=-1
lmet=40
hmet=70

RunYear=2018
preVFP=0
ch=1
iso=4

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


#anatype=2 # validation region 0 < MT < 100
#lmt=0
#hmt=100
#lmet=0
#hmet=-1

#RunYear=2016
#preVFP=0
#
#ch=1
#iso=4
#
#rm BkgPredConfig.txt
#echo 'ichannel' $ch  >> BkgPredConfig.txt
#echo 'anatype'  $anatype >>  BkgPredConfig.txt
#echo 'lowMt'    $lmt >> BkgPredConfig.txt
#echo 'highMt'   $hmt >> BkgPredConfig.txt
#echo 'lowMET'   $lmet >>BkgPredConfig.txt
#echo 'highMET'  $hmet >>BkgPredConfig.txt
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
#RunYear=2016
#preVFP=0
#
#ch=2
#iso=4
#
#rm BkgPredConfig.txt
#echo 'ichannel' $ch  >> BkgPredConfig.txt
#echo 'anatype'  $anatype >>  BkgPredConfig.txt
#echo 'lowMt'    $lmt >> BkgPredConfig.txt
#echo 'highMt'   $hmt >> BkgPredConfig.txt
#echo 'lowMET'   $lmet >>BkgPredConfig.txt
#echo 'highMET'  $hmet >>BkgPredConfig.txt
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


#########################################################






#different lepton pT ranges : 0-50-70-100-1000
#rm BkgPredConfig.txt
#echo 'ichannel' $ch  >> BkgPredConfig.txt
#echo 'anatype'  $anatype >>  BkgPredConfig.txt
#echo 'lowMt'    $lmt >> BkgPredConfig.txt
#echo 'highMt'   $hmt >> BkgPredConfig.txt
#echo 'lowMET'   $lmet >>BkgPredConfig.txt
#echo 'highMET'  $hmet >>BkgPredConfig.txt
#echo 'lowPt'    0 >> BkgPredConfig.txt
#echo 'highPt'   50 >> BkgPredConfig.txt
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
#
#rm BkgPredConfig.txt
#echo 'ichannel' $ch  >> BkgPredConfig.txt
#echo 'anatype'  $anatype >>  BkgPredConfig.txt
#echo 'lowMt'    $lmt >> BkgPredConfig.txt
#echo 'highMt'   $hmt >> BkgPredConfig.txt
#echo 'lowMET'   $lmet >>BkgPredConfig.txt
#echo 'highMET'  $hmet >>BkgPredConfig.txt
#echo 'lowPt'    50 >> BkgPredConfig.txt
#echo 'highPt'   70 >> BkgPredConfig.txt
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
#rm BkgPredConfig.txt
#echo 'ichannel' $ch  >> BkgPredConfig.txt
#echo 'anatype'  $anatype >>  BkgPredConfig.txt
#echo 'lowMt'    $lmt >> BkgPredConfig.txt
#echo 'highMt'   $hmt >> BkgPredConfig.txt
#echo 'lowMET'   $lmet >>BkgPredConfig.txt
#echo 'highMET'  $hmet >>BkgPredConfig.txt
#echo 'lowPt'    70 >> BkgPredConfig.txt
#echo 'highPt'   100 >> BkgPredConfig.txt
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
#rm BkgPredConfig.txt
#echo 'ichannel' $ch  >> BkgPredConfig.txt
#echo 'anatype'  $anatype >>  BkgPredConfig.txt
#echo 'lowMt'    $lmt >> BkgPredConfig.txt
#echo 'highMt'   $hmt >> BkgPredConfig.txt
#echo 'lowMET'   $lmet >>BkgPredConfig.txt
#echo 'highMET'  $hmet >>BkgPredConfig.txt
#echo 'lowPt'    100 >> BkgPredConfig.txt
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
