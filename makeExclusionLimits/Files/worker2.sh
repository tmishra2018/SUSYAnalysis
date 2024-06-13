#!/bin/sh

#analyzeLightBSM 800 127 FastSim T5qqqqHg
gluinoMass=$1
nlspMass=$2
outName=$3
ext=$4
v_12=false
hist1=METv25Bin
currDir=$(pwd)
hist=AllSBins_v7_CD
outRootFile="${outName}_${gluinoMass}_${nlspMass}_${ext}"

echo "========== HiggsAnalysis/CombinedLimit ==============="
echo "${pwd}"


echo "making datacards"
mkdir dataCards
mkdir dataCards/${outRootFile}_${hist1}
# for hist in METvBin1_nj2to4_nbjet0 METvBin2_nj5to6_nbjet0 METvBin2_nj7toinf_nbjet0 METvBin2_nj2to4_nbjetnot0 METvBin2_nj5to6_nbjetnot0 METvBin2_nj7toinf_nbjetnot0
# do
echo "mkdir dataCards/${outRootFile}_${hist}"
mkdir dataCards/${outRootFile}_${hist}
echo "ls dataCards"
ls dataCards
# root -l -q -b 'makeDatacard_Simple_bybhumi.C('${gluinoMass}','${nlspMass}',"'${outRootFile}'.root","'${hist}'",'$v_12')'
root -l -q -b 'makeDatacard_SBins.C('${gluinoMass}','${nlspMass}',"'${outRootFile}'.root")'
echo "ls dataCards"
ls dataCards
# done
# rm ${outRootFile}.root


cd CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit
cmsenv

# for hist in METvBin1_nj2to4_nbjet0 METvBin2_nj5to6_nbjet0 METvBin2_nj7toinf_nbjet0 METvBin2_nj2to4_nbjetnot0 METvBin2_nj5to6_nbjetnot0 METvBin2_nj7toinf_nbjetnot0
# do
# cp ${currDir}/dataCards/${outRootFile}_${hist}/*.txt  ${currDir}/dataCards/${outRootFile}_${hist1}/.
# done


echo "combineCards.py ${currDir}/dataCards/${outRootFile}_${hist}/*.txt > dataCard_${outName}_${gluinoMass}_${nlspMass}_${hist}.txt"
combineCards.py ${currDir}/dataCards/${outRootFile}_${hist}/*.txt > dataCard_${outName}_${gluinoMass}_${nlspMass}_${hist}.txt
echo "combine -M AsymptoticLimits dataCard_${outName}_${gluinoMass}_${nlspMass}_${hist}.txt -t -1 -n ${outName}_${gluinoMass}_${nlspMass}_${hist} -m ${mH}"
combine -M AsymptoticLimits dataCard_${outName}_${gluinoMass}_${nlspMass}_${hist}.txt -t -1 


rm *.txt
cd ${currDir}
rm -rf dataCards
