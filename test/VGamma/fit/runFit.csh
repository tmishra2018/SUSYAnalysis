#!/bin/bash
export RunYear=2017
#export VFP="preVFP"
#export VFP="postVFP"
export VFP=""

rm /eos/uscms/store/user/tmishra/VGamma/VGamma_scalefactor_mg_${RunYear}${VFP}.txt
for i in {0..50}
do
  root -b -q "Fitfractionmg.C+($i,40,70,0,1000,4)"
done

#for i in {0..49}
#do
#  root -b -q "Fitfractionmg.C+($i,40,70,0,50,4)"
#done
#for i in {0..49}
#do
#  root -b -q "Fitfractionmg.C+($i,40,70,50,70,4)"
#done
#for i in {0..49}
#do
#  root -b -q "Fitfractionmg.C+($i,40,70,70,100,4)"
#done
#for i in {0..49}
#do
#  root -b -q "Fitfractionmg.C+($i,40,70,100,1000,4)"
#done
