#!/bin/bash
# detType == 1 for EB and detType == 2 for EE
for (( i = 1 ; i < 2 ; i++ ))
do
	j=$((i+1))
 # root -b -q "newFit.C+($i)"
  root -l -q "fitJetFunc.C+($i)"
  #root -l -q "fitJetFuncMC.C+($i)"
done
