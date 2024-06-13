#!/bin/bash
#ptbins=( 35 50 100 150 200 300 400 600 800 10000)
ptbins=( 35 40 50 60 70 80 100 150 200 250 300 800)
for (( i = 0 ; i < ${#ptbins[@]} ; i++ ))
do
	j=$((i+1))
  #root -b -q "FitJetFake.C+(${ptbins[$i]},${ptbins[$j]},2)"
  root -b -q "FitJetFake.C+(${ptbins[$i]},${ptbins[$j]},2)"
done
