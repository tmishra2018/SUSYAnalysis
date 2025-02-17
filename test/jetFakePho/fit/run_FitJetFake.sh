#!/bin/bash
  # FitJetFake(float lowercut, float uppercut, int detType, float loweta, float higheta) 
#ptbins=( 30 32 34 36 38 40 45 50 55 60 65 70 80 90 120 150 180 300)
#ptbins=( 30 35 40 45 50 55 60 75 90 130 170 210 300)



ptbins=(30 35 40 45 50 55 60 90 170 300)

for (( i = 0 ; i < ${#ptbins[@]}-1 ; i++ ))
do
	j=$((i+1))
  root -b -q "FitJetFake.C+(${ptbins[$i]},${ptbins[$j]},1,0,1.4442)"

done
