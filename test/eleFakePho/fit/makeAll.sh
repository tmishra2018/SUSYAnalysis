#!/bin/bash
set -e

RUN_ALL_BINS=0

RUN_ETA_BINS=0
RUN_PT_BINS=0
RUN_NVTX_BINS=0
RUN_HT_BINS=1

if [[ ${RUN_ALL_BINS} -eq 1 ]]; then
  RUN_ETA_BINS=1
  RUN_PT_BINS=1
  RUN_NVTX_BINS=1
  RUN_HT_BINS=1
fi

#Etabins=(1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2 2.05 2.1 2.15 2.2 2.25 2.3 2.35 2.4 2.45 2.5)

ptbins=(30 35 40 45 50 55 60 65 70 75 80 90 100 120 150 180 10000)
Etabins=(0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2.0 2.05 2.1 2.15 2.2 2.25 2.3 2.35 2.4 2.45 2.5)
#VertexBins=(1 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 10000)
VertexBins=(8 18 28 38 10000)
HTBins=(0 50 100 250 400 500 2000)

rm -f RunnerScript.txt

if [[ ${RUN_ETA_BINS} -eq 1 ]]; then
  for (( i=0; i<${#Etabins[@]}-1; i++ )); do
    j=$((i+1))
    for ch in 1 2 3; do
      for iso in 0 1; do
        echo "${ch} 1 ${iso} ${Etabins[$i]} ${Etabins[$j]} 60 120" >> RunnerScript.txt
      done
    done
  done
fi

if [[ ${RUN_PT_BINS} -eq 1 ]]; then
  for (( i=0; i<${#ptbins[@]}-1; i++ )); do
    j=$((i+1))
    for ch in 1 2 3; do
      for iso in 0 1; do
        echo "${ch} 0 ${iso} ${ptbins[$i]} ${ptbins[$j]} 60 120" >> RunnerScript.txt
      done
    done
  done
fi

if [[ ${RUN_NVTX_BINS} -eq 1 ]]; then
  for (( i=0; i<${#VertexBins[@]}-1; i++ )); do
    j=$((i+1))
    for ch in 1; do # only for Bw-ker
    #for ch in 1 2 3; do
      for iso in 0 1; do
        echo "${ch} 2 ${iso} ${VertexBins[$i]} ${VertexBins[$j]} 60 120" >> RunnerScript.txt
      done
    done
  done
fi

if [[ ${RUN_HT_BINS} -eq 1 ]]; then
  for (( i=0; i<${#HTBins[@]}-1; i++ )); do
    j=$((i+1))
    for ch in 1 2 3; do
      for iso in 0 1; do
        echo "${ch} 3 ${iso} ${HTBins[$i]} ${HTBins[$j]} 60 120" >> RunnerScript.txt
      done
    done
  done
fi


g++ -fPIC -g -Wall -Wno-maybe-uninitialized \
../../../lib/RooDCBShape.o \
../../../lib/RooUserPoly.o \
../../../lib/libAnaClasses.so \
../../../lib/libRooFitClasses.so \
-L$ROOTSYS -lRooFit -lHtml -lMinuit -lRooFitCore -lRooStats -lHistFactory \
`root-config --cflags` FitKer.C -o FitKer.exe \
`root-config --libs` \
-L$ROOTSYS/lib -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d \
-lGpad -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics \
-lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic

cat > configFit_Bw_cmsshape.txt <<EOF
useCMSShape 1
useExpo 0
useKer 0
useDY 0
SaveOutput 1
EOF

cat > configFit_Bw_ker.txt <<EOF
useCMSShape 0
useExpo 0
useKer 1
useDY 0
SaveOutput 1
EOF

cat > configFit_DY_ker.txt <<EOF
useCMSShape 0
useExpo 0
useKer 1
useDY 1
SaveOutput 1
EOF

cat > configFit_Bw_expo.txt <<EOF
useCMSShape 0
useExpo 1
useKer 0
useDY 0
SaveOutput 1
EOF

cat > configFit_cb_poly.txt <<EOF
useCMSShape 0
useExpo 0
useKer 0
useDY 0
usePoly 1
useCrystalBall 1
SaveOutput 1
EOF


#./FitKer.exe 1 0 0 40 45 60 120

cd /uscms/homes/t/tmishra/work/CMSSW_14_0_7/src

tar -czvf files.tar.gz \
SUSYAnalysis/test/eleFakePho/fit/make.sh \
SUSYAnalysis/test/eleFakePho/fit/command_HTcondor.sh \
SUSYAnalysis/test/eleFakePho/fit/command_HTcondor.py \
SUSYAnalysis/test/eleFakePho/fit/FitKer.C \
SUSYAnalysis/lib/* \
SUSYAnalysis/include/*.h \
SUSYAnalysis/src/*.cc \
SUSYAnalysis/include/tdrstyle.C

cp files.tar.gz SUSYAnalysis/test/eleFakePho/fit/
