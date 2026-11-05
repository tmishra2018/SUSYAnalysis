#!/bin/bash
# e-gamma trigger SFs from Trigger/files/plot_egTrigger_{Data,DY}_<period>.root
#   -> Trigger/files/diphoton_pholeg_<period>.root (photon leg) and diphoton_eleg_<period>.root (electron leg)
# Run from SUSYAnalysis/test/Trigger after cmsenv (after the Run-egTrig.sh ntuples).
cd "$(dirname "$0")"
g++ `root-config --cflags` plotTrigger.C -o plotTrigger.exe `root-config --libs` || exit 1
for period in "2016 1" "2016 0" "2017 0" "2018 0"; do
  for leading in 1 0; do
    ./plotTrigger.exe ${period} ${leading} > /dev/null || echo "FAILED plotTrigger ${period} ${leading}"
  done
done
