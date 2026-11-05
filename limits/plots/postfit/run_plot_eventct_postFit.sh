#!/bin/bash
set -e
cd "$(dirname "$0")"

POSTFITDIR="$(cd ../../fits/postfit_unblind_T5WG && pwd)"   # written by limits/fits/run_postfit.sh
NBIN=18
YEAR_ORDER="2016preVFP 2016postVFP 2017 2018 Run2"

hadd -f /uscms_data/d3/tmishra/Output/signalTree_egamma_signal_678.root /uscms_data/d3/tmishra/Output/signalTree_egamma_signal_201*.root
hadd -f /uscms_data/d3/tmishra/Output/signalTree_mg_signal_678.root /uscms_data/d3/tmishra/Output/signalTree_mg_signal_201*.root

hadd -f /uscms/home/tmishra/nobackup/signal_trees/signalTree_T5WG_678.root /uscms/home/tmishra/nobackup/signal_trees/signalTree_T5WG_20*.root
hadd -f /uscms/home/tmishra/nobackup/signal_trees/signalTree_T6WG_678.root /uscms/home/tmishra/nobackup/signal_trees/signalTree_T6WG_20*.root
hadd -f /uscms/home/tmishra/nobackup/signal_trees/signalTree_TChiWG_678.root /uscms/home/tmishra/nobackup/signal_trees/signalTree_TChiWG_20*.root

mkdir -p logs   
for YEAR in ${YEAR_ORDER}; do
  echo ""
  echo "=============================================="
  echo " plot_eventct: ${YEAR}"
  echo "=============================================="
  root -l -b -q "plot_eventct.C(${NBIN}, \"${POSTFITDIR}\", \"${YEAR}\")"
done

echo ""
echo "=== ALL DONE ==="
