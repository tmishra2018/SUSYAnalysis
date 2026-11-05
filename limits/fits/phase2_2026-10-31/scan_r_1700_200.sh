#!/bin/bash
# 1D likelihood scan of r for T5WG 1700/200, Run 2, observed data (scan_r_t5wg_1700_200_Run2.pdf).
set -e
C=/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_Run2_unblind/counting_t5wg_1700_200.txt
W=scan_t5wg_1700_200_Run2; mkdir -p $W; cd $W; cp $C .
text2workspace.py counting_t5wg_1700_200.txt -o ws.root -m 125
combine -M MultiDimFit ws.root -m 125 --algo grid --points 60 --rMin 0 --rMax 3 -n _t5wg_1700_200_Run2 --robustFit 1 --cminDefaultMinimizerStrategy 0
plot1DScan.py higgsCombine_t5wg_1700_200_Run2.MultiDimFit.mH125.root -o scan_r_t5wg_1700_200_Run2 --POI r
