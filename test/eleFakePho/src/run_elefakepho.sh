#!/bin/bash
NPAR=${1:-4}
cd "$(dirname "$0")"
g++ `root-config --cflags` analysis_elefakepho.C -o analysis_elefakepho.exe `root-config --libs` || exit 1
mkdir -p logs_run

JOBS="
2017 '' DY
2016 preVFP DY
2016 postVFP DY
2018 '' DY
2016 B-ver1 Data
2016 C Data
2016 D Data
2016 E Data
2016 F_preVFP Data
2016 F_postVFP Data
2016 G Data
2016 H Data
2017 B Data
2017 C Data
2017 D Data
2017 E Data
2017 F Data
2018 A Data
2018 B Data
2018 C Data
2018 D Data
"
# 2016 B-ver2 Data: waiting for the reprocessed SingleElectron_2016B-ver2 ntuple

echo "$JOBS" | grep -v '^$' | xargs -P ${NPAR} -L 1 bash -c '
  tag="$0_$1_$2"
  ./analysis_elefakepho.exe "$0" "$1" "$2" > logs_run/${tag}.log 2>&1 \
    && echo "done   ${tag}" || echo "FAILED ${tag} (exit $?), see logs_run/${tag}.log"'
echo "All jobs finished"
