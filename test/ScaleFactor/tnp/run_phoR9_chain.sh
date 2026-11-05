#!/bin/bash
# Photon R9 SF production: for each era, photon-probe tree conversion then run_phoR9.sh. Stops the chain at the first failure.
# usage: setsid nohup bash run_phoR9_chain.sh <era> [<era> ...] > $LOGS/phoR9_chain<X>.log 2>&1 &
HERE=$(cd $(dirname $0) && pwd)
LOGS=/uscms_data/d3/tmishra/tnp_work/logs
for TAG in "$@"; do
  echo "[$(date +%H:%M)] convert $TAG"
  bash $HERE/convert_year.sh $TAG pho > $LOGS/convert_pho_$TAG.log 2>&1 || { echo "[$(date +%H:%M)] convert $TAG FAILED"; exit 1; }
  echo "[$(date +%H:%M)] sf $TAG"
  bash $HERE/run_phoR9.sh $TAG > $LOGS/sf_phoR9_$TAG.log 2>&1 || { echo "[$(date +%H:%M)] sf $TAG FAILED"; exit 1; }
  echo "[$(date +%H:%M)] $TAG done"
done
