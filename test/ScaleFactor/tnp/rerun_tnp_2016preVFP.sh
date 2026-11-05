#!/bin/bash
# Step 8: redo the 2016preVFP R9/miniIso electron SFs with the new SingleElectron 2016B-ver2 and 2016C files.
# The old work dir (run_sf.sh would reuse its cached histograms) and SF results are kept as *_old20261019.
# usage: setsid nohup bash rerun_tnp_2016preVFP.sh > /eos/uscms/store/user/tmishra/TnP_ggNtuple/rerun_2016preVFP.log 2>&1 &
set -e
HERE=$(cd $(dirname $0) && pwd)
WORK=/uscms_data/d3/tmishra/tnp_work/2016preVFP
RES=$HERE/results/2016preVFP
O=/eos/uscms/store/user/tmishra/TnP_ggNtuple/2016preVFP

[ -e ${WORK}_old20261019 ] || mv $WORK ${WORK}_old20261019
[ -e ${RES}_old20261019 ]  || cp -a $RES ${RES}_old20261019

echo "[$(date '+%F %T')] convert_year.sh 2016preVFP on $(hostname)"
touch $O/.start
bash $HERE/convert_year.sh 2016preVFP
for f in $O/tnp_*.root; do
  [ $f -nt $O/.start ] && tail -1 $f.log | grep -q "^Wrote" || { echo "conversion failed: $f"; exit 1; }
done

echo "[$(date '+%F %T')] run_sf.sh 2016preVFP"
bash $HERE/run_sf.sh 2016preVFP
echo "[$(date '+%F %T')] done; SFs in $RES (old ones in ${RES}_old20261019)"
