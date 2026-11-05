#!/bin/bash
# Same fix as run_elefakepho.sh: starting all 18 jobs at once runs an LPC interactive node out of memory.
# Run at most NPAR at a time (default 4), one log per job, report failures, then hadd the eras of each
# period into plot_bgtemplate_FullEcal_<period>.root (read by ../fit) if all its eras succeeded.
# Usage: bash run_bgTemplate.sh [NPAR]
NPAR=${1:-4}
cd "$(dirname "$0")"
g++ `root-config --cflags` analysis_bgtemplate.C -o analysis_bgtemplate.exe `root-config --libs` || exit 1
mkdir -p logs_run_bg
rm -f logs_run_bg/*.ok

JOBS="
2017 B
2017 C
2017 D
2017 E
2017 F
2018 A
2018 B
2018 C
2018 D
2016 B-ver1
2016 B-ver2
2016 C
2016 D
2016 E
2016 F_preVFP
2016 F_postVFP
2016 G
2016 H
"

echo "$JOBS" | grep -v '^$' | xargs -P ${NPAR} -L 1 bash -c '
  tag="$0$1"
  ./analysis_bgtemplate.exe "$0" "$1" > logs_run_bg/${tag}.log 2>&1 \
    && grep -q "BGTree events" /eos/uscms/store/user/tmishra/elefakepho/logs/plot_bgtemplate_FullEcal_${tag}.log \
    && { touch logs_run_bg/${tag}.ok; echo "done   ${tag}"; } \
    || echo "FAILED ${tag}, see logs_run_bg/${tag}.log"'

D=/eos/uscms/store/user/tmishra/elefakepho/files
merge() {
  local period=$1; shift
  local inputs=()
  for era in "$@"; do
    [ -f logs_run_bg/${era}.ok ] || { echo "NOT merging ${period}: ${era} failed"; return; }
    inputs+=("$D/plot_bgtemplate_FullEcal_${era}.root")
  done
  hadd -f $D/plot_bgtemplate_FullEcal_${period}.root "${inputs[@]}" > logs_run_bg/hadd_${period}.log 2>&1 \
    && echo "merged ${period}" || echo "FAILED hadd ${period}, see logs_run_bg/hadd_${period}.log"
}
merge 2016preVFP  2016B-ver1 2016B-ver2 2016C 2016D 2016E 2016F_preVFP
merge 2016postVFP 2016F_postVFP 2016G 2016H
merge 2017 2017B 2017C 2017D 2017E 2017F
merge 2018 2018A 2018B 2018C 2018D
echo "All jobs finished"
