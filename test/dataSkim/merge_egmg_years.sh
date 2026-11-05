#!/usr/bin/env bash
set -euo pipefail

DIR=${DIR:-/eos/uscms/store/user/tmishra/eg_mg_treesData}
CH_SEL=${1:-all}
YEAR_SEL=${2:-all}
DRYRUN=${DRYRUN:-0}
BACKUP=${BACKUP:-0}

declare -A ERAS=(
  [2016preVFP]="2016B-ver1 2016B-ver2 2016C 2016D 2016E 2016F_preVFP"
  [2016postVFP]="2016F_postVFP 2016G 2016H"
  [2017]="2017B 2017C 2017D 2017E 2017F"
  [2018]="2018A 2018B 2018C 2018D"
)
declare -A CHANNELS=( [eg]="egsignal_DoubleEG" [mg]="mgsignal_MuonEG" )
command -v hadd >/dev/null || { echo "ERROR: hadd not found, run cmsenv first." >&2; exit 1; }

fail=0
for ch in eg mg; do
  [[ "$CH_SEL" == all || "$CH_SEL" == "$ch" ]] || continue
  for year in 2016preVFP 2016postVFP 2017 2018; do
    [[ "$YEAR_SEL" == all || "$YEAR_SEL" == "$year" ]] || continue
    out="$DIR/resTree_${CHANNELS[$ch]}_${year}.root"
    inputs=(); missing=()
    for era in ${ERAS[$year]}; do
      f="$DIR/resTree_${CHANNELS[$ch]}_${era}.root"
      if [[ -s "$f" ]]; then inputs+=("$f"); else missing+=("$era"); fi
    done
    if (( ${#missing[@]} > 0 )); then
      echo "SKIP $ch $year: missing era files: ${missing[*]}" >&2
      fail=1; continue
    fi
    # jobs may run on another node: refuse files written in the last 10 min or not closed cleanly
    recent=$(find "${inputs[@]}" -mmin -10 2>/dev/null || true)
    if [[ -n "$recent" ]]; then
      echo "SKIP $ch $year: still being written (modified <10 min ago): $(echo $recent | xargs -n1 basename | tr '\n' ' ')" >&2
      fail=1; continue
    fi
    bad=()
    for f in "${inputs[@]}"; do
      ok=$(root -l -b -q -e "TFile *f=TFile::Open(\"$f\"); printf(\"OK=%d\\n\", (f && !f->IsZombie() && !f->TestBit(TFile::kRecovered) && f->Get(\"signalTree\")) ? 1 : 0);" 2>/dev/null | grep -o "OK=[01]" || true)
      [[ "$ok" == "OK=1" ]] || bad+=("$(basename "$f")")
    done
    if (( ${#bad[@]} > 0 )); then
      echo "SKIP $ch $year: not closed cleanly / no signalTree: ${bad[*]}" >&2
      fail=1; continue
    fi
    echo "== $ch $year: ${#inputs[@]} eras -> $(basename "$out")"
    if [[ "$DRYRUN" == 1 ]]; then
      echo "hadd -f $out ${inputs[*]}"; continue
    fi
    if [[ "$BACKUP" == 1 && -f "$out" ]]; then
      mv "$out" "${out%.root}_prev.root"
    fi
    hadd -f "$out" "${inputs[@]}" > "$DIR/logs/hadd_${CHANNELS[$ch]}_${year}.log" 2>&1 \
      || { echo "ERROR: hadd failed for $ch $year, see $DIR/logs/hadd_${CHANNELS[$ch]}_${year}.log" >&2; fail=1; continue; }
    # entries of the signal tree, as a quick check
    root -l -b -q -e "TFile f(\"$out\"); TTree *t=(TTree*)f.Get(\"signalTree\"); printf(\"   signalTree entries: %lld\\n\", t ? t->GetEntries() : -1LL);" 2>/dev/null | grep entries || true
  done
done

if (( fail )); then echo "Finished with skipped or failed merges." >&2; exit 1; fi
echo "All merges done."
