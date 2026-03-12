#!/bin/bash

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MACRO_PATH="${SCRIPT_DIR}/FitJetFake.C"

ptbins=(30 32 34 36 38 40 45 50 55 60 65 70 80 90 120 150 180 220 300)

echo "[INFO] Precompiling FitJetFake.C once..."
use_precompiled=1
if root -b -q "${MACRO_PATH}+(2018, 0, 1, 30, 32, 1, 0, 1.4442)" >/dev/null 2>&1; then
    echo "[INFO] Precompile done. Launching parallel jobs."
else
    echo "[WARN] Precompile failed. Falling back to per-bin ACLiC mode."
    use_precompiled=0
fi

run_one_combo() {
    local RunYear="$1"
    local ISpreVFP="$2"
    local eventType="$3"

    echo "[START] Year=${RunYear}, ISpreVFP=${ISpreVFP}, eventType=${eventType}"

    local failures=0
    for (( i = 0; i < ${#ptbins[@]} - 1; i++ )); do
        local j=$((i + 1))
        local ptmin=${ptbins[$i]}
        local ptmax=${ptbins[$j]}

        echo "Running for Year=${RunYear}, ISpreVFP=${ISpreVFP}, eventType=${eventType}, pt=${ptmin}-${ptmax}"
        if [ "$use_precompiled" -eq 1 ]; then
            run_cmd="${MACRO_PATH}(${RunYear}, ${ISpreVFP}, ${eventType}, ${ptmin}, ${ptmax}, 1, 0, 1.4442)"
        else
            run_cmd="${MACRO_PATH}+(${RunYear}, ${ISpreVFP}, ${eventType}, ${ptmin}, ${ptmax}, 1, 0, 1.4442)"
        fi

        if ! root -b -q "$run_cmd"; then
            echo "[WARN] Failed for Year=${RunYear}, ISpreVFP=${ISpreVFP}, eventType=${eventType}, pt=${ptmin}-${ptmax}"
            failures=$((failures + 1))
        fi
    done

    if [ "$failures" -gt 0 ]; then
        echo "[DONE ] Year=${RunYear}, ISpreVFP=${ISpreVFP}, eventType=${eventType} with ${failures} failed bins"
    else
        echo "[DONE ] Year=${RunYear}, ISpreVFP=${ISpreVFP}, eventType=${eventType} (all bins done)"
    fi
}

job_count=0
for RunYear in 2016 2017 2018; do
    if [ "$RunYear" -eq 2016 ]; then
        ISpreVFP_list=(0 1)
    else
        ISpreVFP_list=(0)
    fi

    for ISpreVFP in "${ISpreVFP_list[@]}"; do
        for eventType in 1 2; do
            run_one_combo "$RunYear" "$ISpreVFP" "$eventType" &
            job_count=$((job_count + 1))
        done
    done
done

echo "Launched ${job_count} parallel jobs (expected: 8)."
wait
