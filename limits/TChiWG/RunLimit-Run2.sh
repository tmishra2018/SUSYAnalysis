#!/bin/bash

base_cards="/eos/uscms/store/user/tmishra/CombinedLimit/TChiWG/cards"
cardpath_exp="${base_cards}/cards_Run2_splusb"
cardpath_obs="${base_cards}/cards_Run2_unblind"
limitpath="/eos/uscms/store/user/tmishra/CombinedLimit/TChiWG/limits/limits_Run2"

mkdir -p ${limitpath}
rm -f ${limitpath}/*
PARALLEL_JOBS=4

# $1 = cardpath, $2 = tag ("expected" or "observed"), $3 = extra combine args
run_batch() {
    local cardpath=$1
    local tag=$2
    local extra_args=$3
    local PIDS=()

    for file in ${cardpath}/counting_tchiwg_*.txt; do
        while [ "${#PIDS[@]}" -ge "${PARALLEL_JOBS}" ]; do
            PIDS=($(for pid in "${PIDS[@]}"; do kill -0 "$pid" 2>/dev/null && echo "$pid"; done))
            sleep 3
        done

        filename=$(basename "$file")
        chargino=$(echo "$filename" | cut -d'_' -f3 | cut -d'.' -f1)

        echo "[${tag}] Processing chargino=$chargino"
        combine -M AsymptoticLimits "$file" ${extra_args} -n "tchiwg${tag}" \
            | tee "${limitpath}/limit_${tag}_${chargino}.txt" &
        PIDS+=($!)
    done
    wait
}

#echo "=== EXPECTED limits (s+b pseudo-data, cards_Run2_splusb, no --run blind) ==="
#run_batch "${cardpath_exp}" "expected" ""

echo "=== OBSERVED limits (real unblinded data, cards_Run2_unblind) ==="
run_batch "${cardpath_obs}" "observed" ""

echo "All jobs completed."
