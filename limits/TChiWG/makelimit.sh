#!/bin/bash

path="/eos/uscms/store/user/tmishra/CombinedLimit/TChiWG/limits/limits_Run2"

expfile="limitvalue_Run2_expected.txt"
obsfile="limitvalue_Run2_observed.txt"

rm -f "${expfile}" "${obsfile}"

n_exp_ok=0
n_exp_missing=0
n_obs_ok=0
n_obs_missing=0

for i in {0..39}
do
    j=$((300 + i * 25))

    expfilepath="${path}/limit_expected_${j}.txt"
    obsfilepath="${path}/limit_observed_${j}.txt"

    # --- Expected ---
    if [[ -f "${expfilepath}" ]]; then
        value=$(awk -F '<' '/50.0%: r </ {print $2}' "${expfilepath}")
        sigma1up=$(awk -F '<' '/84.0%: r </ {print $2}' "${expfilepath}")
        sigma1lo=$(awk -F '<' '/16.0%: r </ {print $2}' "${expfilepath}")
        sigma2up=$(awk -F '<' '/97.5%: r </ {print $2}' "${expfilepath}")
        sigma2lo=$(awk -F '<' '/2.5%: r </ {print $2}' "${expfilepath}")
        if [[ -n "${value}" ]]; then
            echo "${j} ${value} ${sigma1up} ${sigma1lo} ${sigma2up} ${sigma2lo}" >> "${expfile}"
            n_exp_ok=$((n_exp_ok + 1))
        else
            echo "WARNING: no expected quantiles parsed from ${expfilepath} -- skipping mass ${j}" >&2
            n_exp_missing=$((n_exp_missing + 1))
        fi
    else
        echo "WARNING: ${expfilepath} not found -- skipping expected for mass ${j}" >&2
        n_exp_missing=$((n_exp_missing + 1))
    fi

    # --- Observed ---
    if [[ -f "${obsfilepath}" ]]; then
        obsvalue=$(awk -F '<' '/Observed Limit: r </ {print $2}' "${obsfilepath}")
        if [[ -n "${obsvalue}" ]]; then
            echo "${j} ${obsvalue}" >> "${obsfile}"
            n_obs_ok=$((n_obs_ok + 1))
        else
            echo "WARNING: no 'Observed Limit' line parsed from ${obsfilepath} -- was it accidentally run with --run blind?" >&2
            n_obs_missing=$((n_obs_missing + 1))
        fi
    else
        echo "WARNING: ${obsfilepath} not found -- skipping observed for mass ${j}" >&2
        n_obs_missing=$((n_obs_missing + 1))
    fi
done

echo "Expected: wrote ${n_exp_ok} mass point(s) to ${expfile} (${n_exp_missing} missing/skipped)."
echo "Observed: wrote ${n_obs_ok} mass point(s) to ${obsfile} (${n_obs_missing} missing/skipped)."
