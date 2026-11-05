#!/bin/bash
# =============================================================================
# runImpacts.sh
# Generate impact plots for T5WG representative mass points in three modes:
#
#   bonly   -- Asimov dataset from b-only model   (--expectSignal 0)
#              reads datacards from cards_<YEAR>_bonly/
#   splusb  -- Asimov dataset from s+b model      (--expectSignal 1)
#              reads datacards from cards_<YEAR>_splusb/
#   unblind -- Real observed data                 (no --expectSignal flag)
#              reads datacards from cards_<YEAR>_unblind/
#
# Workspace is built once per mass point per mode and reused.
# Output PDFs are named  plot_impacts_t5wg_<mass>_<mode>.pdf
#
# Parallelism:
#   PARALLEL_JOBS  -- number of (mass x mode) jobs to run simultaneously
#   PARALLEL       -- number of NP scan threads inside each combineTool call
# =============================================================================
set +e   # Do not abort on non-zero exit; failed NP fits are non-fatal

CARD_BASE="/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards"
MASS_POINTS=("1700_1100")   # representative points
#MASS_POINTS=("1700_1100" "1300_500" "2100_100")   # representative points
#YEARS=("Run2")
YEARS=("2016preVFP" "2016postVFP" "2017" "2018" "Run2")
MODES=("bonly" "splusb")
#MODES=("bonly" "splusb" "unblind")

# Outer parallelism: simultaneous (mass x mode) jobs
PARALLEL_JOBS=2   # 2 jobs x 4 NP threads = 8 cores fully utilized

# Inner parallelism: NP scan threads per combineTool call
PARALLEL=4

# Fixed mass for workspace/combine (-m flag)
FIXED_MASS=125

# Minimizer options common to all combine calls
# --rMin/rMax widened to prevent "closed range" crossing errors
# --cminDefaultMinimizerStrategy 0  faster/less strict minimization
# --X-rtd FITTER_NEW_CROSSING_ALGO  fixes "Closed range without finding crossing"
# --X-rtd FITTER_NEVER_GIVE_UP      retries failed fits instead of dying
# --X-rtd MINIMIZER_analytic        analytic gradients for faster convergence
MINOPTS="--rMin -20 \
         --rMax 20 \
         --robustFit 1 \
         --cminDefaultMinimizerStrategy 0 \
         --X-rtd FITTER_NEW_CROSSING_ALGO \
         --X-rtd FITTER_NEVER_GIVE_UP \
         --X-rtd MINIMIZER_analytic"

# Directory to collect all logs/status from parallel jobs
LOGDIR="parallel_logs"
mkdir -p "${LOGDIR}"

# =============================================================================
run_impacts() {
    local mass=$1
    local mode=$2
    local workspace=$3

    local tag="t5wg_${mass}_${mode}"
    echo ""
    echo "--- Impacts: ${mode} ---"

    # Build the mode-specific --expectSignal flag
    local expect_signal_flag=""
    local expect_signal_val=0
    if [ "${mode}" = "bonly" ]; then
        expect_signal_flag=""          # no --expectSignal (defaults to 0 / b-only)
        expect_signal_val=0
    elif [ "${mode}" = "splusb" ]; then
        expect_signal_flag="--expectSignal 1"
        expect_signal_val=1
    elif [ "${mode}" = "unblind" ]; then
        expect_signal_flag=""          # real observed data, no Asimov injection
        expect_signal_val=0
    fi

    # Step 1: Initial fit
    echo ">>> Step 1: Initial fit [${tag}]"
    combineTool.py -M Impacts \
        -d ${workspace} \
        -m ${FIXED_MASS} \
        ${expect_signal_flag} \
        --doInitialFit \
        ${MINOPTS} \
        2>&1 | tee log_initialFit_${tag}.txt

    # Check for convergence - warn but do not abort
    if grep -q "ERROR\|Closed range\|not have converged" log_initialFit_${tag}.txt 2>/dev/null; then
        echo "Warning: initial fit issues detected for ${tag} - proceeding anyway"
    else
        echo "Initial fit converged for ${tag}"
    fi

    # Step 2: NP scans
    echo ">>> Step 2: NP scans [${tag}]"
    combineTool.py -M Impacts \
        -d ${workspace} \
        -m ${FIXED_MASS} \
        ${expect_signal_flag} \
        --doFits \
        --job-mode interactive \
        --parallel ${PARALLEL} \
        ${MINOPTS} \
        2>&1 | tee log_fits_${tag}.txt

    # Step 3: Collect results into JSON
    echo ">>> Step 3: Collecting results [${tag}]"
    combineTool.py -M Impacts \
        -d ${workspace} \
        -m ${FIXED_MASS} \
        -o impacts_${tag}.json \
        2>&1 | tee log_collect_${tag}.txt

    if [ ! -f impacts_${tag}.json ]; then
        echo "JSON not created for ${tag}"; return 1
    fi

    # Step 4: Plot PDF named after mode
    echo ">>> Step 4: Plotting [${tag}]"
    plotImpacts.py \
        -i impacts_${tag}.json \
        -o plot_impacts_${tag} \
        2>&1 | tee log_plot_${tag}.txt

    # Verify r-hat is sensible
    python3 -c "
import json, sys
with open('impacts_${tag}.json') as f:
    d = json.load(f)
r_fit = d['POIs'][0]['fit'][1]
r_lo  = d['POIs'][0]['fit'][0]
r_hi  = d['POIs'][0]['fit'][2]
mode  = '${mode}'
expect = ${expect_signal_val}
print(f'  r_hat = {r_fit:.4f}  [{r_lo:.4f}, {r_hi:.4f}]')
if abs(r_fit - expect) < 0.5:
    print(f'consistent with mode={mode} (expected ~{expect})')
else:
    print(f'r={r_fit:.4f} deviates from expected {expect} (mode={mode})')
n_params = len(d['params'])
print(f'  NPs in impact plot: {n_params}')
"

}

# =============================================================================
# process_job: runs the full workspace-build + impact pipeline for one
#              (YEAR, MASS, MODE) combination, called as a background job.
# =============================================================================
process_job() {
    local YEAR=$1
    local MASS=$2
    local MODE=$3

    local JOB_TAG="${YEAR}_${MASS}_${MODE}"
    local JOB_LOG="${LOGDIR}/job_${JOB_TAG}.log"
    local STATUS_FILE="${LOGDIR}/status_${JOB_TAG}"

    echo "Starting job: ${JOB_TAG}" | tee -a "${JOB_LOG}"

    (
        set +e   # Also disable inside subshell

        CARD_DIR="${CARD_BASE}/cards_${YEAR}_${MODE}"
        CARD="${CARD_DIR}/counting_t5wg_${MASS}.txt"

        WORKDIR="impacts_t5wg_${MASS}_${YEAR}_${MODE}"
        mkdir -p "${WORKDIR}"
        cd "${WORKDIR}"

        # Clean previous outputs
        rm -f higgsCombine*.root impacts*.json plot_impacts*.pdf \
              robustHesse*.root workspace*.root log_*.txt

        cp "${CARD}" "counting_t5wg_${MASS}.txt"

        echo ""
        echo ">>> Building workspace: t5wg_${MASS} (${YEAR}, ${MODE})"

        text2workspace.py "counting_t5wg_${MASS}.txt" \
            -m ${FIXED_MASS} \
            -o "workspace_t5wg_${MASS}.root" \
            2>&1 | tee log_t2w.txt

        if [ ! -f "workspace_t5wg_${MASS}.root" ]; then
            echo "Workspace not created for ${JOB_TAG}"
            exit 1
        fi

        WS="workspace_t5wg_${MASS}.root"

        run_impacts "${MASS}" "${MODE}" "${WS}"
        IMPACTS_EXIT=$?

        # Clean intermediate higgsCombine ROOT files
        rm -f higgsCombine*.root

        cd ..
        cp -r "${WORKDIR}" /eos/uscms/store/user/tmishra/Impacts/
        rm -rf "${WORKDIR}"

        exit ${IMPACTS_EXIT}

    ) >> "${JOB_LOG}" 2>&1
    local EXIT_CODE=$?

    if [ ${EXIT_CODE} -eq 0 ]; then
        echo "SUCCESS" > "${STATUS_FILE}"
        echo "Finished: ${JOB_TAG}"
    else
        echo "FAILED"  > "${STATUS_FILE}"
        echo "Failed:   ${JOB_TAG}  (see ${JOB_LOG})"
    fi
}

export -f run_impacts
export -f process_job
export CARD_BASE FIXED_MASS MINOPTS PARALLEL LOGDIR

# =============================================================================
# Main: iterate over all (YEAR, MASS, MODE) combinations and dispatch jobs,
#       capping concurrency at PARALLEL_JOBS.
# =============================================================================
PIDS=()        # track background PIDs
JOB_TAGS=()    # matching tags for wait/report

wait_for_slot() {
    # Block until the number of running background jobs drops below PARALLEL_JOBS
    while [ "${#PIDS[@]}" -ge "${PARALLEL_JOBS}" ]; do
        local new_pids=()
        for pid in "${PIDS[@]}"; do
            if kill -0 "${pid}" 2>/dev/null; then
                new_pids+=("${pid}")
            fi
        done
        PIDS=("${new_pids[@]}")
        [ "${#PIDS[@]}" -ge "${PARALLEL_JOBS}" ] && sleep 5
    done
}

for YEAR in "${YEARS[@]}"; do
    echo "=========================================="
    echo "Year: ${YEAR}"
    echo "=========================================="

    for MASS in "${MASS_POINTS[@]}"; do
        for MODE in "${MODES[@]}"; do

            wait_for_slot

            JOB_TAG="${YEAR}_${MASS}_${MODE}"
            echo "Dispatching: ${JOB_TAG}"

            process_job "${YEAR}" "${MASS}" "${MODE}" &
            PIDS+=($!)
            JOB_TAGS+=("${JOB_TAG}")

        done
    done
done

# Wait for all remaining jobs to finish
echo ""
echo "Waiting for all jobs to complete..."
wait

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "=========================================="
echo "Parallel job summary:"
FAILED=0
for TAG in "${JOB_TAGS[@]}"; do
    STATUS_FILE="${LOGDIR}/status_${TAG}"
    if [ -f "${STATUS_FILE}" ]; then
        STATUS=$(cat "${STATUS_FILE}")
    else
        STATUS="UNKNOWN"
    fi
    printf "  %-40s  %s\n" "${TAG}" "${STATUS}"
    [ "${STATUS}" != "SUCCESS" ] && FAILED=$((FAILED + 1))
done

echo ""
if [ "${FAILED}" -eq 0 ]; then
    echo "All jobs succeeded."
else
    echo "${FAILED} job(s) failed. Check ${LOGDIR}/job_*.log for details."
fi

echo ""
echo "Summary of output PDFs:"
find /eos/uscms/store/user/tmishra/Impacts/ -name "plot_impacts_*.pdf" | sort
echo "=========================================="
