#!/bin/bash
# =============================================================================
# runImpacts_Run2_unblind.sh
# =============================================================================

set +e
# ---------------------------------------------------------------------------
CARD_BASE="/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards"
EOS_OUT="/eos/uscms/store/user/tmishra/Impacts/unblind"

MASS_POINTS=("1700_1100")
# MASS_POINTS=("1700_1100" "1300_500" "2100_100")

#YEARS=("Run2" "2016preVFP" "2016postVFP" "2017" "2018")
YEARS=("Run2" "2016preVFP" "2016postVFP" "2017" "2018")

# Which tasks to run (set to 0 to skip)
RUN_IMPACTS=1        # unblinded impacts
RUN_GOF_SAT=0
RUN_GOF_KS=0

# Number of toys for GOF
# Rule of thumb: >=1000 for a reliable p-value

GOF_TOYS=1
PARALLEL_JOBS=4
PARALLEL=4
FIXED_MASS=125

# Minimizer options (shared by impacts and GOF initial fit)
MINOPTS="--rMin -20 \
         --rMax 20 \
         --robustFit 1 \
         --cminDefaultMinimizerStrategy 0 \
         --X-rtd FITTER_NEW_CROSSING_ALGO \
         --X-rtd FITTER_NEVER_GIVE_UP \
         --X-rtd MINIMIZER_analytic"

# GOF-specific minimizer options:
# GoodnessOfFit does NOT accept --robustFit or --X-rtd MINIMIZER_analytic
# (these are MultiDimFit-only flags and cause "Invalid options" errors)
GOFMINOPTS="--rMin -20 \
            --rMax 20 \
            --cminDefaultMinimizerStrategy 0 \
            --X-rtd FITTER_NEW_CROSSING_ALGO \
            --X-rtd FITTER_NEVER_GIVE_UP"

LOGDIR="parallel_logs_unblind"
mkdir -p "${LOGDIR}"
mkdir -p "${EOS_OUT}"

# =============================================================================
# FUNCTION: run_impacts_unblind
#   Runs unblinded impact plot (floating r, real data).
# =============================================================================
run_impacts_unblind() {
    local mass=$1
    local year=$2
    local workspace=$3

    local tag="t5wg_${mass}_${year}_unblind"
    echo ""
    echo "========== IMPACTS (unblind): ${tag} =========="

    # ------------------------------------------------------------------
    # Step 1: Initial fit  (r floats freely on observed data)
    # ------------------------------------------------------------------
    echo ">>> Step 1: Initial fit [${tag}]"
    combineTool.py -M Impacts \
        -d "${workspace}" \
        -m ${FIXED_MASS} \
        --doInitialFit \
        ${MINOPTS} \
        2>&1 | tee log_initialFit_${tag}.txt

    if grep -q "ERROR\|Closed range\|not have converged" log_initialFit_${tag}.txt 2>/dev/null; then
        echo "WARNING: initial fit issues for ${tag} proceeding anyway"
    else
        echo "Initial fit OK for ${tag}"
    fi

    # ------------------------------------------------------------------
    # Step 2: NP scans (interactive, parallel threads)
    # ------------------------------------------------------------------
    echo ">>> Step 2: NP scans [${tag}]"
    combineTool.py -M Impacts \
        -d "${workspace}" \
        -m ${FIXED_MASS} \
        --doFits \
        --job-mode interactive \
        --parallel ${PARALLEL} \
        ${MINOPTS} \
        2>&1 | tee log_fits_${tag}.txt

    # ------------------------------------------------------------------
    # Step 3: Collect into JSON
    # ------------------------------------------------------------------
    echo ">>> Step 3: Collecting [${tag}]"
    combineTool.py -M Impacts \
        -d "${workspace}" \
        -m ${FIXED_MASS} \
        -o impacts_${tag}.json \
        2>&1 | tee log_collect_${tag}.txt

    if [ ! -f "impacts_${tag}.json" ]; then
        echo "ERROR: JSON not created for ${tag}"; return 1
    fi

    # ------------------------------------------------------------------
    # Step 4: Plot
    # ------------------------------------------------------------------
    echo ">>> Step 4: Plotting [${tag}]"
    plotImpacts.py \
        -i impacts_${tag}.json \
        -o plot_impacts_${tag} \
        --cms-label "Private Work" \
        2>&1 | tee log_plot_${tag}.txt

    # ------------------------------------------------------------------
    # Step 5: Print observed r-hat
    # ------------------------------------------------------------------
    python3 - <<PYEOF
import json
with open('impacts_${tag}.json') as f:
    d = json.load(f)
r_lo, r_fit, r_hi = d['POIs'][0]['fit']
n_np = len(d['params'])
print(f"  [${tag}]  r_hat = {r_fit:.4f}  [{r_lo:.4f}, {r_hi:.4f}]  ({n_np} NPs)")
if r_fit < 0:
    print("  NOTE: r_hat < 0  data is consistent with background-only")
elif r_fit > 2:
    print("  NOTE: r_hat > 2  large signal injection or fit issue, check carefully")
else:
    print("  r_hat looks reasonable for an unblinded fit")
PYEOF

}

# =============================================================================
# FUNCTION: run_gof
#   Goodness-of-fit test at r=0 (background-only hypothesis).
#
#   Physics rationale:
#     We test whether the background model describes the observed data well,
#     WITHOUT injecting any signal (--fixedSignalStrength 0).
#     The test statistic q_GOF is compared to a toy distribution generated
#     under the same b-only hypothesis to get a p-value.
#     A small p-value (< 0.05) indicates tension between data and the model.
#
#   Algorithm choices:
#     saturated  -- recommended for counting experiments (our case).
#                   q = -2 * ln(L_bonly / L_saturated)
#                   where L_saturated is the likelihood at n_obs = n_exp in
#                   every bin, i.e. the best possible fit.
#     KS / AD    -- Kolmogorov-Smirnov / Anderson-Darling, shape-based.
#                   Less appropriate for counting experiments but useful
#                   as a cross-check.
# =============================================================================
run_gof() {
    local mass=$1
    local year=$2
    local workspace=$3
    local algo=$4    # "saturated", "KS", or "AD"

    local tag="t5wg_${mass}_${year}_${algo}"
    echo ""
    echo "========== GOF (${algo}): ${tag} =========="

    # ------------------------------------------------------------------
    # Step 1: Observed test statistic on real data
    #   --fixedSignalStrength 0  → r fixed to 0 (b-only)
    #   No --toys flag           → runs on observed data
    # ------------------------------------------------------------------
    echo ">>> Step 1: Observed GOF statistic [${tag}]"
    combine -M GoodnessOfFit \
        -d "${workspace}" \
        -m ${FIXED_MASS} \
        --algo ${algo} \
        --fixedSignalStrength 0 \
        -n _GOF_obs_${tag} \
        ${GOFMINOPTS} \
        2>&1 | tee log_gof_obs_${tag}.txt

    if [ ! -f "higgsCombine_GOF_obs_${tag}.GoodnessOfFit.mH${FIXED_MASS}.root" ]; then
        echo "ERROR: observed GOF output not found for ${tag}"; return 1
    fi

    # ------------------------------------------------------------------
    # Step 2: Toy distribution (b-only throws)
    #   --toys ${GOF_TOYS}       → generate toys under b-only hypothesis
    #   --fixedSignalStrength 0  → r=0 for toy generation
    #   --toysFrequentist        → frequentist toys (profiled NPs)
    #   --seed 12345             → reproducible
    # ------------------------------------------------------------------
    echo ">>> Step 2: Toy distribution [${tag}] (${GOF_TOYS} toys)"
    combine -M GoodnessOfFit \
        -d "${workspace}" \
        -m ${FIXED_MASS} \
        --algo ${algo} \
        --fixedSignalStrength 0 \
        --toysFrequentist \
        --toys ${GOF_TOYS} \
        --seed 12345 \
        -n _GOF_toys_${tag} \
        ${GOFMINOPTS} \
        2>&1 | tee log_gof_toys_${tag}.txt

    if [ ! -f "higgsCombine_GOF_toys_${tag}.GoodnessOfFit.mH${FIXED_MASS}.123456.root" ] && \
       [ ! -f "higgsCombine_GOF_toys_${tag}.GoodnessOfFit.mH${FIXED_MASS}.12345.root" ]; then
        echo "WARNING: toys output filename not found by expected name — checking..."
        ls higgsCombine_GOF_toys_${tag}.GoodnessOfFit.mH*.root 2>/dev/null || \
            { echo "ERROR: toy output not found for ${tag}"; return 1; }
    fi

    # ------------------------------------------------------------------
    # Step 3: Collect observed + toys into JSON using combineTool
    # ------------------------------------------------------------------
    echo ">>> Step 3: Collecting GOF results [${tag}]"
    local obs_file="higgsCombine_GOF_obs_${tag}.GoodnessOfFit.mH${FIXED_MASS}.root"
    local toy_file
    toy_file=$(ls higgsCombine_GOF_toys_${tag}.GoodnessOfFit.mH${FIXED_MASS}.*.root 2>/dev/null | head -1)

    combineTool.py -M CollectGoodnessOfFit \
        --input "${obs_file}" "${toy_file}" \
        -o gof_${tag}.json \
        2>&1 | tee log_gof_collect_${tag}.txt

    if [ ! -f "gof_${tag}.json" ]; then
        echo "ERROR: GOF JSON not created for ${tag}"; return 1
    fi

    # ------------------------------------------------------------------
    # Step 4: Plot GOF distribution + print p-value
    # plotGof.py does not exist in CombineHarvester v10.3.0, so we use
    # a self-contained Python/matplotlib plot instead.
    # JSON structure (confirmed): d["125.0"] = {"obs": float, "p": float, "toy": [list]}
    # ------------------------------------------------------------------
    echo ">>> Step 4: Plotting GOF + p-value [${tag}]"
    python3 - <<PYEOF
import json, numpy as np, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

tag  = "${tag}"
algo = "${algo}"

try:
    with open(f"gof_{tag}.json") as f:
        d = json.load(f)
except FileNotFoundError:
    print(f"  ERROR: gof_{tag}.json not found")
    sys.exit(1)

# CombineHarvester v10.3.0 structure: keyed by mass as string
mass_key = list(d.keys())[0]   # e.g. "125.0"
entry    = d[mass_key]

# entry["obs"] and entry["p"] can be a list (e.g. [45.3]) or a scalar
# depending on CombineHarvester version and number of toys — always flatten
def to_scalar(x):
    if isinstance(x, (list, np.ndarray)):
        return float(np.atleast_1d(x).flatten()[0])
    return float(x)

obs_val  = to_scalar(entry["obs"])
toy_vals = np.array(entry["toy"]).flatten()
p_val    = to_scalar(entry["p"])
p_manual = float(np.mean(toy_vals >= obs_val))

# ------------------------------------------------------------------
# Print summary
# ------------------------------------------------------------------
print(f"  [{tag}]  GOF ({algo})")
print(f"  Observed q_GOF = {obs_val:.4f}")
print(f"  Toys           = {len(toy_vals)}")
print(f"  Toy median     = {np.median(toy_vals):.4f}")
print(f"  Toy mean       = {np.mean(toy_vals):.4f}")
print(f"  p-value (JSON) = {p_val:.4f}")
print(f"  p-value (check)= {p_manual:.4f}")
if p_val < 0.01:
    print("  *** p < 0.01 : SIGNIFICANT tension between data and model ***")
elif p_val < 0.05:
    print("  **  p < 0.05 : marginal tension — investigate further")
else:
    print("      p >= 0.05 : data consistent with background model — GOF passed")

# ------------------------------------------------------------------
# Plot: toy distribution as histogram, observed as vertical line
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))

# Histogram of toys
n_bins = 40
counts, edges, patches = ax.hist(
    toy_vals, bins=n_bins,
    color="#5b9bd5", edgecolor="white", linewidth=0.4,
    label=f"Toys (N={len(toy_vals)})"
)

# Shade tail beyond observed (toys >= obs)
bin_centers = 0.5 * (edges[:-1] + edges[1:])
for patch, left_edge in zip(patches, edges[:-1]):
    if left_edge >= obs_val:
        patch.set_facecolor("#e05c5c")

# Observed line
ax.axvline(obs_val, color="black", linewidth=2.0,
           label=f"Observed  q = {obs_val:.2f}")

# Median line
ax.axvline(np.median(toy_vals), color="#2e75b6", linewidth=1.5,
           linestyle="--", label=f"Toy median = {np.median(toy_vals):.2f}")

# p-value annotation
ax.text(0.97, 0.95, f"p-value = {p_val:.3f}",
        transform=ax.transAxes, ha="right", va="top",
        fontsize=13, fontweight="bold",
        color="green" if p_val >= 0.05 else "red")

ax.set_xlabel(f"Test statistic  q_{{GOF}}  ({algo})", fontsize=13)
ax.set_ylabel("Toys", fontsize=13)
ax.set_title(f"Goodness-of-Fit   {tag.replace('_', ' ')}", fontsize=13)
ax.legend(fontsize=11)

# Add red patch to legend for tail
tail_patch = mpatches.Patch(color="#e05c5c",
                             label=f"Tail (toys ≥ obs):  {int(np.sum(toy_vals >= obs_val))}")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles + [tail_patch], fontsize=10)

plt.tight_layout()
outname = f"plot_gof_{tag}.pdf"
plt.savefig(outname, bbox_inches="tight")
print(f"  Saved: {outname}")
plt.close()
PYEOF

}

# =============================================================================
# FUNCTION: process_job
#   Orchestrates all tasks for one (YEAR, MASS) combination.
#   Builds workspace once, then calls run_impacts_unblind and/or run_gof.
# =============================================================================
process_job() {
    local YEAR=$1
    local MASS=$2

    local JOB_TAG="${YEAR}_${MASS}_unblind"
    local JOB_LOG="${LOGDIR}/job_${JOB_TAG}.log"
    local STATUS_FILE="${LOGDIR}/status_${JOB_TAG}"

    echo "Starting job: ${JOB_TAG}" | tee -a "${JOB_LOG}"

    (
        set +e

        CARD_DIR="${CARD_BASE}/cards_${YEAR}_unblind"
        CARD="${CARD_DIR}/counting_t5wg_${MASS}.txt"

        if [ ! -f "${CARD}" ]; then
            echo "ERROR: datacard not found: ${CARD}"
            exit 1
        fi

        WORKDIR="unblind_t5wg_${MASS}_${YEAR}"
        mkdir -p "${WORKDIR}"
        cd "${WORKDIR}"

        # Clean previous outputs
        rm -f higgsCombine*.root gof_*.json impacts_*.json \
              plot_impacts_*.pdf plot_gof_*.pdf \
              robustHesse*.root workspace*.root log_*.txt

        cp "${CARD}" "counting_t5wg_${MASS}.txt"

        # ----------------------------------------------------------------
        # Build workspace (shared by impacts and GOF)
        # ----------------------------------------------------------------
        echo ""
        echo ">>> Building workspace: t5wg_${MASS} (${YEAR}, unblind)"
        text2workspace.py "counting_t5wg_${MASS}.txt" \
            -m ${FIXED_MASS} \
            -o "workspace_t5wg_${MASS}.root" \
            2>&1 | tee log_t2w.txt

        if [ ! -f "workspace_t5wg_${MASS}.root" ]; then
            echo "ERROR: workspace not created for ${JOB_TAG}"
            exit 1
        fi

        WS="workspace_t5wg_${MASS}.root"
        local OVERALL_EXIT=0

        # ----------------------------------------------------------------
        # Task 1: Unblinded impacts
        # ----------------------------------------------------------------
        if [ "${RUN_IMPACTS}" -eq 1 ]; then
            run_impacts_unblind "${MASS}" "${YEAR}" "${WS}"
            [ $? -ne 0 ] && OVERALL_EXIT=1
        fi

        # ----------------------------------------------------------------
        # Task 2: GOF with saturated model (recommended)
        # ----------------------------------------------------------------
        if [ "${RUN_GOF_SAT}" -eq 1 ]; then
            run_gof "${MASS}" "${YEAR}" "${WS}" "saturated"
            [ $? -ne 0 ] && OVERALL_EXIT=1
        fi

        # ----------------------------------------------------------------
        # Task 3: GOF with KS (optional cross-check)
        # ----------------------------------------------------------------
        if [ "${RUN_GOF_KS}" -eq 1 ]; then
            run_gof "${MASS}" "${YEAR}" "${WS}" "KS"
            [ $? -ne 0 ] && OVERALL_EXIT=1
        fi

        # ----------------------------------------------------------------
        # Clean large intermediate files, keep only final outputs
        # ----------------------------------------------------------------
        rm -f higgsCombine*.root robustHesse*.root

        # ----------------------------------------------------------------
        # Copy results to EOS
        # ----------------------------------------------------------------
        cd ..
        echo ">>> Copying ${WORKDIR} to ${EOS_OUT}/"
        cp -r "${WORKDIR}" "${EOS_OUT}/"
        rm -rf "${WORKDIR}"

        exit ${OVERALL_EXIT}

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

export -f run_impacts_unblind
export -f run_gof
export -f process_job
export CARD_BASE EOS_OUT FIXED_MASS MINOPTS GOFMINOPTS PARALLEL LOGDIR
export RUN_IMPACTS RUN_GOF_SAT RUN_GOF_KS GOF_TOYS

# =============================================================================
# Main: dispatch (YEAR x MASS) jobs, capped at PARALLEL_JOBS
# =============================================================================
PIDS=()
JOB_TAGS=()

wait_for_slot() {
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

echo "=========================================="
echo " Unblinded Impacts + GOF"
echo " Mass points : ${MASS_POINTS[*]}"
echo " Years       : ${YEARS[*]}"
echo " Tasks       : impacts=${RUN_IMPACTS}  gof_sat=${RUN_GOF_SAT}  gof_ks=${RUN_GOF_KS}"
echo " GOF toys    : ${GOF_TOYS}"
echo " Parallel    : ${PARALLEL_JOBS} outer x ${PARALLEL} NP threads"
echo "=========================================="

for YEAR in "${YEARS[@]}"; do
    echo ""
    echo "--- Year: ${YEAR} ---"
    for MASS in "${MASS_POINTS[@]}"; do

        wait_for_slot

        JOB_TAG="${YEAR}_${MASS}_unblind"
        echo "Dispatching: ${JOB_TAG}"

        process_job "${YEAR}" "${MASS}" &
        PIDS+=($!)
        JOB_TAGS+=("${JOB_TAG}")

    done
done

echo ""
echo "Waiting for all jobs to complete..."
wait

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "=========================================="
echo "Job summary:"
FAILED=0
for TAG in "${JOB_TAGS[@]}"; do
    STATUS_FILE="${LOGDIR}/status_${TAG}"
    STATUS=$( [ -f "${STATUS_FILE}" ] && cat "${STATUS_FILE}" || echo "UNKNOWN" )
    printf "  %-45s  %s\n" "${TAG}" "${STATUS}"
    [ "${STATUS}" != "SUCCESS" ] && FAILED=$((FAILED + 1))
done

echo ""
if [ "${FAILED}" -eq 0 ]; then
    echo "All jobs succeeded."
else
    echo "${FAILED} job(s) failed. Check ${LOGDIR}/job_*.log for details."
fi

echo ""
echo "Output PDFs on EOS:"
find "${EOS_OUT}" -name "plot_impacts_*unblind*.pdf" -o \
                  -name "plot_gof_*.pdf" 2>/dev/null | sort
echo "=========================================="
