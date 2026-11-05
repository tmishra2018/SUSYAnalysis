#!/bin/bash
PARALLEL_JOBS=8   # tune to match available cores (check with: nproc)
LOGDIR="analysis_logs"
mkdir -p "${LOGDIR}"

# -----------------------------------------------------------------------------
# Step 1: Compile both executables (serial, must finish before running)
# -----------------------------------------------------------------------------
echo ">>> Compiling analysis_mgMC.C ..."
g++ `root-config --cflags` ../../lib/libAnaClasses.so analysis_mgMC.C \
    -o analysis_mgMC.exe `root-config --libs`
if [ $? -ne 0 ]; then echo "Compilation of analysis_mgMC failed"; exit 1; fi

echo ">>> Compiling analysis_egMC.C ..."
g++ `root-config --cflags` ../../lib/libAnaClasses.so analysis_egMC.C \
    -o analysis_egMC.exe `root-config --libs`
if [ $? -ne 0 ]; then echo "Compilation of analysis_egMC failed"; exit 1; fi

echo "Compilation done."
echo ""

# -----------------------------------------------------------------------------
# Step 2: Define jobs
# exe, year, preVFP, sample
# -----------------------------------------------------------------------------
JOBS=(
    "analysis_mgMC.exe 2016 1 DYJetsToLL"
    "analysis_mgMC.exe 2016 1 TTJets"
    "analysis_mgMC.exe 2016 1 WW"
    "analysis_mgMC.exe 2016 1 WZ"
    "analysis_mgMC.exe 2016 1 WJetsToLNu"
    "analysis_mgMC.exe 2016 1 WGToLNuG"
    "analysis_mgMC.exe 2016 0 DYJetsToLL"
    "analysis_mgMC.exe 2016 0 TTJets"
    "analysis_mgMC.exe 2016 0 WW"
    "analysis_mgMC.exe 2016 0 WZ"
    "analysis_mgMC.exe 2016 0 WJetsToLNu"
    "analysis_mgMC.exe 2016 0 WGToLNuG"
    "analysis_mgMC.exe 2017 1 DYJetsToLL"
    "analysis_mgMC.exe 2017 1 TTJets"
    "analysis_mgMC.exe 2017 1 WW"
    "analysis_mgMC.exe 2017 1 WZ"
    "analysis_mgMC.exe 2017 1 WJetsToLNu"
    "analysis_mgMC.exe 2017 1 WGToLNuG"
    "analysis_mgMC.exe 2018 1 DYJetsToLL"
    "analysis_mgMC.exe 2018 1 TTJets"
    "analysis_mgMC.exe 2018 1 WW"
    "analysis_mgMC.exe 2018 1 WZ"
    "analysis_mgMC.exe 2018 1 WJetsToLNu"
    "analysis_mgMC.exe 2018 1 WGToLNuG"
    "analysis_egMC.exe 2016 1 DYJetsToLL"
    "analysis_egMC.exe 2016 1 TTJets"
    "analysis_egMC.exe 2016 1 WW"
    "analysis_egMC.exe 2016 1 WZ"
    "analysis_egMC.exe 2016 1 WJetsToLNu"
    "analysis_egMC.exe 2016 1 WGToLNuG"
    "analysis_egMC.exe 2016 0 DYJetsToLL"
    "analysis_egMC.exe 2016 0 TTJets"
    "analysis_egMC.exe 2016 0 WW"
    "analysis_egMC.exe 2016 0 WZ"
    "analysis_egMC.exe 2016 0 WJetsToLNu"
    "analysis_egMC.exe 2016 0 WGToLNuG"
    "analysis_egMC.exe 2017 1 DYJetsToLL"
    "analysis_egMC.exe 2017 1 TTJets"
    "analysis_egMC.exe 2017 1 WW"
    "analysis_egMC.exe 2017 1 WZ"
    "analysis_egMC.exe 2017 1 WJetsToLNu"
    "analysis_egMC.exe 2017 1 WGToLNuG"
    "analysis_egMC.exe 2018 1 DYJetsToLL"
    "analysis_egMC.exe 2018 1 TTJets"
    "analysis_egMC.exe 2018 1 WW"
    "analysis_egMC.exe 2018 1 WZ"
    "analysis_egMC.exe 2018 1 WJetsToLNu"
    "analysis_egMC.exe 2018 1 WGToLNuG"
)

# -----------------------------------------------------------------------------
# Step 3: Dispatch jobs with concurrency cap
# -----------------------------------------------------------------------------
PIDS=()
TAGS=()

wait_for_slot() {
    while [ "${#PIDS[@]}" -ge "${PARALLEL_JOBS}" ]; do
        local new_pids=()
        for pid in "${PIDS[@]}"; do
            if kill -0 "${pid}" 2>/dev/null; then
                new_pids+=("${pid}")
            fi
        done
        PIDS=("${new_pids[@]}")
        [ "${#PIDS[@]}" -ge "${PARALLEL_JOBS}" ] && sleep 2
    done
}

for JOB in "${JOBS[@]}"; do
    wait_for_slot

    # Parse job fields
    EXE=$(echo "${JOB}"    | awk '{print $1}')
    YEAR=$(echo "${JOB}"   | awk '{print $2}')
    PREVFP=$(echo "${JOB}" | awk '{print $3}')
    SAMPLE=$(echo "${JOB}" | awk '{print $4}')

    TAG="${EXE%.exe}_${YEAR}_${PREVFP}_${SAMPLE}"
    LOG="${LOGDIR}/${TAG}.log"

    echo "Dispatching: ./${EXE} ${YEAR} ${PREVFP} ${SAMPLE}"

    ./${EXE} ${YEAR} ${PREVFP} ${SAMPLE} > "${LOG}" 2>&1 &
    PIDS+=($!)
    TAGS+=("${TAG}")
done

echo ""
echo "Waiting for all jobs to finish..."
wait

# -----------------------------------------------------------------------------
# Step 4: Summary
# -----------------------------------------------------------------------------
echo ""
echo "=========================================="
echo "Job summary:"
FAILED=0
for i in "${!TAGS[@]}"; do
    TAG="${TAGS[$i]}"
    LOG="${LOGDIR}/${TAG}.log"
    # Check log for common error indicators
    if grep -qi "error\|segfault\|abort\|terminate" "${LOG}" 2>/dev/null; then
        printf "FAILED   %s  (see %s)\n" "${TAG}" "${LOG}"
        FAILED=$((FAILED + 1))
    else
        printf "  OK       %s\n" "${TAG}"
    fi
done

echo ""
if [ "${FAILED}" -eq 0 ]; then
    echo "All ${#TAGS[@]} jobs completed successfully."
else
    echo "${FAILED}/${#TAGS[@]} job(s) had errors. Check ${LOGDIR}/ for details."
fi
echo "=========================================="
