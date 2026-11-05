#!/bin/bash
# =============================================================================
# split_and_run.sh
#
# Splits a large .txt file list into chunks of N files, runs analysis_pMSSM.exe
# on each chunk in parallel (up to MAX_JOBS simultaneously), then merges all
# output .root files with hadd.
#
# Usage:
#   ./split_and_run.sh <input_list.txt> <Year> <preVFP 0/1> [chunk_size] [max_jobs]
#
# Examples:
#   ./split_and_run.sh PMSSM_set_2_LL_1_UL18.txt 2018 0          # 200 files/chunk, 8 parallel
#   ./split_and_run.sh PMSSM_set_2_LL_1_UL18.txt 2018 0 150 6    # 150 files/chunk, 6 parallel
#
# Output:
#   <base>_chunks/          — working directory with chunk lists + per-chunk outputs
#   <base>_final.root       — merged output from all chunks
# =============================================================================

set -euo pipefail

# ---------- Arguments -------------------------------------------------------
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <input_list.txt> <Year> <preVFP 0/1> [chunk_size=200] [max_jobs=8]"
    exit 1
fi

INPUT_LIST="$1"
YEAR="$2"
PREVFP="$3"
CHUNK_SIZE="${4:-3}"
MAX_JOBS="${5:-8}"

if [ ! -f "$INPUT_LIST" ]; then
    echo "[split_and_run] ERROR: input list not found: $INPUT_LIST"
    exit 1
fi

# Derive base name (strip .txt)
BASE="${INPUT_LIST%.txt}"
CHUNK_DIR="${BASE}_chunks"
FINAL_OUTPUT="${BASE}_final.root"
EXE="./analysis_pMSSM.exe"

if [ ! -x "$EXE" ]; then
    echo "[split_and_run] ERROR: $EXE not found or not executable."
    echo "                Run run_pMSSM.sh first to compile, or pass the path."
    exit 1
fi

mkdir -p "$CHUNK_DIR"

# ---------- Split -----------------------------------------------------------
echo "[split_and_run] Splitting $INPUT_LIST into chunks of $CHUNK_SIZE files..."

# Filter blank/comment lines first into a clean temp list
CLEAN_LIST="${CHUNK_DIR}/clean_input.txt"
grep -v '^\s*$' "$INPUT_LIST" | grep -v '^\s*#' > "$CLEAN_LIST"

TOTAL_FILES=$(wc -l < "$CLEAN_LIST")
N_CHUNKS=$(( (TOTAL_FILES + CHUNK_SIZE - 1) / CHUNK_SIZE ))

echo "[split_and_run] Total files : $TOTAL_FILES"
echo "[split_and_run] Chunk size  : $CHUNK_SIZE"
echo "[split_and_run] N chunks    : $N_CHUNKS"
echo "[split_and_run] Parallel    : $MAX_JOBS"
echo "------------------------------------------------------------"

# Split into numbered chunk files: chunk_001.txt, chunk_002.txt, ...
split -l "$CHUNK_SIZE" -d -a 3 "$CLEAN_LIST" "${CHUNK_DIR}/chunk_"
# Give them .txt extension
for f in "${CHUNK_DIR}"/chunk_[0-9]*; do
    [[ "$f" == *.txt ]] && continue
    mv "$f" "${f}.txt"
done

CHUNK_FILES=("${CHUNK_DIR}"/chunk_*.txt)
echo "[split_and_run] Created ${#CHUNK_FILES[@]} chunk list(s)."

# ---------- Worker function -------------------------------------------------
run_chunk() {
    local CHUNK_TXT="$1"
    local YEAR="$2"
    local PREVFP="$3"
    local EXE="$4"

    local CHUNK_BASE="${CHUNK_TXT%.txt}"
    local CHUNK_OUT="${CHUNK_BASE}_output.root"
    local CHUNK_LOG="${CHUNK_BASE}.log"
    local CHUNK_NAME
    CHUNK_NAME=$(basename "$CHUNK_TXT")

    echo "[chunk] START  $CHUNK_NAME  ($(wc -l < "$CHUNK_TXT") files)"

    if "$EXE" "$YEAR" "$PREVFP" "$CHUNK_TXT" "$CHUNK_OUT" > "$CHUNK_LOG" 2>&1; then
        echo "[chunk] DONE   $CHUNK_NAME  ->  $(basename "$CHUNK_OUT")"
    else
        EXIT_CODE=$?
        echo "[chunk] FAILED $CHUNK_NAME  (exit $EXIT_CODE) — see $CHUNK_LOG"
        # Write a sentinel so the merge step can warn about missing chunks
        echo "FAILED" > "${CHUNK_BASE}.status"
        return 0   # don't abort the whole parallel run
    fi
}

export -f run_chunk

# ---------- Run in parallel with a semaphore --------------------------------
echo "[split_and_run] Launching parallel jobs (max $MAX_JOBS at a time)..."

PIDS=()
RUNNING=0

for CHUNK_TXT in "${CHUNK_FILES[@]}"; do
    # Wait if we're already at MAX_JOBS
    while [ "$RUNNING" -ge "$MAX_JOBS" ]; do
        # Reap any finished jobs
        NEW_PIDS=()
        for pid in "${PIDS[@]}"; do
            if kill -0 "$pid" 2>/dev/null; then
                NEW_PIDS+=("$pid")
            else
                (( RUNNING-- )) || true
            fi
        done
        PIDS=("${NEW_PIDS[@]}")
        sleep 2
    done

    run_chunk "$CHUNK_TXT" "$YEAR" "$PREVFP" "$EXE" &
    PIDS+=($!)
    (( RUNNING++ )) || true
done

# Wait for all remaining jobs
echo "[split_and_run] Waiting for remaining jobs..."
for pid in "${PIDS[@]}"; do
    wait "$pid" || true
done

echo "------------------------------------------------------------"
echo "[split_and_run] All chunks finished."

# ---------- Check for failures ----------------------------------------------
FAILED_CHUNKS=()
for CHUNK_TXT in "${CHUNK_FILES[@]}"; do
    CHUNK_BASE="${CHUNK_TXT%.txt}"
    if [ -f "${CHUNK_BASE}.status" ]; then
        FAILED_CHUNKS+=("$(basename "$CHUNK_TXT")")
    fi
done

if [ "${#FAILED_CHUNKS[@]}" -gt 0 ]; then
    echo "[split_and_run] WARNING: ${#FAILED_CHUNKS[@]} chunk(s) failed:"
    for f in "${FAILED_CHUNKS[@]}"; do
        echo "    $f"
    done
    echo "    Check the corresponding .log files in $CHUNK_DIR"
fi

# ---------- Merge -----------------------------------------------------------
echo "[split_and_run] Collecting output files for merge..."
OUTPUT_FILES=()
for CHUNK_TXT in "${CHUNK_FILES[@]}"; do
    CHUNK_OUT="${CHUNK_TXT%.txt}_output.root"
    CHUNK_BASE="${CHUNK_TXT%.txt}"
    # Skip if output is missing, marked failed, or has no keys (partial write from crash)
    if [ ! -f "$CHUNK_OUT" ]; then
        echo "[split_and_run] WARNING: no output file for $(basename "$CHUNK_TXT") — skipping"
        continue
    fi
    if [ -f "${CHUNK_BASE}.status" ]; then
        echo "[split_and_run] WARNING: $(basename "$CHUNK_TXT") marked failed — skipping its output"
        rm -f "$CHUNK_OUT"   # remove the empty/partial file so hadd never sees it
        continue
    fi
    OUTPUT_FILES+=("$CHUNK_OUT")
done

if [ "${#OUTPUT_FILES[@]}" -eq 0 ]; then
    echo "[split_and_run] ERROR: no output files to merge. All chunks failed?"
    exit 1
fi

echo "[split_and_run] Merging ${#OUTPUT_FILES[@]} / $N_CHUNKS output files..."
hadd -f "$FINAL_OUTPUT" "${OUTPUT_FILES[@]}"

echo "============================================================"
echo "[split_and_run] DONE"
echo "  Chunks run    : ${#OUTPUT_FILES[@]} / $N_CHUNKS succeeded"
echo "  Failed chunks : ${#FAILED_CHUNKS[@]}"
echo "  Final output  : $FINAL_OUTPUT"
echo "============================================================"
