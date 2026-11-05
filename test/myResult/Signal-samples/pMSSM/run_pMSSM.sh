#!/bin/bash
# run_pMSSM.sh
set -e

# --------------------------------------------------------------------------
# Build
# --------------------------------------------------------------------------
echo "[pMSSM] Compiling..."
g++ `root-config --cflags` \
    ../../../lib/libAnaClasses.so \
    analysis_pMSSM.C \
    -o analysis_pMSSM.exe \
    `root-config --libs` \
    -std=c++17
echo "[pMSSM] Build successful."

# --------------------------------------------------------------------------
# Helper: extract Year from filename (preVFP always 0)
# --------------------------------------------------------------------------
parse_year_from_filename() {
    local fname
    fname=$(basename "$1")

    if echo "$fname" | grep -qi "UL18"; then
        YEAR=2018
    elif echo "$fname" | grep -qi "UL17"; then
        YEAR=2017
    elif echo "$fname" | grep -qi "UL16"; then
        YEAR=2016
    else
        echo "[pMSSM] WARNING: cannot determine year from filename: $fname"
        YEAR=0
    fi
    PREVFP=0
}

make_output_tag() {
    local fname
    fname=$(basename "$1")
    fname="${fname%.*}"

    local set_tag sample_tag ul_tag
    set_tag=$(echo "$fname" | grep -oP 'set_\d+')
    
    ul_tag=$(echo "$fname" | grep -oiP 'UL1[678]' | head -1)

    sample_tag=$(echo "$fname" | grep -oP "(?<=set_\d_).*(?=_${ul_tag})")

    if [ -z "$sample_tag" ]; then
        sample_tag="sample"
    fi
    echo "PMSSM_${set_tag}_${sample_tag}_${ul_tag}"
}
# --------------------------------------------------------------------------
# Explicit command-line mode (four arguments passed)
# --------------------------------------------------------------------------
if [ "$#" -eq 4 ]; then
    YEAR=$1
    PREVFP=$2
    INPUT=$3
    OUTPUT=$4
    echo "[pMSSM] Explicit mode: Year=$YEAR preVFP=$PREVFP"
    echo "        Input:  $INPUT"
    echo "        Output: $OUTPUT"
    ./analysis_pMSSM.exe "$YEAR" "$PREVFP" "$INPUT" "$OUTPUT"
    exit 0
fi

# --------------------------------------------------------------------------
# Batch / Automatic Scanning Mode
# --------------------------------------------------------------------------
INPUT_DIR="/eos/uscms/store/user/lpcsusyphotons/SoftPhoton/pMSSM_lepPho/processedSamples"
OUTPUT_DIR="/uscms/home/tmishra/nobackup/pMSSM_output"
mkdir -p "$OUTPUT_DIR"

echo "[pMSSM] Scanning: $INPUT_DIR"
shopt -s nullglob

# Automatically collect both individual .root files AND .txt list files
FILES=("${INPUT_DIR}"/*.root "${INPUT_DIR}"/*.txt)

if [ "${#FILES[@]}" -eq 0 ]; then
    echo "[pMSSM] ERROR: no .root or .txt files found in $INPUT_DIR"
    exit 1
fi

echo "[pMSSM] Found ${#FILES[@]} target input file(s) [ROOT or TXT lists]."
echo "------------------------------------------------------------"

NRUN=0
NSKIP=0

for INPUT in "${FILES[@]}"; do
    FNAME=$(basename "$INPUT")

    # Optional era filter
    if [ -n "${FILTER_UL:-}" ] && ! echo "$FNAME" | grep -qi "$FILTER_UL"; then
        echo "[pMSSM] Skipping (filtered): $FNAME"
        (( NSKIP++ )) || true
        continue
    fi

    parse_year_from_filename "$INPUT"

    if [ "$YEAR" -eq 0 ]; then
        echo "[pMSSM] Skipping (unknown year): $FNAME"
        (( NSKIP++ )) || true
        continue
    fi

    OUT_TAG=$(make_output_tag "$INPUT")
    OUTPUT="${OUTPUT_DIR}/${OUT_TAG}_yields.root"

    echo "[pMSSM] --> Processing target: $FNAME"
    echo "         Year=$YEAR  preVFP=$PREVFP"
    echo "         Output: $(basename "$OUTPUT")"

    # Executive run: C++ code handles whether this string maps to a file list or single ROOT file
    ./analysis_pMSSM.exe "$YEAR" "$PREVFP" "$INPUT" "$OUTPUT"

    echo "         Done."
    echo "------------------------------------------------------------"
    (( NRUN++ )) || true
done

echo "[pMSSM] Finished. Ran: $NRUN  Skipped: $NSKIP"
