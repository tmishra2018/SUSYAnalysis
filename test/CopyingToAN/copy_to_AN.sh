#!/bin/bash

set -euo pipefail

rm -f 2016preVFP_.txt 2016postVFP_.txt 2017_.txt 2018_.txt

echo "Running copy.sh for each year..."
bash listing.sh >> 2016preVFP_.txt
bash listing.sh >> 2016postVFP_.txt
bash listing.sh >> 2017_.txt
bash listing.sh >> 2018_.txt

echo "Cleaning up text files..."

for FILE in 2016preVFP_.txt 2016postVFP_.txt 2017_.txt 2018_.txt; do
    echo "  Cleaning $FILE..."
    case "$FILE" in
        2016preVFP_*.txt)
            sed -i '/2016preVFP/!d' "$FILE"
            ;;
        2016postVFP_*.txt)
            sed -i '/2016postVFP/!d' "$FILE"
            ;;
        2017_*.txt)
            sed -i '/2017/!d' "$FILE"
            ;;
        2018_*.txt)
	    sed -i '/2018\|678/!d' "$FILE"
            ;;
    esac
done

echo "✅ Cleanup done."

echo "Starting file transfers..."

PLOTS_SRC="/eos/uscms/store/user/tmishra/Background/plots"
PDF_SRC="/uscms_data/d3/tmishra/Output"
DEST_BASE="/afs/cern.ch/work/t/trmishra/AnalysisNote/AN-21-007"
USER="trmishra@lxplus.cern.ch"

YEARS=("2016preVFP" "2016postVFP" "2017" "2018")

function green() { echo -e "\033[0;32m$1\033[0m"; }
function yellow() { echo -e "\033[1;33m$1\033[0m"; }
function red() { echo -e "\033[0;31m$1\033[0m"; }

copy_matching_files() {
    local pattern="$1"
    local year="$2"
    local dest="$3"
    cp ${PLOTS_SRC}/201*/*gamma/${pattern}_${year}.png "$dest" 2>/dev/null || true
}

copy_listed_files() {
    local txtfile="$1"
    local dest="$2"
    if [[ -f "$txtfile" ]]; then
        green "  - Copying files listed in $(basename "$txtfile")..."
        while IFS= read -r line || [[ -n "$line" ]]; do
            [[ -z "$line" ]] && continue  # skip empty lines
            if [[ -f "$line" ]]; then
                cp "$line" "$dest"
            else
                yellow "    Warning: Listed file $line not found, skipping."
            fi
        done < "$txtfile"
    else
        yellow "  - Warning: Text file $txtfile not found for this year."
    fi
}

for YEAR in "${YEARS[@]}"; do
    TMPDIR="tmp_${YEAR}"
    DESTDIR="${DEST_BASE}/Figures"
    case "$YEAR" in
        2016postVFP) DESTDIR="${DEST_BASE}/Figures16postVFP" ;;
        2017)        DESTDIR="${DEST_BASE}/Figures17" ;;
        2018)        DESTDIR="${DEST_BASE}/Figures18" ;;
    esac

    green "Processing $YEAR..."
    mkdir -p "$TMPDIR"

    green "  - Copying VALID and BKG .png files..."
    for var in ht met _pt; do
        copy_matching_files "VALID_*${var}*" "$YEAR" "$TMPDIR"
    done
    for var in _mt _pt _ht; do
        copy_matching_files "BKG_*${var}*" "$YEAR" "$TMPDIR"
    done

    green "  - Copying signal count PDF if it exists..."
    PDF_FILE="${PDF_SRC}/signalCount_${YEAR}_NoData.pdf"
    if [[ -f "$PDF_FILE" ]]; then
        cp "$PDF_FILE" "$TMPDIR"
    else
        yellow "  - Warning: Signal count file $PDF_FILE not found."
    fi

    green "  - Copying additional files from text list..."
    TXT_FILE="${YEAR}_.txt"
    copy_listed_files "$TXT_FILE" "$TMPDIR"

    green "  - Syncing to lxplus:${DESTDIR}/..."
    rsync -avz --progress "$TMPDIR/" "${USER}:${DESTDIR}/"

    green "  - Cleaning up local temp folder..."
    rm -rf "$TMPDIR"

    green "Done with $YEAR"
    echo "------------------------------"
done

green "✅ All background plot transfers completed successfully!"
