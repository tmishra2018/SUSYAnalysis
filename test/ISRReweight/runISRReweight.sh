#!/bin/bash
# runISRReweight.sh
# ===================
# Master script to run the improved ISR reweighting chain.
#
# Executes in this order:
#   Step 1: Derive ISR weights for mg channel (MuonEG data, ZGToLLG MC)  [Rec.1]
#   Step 2: Derive ISR weights for eg channel (DoubleEG data, ZGToLLG MC) [Rec.1]
#   Step 3: Derive ISR weights for mg channel, WGToLNuG MC                [Rec.3]
#   Step 4: Derive ISR weights for eg channel, WGToLNuG MC                [Rec.3]
#   Step 5: Rare backgrounds (TTG/WWG/WZG) — NO weight derivation         [Rec.4]
#
# plot_ISRweight_Run2_improved.C arguments:
#   (RunYear, preVFP, channel, procTag, useGenPt)
#   useGenPt=true  -> Rec.2: use gen-level bosonPt as reweight variable
#   useGenPt=false -> use reco ISRJetPt (original behaviour)
#
# Normalization is derived SEPARATELY inside the script (Rec.5):
#   shape-only weights written to ISRweights_*_shapeOnly
#   shape+norm weights written to ISRweights_*_total
#
# Usage:
#   bash runISRReweight.sh 2018 0    # year=2018, preVFP=false
#   bash runISRReweight.sh 2016 1    # year=2016, preVFP=true
#   bash runISRReweight.sh 2016 0    # year=2016, preVFP=false
#   bash runISRReweight.sh 2017 0    # year=2017

set -e  # exit on error

YEAR=${1:-2018}
PREVFP=${2:-0}    # 0=false/postVFP, 1=true/preVFP
USEGEN=true       # Rec.2: use gen-level bosonPt

echo "========================================================"
echo " ISR Reweighting Chain — Year=${YEAR}  preVFP=${PREVFP}"
echo " useGenPt=${USEGEN}"
echo "========================================================"

# ---- Helper to run root macro with error checking -----------------------
run_root(){
    local CMD=$1
    echo ""
    echo ">>> Running: ${CMD}"
    root -b -q "${CMD}"
    if [ $? -ne 0 ]; then
        echo "[ERROR] Failed: ${CMD}"
        exit 1
    fi
}

# ========================================================================
# STEP 1: mg channel — ZGToLLG (Rec.1 + Rec.3)
# ========================================================================
echo ""
echo "=== STEP 1: mg channel, ZGToLLG ==="
run_root "plot_ISRweight_Run2_improved.C+(${YEAR},${PREVFP},\"mg\",\"ZGToLLG\",${USEGEN})"

# ========================================================================
# STEP 2: eg channel — ZGToLLG (Rec.1 + Rec.3)
# Requires resTree_ISR_data_eg_{year}.root from analysis_ISR_eg.C
# ========================================================================
echo ""
echo "=== STEP 2: eg channel, ZGToLLG ==="
run_root "plot_ISRweight_Run2_improved.C+(${YEAR},${PREVFP},\"eg\",\"ZGToLLG\",${USEGEN})"

# ========================================================================
# STEP 3: mg channel — WGToLNuG (Rec.3: separate Wgamma weight)
# Requires resTree_ISR_WGToLNuG_{year}.root from analysis_ISRMC_perprocess.C
# with isWchannel=true, PROC_TAG="WGToLNuG"
# ========================================================================
echo ""
echo "=== STEP 3: mg channel, WGToLNuG ==="
run_root "plot_ISRweight_Run2_improved.C+(${YEAR},${PREVFP},\"mg\",\"WGToLNuG\",${USEGEN})"

# ========================================================================
# STEP 4: eg channel — WGToLNuG (Rec.1 + Rec.3)
# ========================================================================
echo ""
echo "=== STEP 4: eg channel, WGToLNuG ==="
run_root "plot_ISRweight_Run2_improved.C+(${YEAR},${PREVFP},\"eg\",\"WGToLNuG\",${USEGEN})"

# ========================================================================
# STEP 5: Rec.4 — Rare backgrounds get NO ISR weight
# Print reminder instead
# ========================================================================
echo ""
echo "=== STEP 5: Rare backgrounds (TTG, WWG, WZG) ==="
echo "  Rec.4: No ISR weight derivation for rare backgrounds."
echo "  These processes receive ISRWeight=1 with 100% normalization"
echo "  uncertainty in apply_ISRWeight.C."
echo "  Assign a dedicated normalization uncertainty in your"
echo "  background estimation instead."

# ========================================================================
# STEP 6: Print summary of output files
# ========================================================================
echo ""
echo "========================================================"
echo " Output weight files:"
echo "========================================================"

VFP_TAG=""
if [ ${YEAR} -eq 2016 ] && [ ${PREVFP} -eq 1 ]; then VFP_TAG="preVFP"; fi
if [ ${YEAR} -eq 2016 ] && [ ${PREVFP} -eq 0 ]; then VFP_TAG="postVFP"; fi

for CHAN in mg eg; do
    for PROC in ZGToLLG WGToLNuG; do
        FNAME="ISRweights_${CHAN}_${PROC}_${YEAR}${VFP_TAG}.txt"
        if [ -f "${FNAME}" ]; then
            echo "  [OK] ${FNAME}"
            echo "       Contents:"
            grep -v "^#" ${FNAME} | head -10 | sed 's/^/         /'
        else
            echo "  [MISSING] ${FNAME}"
        fi
    done
done

echo ""
echo "To apply these weights in your analysis, add to analysis_VGBkg.C:"
echo '  #include "apply_ISRWeight.C"'
echo '  // Before event loop:'
echo '  initISRReaders("2018");'
echo '  // In event loop (per MC event):'
echo '  float isrUnc = 0;'
echo '  float isrW = getISRWeight("mg", bosonPt_gen, mcType, false, &isrUnc);'
echo '  // isrUnc = 1.0 for rare backgrounds (100% uncertainty flag)'
echo ""
echo "========================================================"
echo " Done."
echo "========================================================"
