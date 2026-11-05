#!/bin/bash

years=("2016preVFP" "2016postVFP" "2017" "2018")

get_card_path() {
    local year=$1
    case $year in
        "2016preVFP")  echo "/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_2016preVFP_new" ;;
        "2016postVFP") echo "/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_2016postVFP_new" ;;
        "2017")        echo "/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_2017_new" ;;
        "2018")        echo "/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_2018_new" ;;
    esac
}

for year in "${years[@]}"; do
    echo "=========================================="
    echo "Processing year: $year"
    echo "=========================================="

    mkdir -p Signal_t5wg_${year}
    cd Signal_t5wg_${year} || { echo "Directory not found!"; exit 1; }

    rm -f higgsCombine*.root impacts*.json robustHesse*.root workspace*.root gmn_params.txt

    CARD_PATH=$(get_card_path ${year})
    cp ${CARD_PATH}/counting_t5wg_1700_1100.txt . || {
        echo "❌ Card not found at ${CARD_PATH}"; cd ..; continue
    }

    # ✅ Extract gmN parameter names to freeze them later
    grep "gmN" counting_t5wg_1700_1100.txt | awk '{print $1}' \
        | tr '\n' ',' | sed 's/,$//' > gmn_params.txt
    GMN_PARAMS=$(cat gmn_params.txt)
    echo "gmN params to freeze: $(echo $GMN_PARAMS | tr ',' '\n' | wc -l) parameters"

    # ✅ Build workspace from ORIGINAL (unmodified) card
    text2workspace.py counting_t5wg_1700_1100.txt \
        -o workspace_t5wg_1700_1100.root

    # ✅ Step 1: Initial fit
    # --freezeParameters freezes all gmN CR counts at their prefit (background-only) values
    # --setParameters r=0 fixes signal strength to 0
    # Together these force a true background-only Asimov fit
    combineTool.py -M Impacts \
        -d workspace_t5wg_1700_1100.root \
        -m 1700 \
        -t -1 --expectSignal 0 \
        --bypassFrequentistFit \
        --freezeParameters ${GMN_PARAMS} \
        --doInitialFit --parallel 10 \
        --cminDefaultMinimizerStrategy 0 \
        --rMin -20 --rMax 20 --robustFit 1 \
        --X-rtd MINIMIZER_analytic --X-rtd FAST_VERTICAL_MORPH

    # ✅ Verify r=0
    INITFILE=$(ls higgsCombine_initialFit*.MultiDimFit*.root 2>/dev/null | head -1)
    echo "=== Verifying r after initial fit ==="
    python3 -c "
import uproot, sys
with uproot.open('${INITFILE}') as f:
    r = float(f['limit']['r'].array()[0])
    print(f'r = {r:.6f}')
    if abs(r) < 0.001:
        print('✅ CORRECTLY BLINDED')
    else:
        print(f'❌ NOT BLINDED r={r:.4f}')
        sys.exit(1)
" || { echo "❌ Blinding failed for $year"; cd ..; continue; }

    # ✅ Step 2: NP fits
    combineTool.py -M Impacts \
        -d workspace_t5wg_1700_1100.root \
        -m 1700 \
        -t -1 --expectSignal 0 \
        --bypassFrequentistFit \
        --freezeParameters ${GMN_PARAMS} \
        --doFits --parallel 10 \
        --cminDefaultMinimizerStrategy 0 \
        --rMin -20 --rMax 20 --robustFit 1 \
        --X-rtd MINIMIZER_analytic --X-rtd FAST_VERTICAL_MORPH

    # ✅ Step 3: Collect JSON
    combineTool.py -M Impacts \
        -d workspace_t5wg_1700_1100.root \
        -m 1700 \
        -t -1 --expectSignal 0 \
        --freezeParameters ${GMN_PARAMS} \
        -o impacts_t5wg_1700_1100.json \
        --parallel 10

    [ -f impacts_t5wg_1700_1100.json ] || {
        echo "❌ JSON not created for $year"; cd ..; continue
    }

    # ✅ Final verification
    python3 -c "
import json
with open('impacts_t5wg_1700_1100.json') as f:
    d = json.load(f)
r_fit = d['POIs'][0]['fit'][1]
pulls = [abs(p['fit'][1]) for p in d['params']]
print(f'r_hat    = {r_fit:.4f}  (expect 0.00)')
print(f'Max pull = {max(pulls):.4f}  (expect ~0.00)')
print('✅ BLINDED' if abs(r_fit) < 0.01 else '❌ NOT BLINDED')
"
    plotImpacts.py -i impacts_t5wg_1700_1100.json \
        -o impacts_t5wg_1700_1100_${year}

    INPUT_PDF=impacts_t5wg_1700_1100_${year}.pdf
    if [ -f "${INPUT_PDF}" ] && command -v gs &>/dev/null; then
        gs -dNOPAUSE -dBATCH -sDEVICE=pdfwrite \
            -dFirstPage=1 -dLastPage=1 \
            -sOutputFile=impacts_t5wg_1700_1100_${year}_page1.pdf ${INPUT_PDF}
        gs -dNOPAUSE -dBATCH -sDEVICE=pdfwrite \
            -dFirstPage=2 -dLastPage=9999 \
            -sOutputFile=impacts_t5wg_1700_1100_${year}_allpages.pdf ${INPUT_PDF}
        echo "✅ PDFs split for $year"
    fi

    cd ..
    echo "✅ Done: $year"
done
