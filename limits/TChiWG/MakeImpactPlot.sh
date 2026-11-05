#!/bin/bash

years=("2016preVFP" "2016postVFP" "2017" "2018")

# ---- TChiWG mass point (single mass) ----
MASS=825

for year in "${years[@]}"; do
    echo "Processing year: $year"

    mkdir -p Signal_tchiwg_${year}
    cd Signal_tchiwg_${year} || { echo "Directory Signal_tchiwg_${year} not found!"; exit 1; }

    rm -f *

    # Copy datacard
    cp /eos/uscms/store/user/tmishra/CombinedLimit/TChiWG/cards/cards_${year}/counting_tchiwg_${MASS}.txt .

    # Create workspace
    text2workspace.py counting_tchiwg_${MASS}.txt \
        -o workspace_tchiwg_${MASS}.root

    # Initial fit
    combineTool.py -M Impacts \
        -d workspace_tchiwg_${MASS}.root \
        -m ${MASS} -t -1 \
        --doInitialFit \
        --parallel 10 \
        --cminDefaultMinimizerStrategy 0 \
        --rMin -20 --rMax 20 \
        --robustFit 1 \
        --X-rtd MINIMIZER_analytic \
        --X-rtd FAST_VERTICAL_MORPH

    # Per-nuisance fits
    combineTool.py -M Impacts \
        -d workspace_tchiwg_${MASS}.root \
        -m ${MASS} -t -1 \
        --doFits \
        --parallel 10 \
        --cminDefaultMinimizerStrategy 0 \
        --rMin -20 --rMax 20 \
        --robustFit 1 \
        --X-rtd MINIMIZER_analytic \
        --X-rtd FAST_VERTICAL_MORPH

    # Collect impacts
    combineTool.py -M Impacts \
        -d workspace_tchiwg_${MASS}.root \
        -m ${MASS} \
        -o impacts_tchiwg_${MASS}.json \
        --parallel 10

    # Plot impacts
    plotImpacts.py \
        -i impacts_tchiwg_${MASS}.json \
        -o impacts_tchiwg_${MASS}_${year}

    cd ..
done
