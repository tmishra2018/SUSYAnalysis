#!/bin/bash

#Blinded:
years=("2016preVFP" "2016postVFP" "2017" "2018")
for year in "${years[@]}"; do
    echo "Processing year: $year"
    mkdir -p Signal_t5wg_${year}
    cd Signal_t5wg_${year} || { echo "Directory Signal_t5wg_${year} not found!"; exit 1; }
    rm -f *

    cp /eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_${year}_new/counting_t5wg_1700_1100.txt .
    text2workspace.py counting_t5wg_1700_1100.txt -o workspace_t5wg_1700_1100.root

    rm -f higgsCombine_*.root

    combineTool.py -M Impacts -d workspace_t5wg_1700_1100.root -m 1700 -t -1 \
        --doInitialFit --parallel 10 --cminDefaultMinimizerStrategy 0 \
        --rMin -20 --rMax 20 --robustFit 1 \
        --X-rtd MINIMIZER_analytic --X-rtd FAST_VERTICAL_MORPH

    combineTool.py -M Impacts -d workspace_t5wg_1700_1100.root -m 1700 -t -1 \
        --doFits --parallel 10 --cminDefaultMinimizerStrategy 0 \
        --rMin -20 --rMax 20 --robustFit 1 \
        --X-rtd MINIMIZER_analytic --X-rtd FAST_VERTICAL_MORPH

    combineTool.py -M Impacts -d workspace_t5wg_1700_1100.root -m 1700 -t -1 \
        -o impacts_t5wg_1700_1100.json --parallel 10

    plotImpacts.py -i impacts_t5wg_1700_1100.json -o impacts_t5wg_1700_1100_${year}

    INPUT_PDF=impacts_t5wg_1700_1100_${year}.pdf

    if command -v pdftk &>/dev/null; then
        echo "Splitting with pdftk..."
        pdftk ${INPUT_PDF} cat 1    output impacts_t5wg_1700_1100_${year}_page1.pdf

    elif command -v gs &>/dev/null; then
        echo "Splitting with ghostscript..."
        gs -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dFirstPage=1 -dLastPage=1 \
            -sOutputFile=impacts_t5wg_1700_1100_${year}_page1.pdf ${INPUT_PDF}
    else
        echo "⚠️  Neither pdftk nor gs found — install one to split PDFs"
    fi
    cd ..
done

#unblinded:

#text2workspace.py Signal_2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_600_MH4_200/xtt_cmb_1_13TeV_201718.txt -o Signal_2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_600_MH4_200/cmb_workspace_600_201718.root
#combineTool.py -M Impacts -d Signal_2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_600_MH4_200/cmb_workspace_600_201718.root -m 125  --doInitialFit --parallel 10 --cminDefaultMinimizerStrategy 0 --rMin -20 --rMax 20 --robustFit 1 --X-rtd MINIMIZER_analytic --X-rtd FAST_VERTICAL_MORPH
#combineTool.py -M Impacts -d Signal_2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_600_MH4_200/cmb_workspace_600_201718.root -m 125  --doFits --parallel 10 --cminDefaultMinimizerStrategy 0 --rMin -20 --rMax 20 --robustFit 1 --X-rtd MINIMIZER_analytic --X-rtd FAST_VERTICAL_MORPH
#combineTool.py -M Impacts -d Signal_2HDMa_gg_sinp_0p35_tanb_1p0_mXd_10_MH3_600_MH4_200/cmb_workspace_600_201718.root -m 125 -o impacts_cmb_100_201718_signal0.json --parallel 10
#plotImpacts.py -i impacts_cmb_100_201718_signal0.json -o impacts_cmb_100_201718_signal0
