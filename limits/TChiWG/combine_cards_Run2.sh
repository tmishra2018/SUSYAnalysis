#!/bin/bash
dir_2016_preVFP="/eos/uscms/store/user/tmishra/CombinedLimit/TChiWG/cards/cards_2016preVFP_splusb"
dir_2016_postVFP="/eos/uscms/store/user/tmishra/CombinedLimit/TChiWG/cards/cards_2016postVFP_splusb"
dir_2017="/eos/uscms/store/user/tmishra/CombinedLimit/TChiWG/cards/cards_2017_splusb"
dir_2018="/eos/uscms/store/user/tmishra/CombinedLimit/TChiWG/cards/cards_2018_splusb"
dir_run2="/eos/uscms/store/user/tmishra/CombinedLimit/TChiWG/cards/cards_Run2_splusb"

# Create Run2 directory
mkdir -p "$dir_run2"
extract_m1() {
    local dir=$1
    local output_file=$2

    find "$dir" -type f -name 'counting_tchiwg_*.txt' | \
    sed -E 's|.*counting_tchiwg_([0-9]+)\.txt|\1|' | \
    sort -u > "$output_file"
}
extract_m1 "$dir_2016_preVFP" "m1_2016_preVFP.txt"
extract_m1 "$dir_2016_postVFP" "m1_2016_postVFP.txt"
extract_m1 "$dir_2017" "m1_2017.txt"
extract_m1 "$dir_2018" "m1_2018.txt"
combine_cards() {
    local m1=$1

    if [[ -f "$dir_2016_preVFP/counting_tchiwg_${m1}.txt" && \
          -f "$dir_2016_postVFP/counting_tchiwg_${m1}.txt" && \
          -f "$dir_2017/counting_tchiwg_${m1}.txt" && \
          -f "$dir_2018/counting_tchiwg_${m1}.txt" ]]; then
        combineCards.py \
            "$dir_2016_preVFP/counting_tchiwg_${m1}.txt" \
            "$dir_2016_postVFP/counting_tchiwg_${m1}.txt" \
            "$dir_2017/counting_tchiwg_${m1}.txt" \
            "$dir_2018/counting_tchiwg_${m1}.txt" \
            > "$dir_run2/counting_tchiwg_${m1}.txt"
    fi
}
cat m1_2016_preVFP.txt m1_2016_postVFP.txt m1_2017.txt m1_2018.txt | sort -u | while read -r m1; do
    combine_cards "$m1"
done
rm m1_2016_preVFP.txt m1_2016_postVFP.txt m1_2017.txt m1_2018.txt
echo "Process completed."
