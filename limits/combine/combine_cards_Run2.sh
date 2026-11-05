#!/bin/bash
# Check for model argument
if [[ $# -ne 1 ]]; then
    echo "Usage: ./combine_cards_Run2.sh [T5WG|T6WG]"
    exit 1
fi
model=$1
model_lower=$(echo "$model" | tr '[:upper:]' '[:lower:]')
# Define directories based on model
dir_2016_preVFP="/eos/uscms/store/user/tmishra/CombinedLimit/${model}/cards/cards_2016preVFP_unblind"
dir_2016_postVFP="/eos/uscms/store/user/tmishra/CombinedLimit/${model}/cards/cards_2016postVFP_unblind"
dir_2017="/eos/uscms/store/user/tmishra/CombinedLimit/${model}/cards/cards_2017_unblind"
dir_2018="/eos/uscms/store/user/tmishra/CombinedLimit/${model}/cards/cards_2018_unblind"
dir_run2="/eos/uscms/store/user/tmishra/CombinedLimit/${model}/cards/cards_Run2_unblind"
# Create Run2 directory
mkdir -p "$dir_run2"
# Function to extract m1 and m2 values from filenames
extract_m1_m2() {
    local dir=$1
    local output_file=$2

    find "$dir" -type f -name "counting_${model_lower}_*.txt" | \
    sed -E "s|.*counting_${model_lower}_([0-9]+)_([0-9]+)\.txt|\1 \2|" | \
    sort -u > "$output_file"
}
# Extract m1 and m2 values for each year
extract_m1_m2 "$dir_2016_preVFP" "m1_m2_2016_preVFP.txt"
extract_m1_m2 "$dir_2016_postVFP" "m1_m2_2016_postVFP.txt"
extract_m1_m2 "$dir_2017" "m1_m2_2017.txt"
extract_m1_m2 "$dir_2018" "m1_m2_2018.txt"

# Function to combine cards for a given m1, m2
combine_cards() {
    local m1=$1
    local m2=$2

    if [[ -f "$dir_2016_preVFP/counting_${model_lower}_${m1}_${m2}.txt" && \
          -f "$dir_2016_postVFP/counting_${model_lower}_${m1}_${m2}.txt" && \
          -f "$dir_2017/counting_${model_lower}_${m1}_${m2}.txt" && \
          -f "$dir_2018/counting_${model_lower}_${m1}_${m2}.txt" ]]; then

        combineCards.py \
            "$dir_2016_preVFP/counting_${model_lower}_${m1}_${m2}.txt" \
            "$dir_2016_postVFP/counting_${model_lower}_${m1}_${m2}.txt" \
            "$dir_2017/counting_${model_lower}_${m1}_${m2}.txt" \
            "$dir_2018/counting_${model_lower}_${m1}_${m2}.txt" \
            > "$dir_run2/counting_${model_lower}_${m1}_${m2}.txt"
    fi
}

# Read each m1 and m2 pair and combine cards if files exist in all directories
cat m1_m2_2016_preVFP.txt m1_m2_2016_postVFP.txt m1_m2_2017.txt m1_m2_2018.txt | sort -u | while read -r m1 m2; do
    combine_cards "$m1" "$m2"
done

# Clean up intermediate files
rm m1_m2_2016_preVFP.txt m1_m2_2016_postVFP.txt m1_m2_2017.txt m1_m2_2018.txt

echo "[$model] Combine cards process completed."
