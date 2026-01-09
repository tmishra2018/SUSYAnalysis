#!/bin/bash
years=("2016preVFP" "2016postVFP" "2017" "2018")
base_dir="/eos/uscms/store/user/tmishra/elefakepho/DATAResult"

for year in "${years[@]}"; do
    infile="${base_dir}${year}/result_eta_dependence_Data.txt"
    outfile="${base_dir}${year}/result_eta_dependence_Data.txt"
    ee_outfile="${base_dir}${year}/result_eta_dependence_Data_EE.txt"

    # Extract everything between the braces { ... }
    content=$(grep -oP '\{.*\}' "$infile" | tr -d '{}')
    
    # Convert to array
    IFS=',' read -ra values <<< "$content"
    
    # Total values check
    total=${#values[@]}
    if [ "$total" -ne 48 ]; then
        echo "Warning: $infile does not contain 48 values. Found $total. Skipping."
        continue
    fi

    # Write new EE file with last 19 values
    echo -n "double etaRatesEE_${year}[] = {" > "$ee_outfile"
    for ((i=29; i<48; i++)); do
        if [ $i -ne 47 ]; then
            echo -n "${values[$i]}," >> "$ee_outfile"
        else
            echo -n "${values[$i]}" >> "$ee_outfile"
        fi
    done
    echo "};" >> "$ee_outfile"
    echo "Written: $ee_outfile"

    # Rewrite original file with only first 29 values
    echo -n "double etaRatesEB_${year}[] = {" > "$outfile"
    for ((i=0; i<29; i++)); do
        if [ $i -ne 28 ]; then
            echo -n "${values[$i]}," >> "$outfile"
        else
            echo -n "${values[$i]}" >> "$outfile"
        fi
    done
    echo "};" >> "$outfile"
    echo "Updated: $outfile"
done
