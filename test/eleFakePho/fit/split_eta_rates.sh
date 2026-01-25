#!/bin/bash
years=("2016preVFP" "2016postVFP" "2017" "2018")
base_dir="/eos/uscms/store/user/tmishra/elefakepho/DATAResult"

for year in "${years[@]}"; do
    infile="${base_dir}${year}/result_eta_dependence_Data.txt"
    eb_outfile="${base_dir}${year}/result_eta_dependence_Data.txt"
    ee_outfile="${base_dir}${year}/result_eta_dependence_Data_EE.txt"

    # Extract everything between the braces { ... }
    content=$(grep -oP '\{.*\}' "$infile" | tr -d '{}')
    
    # Convert to array
    IFS=',' read -ra values <<< "$content"
    
    # Total number of points
    total=${#values[@]}
    echo "Year: $year | Total points in original file: $total"

    # If fewer than 48 points, skip this file
    if [ "$total" -lt 48 ]; then
        echo "Warning: $infile has less than 48 points. Skipping this file."
        echo "-----------------------------------------"
        continue
    fi

    # Enforce 29 points for EB, 19 for EE
    eb_count=29
    ee_count=$((total - eb_count))
    echo "Points assigned: EB = $eb_count | EE = $ee_count"

    # Write EB file (first eb_count points)
    echo -n "double etaRatesEB_${year}[] = {" > "$eb_outfile"
    for ((i=0; i<eb_count; i++)); do
        if [ $i -ne $((eb_count-1)) ]; then
            echo -n "${values[$i]}," >> "$eb_outfile"
        else
            echo -n "${values[$i]}" >> "$eb_outfile"
        fi
    done
    echo "};" >> "$eb_outfile"
    echo "Written EB: $eb_outfile"

    # Write EE file (remaining points)
    echo -n "double etaRatesEE_${year}[] = {" > "$ee_outfile"
    for ((i=eb_count; i<total; i++)); do
        if [ $i -ne $((total-1)) ]; then
            echo -n "${values[$i]}," >> "$ee_outfile"
        else
            echo -n "${values[$i]}" >> "$ee_outfile"
        fi
    done
    echo "};" >> "$ee_outfile"
    echo "Written EE: $ee_outfile"

    echo "-----------------------------------------"
done
