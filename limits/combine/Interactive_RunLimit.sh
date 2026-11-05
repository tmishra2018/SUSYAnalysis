#!/bin/bash

cardpath="/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_Run2"
limitpath="/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/limits/limits_Run2"

rm /eos/uscms/store/user/tmishra/CombinedLimit/T5WG/limits/limits_Run2/*
for file in ${cardpath}/counting_t5Wg_*.txt; do
    filename=$(basename "$file")  # Get the filename only
    gluino=$(echo "$filename" | cut -d'_' -f3)  # Extract 3rd field (gluino)
    chargino=$(echo "$filename" | cut -d'_' -f4 | cut -d'.' -f1) 
    
    echo "Processing gluino=$gluino, chargino=$chargino"

    combine -M AsymptoticLimits ${cardpath}/counting_t5Wg_${gluino}_${chargino}.txt \
        --run blind --expectSignal=1 -n t5wg | tee ${limitpath}/limit_${gluino}_${chargino}.txt
done

echo "All jobs completed."
