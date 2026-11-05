#!/bin/bash
# Writes the datacards to EOS (CombinedLimit/<model>/cards). The python scripts work inside test/myResult.
cd "$(dirname "$0")"
years=(2016 2017 2018)
bins=(0 1)
modes=(splusb)
#modes=(bonly splusb unblind)
cards=(writeT5WGcard.py writeT6WGcard.py writeTChiWGcard.py)

mass=18
max_jobs=6   # control parallel jobs

# -----------------------------
# Step 1: createDataCard (serial)
# -----------------------------
echo "=== Running createDataCard (sequential) ==="

for year in "${years[@]}"; do
  for bin in "${bins[@]}"; do

    # Skip bin=1 for non-2016
    if [[ "$year" != "2016" && "$bin" == "1" ]]; then
      continue
    fi

    for mode in "${modes[@]}"; do
      echo "createDataCard: $year $bin $mode"
      python3 createDataCard.py $mass $year $bin $mode
    done
  done
done

# -----------------------------
# Step 2: write*.py (parallel)
# -----------------------------
echo "=== Running write* scripts (parallel) ==="

job_count=0

for year in "${years[@]}"; do
  for bin in "${bins[@]}"; do

    if [[ "$year" != "2016" && "$bin" == "1" ]]; then
      continue
    fi

    for mode in "${modes[@]}"; do
      for card in "${cards[@]}"; do

        echo "Running: $card $year $bin $mode"

        python3 $card $mass $year $bin $mode &

        ((job_count++))

        # Limit parallel jobs
        if (( job_count % max_jobs == 0 )); then
          wait
        fi

      done
    done
  done
done

# Wait for remaining jobs
wait

echo "=== All jobs completed ==="
