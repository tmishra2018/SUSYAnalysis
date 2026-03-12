#!/usr/bin/env bash
set -euo pipefail

if [[ -f "analysis_eg.C" && -f "analysis_mg.C" ]]; then
  WORKDIR="."
else
  echo "ERROR: cannot find analysis_eg.C and analysis_mg.C" >&2
  exit 1
fi

cd "$WORKDIR"

echo 'start'

MAX_JOBS=${MAX_JOBS:-$(nproc)}

g++ $(root-config --cflags) ../../lib/libAnaClasses.so analysis_eg.C -o analysis_eg.exe $(root-config --libs)
g++ $(root-config --cflags) ../../lib/libAnaClasses.so analysis_mg.C -o analysis_mg.exe $(root-config --libs)

TASKS=(

#  "analysis_eg.exe 2018 A"
#  "analysis_eg.exe 2018 B"
#  "analysis_eg.exe 2018 C"
#  "analysis_eg.exe 2018 D"

#  "analysis_eg.exe 2017 B"
#  "analysis_eg.exe 2017 C"
#  "analysis_eg.exe 2017 D"
#  "analysis_eg.exe 2017 E"
#  "analysis_eg.exe 2017 F"

#  "analysis_eg.exe 2016 B-ver1"
#  "analysis_eg.exe 2016 B-ver2"
#  "analysis_eg.exe 2016 C"
#  "analysis_eg.exe 2016 D"
#  "analysis_eg.exe 2016 E"
#  "analysis_eg.exe 2016 F_preVFP"
#  "analysis_eg.exe 2016 F_postVFP"
#  "analysis_eg.exe 2016 G"
#  "analysis_eg.exe 2016 H"

  "analysis_mg.exe 2017 D"
  "analysis_mg.exe 2017 E"
  "analysis_mg.exe 2017 F"
  "analysis_mg.exe 2017 B"
  "analysis_mg.exe 2017 C"

  "analysis_mg.exe 2016 B-ver1"
  "analysis_mg.exe 2016 B-ver2"
  "analysis_mg.exe 2016 C"
  "analysis_mg.exe 2016 D"
  "analysis_mg.exe 2016 E"
  "analysis_mg.exe 2016 F_preVFP"
  "analysis_mg.exe 2016 F_postVFP"
  "analysis_mg.exe 2016 G"
  "analysis_mg.exe 2016 H"

  "analysis_mg.exe 2018 A"
  "analysis_mg.exe 2018 B"
  "analysis_mg.exe 2018 C"
  "analysis_mg.exe 2018 D"
)

run_task() {
  local cmd="$1"
  echo "[RUN] ./${cmd}"
  ./${cmd}
}

active_jobs=0
fail_count=0

for task in "${TASKS[@]}"; do
  run_task "$task" &
  ((active_jobs+=1))

  if (( active_jobs >= MAX_JOBS )); then
    if ! wait -n; then
      ((fail_count+=1))
    fi
    ((active_jobs-=1))
  fi
done

# Wait remaining jobs
while (( active_jobs > 0 )); do
  if ! wait -n; then
    ((fail_count+=1))
  fi
  ((active_jobs-=1))
done

if (( fail_count > 0 )); then
  echo "Completed with ${fail_count} failed job(s)." >&2
  exit 1
fi

echo "All jobs completed successfully using MAX_JOBS=${MAX_JOBS}."
