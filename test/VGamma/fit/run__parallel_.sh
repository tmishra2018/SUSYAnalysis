#!/bin/bash

export eg=true
export mg=true

output_dir="/eos/uscms/store/user/tmishra/VGamma"

MAX_PARALLEL_JOBS=30

RunYears=("2016" "2017" "2018")
VFPs_2016=("preVFP" "postVFP")
VFPs_other=("")

FULL_MET_LOW=40
FULL_MET_HIGH=70
FULL_PT_LOW=0
FULL_PT_HIGH=1000


wait_for_slot() {
  while [ $(jobs -r | wc -l) -ge $MAX_PARALLEL_JOBS ]; do
    sleep 2
  done
}

run_toys() {

  macro=$1
  lmet=$2
  hmet=$3
  lpt=$4
  hpt=$5
  RunYear=$6
  VFP_flag=$7

  echo "Running $macro | Year=$RunYear VFP=$VFP_flag MET=$lmet-$hmet PT=$lpt-$hpt"

  for i in {0..999}; do 
    root -b -q "${macro}.C+($i,$lmet,$hmet,$lpt,$hpt,4,$RunYear,$VFP_flag)"
  done
}

launch_jobs() {

  RunYear=$1
  VFP=$2

  if [ "$VFP" = "preVFP" ]; then
    VFP_flag=1
  else
    VFP_flag=0
  fi

  if [ "$eg" = true ]; then
    wait_for_slot
    run_toys Fitfractioneg $FULL_MET_LOW $FULL_MET_HIGH \
             $FULL_PT_LOW $FULL_PT_HIGH $RunYear $VFP_flag &
  fi

  if [ "$mg" = true ]; then
    wait_for_slot
    run_toys Fitfractionmg $FULL_MET_LOW $FULL_MET_HIGH \
             $FULL_PT_LOW $FULL_PT_HIGH $RunYear $VFP_flag &
  fi

}

############################################
# Main Loop
############################################

for RunYear in "${RunYears[@]}"; do

  if [ "$RunYear" = "2016" ]; then
    VFP_LIST=("${VFPs_2016[@]}")
  else
    VFP_LIST=("${VFPs_other[@]}")
  fi

  for VFP in "${VFP_LIST[@]}"; do
    launch_jobs $RunYear $VFP
  done

done

############################################
# Wait for all jobs
############################################

wait
echo "All parallel jobs finished successfully."
