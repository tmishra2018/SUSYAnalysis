#!/bin/bash

############################################
# Configuration
############################################

export eg=true
export mg=true

output_dir="/eos/uscms/store/user/tmishra/VGamma"

MAX_PARALLEL_JOBS=30

############################################
# Years and VFP
############################################

RunYears=("2016" "2017" "2018")
VFPs_2016=("preVFP" "postVFP")
VFPs_other=("")

############################################
# Binning Definitions
############################################

# Full range
FULL_MET_LOW=40
FULL_MET_HIGH=70
FULL_PT_LOW=0
FULL_PT_HIGH=1000

# MET scan (inside 40–70)
MET_SCAN_BINS=(
"40 50"
"50 60"
"60 70"
)

# PT scan
PT_SCAN_BINS=(
"0 50"
"50 70"
"70 100"
"100 1000"
)

############################################
# Parallel job controller
############################################
wait_for_slot() {
  while [ $(jobs -r | wc -l) -ge $MAX_PARALLEL_JOBS ]; do
    sleep 2
  done
}

############################################
# Toy runner (Sequential)
############################################
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

############################################
# Launch jobs for a given config
############################################
launch_jobs() {

  RunYear=$1
  VFP=$2

  if [ "$VFP" = "preVFP" ]; then
    VFP_flag=1
  else
    VFP_flag=0
  fi

  if [ "$eg" = true ]; then
  	if [ "$RunYear" -eq 2017 ] || [ "$RunYear" -eq 2018 ]; then
    		eg_output_file="${output_dir}/VGamma_scalefactor_eg_${RunYear}.txt"
    		rm -f "$eg_output_file"
    		rm -f ${output_dir}/VGamma_scalefactor_eg_${RunYear}_MET*_lep*.txt 2>/dev/null

  	else
    		eg_output_file="${output_dir}/VGamma_scalefactor_eg_${RunYear}${VFP}.txt"
    		rm -f "$eg_output_file"
    		rm -f ${output_dir}/VGamma_scalefactor_eg_${RunYear}${VFP}_MET*_lep*.txt 2>/dev/null
  	fi
  fi


  if [ "$mg" = true ]; then

  	if [ "$RunYear" -eq 2017 ] || [ "$RunYear" -eq 2018 ]; then
    		mg_output_file="${output_dir}/VGamma_scalefactor_mg_${RunYear}.txt"
    		rm -f "$mg_output_file"
    		rm -f ${output_dir}/VGamma_scalefactor_mg_${RunYear}_MET*_lep*.txt 2>/dev/null
  	else
    		mg_output_file="${output_dir}/VGamma_scalefactor_mg_${RunYear}${VFP}.txt"
    		rm -f "$mg_output_file"
    		rm -f ${output_dir}/VGamma_scalefactor_mg_${RunYear}${VFP}_MET*_lep*.txt 2>/dev/null
  	fi
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

  ##########################################
  # MET SCAN (PT FULL)
  ##########################################

  for metbin in "${MET_SCAN_BINS[@]}"; do

    lmet=$(echo $metbin | cut -d' ' -f1)
    hmet=$(echo $metbin | cut -d' ' -f2)

    if [ "$eg" = true ]; then
      wait_for_slot
      run_toys Fitfractioneg $lmet $hmet 0 1000 $RunYear $VFP_flag &
    fi

    if [ "$mg" = true ]; then
      wait_for_slot
      run_toys Fitfractionmg $lmet $hmet 0 1000 $RunYear $VFP_flag &
    fi

  done

  ##########################################
  # PT SCAN (MET FULL)
  ##########################################

  for ptbin in "${PT_SCAN_BINS[@]}"; do

    lpt=$(echo $ptbin | cut -d' ' -f1)
    hpt=$(echo $ptbin | cut -d' ' -f2)

    if [ "$eg" = true ]; then
      wait_for_slot
      run_toys Fitfractioneg 40 70 $lpt $hpt $RunYear $VFP_flag &
    fi

    if [ "$mg" = true ]; then
      wait_for_slot
      run_toys Fitfractionmg 40 70 $lpt $hpt $RunYear $VFP_flag &
    fi

  done
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
