#!/bin/bash

############################################
# Configuration
############################################

export eg=true
export mg=true

doFullRange=true
doPhoScan=false   # photon-pT bins (the default binning of the V+gamma / fake-lepton factors)
doHTScan=false

doMETScan=false
doPTScan=false

# FIT_MODE=fast (default): run_toys.C, 100 toys per ROOT session (minutes)
# FIT_MODE=slow           : one ROOT process per toy (hours), kept for checks:  FIT_MODE=slow bash run_fit_ALL_HT.sh
FIT_MODE=${FIT_MODE:-fast}
NTOYS=${NTOYS:-1000}      # toys 0..NTOYS-1; toy 0 is the nominal fit

output_dir="/eos/uscms/store/user/tmishra/VGamma"

MAX_PARALLEL_JOBS=30

############################################
# Years and VFP
############################################

RunYears=("678")
#RunYears=("2016" "2017" "2018")
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

# MET scan (inside 40-70)
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

# photon-pT scan (MET, PT and HT fixed to full range)
PHO_SCAN_BINS=(
"0 100"
"100 200"
"200 5000"
)

# HT scan (MET and PT fixed to full range)
HT_SCAN_BINS=(
"0 100"
"100 400"
"400 3000"
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
# Toy runner (sequential per fit)
############################################
run_toys() {

  macro=$1
  lmet=$2
  hmet=$3
  lpt=$4
  hpt=$5
  RunYear=$6
  VFP_flag=$7
  lht=$8
  hht=$9
  lpho=${10:--1}
  hpho=${11:--1}
  ch=${macro#Fitfraction}

  echo "Running $macro ($FIT_MODE) | Year=$RunYear VFP=$VFP_flag MET=$lmet-$hmet PT=$lpt-$hpt HT=$lht-$hht PHO=$lpho-$hpho"

  if [ "$FIT_MODE" = slow ]; then
    for ((i=0; i<NTOYS; i++)); do
      root -b -q "${macro}.C+($i,$lmet,$hmet,$lpt,$hpt,4,$RunYear,$VFP_flag,$lht,$hht,1,$lpho,$hpho)"
    done
  else
    for ((b=0; b<NTOYS; b+=100)); do
      last=$((b+99 < NTOYS-1 ? b+99 : NTOYS-1))
      root -b -q "run_toys.C(\"$ch\",$b,$last,$RunYear,$VFP_flag,$lpho,$hpho,$lmet,$hmet,$lpt,$hpt,$lht,$hht)"
    done
  fi
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
    		if [ "$doFullRange" = true ]; then rm -f "$eg_output_file"; fi
    		if [ "$doMETScan" = true ] || [ "$doPTScan" = true ]; then rm -f ${output_dir}/VGamma_scalefactor_eg_${RunYear}_MET*_lep*.txt 2>/dev/null; fi
    		if [ "$doHTScan" = true ]; then rm -f ${output_dir}/VGamma_scalefactor_eg_${RunYear}_HT*.txt 2>/dev/null; fi

  	else
    		eg_output_file="${output_dir}/VGamma_scalefactor_eg_${RunYear}${VFP}.txt"
    		if [ "$doFullRange" = true ]; then rm -f "$eg_output_file"; fi
    		if [ "$doMETScan" = true ] || [ "$doPTScan" = true ]; then rm -f ${output_dir}/VGamma_scalefactor_eg_${RunYear}${VFP}_MET*_lep*.txt 2>/dev/null; fi
    		if [ "$doHTScan" = true ]; then rm -f ${output_dir}/VGamma_scalefactor_eg_${RunYear}${VFP}_HT*.txt 2>/dev/null; fi
  	fi
  fi


  if [ "$mg" = true ]; then

  	if [ "$RunYear" -eq 2017 ] || [ "$RunYear" -eq 2018 ]; then
    		mg_output_file="${output_dir}/VGamma_scalefactor_mg_${RunYear}.txt"
    		if [ "$doFullRange" = true ]; then rm -f "$mg_output_file"; fi
    		if [ "$doMETScan" = true ] || [ "$doPTScan" = true ]; then rm -f ${output_dir}/VGamma_scalefactor_mg_${RunYear}_MET*_lep*.txt 2>/dev/null; fi
    		if [ "$doHTScan" = true ]; then rm -f ${output_dir}/VGamma_scalefactor_mg_${RunYear}_HT*.txt 2>/dev/null; fi
  	else
    		mg_output_file="${output_dir}/VGamma_scalefactor_mg_${RunYear}${VFP}.txt"
    		if [ "$doFullRange" = true ]; then rm -f "$mg_output_file"; fi
    		if [ "$doMETScan" = true ] || [ "$doPTScan" = true ]; then rm -f ${output_dir}/VGamma_scalefactor_mg_${RunYear}${VFP}_MET*_lep*.txt 2>/dev/null; fi
    		if [ "$doHTScan" = true ]; then rm -f ${output_dir}/VGamma_scalefactor_mg_${RunYear}${VFP}_HT*.txt 2>/dev/null; fi
  	fi
  fi

  if [ "$doFullRange" = true ]; then
    if [ "$eg" = true ]; then
      wait_for_slot
      run_toys Fitfractioneg $FULL_MET_LOW $FULL_MET_HIGH \
               $FULL_PT_LOW $FULL_PT_HIGH $RunYear $VFP_flag -1 -1 &
    fi

    if [ "$mg" = true ]; then
      wait_for_slot
      run_toys Fitfractionmg $FULL_MET_LOW $FULL_MET_HIGH \
               $FULL_PT_LOW $FULL_PT_HIGH $RunYear $VFP_flag -1 -1 &
    fi
  fi

  ##########################################
  # PHOTON-PT SCAN (MET, PT, HT FULL)
  ##########################################

  if [ "$doPhoScan" = true ]; then
  rm -f ${output_dir}/VGamma_scalefactor_{eg,mg}_${RunYear}${VFP}_pho*.txt 2>/dev/null
  for phobin in "${PHO_SCAN_BINS[@]}"; do

    lpho=$(echo $phobin | cut -d' ' -f1)
    hpho=$(echo $phobin | cut -d' ' -f2)

    if [ "$eg" = true ]; then
      wait_for_slot
      run_toys Fitfractioneg 40 70 0 1000 $RunYear $VFP_flag -1 -1 $lpho $hpho &
    fi

    if [ "$mg" = true ]; then
      wait_for_slot
      run_toys Fitfractionmg 40 70 0 1000 $RunYear $VFP_flag -1 -1 $lpho $hpho &
    fi

  done
  fi

  ##########################################
  # MET SCAN (PT FULL)
  ##########################################

  if [ "$doMETScan" = true ]; then
  for metbin in "${MET_SCAN_BINS[@]}"; do

    lmet=$(echo $metbin | cut -d' ' -f1)
    hmet=$(echo $metbin | cut -d' ' -f2)

    if [ "$eg" = true ]; then
      wait_for_slot
      run_toys Fitfractioneg $lmet $hmet 0 1000 $RunYear $VFP_flag -1 -1 &
    fi

    if [ "$mg" = true ]; then
      wait_for_slot
      run_toys Fitfractionmg $lmet $hmet 0 1000 $RunYear $VFP_flag -1 -1 &
    fi

  done
  fi

  ##########################################
  # PT SCAN (MET FULL)
  ##########################################

  if [ "$doPTScan" = true ]; then
  for ptbin in "${PT_SCAN_BINS[@]}"; do

    lpt=$(echo $ptbin | cut -d' ' -f1)
    hpt=$(echo $ptbin | cut -d' ' -f2)

    if [ "$eg" = true ]; then
      wait_for_slot
      run_toys Fitfractioneg 40 70 $lpt $hpt $RunYear $VFP_flag -1 -1 &
    fi

    if [ "$mg" = true ]; then
      wait_for_slot
      run_toys Fitfractionmg 40 70 $lpt $hpt $RunYear $VFP_flag -1 -1 &
    fi

  done
  fi

  ##########################################
  # HT SCAN (MET AND PT FULL)
  ##########################################

  if [ "$doHTScan" = true ]; then
  for htbin in "${HT_SCAN_BINS[@]}"; do

    lht=$(echo $htbin | cut -d' ' -f1)
    hht=$(echo $htbin | cut -d' ' -f2)

    if [ "$eg" = true ]; then
      wait_for_slot
      run_toys Fitfractioneg $FULL_MET_LOW $FULL_MET_HIGH \
               $FULL_PT_LOW $FULL_PT_HIGH $RunYear $VFP_flag $lht $hht &
    fi

    if [ "$mg" = true ]; then
      wait_for_slot
      run_toys Fitfractionmg $FULL_MET_LOW $FULL_MET_HIGH \
               $FULL_PT_LOW $FULL_PT_HIGH $RunYear $VFP_flag $lht $hht &
    fi

  done
  fi
}

############################################
# Main Loop
############################################

# compile once, so that the parallel ROOT sessions do not all try to rebuild the macros
for m in Fitfractioneg Fitfractionmg; do root -l -b -q -e ".L ${m}.C+" > /dev/null 2>&1; done

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
