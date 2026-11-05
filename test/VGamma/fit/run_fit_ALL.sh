#!/bin/bash

export eg=true
export mg=true

# FIT_MODE=fast (default): run_toys.C, 100 toys per ROOT session (minutes)
# FIT_MODE=slow           : one ROOT process per toy (hours), kept for checks:  FIT_MODE=slow bash run_fit_ALL.sh
FIT_MODE=${FIT_MODE:-fast}
# inclusive fit (MET 40-70, lepton pT 0-1000): toys 0..999, toy 0 = nominal
fit_toys(){ # ch RunYear VFP_flag
  if [ "$FIT_MODE" = slow ]; then
    for i in {0..999}; do root -b -q "Fitfraction$1.C+($i,40,70,0,1000,4,$2,$3)"; done
  else
    for b in {0..9}; do root -b -q "run_toys.C(\"$1\",$((b*100)),$((b*100+99)),$2,$3,-1,-1)"; done
  fi
}

output_dir="/eos/uscms/store/user/tmishra/VGamma"
#/eos/uscms/store/user/tmishra/VGamma/VGamma_scalefactor_*.txt

declare -a RunYears=("2016" "2017" "2018")
declare -A VFPs=( ["2016"]="preVFP postVFP" ["2017"]="noVFP" ["2018"]="noVFP" )

if [ -z "$eg" ] && [ -z "$mg" ]; then
  echo "Please set either eg=true or mg=true before running this script."
  return 1
fi

for RunYear in "${RunYears[@]}"; do
  for VFP in ${VFPs[$RunYear]}; do
    if [ "$VFP" = "preVFP" ]; then
      VFP_flag=1
    else
      VFP_flag=0
    fi

    if [ "$eg" = true ]; then
      if [ "$RunYear" -eq 2017 ] || [ "$RunYear" -eq 2018 ]; then
    	eg_output_file="${output_dir}/VGamma_scalefactor_eg_${RunYear}.txt"
      else
    	eg_output_file="${output_dir}/VGamma_scalefactor_eg_${RunYear}${VFP}.txt"
      fi	
      echo "Removing old file: $eg_output_file"
      rm -f "$eg_output_file"

      echo "Running Fitfractioneg for RunYear=$RunYear, VFP=$VFP ($FIT_MODE)"
      fit_toys eg $RunYear $VFP_flag
    fi

    if [ "$mg" = true ]; then
      if [ "$RunYear" -eq 2017 ] || [ "$RunYear" -eq 2018 ]; then
    	mg_output_file="${output_dir}/VGamma_scalefactor_mg_${RunYear}.txt"
      else
    	mg_output_file="${output_dir}/VGamma_scalefactor_mg_${RunYear}${VFP}.txt"
      fi	
      echo "Removing old file: $mg_output_file"
      rm -f "$mg_output_file"

      echo "Running Fitfractionmg for RunYear=$RunYear, VFP=$VFP ($FIT_MODE)"
      fit_toys mg $RunYear $VFP_flag
    fi
  done
done
echo "All tasks completed successfully."
