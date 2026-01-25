#!/bin/bash

export eg=true   
export mg=true  

output_dir="/eos/uscms/store/user/tmishra/VGamma"

declare -a RunYears=("678")
declare -A VFPs=( ["678"]="noVFP" )

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
      eg_output_file="${output_dir}/VGamma_scalefactor_eg_${RunYear}.txt"
      echo "Removing old file: $eg_output_file"
      rm -f "$eg_output_file"

      echo "Running Fitfractioneg for RunYear=$RunYear, VFP=$VFP"
      for i in {0..0}; do
        root_command="root -b -q \"Fitfractioneg.C+($i,40,70,0,1000,4, $RunYear, $VFP_flag)\""
        echo "Running: $root_command"
        eval $root_command
      done
    fi

    if [ "$mg" = true ]; then
      mg_output_file="${output_dir}/VGamma_scalefactor_mg_${RunYear}.txt"
      echo "Removing old file: $mg_output_file"
      rm -f "$mg_output_file"

      echo "Running Fitfractionmg for RunYear=$RunYear, VFP=$VFP"
      for i in {0..0}; do
        root_command="root -b -q \"Fitfractionmg.C+($i,40,70,0,1000,4, $RunYear, $VFP_flag)\""
        echo "Running: $root_command"
        eval $root_command
      done
    fi
  done
done
echo "All tasks completed successfully."
