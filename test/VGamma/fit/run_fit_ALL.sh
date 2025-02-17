#!/bin/bash

export eg=true   
export mg=true  

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

      echo "Running Fitfractioneg for RunYear=$RunYear, VFP=$VFP"
      for i in {0..999}; do
        root_command="root -b -q \"Fitfractioneg.C+($i,40,70,0,1000,4, $RunYear, $VFP_flag)\""
        echo "Running: $root_command"
        eval $root_command
      done
    fi

    if [ "$mg" = true ]; then
      if [ "$RunYear" -eq 2017 ] || [ "$RunYear" -eq 2018 ]; then
    	mg_output_file="${output_dir}/VGamma_scalefactor_mg_${RunYear}.txt"
      else
    	mg_output_file="${output_dir}/VGamma_scalefactor_mg_${RunYear}${VFP}.txt"
      fi	
      echo "Removing old file: $mg_output_file"
      rm -f "$mg_output_file"

      echo "Running Fitfractionmg for RunYear=$RunYear, VFP=$VFP"
      for i in {0..999}; do
        root_command="root -b -q \"Fitfractionmg.C+($i,40,70,0,1000,4, $RunYear, $VFP_flag)\""
        echo "Running: $root_command"
        eval $root_command
      done
    fi
  done
done

echo "All tasks completed successfully."
