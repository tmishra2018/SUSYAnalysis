#!/bin/bash

# Signal region plots. 

dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

NBIN=18
METbin1=200
METbin2=400
HTbin1=100
HTbin2=400
PHOETbin=200

rm -f binConfig.txt
echo 'NBIN' $NBIN >> binConfig.txt
echo 'METbin1' $METbin1 >> binConfig.txt
echo 'METbin2' $METbin2 >> binConfig.txt
echo 'HTbin1' $HTbin1 >> binConfig.txt
echo 'HTbin2' $HTbin2 >> binConfig.txt
echo 'PHOETbin' $PHOETbin >> binConfig.txt

RunYears=("2016" "2017" "2018")

for RunYear in "${RunYears[@]}"; do
    if [ "$RunYear" -eq 2016 ]; then
        preVFPs=("0" "1") 
    else
        preVFPs=("0")  
    fi

    for preVFP in "${preVFPs[@]}"; do
        if [ "$RunYear" -eq 2016 ]; then
            if [ "$preVFP" -eq 1 ]; then
                VFP_string="preVFP"
            else
                VFP_string="postVFP"
            fi
        else
            VFP_string=""
        fi

        echo "Running for RunYear=$RunYear and preVFP=$preVFP"

        # ---- First Signal Configuration (eg) ----
        ch=1
        anatype=3
        lmt=100
        hmt=-1
        lmet=120
        hmet=-1
        iso=4
        lpt=0
        hpt=-1

        rm -f SigConfig.txt
        echo 'ichannel' $ch >> SigConfig.txt
        echo 'anatype' $anatype >> SigConfig.txt
        echo 'lowMt' $lmt >> SigConfig.txt
        echo 'highMt' $hmt >> SigConfig.txt
        echo 'lowMET' $lmet >> SigConfig.txt
        echo 'highMET' $hmet >> SigConfig.txt
        echo 'lowPt' $lpt >> SigConfig.txt
        echo 'highPt' $hpt >> SigConfig.txt
        echo 'lepIso' $iso >> SigConfig.txt
        echo 'RunYear' $RunYear >> SigConfig.txt
	echo 'preVFP' $preVFP >> SigConfig.txt

        # Run the first set of background predictions (eg)
        if [ "$RunYear" -eq 2016 ]; then
            root -l -q pred_VGBkg.C++ > "logs/VG_eg_${RunYear}${VFP_string}.log"
            root -l -q pred_eleBkg.C++ > "logs/ele_eg_${RunYear}${VFP_string}.log"
            root -l -q pred_jetBkg.C++ > "logs/jet_eg_${RunYear}${VFP_string}.log"
            root -l -q pred_qcdBkg.C++ > "logs/qcd_eg_${RunYear}${VFP_string}.log"
            root -l -q pred_rareBkg.C++ > "logs/rare_eg_${RunYear}${VFP_string}.log"
      	    root -l -q pred_sig.C++ > "logs/data_eg_${RunYear}${VFP_string}.log"
        else
            root -l -q pred_VGBkg.C++ > "logs/VG_eg_${RunYear}.log"
            root -l -q pred_eleBkg.C++ > "logs/ele_eg_${RunYear}.log"
            root -l -q pred_jetBkg.C++ > "logs/jet_eg_${RunYear}.log"
            root -l -q pred_qcdBkg.C++ > "logs/qcd_eg_${RunYear}.log"
            root -l -q pred_rareBkg.C++ > "logs/rare_eg_${RunYear}.log"
            root -l -q pred_sig.C++ > "logs/data_eg_${RunYear}.log"
        fi

        # ---- Second Signal Configuration (mg) ----
        ch=2
        anatype=3
        lmt=100
        hmt=-1
        lmet=120
        hmet=-1
        iso=4
        lpt=0
        hpt=-1

        rm -f SigConfig.txt
        echo 'ichannel' $ch >> SigConfig.txt
        echo 'anatype' $anatype >> SigConfig.txt
        echo 'lowMt' $lmt >> SigConfig.txt
        echo 'highMt' $hmt >> SigConfig.txt
        echo 'lowMET' $lmet >> SigConfig.txt
        echo 'highMET' $hmet >> SigConfig.txt
        echo 'lowPt' $lpt >> SigConfig.txt
        echo 'highPt' $hpt >> SigConfig.txt
        echo 'lepIso' $iso >> SigConfig.txt
        echo 'RunYear' $RunYear >> SigConfig.txt
	echo 'preVFP' $preVFP >> SigConfig.txt


        if [ "$RunYear" -eq 2016 ]; then
            root -l -q pred_VGBkg.C++ > "logs/VG_mg_${RunYear}${VFP_string}.log"
            root -l -q pred_eleBkg.C++ > "logs/ele_mg_${RunYear}${VFP_string}.log"
            root -l -q pred_jetBkg.C++ > "logs/jet_mg_${RunYear}${VFP_string}.log"
            root -l -q pred_qcdBkg.C++ > "logs/qcd_mg_${RunYear}${VFP_string}.log"
            root -l -q pred_rareBkg.C++ > "logs/rare_mg_${RunYear}${VFP_string}.log"
      	    root -l -q pred_sig.C++ > "logs/data_mg_${RunYear}${VFP_string}.log"
        else
            root -l -q pred_VGBkg.C++ > "logs/VG_mg_${RunYear}.log"
            root -l -q pred_eleBkg.C++ > "logs/ele_mg_${RunYear}.log"
            root -l -q pred_jetBkg.C++ > "logs/jet_mg_${RunYear}.log"
            root -l -q pred_qcdBkg.C++ > "logs/qcd_mg_${RunYear}.log"
            root -l -q pred_rareBkg.C++ > "logs/rare_mg_${RunYear}.log"
            root -l -q pred_sig.C++ > "logs/data_mg_${RunYear}.log"
        fi
 	
	root -l -q Signal-samples/analysis_TChiWG.C++  # Run it here as it takes the values from SigConfig.txt

        if [ "$RunYear" -eq 2016 ]; then
            rm "logs/eventcount_${RunYear}${VFP_string}.txt"
            root -l -q "plot_eventct.C+($NBIN)" >> "logs/eventcount_${RunYear}${VFP_string}.txt"
        else
            rm "logs/eventcount_${RunYear}.txt"
            root -l -q "plot_eventct.C+($NBIN)" >> "logs/eventcount_${RunYear}.txt"
        fi

# ======================================================================================================	
	
        python3 createDataCard.py ${NBIN} ${RunYear} ${preVFP}                  # takes input from plot_eventct.C
        
	if [ "$RunYear" -eq 2016 ]; then
		rm /eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_${RunYear}${VFP_string}/*
		rm /eos/uscms/store/user/tmishra/CombinedLimit/T6WG/cards/cards_${RunYear}${VFP_string}/*
		rm /eos/uscms/store/user/tmishra/CombinedLimit/TChiWG/cards/cards_${RunYear}${VFP_string}/*
		mkdir -p /eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_${RunYear}${VFP_string}
		mkdir -p /eos/uscms/store/user/tmishra/CombinedLimit/T6WG/cards/cards_${RunYear}${VFP_string}
		mkdir -p /eos/uscms/store/user/tmishra/CombinedLimit/TChiWG/cards/cards_${RunYear}${VFP_string}
	else
		rm /eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_${RunYear}/*
		rm /eos/uscms/store/user/tmishra/CombinedLimit/T6WG/cards/cards_${RunYear}/*
		rm /eos/uscms/store/user/tmishra/CombinedLimit/TChiWG/cards/cards_${RunYear}/*
		mkdir -p /eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_${RunYear}             
		mkdir -p /eos/uscms/store/user/tmishra/CombinedLimit/T6WG/cards/cards_${RunYear}             
		mkdir -p /eos/uscms/store/user/tmishra/CombinedLimit/TChiWG/cards/cards_${RunYear}             
	fi
        python3 writeT5WGcard.py ${NBIN} ${RunYear} ${preVFP}                   # takes input from Signal-samples/analysis_TChiWG.C      
        python3 writeT6WGcard.py ${NBIN} ${RunYear} ${preVFP}                   
        python3 writeTChiWGcard.py ${NBIN} ${RunYear} ${preVFP}                     

    done
done

echo "All tasks completed successfully. do the limit plots.. "
cd  /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/HiggsAnalysis/CombinedLimit/tools
#cmsenv
#bash MakeImpactPlot.sh
