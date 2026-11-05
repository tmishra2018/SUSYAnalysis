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

	ch=1
        anatype=3
        lmet=0
        hmet=70
        lmt=0
        hmt=-1
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

        #if [ "$RunYear" -eq 2016 ]; then
        #    root -l -q pred_VGBkg.C++ > "logs/VG_eg_${RunYear}${VFP_string}.log" 
        #    root -l -q pred_eleBkg_pT_nVtx.C++ > "logs/ele_eg_${RunYear}${VFP_string}.log" 
        #    root -l -q pred_jetBkg_pT_nVtx.C++ > "logs/jet_eg_${RunYear}${VFP_string}.log" 
        #    root -l -q pred_qcdBkg.C++ > "logs/qcd_eg_${RunYear}${VFP_string}.log" 
        #    root -l -q pred_rareBkg.C++ > "logs/rare_eg_${RunYear}${VFP_string}.log" 
      	#    root -l -q pred_sig.C++ > "logs/data_eg_${RunYear}${VFP_string}.log" 
        #else
        #    root -l -q pred_VGBkg.C++ > "logs/VG_eg_${RunYear}.log" 
        #    root -l -q pred_eleBkg_pT_nVtx.C++ > "logs/ele_eg_${RunYear}.log" 
        #    root -l -q pred_jetBkg_pT_nVtx.C++ > "logs/jet_eg_${RunYear}.log" 
        #    root -l -q pred_qcdBkg.C++ > "logs/qcd_eg_${RunYear}.log" 
        #    root -l -q pred_rareBkg.C++ > "logs/rare_eg_${RunYear}.log" 
        #    root -l -q pred_sig.C++ > "logs/data_eg_${RunYear}.log" 
        #fi
	root -l -b -q 'plot_invmass.C(1)'
    done
done

#bash makeCards.sh
#bash makePlot.sh
#root -l -b -q 'plotSyst_VGBkg.C("/uscms_data/d3/tmishra/Output/signalTree_egamma_VGBkg_678.root", "/uscms_data/d3/tmishra/Output/signalTree_mg_VGBkg_678.root","2018")'

#echo "All tasks completed successfully. do the limit plots.. "
#cd  /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/SUSYAnalysis/limits/combine
#cmsenv
#bash MakeImpactPlot.sh
