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

	root -l -q "plot_newbkg_NoData.C"
	root -l -q "plot_pie_chart.C(\"$log_file\",\"$RunYear\",\"$VFP_string\")"
	
	log_file="../Background/logs/BKG_${RunYear}${VFP_string}.log"
	root -l -q "plot_pie_chart.C(\"$log_file\",\"$RunYear\",\"$VFP_string\")"

	log_file="../Background/logs/VALID_${RunYear}${VFP_string}.log"
	root -l -q "plot_pie_chart.C(\"$log_file\",\"$RunYear\",\"$VFP_string\")"
        
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
	
	root -l -q "plot_newbkg_NoData.C"
	root -l -q "plot_pie_chart.C(\"$log_file\",\"$RunYear\",\"$VFP_string\")"
	
	log_file="../Background/logs/BKG_${RunYear}${VFP_string}.log"
	root -l -q "plot_pie_chart.C(\"$log_file\",\"$RunYear\",\"$VFP_string\")"

	log_file="../Background/logs/VALID_${RunYear}${VFP_string}.log"
	root -l -q "plot_pie_chart.C(\"$log_file\",\"$RunYear\",\"$VFP_string\")"
	
	log_file="logs/eventcount_${RunYear}${VFP_string}.txt"
	rm -f "$log_file"
	root -l -q "plot_eventct.C+($NBIN)" >> "logs/eventcount_${RunYear}${VFP_string}.txt"
	root -l -q "plot_eventct_NoData.C+($NBIN)"
    done
done

############################################
# Combine all years (Run2 = 678)
############################################

rm -f binConfig.txt
echo 'NBIN' $NBIN >> binConfig.txt
echo 'METbin1' $METbin1 >> binConfig.txt
echo 'METbin2' $METbin2 >> binConfig.txt
echo 'HTbin1' $HTbin1 >> binConfig.txt
echo 'HTbin2' $HTbin2 >> binConfig.txt
echo 'PHOETbin' $PHOETbin >> binConfig.txt

OUTDIR="/uscms_data/d3/tmishra/Output"
SIGDIR="/uscms/home/tmishra/nobackup/signal_trees"

BKG_LIST=("eleBkg" "jetbkg" "qcd" "VGBkg" "rareBkg")
YEARS=("2016preVFP" "2016postVFP" "2017" "2018")

echo "=== Hadd Run2 combined ROOT files (678) ==="

# mgamma backgrounds
for bkg in "${BKG_LIST[@]}"; do
    hadd -f ${OUTDIR}/signalTree_mg_${bkg}_678.root \
        ${OUTDIR}/signalTree_mg_${bkg}_2016preVFP.root \
        ${OUTDIR}/signalTree_mg_${bkg}_2016postVFP.root \
        ${OUTDIR}/signalTree_mg_${bkg}_2017.root \
        ${OUTDIR}/signalTree_mg_${bkg}_2018.root
done

# egamma backgrounds
for bkg in "${BKG_LIST[@]}"; do
    hadd -f ${OUTDIR}/signalTree_egamma_${bkg}_678.root \
        ${OUTDIR}/signalTree_egamma_${bkg}_2016preVFP.root \
        ${OUTDIR}/signalTree_egamma_${bkg}_2016postVFP.root \
        ${OUTDIR}/signalTree_egamma_${bkg}_2017.root \
        ${OUTDIR}/signalTree_egamma_${bkg}_2018.root
done

# signals
for chan in mg egamma; do
    hadd -f ${OUTDIR}/signalTree_${chan}_signal_678.root \
        ${OUTDIR}/signalTree_${chan}_signal_2016preVFP.root \
        ${OUTDIR}/signalTree_${chan}_signal_2016postVFP.root \
        ${OUTDIR}/signalTree_${chan}_signal_2017.root \
        ${OUTDIR}/signalTree_${chan}_signal_2018.root
done

hadd -f ${SIGDIR}/signalTree_T6WG_678.root \
    ${SIGDIR}/signalTree_T6WG_2016preVFP.root \
    ${SIGDIR}/signalTree_T6WG_2016postVFP.root \
    ${SIGDIR}/signalTree_T6WG_2017.root \
    ${SIGDIR}/signalTree_T6WG_2018.root

hadd -f ${SIGDIR}/signalTree_T5WG_678.root \
    ${SIGDIR}/signalTree_T5WG_2016preVFP.root \
    ${SIGDIR}/signalTree_T5WG_2016postVFP.root \
    ${SIGDIR}/signalTree_T5WG_2017.root \
    ${SIGDIR}/signalTree_T5WG_2018.root

hadd -f ${SIGDIR}/signalTree_TChiWG_678.root \
    ${SIGDIR}/signalTree_TChiWG_2016preVFP.root \
    ${SIGDIR}/signalTree_TChiWG_2016postVFP.root \
    ${SIGDIR}/signalTree_TChiWG_2017.root \
    ${SIGDIR}/signalTree_TChiWG_2018.root
############################################
# Run plots for combined Run2 (678)
############################################

RunYear=678
preVFP=0
VFP_string=""

echo "Running combined Run2 plots (RunYear=678)"

rm -f SigConfig.txt
echo 'ichannel' 1 >> SigConfig.txt
echo 'anatype' 3 >> SigConfig.txt
echo 'lowMt' 100 >> SigConfig.txt
echo 'highMt' -1 >> SigConfig.txt
echo 'lowMET' 120 >> SigConfig.txt
echo 'highMET' -1 >> SigConfig.txt
echo 'lowPt' 0 >> SigConfig.txt
echo 'highPt' -1 >> SigConfig.txt
echo 'lepIso' 4 >> SigConfig.txt
echo 'RunYear' 678 >> SigConfig.txt
echo 'preVFP' 0 >> SigConfig.txt

log_file="logs/eventcount_678.txt"
rm -f "$log_file"

root -l -q "plot_newbkg_NoData.C"
root -l -q "plot_pie_chart.C(\"$log_file\",\"678\",\"\")"

rm -f SigConfig.txt
echo 'ichannel' 2 >> SigConfig.txt
echo 'anatype' 3 >> SigConfig.txt
echo 'lowMt' 100 >> SigConfig.txt
echo 'highMt' -1 >> SigConfig.txt
echo 'lowMET' 120 >> SigConfig.txt
echo 'highMET' -1 >> SigConfig.txt
echo 'lowPt' 0 >> SigConfig.txt
echo 'highPt' -1 >> SigConfig.txt
echo 'lepIso' 4 >> SigConfig.txt
echo 'RunYear' 678 >> SigConfig.txt
echo 'preVFP' 0 >> SigConfig.txt

log_file="logs/eventcount_678.txt"
rm -f "$log_file"

root -l -q "plot_newbkg_NoData.C"
root -l -q "plot_pie_chart.C(\"$log_file\",\"678\",\"\")"

#########################################################################

root -l -q "plot_eventct.C+($NBIN)" >> "$log_file"
root -l -q "plot_eventct_NoData.C+($NBIN)"
