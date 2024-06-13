#!/bin/sh
input_Scan=$1
anaExe="multijet"
anaArg="FastSim_T5ttttZg"
exeAtWorker="worker2.sh"
dataSetType="T5ttttZg"
#filesToTransfer="GJets.root,QCD.root,TTGJets.root,TTJetsHT.root,WGJetsToLNuG.root,WJetsToLNu.root,ZGJetsToNuNuG.root,ZJetsToNuNu.root,makeDatacard_SBins.C,${anaExe},combineCards.py,higgsCombine.tar"
filesToTransfer="makeDatacard_SBinsV7_v6.C,SBinHists.root,mGl_Xsecpb_absUnc.root,${anaExe},combineCards.py,higgsCombine.tar,PileupHistograms_0121_69p2mb_pm4p6.root,T5bbbbZg_MassScan.root,T5qqqqHg_MassScan.root,T5ttttZg_MassScan.root,T6ttZg_MassScan.root,SignalCont_T6ttZg_LostMu.root,SignalCont_T5ttttZg_LostMu.root,btagFiles.tar"

while read -a massP
do 
    echo ${massP[0]} ${massP[1]}
    jdl_file="condor_${dataSetType}_${anaArg}_${massP[0]}_${massP[1]}_job.jdl"
    log_prefix="condor_${dataSetType}_${anaArg}_${massP[0]}_${massP[1]}_job"
    echo "universe = vanilla">$jdl_file
    echo "Executable = $exeAtWorker">>$jdl_file
    echo "Should_Transfer_Files = YES">>$jdl_file
    echo "WhenToTransferOutput = ON_EXIT_OR_EVICT">>$jdl_file
    echo "Transfer_Input_Files = ${filesToTransfer}">>$jdl_file
    echo "Output = ${log_prefix}.stdout">>$jdl_file
    echo "Error = ${log_prefix}.stderr">>$jdl_file
    echo "Log = ${log_prefix}.condor">>$jdl_file
    echo "notification = never">>$jdl_file
    echo "Arguments = ${anaExe} ${massP[0]} ${massP[1]} ${anaArg} ${dataSetType}">>$jdl_file
    echo "Queue">>$jdl_file
    condor_submit $jdl_file
done < $input_Scan

