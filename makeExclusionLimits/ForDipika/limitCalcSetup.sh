#!/bin/sh    
for i in makeDatacard_SBinsV7_v4.C SBinHists.root mGl_Xsecpb_absUnc.root multijet combineCards.py higgsCombine.tar PileupHistograms_0121_69p2mb_pm4p6.root T5bbbbZg_MassScan.root T5qqqqHg_MassScan.root T5ttttZg_MassScan.root T6ttZg_MassScan.root btagFiles.tar genRecoMET_T5bbbbZg_1600_150.root bTagSFupFile_T5bbbbZg_1600_150.root worker2.sh *_MassScan_subset.txt calcLimit.sh FastSim_T5bbbbZg_1600_150.root
do
    echo $i
#    cp /uscms/home/vhegde/nobackup/Physics/skims_Run2ProductionV12/SignalRegion_getLimitUnblinded/$i .
done
#mv genRecoMET_T5bbbbZg_1600_150.root genRecoMET.root
#mv bTagSFupFile_T5bbbbZg_1600_150.root bTagSFupFile.root
echo "Copied reqiured files to your area from /uscms/home/vhegde/nobackup/Physics/skims_Run2ProductionV12/SignalRegion_getLimitUnblinded/ ......"
echo "Setup  CMSSW_7_4_7 . Then if you want to calculate limit for a sample mass point (T5bbbbZg_1600_150) run the following commands:"
echo "\$\$\$\$\$\$\$\$\$\$\$"
echo "mkdir dataCards"
echo "root -l 'makeDatacard_SBinsV7_v4.C(1600,150,\"FastSim_T5bbbbZg_1600_150.root\")' -q"
echo "combineCards.py dataCards/*.txt > Excard.txt"
echo "combine -M Asymptotic Excard.txt"

echo "----------------------------------"
echo "For submitting jobs to condor for different SMS scans, follow these instructions."
echo "1. IMP: Modify these lines in calcLimit.sh depending on the scan for which you want to submit job."
echo "anaArg=\"FastSim_T5ttttZg\"#or FastSim_T5bbbbZg or FastSim_T5qqqqHg or FastSim_T6ttZg"
echo "dataSetType=\"T5ttttZg\"#or T5bbbbZg or T5qqqqHg or T6ttZg"
echo "2. Submit the jobs to condor using following command. Note you need to give correct text file containing Gluion NLSP mass points. This text file should be compatible with what you mentioned in calcLimit.sh"
echo "\$\$\$\$\$\$\$\$\$\$\$"
echo "./calcLimit.sh T5ttttZg_MassScan_subset.txt"
