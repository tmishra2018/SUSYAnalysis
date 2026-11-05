cd /uscms_data/d3/tmishra/datacards
rm -rf input
mkdir -p input

for TOPO in T5WG T6WG TChiWG; do
  for YEAR in 2016preVFP 2016postVFP 2017 2018 Run2; do
    mkdir -p input/${TOPO}/${YEAR}
    cp /eos/uscms/store/user/tmishra/CombinedLimit/${TOPO}/cards/cards_${YEAR}_splusb/counting_*.txt input/${TOPO}/${YEAR}/
  done
done

cp /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/SUSYAnalysis/limits/datacards/syst_SUS2411.yaml systematics.yml
cp /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/SUSYAnalysis/limits/datacards/systematics/check_names.py systematics/
cp /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/SUSYAnalysis/limits/datacards/convert_all_cards.py .
cp /uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/SUSYAnalysis/limits/datacards/convert_dict.yaml .

find input/ -name "counting_*.txt" | wc -l 
du -sh input/

#git add input/ convert_all_cards.py convert_dict.yaml systematics.yml  systematics/
#git commit -m "Updated disctionary and TChiwg datacards"
#git push origin datacards/Run2-SUS2411-renamed-systematics
