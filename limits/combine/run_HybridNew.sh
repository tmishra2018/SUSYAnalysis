for POINT in "1300_1200" "1700_1100" "2100_100"; do
  CARD="/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_2018_new/counting_t5wg_${POINT}.txt"
  combine -M AsymptoticLimits ${CARD} --expectSignal=1  --rMin -20 --rMax 20 | tee limit_asym_${POINT}_blind.txt
  combine -M HybridNew ${CARD} --LHCmode LHC-limits --expectSignal=1  --rMin -20 --rMax 20 | tee limit_hyb_${POINT}_blind.txt
done

echo "Point        Asym_exp   HybridNew_exp"
for POINT in "1300_1200" "1700_1100" "2100_100"; do
  ASYM=$(grep "50.0%" limit_asym_${POINT}_blind.txt | awk '{print $NF}')
  HYB=$(grep "50.0%" limit_hyb_${POINT}_blind.txt   | awk '{print $NF}')
  echo "t5wg_${POINT}   ${ASYM}   ${HYB}"
done
