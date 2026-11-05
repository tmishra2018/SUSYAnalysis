#!/bin/bash
#for METbin1 in {200,220,240,260,300}
#do
#	for METbin2 in {320,340,360,380,400,420,440,460,500}
#	do
#		for i in {0..40}
#		do
#			j=$((300+i*25))
#			combine -M Asymptotic tchiwg/counting_tchiwg_${j}_${METbin1}_${METbin2}.txt | tee tchiwg/limit_${j}_${METbin1}_${METbin2}.txt 
#		done
#	done
#done


#for i in {0..40}
#do
#	j=$((300+i*25))
#	combine -M Asymptotic tchiwgdata/counting_tchiwg_${j}.txt | tee TChiWG_Unblind_Data/limit_${j}.txt 
#done

#for ich in {1..36}
#do
	ich=7
  for i in {0..40}
  do
	  j=$((300+i*25))
	  combine -M Asymptotic ../DataCard/Sensitive/cards/counting_tchiwg_${j}_${ich}.txt | tee ../DataCard/Sensitive/limits/limit_${j}_${ich}.txt 
  done
#done
