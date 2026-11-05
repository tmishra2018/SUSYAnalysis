#!/bin/bash

#for filename in /eos/uscms/store/user/tmishra/CombinedLimit/Sensitive/cards/counting_t5Wg_*_*_*.txt
#		do
#			rm filelist
#			echo $filename > filelist
#			gluino=`awk -F '[^0-9]*' '{print $3}' filelist`
#			chargino=`awk -F '[^0-9]*' '{print $4}' filelist`
#      channel=`awk -F '[^0-9]*' '{print $5}' filelist`
#			echo $gluino ' ' $chargino ' ' $channel 
#			combine -M AsymptoticLimits $filename | tee /eos/uscms/store/user/tmishra/CombinedLimit/Sensitive/limits/limit_${gluino}_${chargino}_${channel}.txt 
#		done

for filename in /eos/uscms/store/user/tmishra/CombinedLimit/Sensitive_myFile/cardsPlan5/counting_t5Wg_*_*.txt
do
	rm t5filelist1
	echo $filename > t5filelist1
	gluino=`awk -F '[^0-9]*' '{print $4}' t5filelist1`
	chargino=`awk -F '[^0-9]*' '{print $5}' t5filelist1`
	echo  $gluino ' ' $chargino  >> t5filelistt
	#combine -M AsymptoticLimits $filename | tee /eos/uscms/store/user/tmishra/CombinedLimit/Sensitive_myFile/limitsPlan5/limit_${gluino}_${chargino}.txt 
	#combine -M Significance $filename -t -1 --expectSignal=1 | tee /eos/uscms/store/user/tmishra/CombinedLimit/Sensitive/signif/signif_${gluino}_${chargino}.txt 
done
