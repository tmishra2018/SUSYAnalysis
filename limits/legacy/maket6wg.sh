#!/bin/bash
	#	for filename in ./t6wgdata/counting_t6Wg_*_*.txt
	#	do
	#		rm t6filelist
	#		echo $filename > t6filelist
	#		gluino=`awk -F '[^0-9]*' '{print $4}' t6filelist`
	#		chargino=`awk -F '[^0-9]*' '{print $5}' t6filelist`
	#		echo $gluino ' ' $chargino 
	#		combine -M Asymptotic $filename | tee T6WG_Unblind_Data/limit_${gluino}_${chargino}.txt 
	#	done

		for filename in ./t6wgexp/counting_t6Wg_*_*.txt
		do
			rm t6filelist
			echo $filename > t6filelist
			gluino=`awk -F '[^0-9]*' '{print $4}' t6filelist`
			chargino=`awk -F '[^0-9]*' '{print $5}' t6filelist`
			echo $gluino ' ' $chargino 
			combine -M Asymptotic $filename | tee T6WG_Exp/limit_${gluino}_${chargino}.txt 
		done
