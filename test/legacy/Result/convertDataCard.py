import re
import os
import ROOT
import argparse
import math

limdir = './combinecard/'

n_channels = 36
n_processes= 6
file_list = ['signalTree_mg_eleBkg.dat',
						'signalTree_mg_jetBkg.dat',
						'signalTree_mg_qcdBkg.dat',
						'signalTree_mg_VGBkg.dat',
						'signalTree_mg_rareBkg.dat',
						'signalTree_egamma_eleBkg.dat',
						'signalTree_egamma_jetBkg.dat',
						'signalTree_egamma_qcdBkg.dat',
						'signalTree_egamma_VGBkg.dat',
						'signalTree_egamma_rareBkg.dat']
pro_names = ['SUSY','elefakepho', 'jetfakepho', 'qcdfakelep','VGamma','rare']
syst_names = ['jes','phoesf','lepesf','scale','e_to_pho_syst','j_to_pho_syst','fakelep_shape','xs','lumi','isr','PU','PDFSCALE','FASTMET']
data_types = ['nom','control','fakerate','stat','jes','phoesf','lepesf','scale','e_to_pho_syst','j_to_pho_syst','fakelep_shape','xs','lumi','isr']
evtyield = {}
for i in range(1,36+1):
	for process in pro_names:
		evtyield[process+str(i)] = {}

for ip in range(1,6):
	ifile = ip - 1
	theFile = open(file_list[ifile], "r")
	for i in range(1,18+1):
		line = theFile.readline()
		theInts = line.split()
		isys = 0
		for datatype in data_types:
			evtyield[pro_names[ip]+str(i)][datatype] = theInts[isys]
			isys += 1
	theFile.close()

for ip in range(1,6):
	ifile = ip - 1 + 5
	theFile = open(file_list[ifile], "r")
	for i in range(19,36+1):
		line = theFile.readline()
		theInts = line.split()
		isys = 0
		for datatype in data_types:
			evtyield[pro_names[ip]+str(i)][datatype] = theInts[isys]
			isys += 1
	theFile.close()

file_out = open( limdir+'newcounting_exp_XXX_YYY.txt', 'w')

file_out.write("imax 36 number of channels\n")
file_out.write("jmax 5  number of backgrounds\n")
file_out.write("kmax *  number of nuisance parameters\n")
file_out.write("------------\n")
#file_out.write("bin             lg_1  lg_2  lg_3  lg_4  lg_5  lg_6  lg_7  lg_8  lg_9  lg_10  lg_11  lg_12  lg_13  lg_14  lg_15  lg_16  lg_17  lg_18  lg_19  lg_20  lg_21  lg_22  lg_23  lg_24  lg_25  lg_26  lg_27  lg_28  lg_29  lg_30  lg_31  lg_32  lg_33  lg_34  lg_35  lg_36   \n")
file_out.write("bin             1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36   \n")
file_out.write("observation    309  494  85  32  64  45  1  1  5  12  23  20  4  12  7  1  1  0  153  276  67   32   46   32   1  1  4  10  21  14  6  9  4  0  1  3  \n")
## file_out.write("observation    153  276  67   32   46   32   1  1  4  10  21  14  6  9  4  0  1  3 309  494  85  32  64  45  1  1  5  12  23  20  4  12  7  1  1  0   \n")
## ## COmbine ##
## ##file_out.write("observation    309  501  86  33  65  45  1  1  5  12  24  20  4  13  7  1  1  0  154  276  67   32   46   32   1  1  4  10  22  15  7  9  4  0  1  3  \n")
 
file_out.write("------------\n")

file_out.write('{:26s}'.format('bin'))
for i in range(1,n_channels+1):
    for j in range(0,n_processes):
        file_out.write('{:12d} '.format(i))
file_out.write('\n')

file_out.write('{:26s}'.format('process'))
for i in range(0,n_channels):
    for p in pro_names: 
        file_out.write('{:>12s} '.format(p))
file_out.write('\n')

file_out.write('{:26s}'.format('process'))
for i in range(0,n_channels):
    for j in range(0,n_processes):
        file_out.write('{:12d} '.format(j))
file_out.write('\n')



file_out.write('{:26s}'.format('rate'))
for ich in range(1,n_channels+1):
    signal_rate = 0
    for p in pro_names:
        evtrate = 0
        if(p == 'SUSY'):
            file_out.write('{:>12s} '.format('NSC'+str(ich)))
        else:
            evtrate = float(evtyield[p+str(ich)]['nom'])
            file_out.write('{:12.3f} '.format(evtrate))
file_out.write('\n')

file_out.write('------------\n')

for syst in syst_names:
    file_out.write('{:15s} {:3s} {:6s}'.format(syst,'lnN',''))
    for ich in range(1,n_channels+1):
        for p in pro_names:
				if(p == 'SUSY'):
					if(syst == 'jes' or syst == 'phoesf' or syst == 'lepesf' or syst == 'isr' or syst == 'PU' or syst == 'PDFSCALE' or syst == 'FASTMET'):
					    file_out.write('{:>12s} '.format(syst+ str(ich)))
					elif(syst == 'lumi'):
					    file_out.write('{:12.3f} '.format(1.025))
					else: 
					    file_out.write('{:>12s} '.format('-'))
				else:
					if(syst == 'PU' or syst == 'PDFSCALE' or syst == 'FASTMET'):
							file_out.write('{:>12s} '.format('-'))
					else:
						syserr = float(evtyield[p+str(ich)][syst])
						if(syserr < 0.0001):
								file_out.write('{:>12s} '.format('-'))
						else:
							if(syst == 'scale' and p == 'VGamma'):
									file_out.write('{:12.3f} '.format(1+syserr))
							elif(syst == 'scale' and p == 'qcdfakelep'):
									file_out.write('{:12.3f} '.format(1-syserr))
							else:
									file_out.write('{:12.3f} '.format(1+syserr))
    file_out.write('\n')

for ich in range(1,n_channels+1):
    file_out.write('{:15s} {:3s} {:6s}'.format('SUSY_stat'+str(ich),'lnN',''))
    for k in range(1,n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('STSC'+ str(ich),'-', '-', '-', '-', '-'))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')

for ich in range(1, n_channels+1):
    nevt = float(evtyield['elefakepho'+str(ich)]['control']) 
    fakerate = float(evtyield['elefakepho'+str(ich)]['fakerate']) 
    file_out.write('{:15s} {:3s} {:6d}'.format('e_to_pho_stat'+str(ich),'gmN',int(nevt)))
    for k in range(1, n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:12.5f} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', fakerate, '-', '-', '-', '-'))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')
for ich in range(1, n_channels+1):
    nevt = float(evtyield['jetfakepho'+str(ich)]['control']) 
    fakerate = float(evtyield['jetfakepho'+str(ich)]['fakerate']) 
    file_out.write('{:15s} {:3s} {:6d}'.format('j_to_pho_stat'+str(ich),'gmN',int(nevt)))
    for k in range(1,n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:>12s} {:12.5f} {:>12s} {:>12s} {:>12s} '.format('-', '-', fakerate, '-', '-', '-'))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')
for ich in range(1,n_channels+1):
    nevt = float(evtyield['qcdfakelep'+str(ich)]['control']) 
    fakerate = float(evtyield['qcdfakelep'+str(ich)]['fakerate']) 
    file_out.write('{:15s} {:3s} {:6d}'.format('j_to_lep_stat'+str(ich),'gmN',int(nevt)))
    for k in range(1, n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:>12s} {:>12s} {:12.5f} {:>12s} {:>12s} '.format('-', '-', '-', fakerate, '-', '-'))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')
for ich in range(1, n_channels+1):
    file_out.write('{:15s} {:3s} {:6s}'.format('rare_stat'+str(ich),'lnN',''))
    staterror= 1.0 + float(evtyield['rare'+str(ich)]['stat'])
    for k in range(1, n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:12.3f} '.format('-', '-', '-', '-', '-',staterror))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')

for ich in range(1, n_channels+1):
    file_out.write('{:15s} {:3s} {:6s}'.format('VG_stat'+str(ich),'lnN',''))
    staterror= 1.0 + float(evtyield['VGamma'+str(ich)]['stat']) 
    for k in range(1, n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:12.3f} {:>12s} '.format('-', '-', '-', '-',staterror, '-'))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')
file_out.close()

