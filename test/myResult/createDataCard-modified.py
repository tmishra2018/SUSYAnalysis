####    script to produce a format for datacard, bin wise entries for each background is stored in SignalSystematic.root
####    lnN stands for Log-normal recommended choice for multiplicative corrections (efficiencies, cross sections, ...).
####    gmN stands for Gamma, and is the recommended choice for the statistical uncertainty on a background coming from the number of events in a control region 
import re
import os
import ROOT
import argparse
from os import system
import sys

limdir = './'

n_channels = int(sys.argv[1])*2
RunYear = sys.argv[2]
#n_channels = 36
n_processes= 6
# 5 bkgrounds, 1 susy process
pro_names = ['SUSY','elefakepho', 'jetfakepho', 'qcdfakelep','VGamma','rare']
syst_names = ['jes','esf','scale','e_to_pho_syst','j_to_pho_syst','fakelep_shape','xs','lumi','isr']

file_in = ROOT.TFile('/uscms_data/d3/tmishra/Output/SignalSystematic_'+RunYear+'.root', 'read')
h_rates = {}
for hname in pro_names:
    h_rates['h_' + hname + '_norm'] = file_in.Get('h_' + hname + '_norm')
    for sys in syst_names:
        h_rates['h_' + hname + '_syserr_' + sys] = file_in.Get('h_' + hname + '_syserr_' + sys)
h_rates['h_elefakepho_controlsample'] = file_in.Get('h_elefakepho_controlsample')
h_rates['h_elefakepho_transferfactor'] = file_in.Get('h_elefakepho_transferfactor')
h_rates['h_jetfakepho_controlsample'] = file_in.Get('h_jetfakepho_controlsample')
h_rates['h_jetfakepho_transferfactor'] = file_in.Get('h_jetfakepho_transferfactor')
h_rates['h_qcdfakelep_controlsample'] = file_in.Get('h_qcdfakelep_controlsample')
h_rates['h_qcdfakelep_transferfactor'] = file_in.Get('h_qcdfakelep_transferfactor')

logfile_out = open( 't5wg.log','w')
file_out = open( limdir+'counting_exp_XXX_YYY_{}_{}.txt'.format(n_channels,RunYear), 'w')

#file_out.write("imax 36 number of channels\n")
file_out.write('imax {:2d} number of channels\n'.format(n_channels))
file_out.write("jmax 5  number of backgrounds\n")
file_out.write("kmax *  number of nuisance parameters\n")
file_out.write("------------\n")

file_out.write('{:16s}'.format('bin'))
for i in range(1,n_channels+1):
        file_out.write('{:>8s} '.format('bin'+ str(i)))
file_out.write('\n')

file_out.write('{:16s}'.format('observation    '))
with open('data_mg.log') as fileMG:
	for line in fileMG:
	          if re.search('bin', line):
                        l = line.split('bin')
                        v = int(l[1].strip())
			file_out.write('{:>4d} '.format(v))

with open('data_eg.log') as fileEG:
	for line in fileEG:
	          if re.search('bin', line):
                        l = line.split('bin')
                        v = int(l[1].strip())
			file_out.write('{:>4d} '.format(v))
file_out.write('\n')

#### from file pred_sig.C
#if n_channels == 22 :
#	file_out.write("observation    321 517 105 36 76 52 1 2 3 1 2 163 297 81 38 55 36 1 1 6 1 1  \n")

#if n_channels == 36 :
	# MET cut starts at 100
	#file_out.write("747  862  138  32  64  45   1   1    5   15  44  28 4 12 7 1 1 0 345  483  109  32  46  32  1  1  4  11  31  17  6  9  4  0  1  3 \n")
	# photon pT boundary 200
	#file_out.write("observation    309  494  85  32  64  45  1  1  5  12  23  20  4  12  7  1  1  0  153  276  67   32   46   32   1  1  4  10  21  14  6  9  4  0  1  3  \n")
	# photon pT boundary 100
#	file_out.write("observation     262  436  60  20  52  37  1  1  3  59  81  45  16  24  15  1  1  2  107  225  54  21  35  24  1 0 2 56 72  27  17  20  12  0  2  5  \n")

#if n_channels == 18 :
#	file_out.write("observation    321  517  105  36  76  52  2  2  5  163  297  81  38  55  36  1  2  7  \n")

#if n_channels == 20 :
#	file_out.write("observation    321  517  105  36  76  52  1  2  3  3  163  297  81  38 55  36  1 1 6 2    \n")

file_out.write("------------\n")

file_out.write('{:26s}'.format('bin'))
for i in range(1,n_channels+1):
    for j in range(0,n_processes):
	file_out.write('{:>12s} '.format('bin'+ str(i)))
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
    bkg_rate = 0
    ratio = 0
    for p in pro_names:
        evtrate = 0
        if(p == 'SUSY'):
            file_out.write('{:>12s} '.format('NSC'+str(ich)))
        else:
	    # contribution from each backgrounds bin wise
            evtrate = h_rates['h_' + p + '_norm'].GetBinContent(ich)
            bkg_rate += evtrate
            file_out.write('{:12.3f} '.format(evtrate))
            logfile_out.write('{:15} {:2d} {:10.5f}\n'.format(p, ich, evtrate))
    ratio = signal_rate/bkg_rate
file_out.write('\n')

file_out.write('------------\n')


# Add it later
#file_out.write('{:15s} {:3s} {:6s}'.format('PDFSCALE','lnN',''))
#for ich in range(1,n_channels+1):
#    for p in pro_names:
#        if(p == 'SUSY'):
#            file_out.write('{:>12s} '.format('PDFSCALES'+str(ich)))
#        else:
#            file_out.write('{:>12s} '.format('-'))
#file_out.write('\n')
#file_out.write('{:15s} {:3s} {:6s}'.format('FASTMET','lnN',''))
#for ich in range(1,n_channels+1):
#    for p in pro_names:
#        if(p == 'SUSY'):
#            file_out.write('{:>12s} '.format('FASTMETS'+str(ich)))
#        else:
#            file_out.write('{:>12s} '.format('-'))
#file_out.write('\n')

file_out.close()
logfile_out.close()

