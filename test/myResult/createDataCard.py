####    script to produce a format for datacard, bin wise entries for each background is stored in SignalSystematic.root
####    lnN stands for Log-normal recommended choice for multiplicative corrections (efficiencies, cross sections, ...).
####    gmN stands for Gamma, and is the recommended choice for the statistical uncertainty on a background coming from the number of events in a control region 
import re
import os
import ROOT
import argparse
from os import system
import sys

limdir = 'logs/'

n_channels = int(sys.argv[1])*2
RunYear = sys.argv[2]
preVFP = sys.argv[3]
#n_channels = 36
n_processes= 6
# 5 bkgrounds, 1 susy process
pro_names = ['SUSY','elefakepho', 'jetfakepho', 'qcdfakelep','VGamma','rare']
syst_names = ['jes','jer','esf','scale','e_to_pho_syst','j_to_pho_syst','xs','lumi']

if RunYear == '2016' and preVFP == '1':
    whichVFP = 'preVFP'
elif RunYear == '2016' and preVFP == '0':
    whichVFP = 'postVFP'
else:
    whichVFP = ''

file_in = ROOT.TFile(f'/uscms_data/d3/tmishra/Output/SignalSystematic_{RunYear}{whichVFP}.root', 'read')
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

logfile_out = open( 'logs/t5wg.log','w')
filename = limdir + 'counting_exp_XXX_YYY_{}_{}{}.txt'.format(n_channels, RunYear, whichVFP)
file_out = open(filename, 'w')

#file_out.write("imax 36 number of channels\n")
file_out.write('imax {:2d} number of channels\n'.format(n_channels))
file_out.write("jmax 5  number of backgrounds\n")
file_out.write("kmax *  number of nuisance parameters\n")
file_out.write("------------\n")
file_out.write("shapes * * FAKE\n")

file_out.write('{:16s}'.format('bin'))
for i in range(1,n_channels+1):
        file_out.write('{:>8s} '.format('bin'+ str(i)))
file_out.write('\n')

file_out.write('{:16s}'.format('observation    '))

with open('logs/data_mg_{}{}.log'.format(RunYear, whichVFP)) as fileMG:
    for line in fileMG:
        if re.search('bin', line):
            l = line.split('bin')
            v = int(l[1].strip())
            file_out.write('{:>4d} '.format(v))

with open('logs/data_eg_{}{}.log'.format(RunYear, whichVFP)) as fileEG:
    for line in fileEG:
        if re.search('bin', line):
            l = line.split('bin')
            v = int(l[1].strip())
            file_out.write('{:>4d} '.format(v))
file_out.write('\n')

file_out.write("------------\n")

file_out.write('{:26s}'.format('bin'))
for i in range(1, n_channels + 1):
    for j in range(0, n_processes):
        file_out.write('{:>12s} '.format('bin' + str(i)))
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

low_p = {}
high_p = {}
for p in pro_names:
  low_p[p] = 2
  high_p[p] = 0

for syst in syst_names:
    file_out.write('{:15s} {:3s} {:6s}'.format(syst,'lnN',''))
    for ich in range(1,n_channels+1):
        for p in pro_names:
            if(p == 'SUSY'):
		# Just a string is saved for SUSY
                if(syst == 'jes'):
                    file_out.write('{:>12s} '.format('JESS'+ str(ich)))
                elif(syst == 'esf'):
                    file_out.write('{:>12s} '.format('ESF'+ str(ich)))
#                elif(syst == 'isr'):
#                    file_out.write('{:>12s} '.format('ISR'+ str(ich)))
                elif(syst == 'lumi'):
                    file_out.write('{:12.3f} '.format(1.026))
                else: 
                    file_out.write('{:>12s} '.format('-'))
            else:
                syserr = 0
                if(h_rates['h_' + p + '_norm'].GetBinContent(ich)==0):
		    # errors from transfer factors taken for bkgs, when central value is zero
                    if(p == 'elefakepho' and syst == 'e_to_pho_syst'):
                        syserr = h_rates['h_elefakepho_transferfactor'].GetBinError( ich) + 1
                        file_out.write('{:12.3f} '.format(syserr))
                    elif(p == 'jetfakepho' and syst == 'j_to_pho_syst'):
                        syserr = h_rates['h_jetfakepho_transferfactor'].GetBinError( ich) + 1
                        file_out.write('{:12.3f} '.format(syserr))
                    elif(p == 'qcdfakelep' and syst == 'scale'):
                        syserr = 1 - h_rates['h_qcdfakelep_transferfactor'].GetBinError(ich)
                        file_out.write('{:12.3f} '.format(syserr))
                    else:
                        file_out.write('{:>12s} '.format('-'))
                else:
		    # when bkg events in bins is non zero
                    syserr = h_rates['h_' + p + '_syserr_' + syst].GetBinContent(ich)/h_rates['h_' + p + '_norm'].GetBinContent(ich)
                    if(syserr < 0):
                         file_out.write('{:>12s} '.format('-'))
                    else:
                         if(p == "qcdfakelep" and syst == "scale"):
                             syserr = 1 - syserr
                         else: 
                             syserr = syserr + 1
                         file_out.write('{:12.3f} '.format(syserr))
    file_out.write('\n')

file_out.write('{:15s} {:3s} {:6s}'.format('PU','lnN',''))
for ich in range(1,n_channels+1):
    for p in pro_names:
        if(p == 'SUSY'):
            file_out.write('{:>12s} '.format('PUS'))
        else:
            file_out.write('{:>12s} '.format('-'))
file_out.write('\n')

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

for ich in range(1,n_channels+1):
    file_out.write('{:15s} {:3s} {:6s}'.format('SUSY_stat'+str(ich),'lnN',''))
    #file_out.write('{:22s} {:3s} {:6s}'.format('SUSY_' + RunYear + whichVFP + '_stat' + str(ich), 'lnN', ''))
    for k in range(1,n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('STSC'+ str(ich),'-', '-', '-', '-', '-'))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')


for ich in range(1,n_channels+1):
    nevt = h_rates['h_elefakepho_controlsample'].GetBinContent(ich)
    fakerate = h_rates['h_elefakepho_transferfactor'].GetBinContent(ich)
    # number of events in control sample and fake rate
    file_out.write('{:22s} {:3s} {:6d}'.format('e_to_pho_'  + RunYear + whichVFP + '_stat'+str(ich),'gmN',int(nevt)))
    for k in range(1,n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:12.5f} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', fakerate, '-', '-', '-', '-'))
            logfile_out.write('{:15} {:2d} {:10.5f}\n'.format('elefakepho', ich, nevt*fakerate))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')
for ich in range(1,n_channels+1):
    nevt = h_rates['h_jetfakepho_controlsample'].GetBinContent(ich)
    fakerate = h_rates['h_jetfakepho_transferfactor'].GetBinContent(ich)
    # number of events in control sample and fake rate
    file_out.write('{:22s} {:3s} {:6d}'.format('j_to_pho_'  + RunYear + whichVFP + '_stat'+str(ich),'gmN',int(nevt)))
    for k in range(1,n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:>12s} {:12.5f} {:>12s} {:>12s} {:>12s} '.format('-', '-', fakerate, '-', '-', '-'))
            logfile_out.write('{:15} {:2d} {:10.5f}\n'.format('jetfakepho', ich, nevt*fakerate))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')
for ich in range(1,n_channels+1):
    nevt = h_rates['h_qcdfakelep_controlsample'].GetBinContent(ich)
    fakerate = h_rates['h_qcdfakelep_transferfactor'].GetBinContent(ich)
    # print(ich, nevt, fakerate)
    # number of events in control sample and fake rate
    file_out.write('{:22s} {:3s} {:6d}'.format('j_to_lep_'  + RunYear + whichVFP + '_stat'+str(ich),'gmN',int(nevt)))
    for k in range(1,n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:>12s} {:>12s} {:12.5f} {:>12s} {:>12s} '.format('-', '-', '-', fakerate, '-', '-'))
            logfile_out.write('{:15} {:2d} {:10.5f}\n'.format('jetfakelep', ich, nevt*fakerate))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')

for ich in range(1,n_channels+1):
    file_out.write('{:15s} {:3s} {:6s}'.format('rare_stat'+str(ich),'lnN',''))
    #file_out.write('{:22s} {:3s} {:6s}'.format('rare_'  + RunYear + whichVFP + '_stat'+str(ich),'lnN',''))
    if(h_rates['h_rare_norm'].GetBinContent(ich) != 0):
    	staterror= 1.0 + h_rates['h_rare_norm'].GetBinError(ich)/h_rates['h_rare_norm'].GetBinContent(ich)
    else:
	    staterror= 1.0	
    # stat error for rare bkgs
    for k in range(1,n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:12.3f} '.format('-', '-', '-', '-', '-',staterror))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')

for ich in range(1,n_channels+1):
    file_out.write('{:15s} {:3s} {:6s}'.format('VG_stat'+str(ich),'lnN',''))
    #file_out.write('{:22s} {:3s} {:6s}'.format('VG_'  + RunYear + whichVFP + '_stat'+str(ich),'lnN',''))
    if(h_rates['h_VGamma_norm'].GetBinContent(ich) != 0):
    	staterror= 1.0 + h_rates['h_VGamma_norm'].GetBinError(ich)/h_rates['h_VGamma_norm'].GetBinContent(ich)
    else:
	    staterror= 1.0	
    # stat error for VG bkgs
    for k in range(1,n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:12.3f} {:>12s} '.format('-', '-', '-', '-',staterror, '-'))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')
file_out.close()
logfile_out.close()
