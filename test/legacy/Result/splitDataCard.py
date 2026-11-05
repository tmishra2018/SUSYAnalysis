import re
import os
import ROOT
import argparse

limdir = './t5wg_test/'

n_channels = 36
n_processes= 6
pro_names = ['SUSY','elefakepho', 'jetfakepho', 'qcdfakelep','VGamma','rare']
syst_names = ['jes','jer','esf','scale','e_to_pho_syst','j_to_pho_syst','fakelep_shape','xs','lumi','isr']

file_in = ROOT.TFile('./SignalSystematic.root', 'read')
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

for ich in range(1,n_channels+1):
  logfile_out= open( 't5wg.log','w')
  file_out = open( limdir+'counting_exp_XXX_YYY_'+str(ich)+'.txt', 'w')
  
  file_out.write("imax 1 number of channels\n")
  file_out.write("jmax 5  number of backgrounds\n")
  file_out.write("kmax *  number of nuisance parameters\n")
  file_out.write("------------\n")
  file_out.write("bin  1  \n")
  file_out.write("observation  324.214  467.461   98.865   27.271   61.307   45.150    1.236    2.914    5.311    6.538   21.444   15.580    5.012    8.368    5.459    0.766    0.634    0.539  163.626  256.955   78.950   17.309   49.956   28.164    1.302    1.160    2.804    5.927   21.614   11.677    4.135    8.876    5.123    0.394    0.532    0.848\n") 
  file_out.write("------------\n")
  
  file_out.write('{:26s}'.format('bin'))
  for j in range(0,n_processes):
    file_out.write('{:12d} '.format(1))
  file_out.write('\n')
  
  file_out.write('{:26s}'.format('process'))
  for p in pro_names: 
      file_out.write('{:>12s} '.format(p))
  file_out.write('\n')
  
  file_out.write('{:26s}'.format('process'))
  for j in range(0,n_processes):
      file_out.write('{:12d} '.format(j))
  file_out.write('\n')
  
  file_out.write('{:26s}'.format('rate'))
  signal_rate = 0
  bkg_rate = 0
  ratio = 0
  for p in pro_names:
      evtrate = 0
      if(p == 'SUSY'):
          file_out.write('{:>12s} '.format('NSC'+str(ich)))
      else:
          evtrate = h_rates['h_' + p + '_norm'].GetBinContent(ich)
          bkg_rate += evtrate
          file_out.write('{:12.3f} '.format(evtrate))
          logfile_out.write('{:15} {:2d} {:10.5f}\n'.format(p, ich, evtrate))
  ratio = signal_rate/bkg_rate
  file_out.write('\n')
  
  file_out.write('------------\n')
  
  for syst in syst_names:
    file_out.write('{:15s} {:3s} {:6s}'.format(syst,'lnN',''))
    for p in pro_names:
        if(p == 'SUSY'):
            if(syst == 'jes'):
                file_out.write('{:>12s} '.format('JESS'+ str(ich)))
            elif(syst == 'jer'):
                file_out.write('{:>12s} '.format('JERS'+ str(ich)))
            elif(syst == 'esf'):
                file_out.write('{:>12s} '.format('ESF'+ str(ich)))
            elif(syst == 'lumi'):
                file_out.write('{:12.3f} '.format(1.026))
            else: 
                file_out.write('{:>12s} '.format('-'))
        else:
            syserr = 0
            if(h_rates['h_' + p + '_norm'].GetBinContent(ich)==0):
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
  
  file_out.write('{:15s} {:3s} {:6s}'.format('SUSY_stat'+str(ich),'lnN',''))
  for k in range(1,n_channels+1):
      if( k == ich):
          file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('STSC'+ str(ich),'-', '-', '-', '-', '-'))
  file_out.write('\n')
  
  nevt = h_rates['h_elefakepho_controlsample'].GetBinContent(ich)
  fakerate = h_rates['h_elefakepho_transferfactor'].GetBinContent(ich)
  file_out.write('{:15s} {:3s} {:6d}'.format('e_to_pho_stat'+str(ich),'gmN',int(nevt)))
  for k in range(1,n_channels+1):
      if( k == ich):
          file_out.write('{:>12s} {:12.5f} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', fakerate, '-', '-', '-', '-'))
          logfile_out.write('{:15} {:2d} {:10.5f}\n'.format('elefakepho', ich, nevt*fakerate))
  file_out.write('\n')

  nevt = h_rates['h_jetfakepho_controlsample'].GetBinContent(ich)
  fakerate = h_rates['h_jetfakepho_transferfactor'].GetBinContent(ich)
  file_out.write('{:15s} {:3s} {:6d}'.format('j_to_pho_stat'+str(ich),'gmN',int(nevt)))
  for k in range(1,n_channels+1):
      if( k == ich):
          file_out.write('{:>12s} {:>12s} {:12.5f} {:>12s} {:>12s} {:>12s} '.format('-', '-', fakerate, '-', '-', '-'))
          logfile_out.write('{:15} {:2d} {:10.5f}\n'.format('jetfakepho', ich, nevt*fakerate))
  file_out.write('\n')

  nevt = h_rates['h_qcdfakelep_controlsample'].GetBinContent(ich)
  fakerate = h_rates['h_qcdfakelep_transferfactor'].GetBinContent(ich)
  file_out.write('{:15s} {:3s} {:6d}'.format('j_to_lep_stat'+str(ich),'gmN',int(nevt)))
  for k in range(1,n_channels+1):
      if( k == ich):
          file_out.write('{:>12s} {:>12s} {:>12s} {:12.5f} {:>12s} {:>12s} '.format('-', '-', '-', fakerate, '-', '-'))
          logfile_out.write('{:15} {:2d} {:10.5f}\n'.format('jetfakelep', ich, nevt*fakerate))
  file_out.write('\n')
  
  file_out.write('{:15s} {:3s} {:6s}'.format('rare_stat'+str(ich),'lnN',''))
  staterror= 1.0 + h_rates['h_rare_norm'].GetBinError(ich)/h_rates['h_rare_norm'].GetBinContent(ich)
  for k in range(1,n_channels+1):
      if( k == ich):
          file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:12.3f} '.format('-', '-', '-', '-', '-',staterror))
  file_out.write('\n')

  file_out.close()
  logfile_out.close()

