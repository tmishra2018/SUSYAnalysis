import re
import os
import ROOT
import argparse
import math

limdir = './'

n_channels = 36
pro_names = ['elefakepho', 'jetfakepho', 'qcdfakelep','VGamma','rare']
syst_names = ['stat','jes','jer','esf','scale','e_to_pho_syst','j_to_pho_syst','fakelep_shape','xs','lumi','isr']

file_in = ROOT.TFile('./SignalSystematic.root', 'read')
h_rates = {}
for hname in pro_names:
    h_rates['h_' + hname + '_norm'] = file_in.Get('h_' + hname + '_norm')
    for sys in syst_names:
        if(sys != 'stat'):
            h_rates['h_' + hname + '_syserr_' + sys] = file_in.Get('h_' + hname + '_syserr_' + sys)
h_rates['h_elefakepho_controlsample'] = file_in.Get('h_elefakepho_controlsample')
h_rates['h_elefakepho_transferfactor'] = file_in.Get('h_elefakepho_transferfactor')
h_rates['h_jetfakepho_controlsample'] = file_in.Get('h_jetfakepho_controlsample')
h_rates['h_jetfakepho_transferfactor'] = file_in.Get('h_jetfakepho_transferfactor')
h_rates['h_qcdfakelep_controlsample'] = file_in.Get('h_qcdfakelep_controlsample')
h_rates['h_qcdfakelep_transferfactor'] = file_in.Get('h_qcdfakelep_transferfactor')
h_rates['h_sig'] = file_in.Get('h_sig')

file_out = open( limdir+'unblinding_XXX_YYY.txt', 'w')

for ich in range(1,n_channels+1):
  channelrate = 0
  channelerror = 0
  for p in pro_names:
      file_out.write('{:2s} {:2d} '.format('ch',ich))
      file_out.write('{:12s} '.format(p))
      evtrate = h_rates['h_' + p + '_norm'].GetBinContent(ich)
      channelrate += evtrate
      file_out.write('{:7.3f} '.format(evtrate))
      totalerror = 0
      for syst in syst_names:
          syserr = 0
          if(h_rates['h_' + p + '_norm'].GetBinContent(ich)==0):
              if(syst == 'stat'):
                transferF = 0
                normN = h_rates['h_' + p + '_controlsample'].GetBinContent(ich)
                alpha = 1.0-0.6827
                if(p == 'elefakepho'):
                    transferF = h_rates['h_elefakepho_transferfactor'].GetBinContent(ich)
                elif(p == 'jetfakepho'):
                    transferF = h_rates['h_jetfakepho_transferfactor'].GetBinContent(ich)
                elif(p == 'qcdfakelep'):
                    transferF = h_rates['h_qcdfakelep_transferfactor'].GetBinContent(ich)
               # upperN = ROOT.Math.gamma_quantile(alpha/2, normN+1, transferF)
               # syserr = upperN
                upperN = 1.29 
                syserr = upperN*transferF
                print upperN, ' ', transferF
                file_out.write('{:7.3f} '.format(syserr))
              else:
                syserr = 0
                file_out.write('{:7d} '.format(-1))
                continue
              if(syserr > 0):
                  totalerror += syserr*syserr 
          else:
              if(syst != 'stat'):
                   syserr = h_rates['h_' + p + '_syserr_' + syst].GetBinContent(ich)/h_rates['h_' + p + '_norm'].GetBinContent(ich)
              else:
                   syserr = h_rates['h_' + p + '_norm'].GetBinError(ich)/h_rates['h_' + p + '_norm'].GetBinContent(ich)
              if(syserr < 0):
                   file_out.write('{:7d} '.format(-1))
              else:
                   file_out.write('{:7.3f} '.format(syserr))
              if(syserr > 0):
                   totalerror += evtrate*syserr*evtrate*syserr
      channelerror += totalerror 
      file_out.write('  {:12.3f} '.format(math.sqrt(totalerror)))
      file_out.write('\n')
  file_out.write('chan {:2d} {:4d} {:10.3f} +- {:7.3f} {:10.3f}\n'.format(ich,int(h_rates['h_sig'].GetBinContent(ich)),channelrate, math.sqrt(channelerror), (h_rates['h_sig'].GetBinContent(ich) - channelrate)/math.sqrt(channelerror + h_rates['h_sig'].GetBinContent(ich))))  

file_out.close()

