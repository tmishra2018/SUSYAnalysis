import re
import os
import ROOT
import math

outdir = '/uscmst1b_scratch/lpc1/3DayLifetime/mengleis/M1M3/'

susy_in = ROOT.TFile('../test_M1M3.root','read')
n_channels = 36
syst_names = ['jes','phoesf','lepesf','lumi','isr','pdf','met']

h_SUSYmass = susy_in.Get('SUSYMass')
h_rates = {}
for i in range(1, n_channels + 1):
    h_rates['ALL_chan' + str(i) + '_rate_nom'] = susy_in.Get('ALL_chan' + str(i) + '_rate_nom')
    h_rates['ALL_chan' + str(i) + '_rate_control'] = susy_in.Get('ALL_chan' + str(i) + '_rate_control')
    for j in syst_names:
        h_rates['susy_chan' + str(i) + '_syserr_' + j] = susy_in.Get('susy_chan' + str(i) + '_syserr_' + j)
h_rates['susy_syserr_PU'] = susy_in.Get('susy_syserr_PU')

file_template = open('./newcounting_exp_XXX_YYY.txt', 'r')
lines = [line for line in file_template.readlines()]

if not os.path.exists(outdir):
    os.makedirs(outdir)
for i in range(1, h_SUSYmass.GetXaxis().GetNbins() + 1):
    for j in range(1, h_SUSYmass.GetYaxis().GetNbins() + 1):
        if(h_SUSYmass.GetBinContent(i,j) <= 0):
            continue
        file_out = open(
            outdir + 'counting_GMSB_'          
            + str(int(h_SUSYmass.GetXaxis().GetBinCenter(i))) + '_'
            + str(int(h_SUSYmass.GetYaxis().GetBinCenter(j))) + '.txt', 'w'
            )

        for l in lines:
            if (re.search('NSC', l)
                    or re.search('STSC', l)
                    or re.search('JES', l)
                    or re.search('ISR', l)
                    or re.search('PhotonSF', l)
                    or re.search('ElectronSF', l)
                    or re.search('PU', l)
                    or re.search('PDFSCALE', l)
                    or re.search('GenMet', l)):
                for k in range(1, n_channels+1):
         
                    binn = k
                    if( k <= 18): binn = k + 18
                    elif( k > 18): binn = k - 18
                    n_nom = h_rates['ALL_chan' + str(k) + '_rate_nom'].GetBinContent(i, j)
                    e_nom = math.sqrt(h_rates['ALL_chan' + str(k) + '_rate_control'].GetBinContent(i, j))
                    unc_phoesf = 1+ h_rates['susy_chan' + str(k) + '_syserr_phoesf'].GetBinContent(i, j)
                    unc_lepesf = 1+ h_rates['susy_chan' + str(k) + '_syserr_lepesf'].GetBinContent(i, j)
                    unc_met =  1+ h_rates['susy_chan' + str(k) + '_syserr_met'].GetBinContent(i, j) 
                    unc_isr =  1+ h_rates['susy_chan' + str(k) + '_syserr_isr'].GetBinContent(i, j) 
                    unc_PU = 1+ h_rates['susy_syserr_PU'].GetBinContent(i, j)
                    unc_stat = 1.0
                    unc_jes = 1.0 + h_rates['susy_chan' + str(k) + '_syserr_jes'].GetBinContent(i, j) 
                    if n_nom != 0.0:
                        unc_stat = 1 + e_nom/h_rates['ALL_chan' + str(k) + '_rate_control'].GetBinContent(i, j) 
                    else:
                        unc_stat = 2.0

                    if re.search('NSC', l):
                        l = l.replace('NSC' + str(binn) + ' ', str(round(n_nom, 3)))
                    elif re.search('STSC', l):
                        l = l.replace('STSC' + str(binn)+ ' ', str(round(unc_stat, 3)))
                    elif re.search('JES', l):
                        l = l.replace('JES' + str(binn)+ ' ', str(round(unc_jes, 3)))
                    elif re.search('PhotonSF', l):
                        l = l.replace('PhotonSF' + str(binn)+ ' ', str(round(unc_phoesf, 3)))
                    elif re.search('ElectronSF', l):
                        l = l.replace('ElectronSF' + str(binn)+ ' ', str(round(unc_lepesf, 3)))
                    elif re.search('PU', l):
                        l = l.replace('PU' +  str(binn)+ ' ', str(round(unc_PU , 3)))
                    elif re.search('PDFSCALE', l):
                        l = l.replace('PDFSCALE' + str(binn)+ ' ', '-') 
                    elif re.search('GenMet', l):
                        l = l.replace('GenMet' + str(binn)+ ' ', str(round(unc_met , 3)))
                    elif re.search('ISR', l):
                        l = l.replace('ISR' + str(binn)+ ' ', str(round(unc_isr , 3)))

            file_out.write(l)

        file_out.close()

