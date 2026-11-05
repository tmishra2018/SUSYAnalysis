import re
import os
import ROOT
from os import system
import sys

br_neu = 0.5
br_cha = 1-br_neu

br_susy = 2*br_neu*br_cha
br_susy = br_susy/0.5
#print br_susy
#br_susy = 1

susy_in = ROOT.TFile('/uscms/home/tmishra/nobackup/signal_trees/signalTree_T5WG_2016.root','read')
n_channels = int(sys.argv[1])*2
#n_channels = 36
syst_names = ['jes','jer','esf','scale','eleshape','jetshape','qcdshape','xs','lumi','isr']

h_SUSYmass = susy_in.Get('SUSYMass')
h_rates = {}
for i in range(1, n_channels + 1):
    h_rates['h_chan' + str(i) + '_rate_nom'] = susy_in.Get('h_chan' + str(i) + '_rate_nom')
    for j in syst_names:
        h_rates['h_chan' + str(i) + '_syserr_' + j] = susy_in.Get('h_chan' + str(i) + '_syserr_' + j)
h_rates['t5wg_h_syserr_PU'] = susy_in.Get('t5wg_h_syserr_PU')

file_template = open('counting_exp_XXX_YYY_{}.txt'.format(n_channels), 'r')
lines = [line for line in file_template.readlines()]

low_p = 2 
high_p = 0

for i in range(1, h_SUSYmass.GetXaxis().GetNbins() + 1):
    for j in range(1, h_SUSYmass.GetYaxis().GetNbins() + 1):
        if(h_SUSYmass.GetBinContent(i,j) <= 0):
            continue
        file_out = open(
            '/tmp/tribeni/counting_t5Wg_'
            + str(int(h_SUSYmass.GetXaxis().GetBinCenter(i))) + '_'
            + str(int(h_SUSYmass.GetYaxis().GetBinCenter(j))) + '.txt', 'w'
            )
        avg_jes = 0
        avg_jer = 0
        avg_esf = 0
        n_nonzero = 0
# Get Average
        for k in range(1, n_channels+1):
            n_nom = h_rates['h_chan' + str(k) + '_rate_nom'].GetBinContent(i, j)
	    if int(h_SUSYmass.GetXaxis().GetBinCenter(i)) == 1750 and int(h_SUSYmass.GetYaxis().GetBinCenter(j)) == 100:
	    	print n_nom
            e_jes = h_rates['h_chan' + str(k) + '_syserr_jes'].GetBinContent(i, j)
            e_jer = h_rates['h_chan' + str(k) + '_syserr_jer'].GetBinContent(i, j)
            e_esf = h_rates['h_chan' + str(k) + '_syserr_esf'].GetBinContent(i, j)
            if n_nom != 0.0:
                avg_jes +=  e_jes / n_nom
                avg_jer +=  e_jer / n_nom 
                avg_esf +=  e_esf / n_nom
                n_nonzero += 1
        avg_jes = avg_jes/n_nonzero
        avg_jer = avg_jer/n_nonzero
        avg_esf = avg_esf/n_nonzero

        for l in lines:
            if (re.search('NSC', l)
                    or re.search('NEVT', l)
                    or re.search('STSC', l)
                    or re.search('STAT', l)
                    or re.search('JESS', l)
                    or re.search('JERS', l)
                    or re.search('ESF', l)
                    or re.search('XSS', l)
                    or re.search('PUS', l)):
                for k in range(1, n_channels+1):
                    n_nom = h_rates['h_chan' + str(k) + '_rate_nom'].GetBinContent(i, j)
                    e_nom = h_rates['h_chan' + str(k) + '_rate_nom'].GetBinError(i, j)
                    e_jes = h_rates['h_chan' + str(k) + '_syserr_jes'].GetBinContent(i, j)
                    e_jer = h_rates['h_chan' + str(k) + '_syserr_jer'].GetBinContent(i, j)
                    e_esf = h_rates['h_chan' + str(k) + '_syserr_esf'].GetBinContent(i, j)
                    unc_PU = h_rates['t5wg_h_syserr_PU'].GetBinContent(i, j) + 1.0
                    unc_stat = 1.0
                    unc_jes = 1.0 
                    unc_jer = 1.0
                    unc_esf = 1.0
                    if n_nom != 0.0:
                        unc_stat = 1.0 +  e_nom/n_nom 
                        unc_jes = 1.0 +  e_jes/n_nom
                        unc_jer = 1.0 +  e_jer/n_nom 
                        unc_esf = 1.0 +  e_esf/n_nom
                    else:
                        unc_stat = 2.0
                        unc_jes = 1.0 + avg_jes 
                        unc_jer = 1.0 + avg_jer 
                        unc_esf = 1.0 + avg_esf 

                    if re.search('NSC', l):
                        l = l.replace('NSC' + str(k) + ' ', str(round(n_nom*br_susy, 3)))
                    elif re.search('STSC', l):
                        l = l.replace('STSC' + str(k)+ ' ', str(round(unc_stat, 3)))
                   # elif re.search('NEVT', l):
                   #     l = l.replace('NEVT' + str(k)+ ' ', str(int(n_nom/n_xs))+' ')
                   #     l = l.replace('STAT' + str(k)+ ' ', str(round(n_xs, 4)))
                    elif re.search('JESS', l):
                        l = l.replace('JESS' + str(k)+ ' ', str(round(unc_jes, 3)))
                    elif re.search('JERS', l):
                        l = l.replace('JERS' + str(k)+ ' ', str(round(unc_jer, 3)))
                    elif re.search('ESF', l):
                        l = l.replace('ESF' + str(k)+ ' ', str(round(unc_esf, 3)))
                    elif re.search('XSS', l):
                        l = l.replace('XSS' + ' ', '-')
                    elif re.search('PUS', l):
                        l = l.replace('PUS' + ' ', str(round(unc_PU , 3)))
                        if(n_nom > 0):
                          if(unc_PU < low_p):
                              low_p = unc_PU
                          if(unc_PU > high_p):
                              high_p = unc_PU

            file_out.write(l)

        file_out.close()
print low_p, high_p

