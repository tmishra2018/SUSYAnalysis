import re
import os
import ROOT

#susy_in = ROOT.TFile('/uscms/home/tmishra/nobackup/signal_trees/signalTree_T5WG_2016.root','read')
susy_in = ROOT.TFile('./signalTree_T5WG.root','read')
n_channels = 36
syst_names = ['stat','jes','jer','esf','scale','eleshape','jetshape','xs','lumi']

h_SUSYmass = susy_in.Get('SUSYMass')
h_rates = {}
for i in range(1, n_channels + 1):
    h_rates['h_chan' + str(i) + '_rate_nom'] = susy_in.Get('h_chan' + str(i) + '_rate_nom')
    for j in syst_names:
        h_rates['h_chan' + str(i) + '_syserr_' + j] = susy_in.Get('h_chan' + str(i) + '_syserr_' + j)

file_template = open('./counting_XXX_YYY.txt', 'r')
lines = [line for line in file_template.readlines()]

for i in range(1, h_SUSYmass.GetXaxis().GetNbins() + 1):
    for j in range(1, h_SUSYmass.GetYaxis().GetNbins() + 1):
        if(h_SUSYmass.GetBinContent(i,j) <= 0):
            continue
        file_out = open(
            './counting/counting_t5Wg_'
            + str(int(h_SUSYmass.GetXaxis().GetBinCenter(i))) + '_'
            + str(int(h_SUSYmass.GetYaxis().GetBinCenter(j))) + '.txt', 'w'
            )

        for l in lines:
            if (re.search('NSC', l)
                    or re.search('STSC', l)
                    or re.search('JESS', l)
                    or re.search('JERS', l)
                    or re.search('ESFS', l)
                    or re.search('XSS', l)):
                for k in range(1, n_channels+1):
                    e_nom = h_rates['h_chan' + str(k) + '_rate_nom'].GetBinError(i, j)
                    n_nom = h_rates['h_chan' + str(k) + '_rate_nom'].GetBinContent(i, j)
                    e_jes = h_rates['h_chan' + str(k) + '_syserr_jes'].GetBinContent(i, j)
                    e_jer = h_rates['h_chan' + str(k) + '_syserr_jer'].GetBinContent(i, j)
                    e_esf = h_rates['h_chan' + str(k) + '_syserr_esf'].GetBinContent(i, j)
                    e_xs  = h_rates['h_chan' + str(k) + '_syserr_xs'].GetBinContent(i, j)
                    if n_nom != 0.0:
                        unc_stat = 1.0 + e_nom / n_nom
                        unc_jes = 1.0 +  e_jes / n_nom
                        unc_jer = 1.0 +  e_jer / n_nom 
                        unc_esf = 1.0 +  e_esf / n_nom
                        unc_xs  = 1.0 +  e_xs  / n_nom
                    else:
                        unc_stat = 1.0
                        unc_jes = 1.0
                        unc_jer = 1.0
                        unc_esf = 1.0
                        unc_xs  = 1.0

                    if re.search('NSC', l):
                        l = l.replace('NSC' + str(k), str(round(n_nom, 3)))
                        #print 'replace NSC',str(k),' ',n_nom
                    elif re.search('STSC', l):
                        l = l.replace('STSC' + str(k), str(round(unc_stat, 3)))
                        #print 'replace STSC',str(k),' ',unc_stat
                    elif re.search('JESS', l):
                        l = l.replace('JESS' + str(k), str(round(unc_jes, 3)))
                    elif re.search('JERS', l):
                        l = l.replace('JERS' + str(k), str(round(unc_jer, 3)))
                    elif re.search('ESFS', l):
                        l = l.replace('ESFS' + str(k), str(round(unc_esf, 3)))
                    elif re.search('XSS', l):
                        l = l.replace('XSS' + str(k), str(round(unc_xs , 3)))

            file_out.write(l)

        file_out.close()

