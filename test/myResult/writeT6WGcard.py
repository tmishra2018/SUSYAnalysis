import re
import os
import ROOT
import sys

n_channels = int(sys.argv[1])*2
RunYear = sys.argv[2]
preVFP = sys.argv[3]

if RunYear == '2016' and preVFP == '1':
    whichVFP = 'preVFP'
elif RunYear == '2016' and preVFP == '0':
    whichVFP = 'postVFP'
else:
    whichVFP = ''

if RunYear == '2016':
    output_file_path = f'/eos/uscms/store/user/tmishra/CombinedLimit/T6WG/cards/cards_{RunYear}{whichVFP}/'
else:
    output_file_path = f'/eos/uscms/store/user/tmishra/CombinedLimit/T6WG/cards/cards_{RunYear}/'

susy_in = ROOT.TFile(f'/uscms/home/tmishra/nobackup/signal_trees/signalTree_T6WG_{RunYear}{whichVFP}.root', 'read')
syst_names = ['jes','jer','esf','scale','eleshape','jetshape','qcdshape','xs','lumi']

h_SUSYmass = susy_in.Get('T6WGSUSYMass')
h_rates = {}
for i in range(1, n_channels + 1):
    h_rates['h_chan' + str(i) + '_rate_nom'] = susy_in.Get('h_chan' + str(i) + '_rate_nom')
    for j in syst_names:
        h_rates['h_chan' + str(i) + '_syserr_' + j] = susy_in.Get('h_chan' + str(i) + '_syserr_' + j)

file_template = open('logs/counting_exp_XXX_YYY_{}_{}{}.txt'.format(n_channels,RunYear,whichVFP), 'r')
lines = [line for line in file_template.readlines()]

for i in range(1, h_SUSYmass.GetXaxis().GetNbins() + 1):
    for j in range(1, h_SUSYmass.GetYaxis().GetNbins() + 1):
        if(h_SUSYmass.GetBinContent(i,j) <= 0):
            continue
        squark_mass = int(h_SUSYmass.GetXaxis().GetBinCenter(i))
        neutralino_mass = int(h_SUSYmass.GetYaxis().GetBinCenter(j))
        file_out = open(f"{output_file_path}counting_t6wg_{squark_mass}_{neutralino_mass}.txt", 'w')

        avg_jes = 0
        avg_jer = 0
        avg_esf = 0
        n_nonzero = 0
# Get Average
        for k in range(1, n_channels+1):
            n_nom = h_rates['h_chan' + str(k) + '_rate_nom'].GetBinContent(i, j)
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
                    unc_stat = 1.0
                    unc_jes = 1.0 
                    unc_jer = 1.0
                    unc_esf = 1.0
                    if n_nom != 0.0:
                        unc_stat = 1.0 + e_nom/n_nom 
                        unc_jes = 1.0 +  e_jes/n_nom
                        unc_jer = 1.0 +  e_jer/n_nom 
                        unc_esf = 1.0 +  e_esf/n_nom
                    else:
                        unc_stat = 2.0
                        unc_jes = 1.0 + avg_jes 
                        unc_jer = 1.0 + avg_jer 
                        unc_esf = 1.0 + avg_esf 

                    if re.search('NSC', l):
                        l = l.replace('NSC' + str(k) + ' ', str(round(n_nom, 3)))
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
                    elif re.search('PUS', l):
                        l = l.replace('PUS' + ' ', str(round(1.05 , 3)))

            file_out.write(l)

        file_out.close()

