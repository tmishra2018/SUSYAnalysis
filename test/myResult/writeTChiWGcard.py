import re
import os
import sys
import ROOT
import argparse

n_channels = int(sys.argv[1]) * 2
RunYear = sys.argv[2]
preVFP = sys.argv[3]

if RunYear == '2016' and preVFP == '1':
    whichVFP = 'preVFP'
elif RunYear == '2016' and preVFP == '0':
    whichVFP = 'postVFP'
else:
    whichVFP = ''

if RunYear == '2016':
    output_file_path = f'/eos/uscms/store/user/tmishra/CombinedLimit/TChiWG/cards/cards_{RunYear}{whichVFP}/'
else:
    output_file_path = f'/eos/uscms/store/user/tmishra/CombinedLimit/TChiWG/cards/cards_{RunYear}/'

susy_in = ROOT.TFile(f'/uscms/home/tmishra/nobackup/signal_trees/signalTree_TChiWG_{RunYear}{whichVFP}.root', 'read')

syst_names = ['jes', 'jer', 'esf', 'scale', 'eleshape', 'jetshape', 'qcdshape', 'xs', 'lumi']
h_SUSYmass = susy_in.Get('h_chan1_rate_nom')
h_rates = {}

for ich in range(1, n_channels + 1):
    h_rates[f'h_chan{ich}_rate_nom'] = susy_in.Get(f'h_chan{ich}_rate_nom')
    for sys in syst_names:
        h_rates[f'h_chan{ich}_syserr_{sys}'] = susy_in.Get(f'h_chan{ich}_syserr_{sys}')

h_rates['tchiwg_h_syserr_PU'] = susy_in.Get('tchiwg_h_syserr_PU')

file_template = open(f'logs/counting_exp_XXX_YYY_{n_channels}_{RunYear}{whichVFP}.txt', 'r')
lines = [line for line in file_template.readlines()]

low_p = 2
high_p = 0

for xbin in range(1, h_SUSYmass.GetXaxis().GetNbins() + 1):
    if h_SUSYmass.GetBinContent(xbin) <= 0:
        continue
    mass1 = h_SUSYmass.GetXaxis().GetBinCenter(xbin)

    file_out = open(f'{output_file_path}counting_tchiwg_{int(mass1)}.txt', 'w')
    
    avg_jes = 0
    avg_jer = 0
    avg_esf = 0
    n_nonzero = 0

    for k in range(1, n_channels + 1):
        n_nom = h_rates[f'h_chan{k}_rate_nom'].GetBinContent(xbin)
        e_jes = h_rates[f'h_chan{k}_syserr_jes'].GetBinContent(xbin)
        e_jer = h_rates[f'h_chan{k}_syserr_jer'].GetBinContent(xbin)
        e_esf = h_rates[f'h_chan{k}_syserr_esf'].GetBinContent(xbin)
        if n_nom != 0.0:
            avg_jes += e_jes / n_nom
            avg_jer += e_jer / n_nom
            avg_esf += e_esf / n_nom
            n_nonzero += 1

    if n_nonzero > 0:
        avg_jes /= n_nonzero
        avg_jer /= n_nonzero
        avg_esf /= n_nonzero

    for l in lines:
        if re.search(r'NSC|NEVT|STSC|STAT|JESS|JERS|ESF|XSS|PUS', l):
            for k in range(1, n_channels + 1):
                n_nom = h_rates[f'h_chan{k}_rate_nom'].GetBinContent(xbin)
                e_nom = h_rates[f'h_chan{k}_rate_nom'].GetBinError(xbin)
                e_jes = h_rates[f'h_chan{k}_syserr_jes'].GetBinContent(xbin)
                e_jer = h_rates[f'h_chan{k}_syserr_jer'].GetBinContent(xbin)
                e_esf = h_rates[f'h_chan{k}_syserr_esf'].GetBinContent(xbin)
                unc_PU = h_rates['tchiwg_h_syserr_PU'].GetBinContent(xbin) + 1.0
                unc_stat = 1.0
                unc_jes = 1.0
                unc_jer = 1.0
                unc_esf = 1.0
                unc_xs = 1.0

                if n_nom != 0.0:
                    unc_stat = 1.0 + e_nom / n_nom
                    unc_jes = 1.0 + e_jes / n_nom
                    unc_jer = 1.0 + e_jer / n_nom
                    unc_esf = 1.0 + e_esf / n_nom
                else:
                    unc_stat = 2.0
                    unc_jes = 1.0 + avg_jes
                    unc_jer = 1.0 + avg_jer
                    unc_esf = 1.0 + avg_esf

                replacements = {
                    f'NSC{k} ': str(round(n_nom, 3)),
                    f'STSC{k} ': str(round(unc_stat, 3)),
                    f'JESS{k} ': str(round(unc_jes, 3)),
                    f'JERS{k} ': str(round(unc_jer, 3)),
                    f'ESF{k} ': str(round(unc_esf, 3)),
                    f'XSS{k} ': str(round(unc_xs, 3)),
                    'PUS ': str(round(unc_PU, 3))
                }

                for key, value in replacements.items():
                    if re.search(key, l):
                        l = l.replace(key, value)

                if n_nom > 0:
                    low_p = min(low_p, unc_PU)
                    high_p = max(high_p, unc_PU)

        file_out.write(l)

    file_out.close()

print(low_p, high_p)
