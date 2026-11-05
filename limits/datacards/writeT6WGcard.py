import re
# QCD scale and PDF nuisance rows are written by createDataCard.py with real
# numeric kappas for VGamma and rare backgrounds (bin-wise, from histograms).
# The SUSY column in those rows is written as a placeholder (QCDS{k}, PDFS{k})
# and substituted here with a flat 5% uncertainty (1.05) per bin, so that the
# signal uncertainty is correctly applied per mass point rather than being a
# constant '-' in the template.
# Two combined nuisance rows in the template:
#   QCDscale  lnN  QCDS{k} ... kappa_VG ... kappa_rare  (SUSY placeholder; bkg bin-wise)

import os
# run inside test/myResult (reads/writes its logs/: data logs, counting_exp templates)
os.chdir(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../test/myResult')))
import ROOT
import sys

n_channels = int(sys.argv[1])*2
n_processes = 6   # SUSY + 5 backgrounds (same as createDataCard.py)
RunYear = sys.argv[2]
preVFP = sys.argv[3]
obs_mode = sys.argv[4] if len(sys.argv) > 4 else 'bonly'   # NEW: unblind | bonly | splusb

if obs_mode not in ('unblind', 'bonly', 'splusb'):
    print(f"ERROR: obs_mode must be one of: unblind, bonly, splusb  (got '{obs_mode}')")
    sys.exit(1)

if RunYear == '2016' and preVFP == '1':
    whichVFP = 'preVFP'
elif RunYear == '2016' and preVFP == '0':
    whichVFP = 'postVFP'
else:
    whichVFP = ''

# NEW: output folder is mode-specific; datacard filename stays the same
if RunYear == '2016':
    base_path = f'/eos/uscms/store/user/tmishra/CombinedLimit/T6WG/cards/cards_{RunYear}{whichVFP}'
else:
    base_path = f'/eos/uscms/store/user/tmishra/CombinedLimit/T6WG/cards/cards_{RunYear}'

output_file_path = f'{base_path}_{obs_mode}/'
os.makedirs(output_file_path, exist_ok=True)

susy_in = ROOT.TFile(f'/uscms/home/tmishra/nobackup/signal_trees/signalTree_T6WG_{RunYear}{whichVFP}.root', 'read')
syst_names = ['jes','jer','esf','esfpho','scale','eleshape','jetshape','qcdshape','xs','lumi','isr','qcdscale']

SIG_PDF_UNC = 1.05

qcdscale_dc_min, qcdscale_dc_max = 1e9, 0.0

jes_dc_min, jes_dc_max = 1e9, 0.0
jer_dc_min, jer_dc_max = 1e9, 0.0
esf_dc_min, esf_dc_max = 1e9, 0.0
esfpho_dc_min, esfpho_dc_max = 1e9, 0.0
pu_dc_min,  pu_dc_max  = 1e9, 0.0
xs_dc_min,  xs_dc_max  = 1e9, 0.0


h_SUSYmass = susy_in.Get('T6WGSUSYMass')
h_rates = {}
for i in range(1, n_channels + 1):
    h_rates['h_chan' + str(i) + '_rate_nom'] = susy_in.Get('h_chan' + str(i) + '_rate_nom')
    for j in syst_names:
        h_rates['h_chan' + str(i) + '_syserr_' + j] = susy_in.Get('h_chan' + str(i) + '_syserr_' + j)
h_rates['t6wg_h_syserr_PU'] = susy_in.Get('t6wg_h_syserr_PU')

for i in range(1, n_channels + 1):
    h_rates['h_chan' + str(i) + '_syserr_pu'] = \
        susy_in.Get('h_chan' + str(i) + '_syserr_pu')

# NEW: template filename includes obs_mode to match createDataCard.py output
file_template = open('logs/counting_exp_XXX_YYY_{}_{}_{}{}.txt'.format(obs_mode, n_channels, RunYear, whichVFP), 'r')
lines = [line for line in file_template.readlines()]

# NEW: for splusb mode, parse the bonly template to get per-bin background sums.
# T6WG is 2D like T5WG; background sums are the same regardless of signal model.
bonly_obs_per_bin = {}
if obs_mode == 'splusb':
    bonly_tmpl_path = 'logs/counting_exp_XXX_YYY_bonly_{}_{}{}.txt'.format(n_channels, RunYear, whichVFP)
    with open(bonly_tmpl_path, 'r') as bonly_tmpl:
        for bline in bonly_tmpl:
            if bline.strip().startswith('observation'):
                vals = bline.split()[1:]   # drop the 'observation' token
                for idx, v in enumerate(vals, start=1):
                    bonly_obs_per_bin[idx] = float(v)
                break

for i in range(1, h_SUSYmass.GetXaxis().GetNbins() + 1):
    for j in range(1, h_SUSYmass.GetYaxis().GetNbins() + 1):
        if(h_SUSYmass.GetBinContent(i,j) <= 0):
            continue
        squark_mass = int(h_SUSYmass.GetXaxis().GetBinCenter(i))
        neutralino_mass = int(h_SUSYmass.GetYaxis().GetBinCenter(j))
        signal_proc = f't6wg_{squark_mass}_{neutralino_mass}'
        # Datacard filename unchanged; only the folder differs
        file_out = open(f"{output_file_path}counting_t6wg_{squark_mass}_{neutralino_mass}.txt", 'w')

        qcd_by_k = {}   # per-bin signal QCD-scale kappa, filled in the channel loop below
        avg_jes = 0
        avg_jer = 0
        avg_esf = 0
        avg_esfpho = 0
        n_nonzero = 0
        # Get Average
        for k in range(1, n_channels+1):
            n_nom = h_rates['h_chan' + str(k) + '_rate_nom'].GetBinContent(i, j)
            e_jes = h_rates['h_chan' + str(k) + '_syserr_jes'].GetBinContent(i, j)
            e_jer = h_rates['h_chan' + str(k) + '_syserr_jer'].GetBinContent(i, j)
            e_esf = h_rates['h_chan' + str(k) + '_syserr_esf'].GetBinContent(i, j)
            e_esfpho = h_rates['h_chan' + str(k) + '_syserr_esfpho'].GetBinContent(i, j)
            if n_nom != 0.0:
                avg_jes +=  e_jes / n_nom
                avg_jer +=  e_jer / n_nom 
                avg_esf +=  e_esf / n_nom
                avg_esfpho +=  e_esfpho / n_nom
                n_nonzero += 1

            if n_nonzero > 0:
                avg_jes = avg_jes / n_nonzero
                avg_jer = avg_jer / n_nonzero
                avg_esf = avg_esf / n_nonzero
                avg_esfpho = avg_esfpho / n_nonzero
            else:
                avg_jes = 0
                avg_jer = 0
                avg_esf = 0
                avg_esfpho = 0

        for l in lines:
            l = l.replace('SUSY',   signal_proc)
            if (re.search('NSC', l)
                    or re.search('NEVT', l)
                    or re.search('STSC', l)
                    or re.search('STAT', l)
                    or re.search('JESS', l)
                    or re.search('JERS', l)
                    or re.search('ESF', l)
                    or re.search('PHOSF', l)
                    or re.search('XSS', l)
                    or re.search('ISR', l)
                    or re.search('PUS', l)
                    or re.search('PDFS', l)
                    or re.search('SPLUSB', l)):
                for k in range(1, n_channels+1):
                    n_nom = h_rates['h_chan' + str(k) + '_rate_nom'].GetBinContent(i, j)
                    e_nom = h_rates['h_chan' + str(k) + '_rate_nom'].GetBinError(i, j)
                    e_jes = h_rates['h_chan' + str(k) + '_syserr_jes'].GetBinContent(i, j)
                    e_jer = h_rates['h_chan' + str(k) + '_syserr_jer'].GetBinContent(i, j)
                    e_esf = h_rates['h_chan' + str(k) + '_syserr_esf'].GetBinContent(i, j)
                    e_esfpho = h_rates['h_chan' + str(k) + '_syserr_esfpho'].GetBinContent(i, j)

                    unc_PU = h_rates['t6wg_h_syserr_PU'].GetBinContent(i, j) + 1.0
                    unc_stat = 1.0
                    unc_jes = 1.0 
                    unc_jer = 1.0
                    unc_esf = 1.0
                    unc_esfpho = 1.0
                    unc_isr = 1.0
                    unc_qcdscale = 1.0

                    if n_nom != 0.0:
                        unc_stat = 1.0 + e_nom/n_nom 
                        unc_jes = 1.0 +  e_jes/n_nom
                        unc_jer = 1.0 +  e_jer/n_nom 
                        unc_esf = 1.0 +  e_esf/n_nom
                        unc_esfpho = 1.0 +  e_esfpho/n_nom
                        unc_isr = 1.0 + h_rates['h_chan' + str(k) + '_syserr_isr'].GetBinContent(i, j)
                        qcdscale_val = h_rates['h_chan' + str(k) + '_syserr_qcdscale'].GetBinContent(i, j)
                        unc_qcdscale = 1.0 + qcdscale_val if qcdscale_val >= 0 else 1.0
                    else:
                        unc_stat = 2.0
                        unc_jes = 1.0 + avg_jes 
                        unc_jer = 1.0 + avg_jer 
                        unc_esf = 1.0 + avg_esf 
                        unc_esfpho = 1.0 + avg_esfpho 
                        unc_isr = -1
                        unc_qcdscale = 1.10

                    if re.search('NSC', l):
                        l = l.replace('NSC' + str(k) + ' ', str(round(n_nom, 3)))
                    elif re.search('STSC', l):
                        l = l.replace('STSC' + str(k)+ ' ', str(round(unc_stat, 3)))
                    elif re.search('JESS', l):
                        l = l.replace('JESS' + str(k)+ ' ', str(round(unc_jes, 3)))
                    elif re.search('JERS', l):
                        l = l.replace('JERS' + str(k)+ ' ', str(round(unc_jer, 3)))
                    elif re.search('ESF', l):
                        l = l.replace('ESF' + str(k)+ ' ', str(round(unc_esf, 3)))
                    elif re.search('PHOSF', l):
                        l = l.replace('PHOSF' + str(k)+ ' ', str(round(unc_esfpho, 3)))
                    elif re.search('ISR', l):
                        if(unc_isr == 0 or unc_isr < 0):
                          l = l.replace('ISR' + str(k)+ ' ', '- ')
                        else:
                          l = l.replace('ISR' + str(k)+ ' ', str(round(unc_isr,3)))
                    elif re.search('PUS', l):
                        n_nom  = h_rates['h_chan' + str(k) + '_rate_nom'].GetBinContent(i, j)
                        pu_err = h_rates['h_chan' + str(k) + '_syserr_pu'].GetBinContent(i, j)

                        if n_nom > 0 and pu_err >= 0:
                            unc_PU = round(pu_err / n_nom + 1.0, 3)
                        elif pu_err == -1:
                            unc_PU = 1.0
                        else:
                            unc_PU = 1.0
                        l = l.replace('PUS' + str(k) + ' ', str(unc_PU) + ' ')
                    elif re.search('SPLUSB', l):
                        # observation = signal_yield + total_background for this bin
                        n_nom_sig = h_rates['h_chan' + str(k) + '_rate_nom'].GetBinContent(i, j)
                        bkg_sum = bonly_obs_per_bin.get(k, 0.0)
                        splusb_val = round(n_nom_sig + bkg_sum, 6)
                        l = l.replace('SPLUSB' + str(k) + ' ', str(splusb_val) + ' ')
                    elif re.search('PDFS', l):
                        l = l.replace('PDFS' + str(k)+ ' ', str(round(SIG_PDF_UNC, 3)) + ' ')

                    jes_dc_min = min(jes_dc_min, unc_jes)
                    jes_dc_max = max(jes_dc_max, unc_jes)

                    qcd_by_k[k] = unc_qcdscale
                    jer_dc_min = min(jer_dc_min, unc_jer)
                    jer_dc_max = max(jer_dc_max, unc_jer)

                    esf_dc_min = min(esf_dc_min, unc_esf)
                    esfpho_dc_min = min(esfpho_dc_min, unc_esfpho)
                    esf_dc_max = max(esf_dc_max, unc_esf)
                    esfpho_dc_max = max(esfpho_dc_max, unc_esfpho)

                    pu_dc_min = min(pu_dc_min, unc_PU)
                    pu_dc_max = max(pu_dc_max, unc_PU)
            file_out.write(l)

        # --- Signal-specific QCD scale and PDF nuisance rows ---
        bkg_dashes = ('{:>12s} ' * (n_processes - 1)).format(*(['-'] * (n_processes - 1)))
        sig_qcd_list = [f'QCDscale_{signal_proc}']
        for row_name in sig_qcd_list:
            file_out.write('{:30s} {:3s} {:6s}'.format(row_name, 'lnN', ''))
            for ich in range(1, n_channels + 1):
                # per-bin value of channel ich (was: unc_qcdscale left over from the last channel, bin 36, for all bins)
                file_out.write('{:12.3f} '.format(round(qcd_by_k.get(ich, 1.0), 3)))
                file_out.write(bkg_dashes)
            file_out.write('\n') # End the row
        file_out.close()

print("\nSyst. unc. ranges (T6WG):")
print(f"  JES : {(jes_dc_min-1)*100:.1f} - {(jes_dc_max-1)*100:.1f} %")
print(f"  JER : {(jer_dc_min-1)*100:.1f} - {(jer_dc_max-1)*100:.1f} %")
print(f"  ESF : {(esf_dc_min-1)*100:.1f} - {(esf_dc_max-1)*100:.1f} %")
print(f"  PHOSF : {(esfpho_dc_min-1)*100:.1f} - {(esfpho_dc_max-1)*100:.1f} %")
print(f"  PU  : {(pu_dc_min -1)*100:.1f} - {(pu_dc_max -1)*100:.1f} %")
