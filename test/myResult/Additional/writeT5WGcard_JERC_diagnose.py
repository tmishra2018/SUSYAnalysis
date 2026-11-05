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
import ROOT
from os import system
import sys

# --- Diagnostics for large signal systematics ---
JEC_WARN = 1.99   # ~100%
JER_WARN = 1.99

SIG_PDF_UNC = 1.05

JEC_100_cases = {}
JER_100_cases = {}

FLAT_PU_UNC = {
    '2016':        1.022,
    '2016preVFP':  1.022,
    '2016postVFP': 1.022,
    '2017':        1.023,
    '2018':        1.025,
}

def classify_instability(n_nom):
    if n_nom < 1e-4:
        return "ultra-low yield"
    elif n_nom < 1e-2:
        return "low yield"
    else:
        return "shape-driven"


br_susy = 1
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
    base_path = f'/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_{RunYear}{whichVFP}'
else:
    base_path = f'/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_{RunYear}'

output_file_path = f'{base_path}_{obs_mode}/'
os.makedirs(output_file_path, exist_ok=True)

susy_in = ROOT.TFile(f'/uscms/home/tmishra/nobackup/signal_trees/signalTree_T5WG_{RunYear}{whichVFP}.root', 'read')
syst_names = ['jes','jer','esf','scale','eleshape','jetshape','qcdshape','xs','lumi','isr','qcdscale']

# Cross section uncertainty for rare process
XSEC_RARE_UNC = 1.50

qcdscale_dc_min, qcdscale_dc_max = 1e9, 0.0

jes_dc_min, jes_dc_max = 1e9, 0.0
jer_dc_min, jer_dc_max = 1e9, 0.0
# "Filtered" versions exclude bins flagged as >= JEC_WARN/JER_WARN (the ~100% low-stat
# instability cases) - used for the quoted range in the final summary.
jes_dc_min_filt, jes_dc_max_filt = 1e9, 0.0
jer_dc_min_filt, jer_dc_max_filt = 1e9, 0.0
esf_dc_min, esf_dc_max = 1e9, 0.0
pu_dc_min,  pu_dc_max  = 1e9, 0.0
xs_dc_min,  xs_dc_max  = 1e9, 0.0

h_SUSYmass = susy_in.Get('SUSYMass')
h_rates = {}
for i in range(1, n_channels + 1):
    h_rates['h_chan' + str(i) + '_rate_nom'] = susy_in.Get('h_chan' + str(i) + '_rate_nom')
    for j in syst_names:
        h_rates['h_chan' + str(i) + '_syserr_' + j] = susy_in.Get('h_chan' + str(i) + '_syserr_' + j)
h_rates['t5wg_h_syserr_PU'] = susy_in.Get('t5wg_h_syserr_PU')

for i in range(1, n_channels + 1):
    h_rates['h_chan' + str(i) + '_syserr_pu'] = \
        susy_in.Get('h_chan' + str(i) + '_syserr_pu')

# NEW: template filename includes obs_mode to match createDataCard.py output
file_template = open('logs/counting_exp_XXX_YYY_{}_{}_{}{}.txt'.format(obs_mode, n_channels, RunYear, whichVFP), 'r')
lines = [line for line in file_template.readlines()]

# NEW: for splusb mode, parse the bonly template to get per-bin background sums.
# These are used below to compute observation = signal + background per mass point.
# Requires the bonly template to have been produced first by createDataCard.py bonly.
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

low_p = 2 
high_p = 0
count = 0
        
for i in range(1, h_SUSYmass.GetXaxis().GetNbins() + 1):
    for j in range(1, h_SUSYmass.GetYaxis().GetNbins() + 1):
        if h_SUSYmass.GetBinContent(i, j) <= 0:
            continue
        count += 1
        gluino_mass = int(h_SUSYmass.GetXaxis().GetBinCenter(i))
        neutralino_mass = int(h_SUSYmass.GetYaxis().GetBinCenter(j))
        signal_proc = f't5wg_{gluino_mass}_{neutralino_mass}'
        # Datacard filename unchanged; only the folder differs
        file_out = open(f"{output_file_path}counting_t5wg_{gluino_mass}_{neutralino_mass}.txt", 'w')
        avg_jes = 0
        avg_jer = 0
        avg_esf = 0
        avg_xs = 0
        n_nonzero = 0

        for k in range(1, n_channels+1):
            n_nom = h_rates['h_chan' + str(k) + '_rate_nom'].GetBinContent(i, j)
            e_jes = h_rates['h_chan' + str(k) + '_syserr_jes'].GetBinContent(i, j)
            e_jer = h_rates['h_chan' + str(k) + '_syserr_jer'].GetBinContent(i, j)
            e_esf = h_rates['h_chan' + str(k) + '_syserr_esf'].GetBinContent(i, j)
            e_xs = h_rates['h_chan' + str(k) + '_syserr_xs'].GetBinContent(i, j)

            if n_nom != 0.0:
                avg_jes +=  e_jes / n_nom
                avg_jer +=  e_jer / n_nom 
                avg_esf +=  e_esf / n_nom
                avg_xs +=  e_xs / n_nom
                n_nonzero += 1
                
        if n_nonzero > 0:
            avg_jes = avg_jes/n_nonzero
            avg_jer = avg_jer/n_nonzero
            avg_esf = avg_esf/n_nonzero
            avg_xs = avg_xs/n_nonzero
        
        else:
            avg_jes = 0
            avg_jer = 0
            avg_esf = 0
            avg_xs = 0

        for l in lines:
            l = l.replace('SUSY',   signal_proc)
            if (re.search('NSC', l)
                    or re.search('NEVT', l)
                    or re.search('STSC', l)
                    or re.search('STAT', l)
                    or re.search('JESS', l)
                    or re.search('JERS', l)
                    or re.search('ESF', l)
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
                    e_xs  = h_rates['h_chan' + str(k) + '_syserr_xs'].GetBinContent(i, j)

                    unc_PU = h_rates['t5wg_h_syserr_PU'].GetBinContent(i, j) + 1.0
                    unc_stat = 1.0
                    unc_jes = 1.0 
                    unc_jer = 1.0
                    unc_esf = 1.0
                    unc_isr = 1.0
                    unc_xs = 1.0
                    unc_qcdscale = 1.0

                    key = (RunYear, gluino_mass, neutralino_mass, k)
                    if n_nom != 0.0:
                        unc_stat = 1.0 +  e_nom/n_nom 
                        unc_jes = 1.0 +  e_jes/n_nom
                        unc_jer = 1.0 +  e_jer/n_nom 
                        unc_esf = 1.0 +  e_esf/n_nom
                        unc_xs  = 1.0 +  e_xs/n_nom
                        unc_isr = 1.0 + h_rates['h_chan' + str(k) + '_syserr_isr'].GetBinContent(i, j)
                        qcdscale_val = h_rates['h_chan' + str(k) + '_syserr_qcdscale'].GetBinContent(i, j)
                        unc_qcdscale = 1.0 + qcdscale_val if qcdscale_val >= 0 else 1.0
                        if unc_jes >= JEC_WARN:
                            JEC_100_cases[key] = (n_nom, unc_jes, classify_instability(n_nom))
                        if unc_jer >= JER_WARN:
                            JER_100_cases[key] = (n_nom, unc_jer, classify_instability(n_nom))

                    else:
                        unc_stat = 2.0
                        unc_jes = 1.0 + avg_jes 
                        unc_jer = 1.0 + avg_jer 
                        unc_esf = 1.0 + avg_esf 
                        unc_xs  = 1.0 + avg_xs
                        unc_isr = -1
                        unc_qcdscale = 1.10
                        if unc_jes >= JEC_WARN:
                            JEC_100_cases[key] = (n_nom, unc_jes, classify_instability(n_nom))
                        if unc_jer >= JER_WARN:
                            JER_100_cases[key] = (n_nom, unc_jer, classify_instability(n_nom))

                    if re.search('NSC', l):
                        l = l.replace('NSC' + str(k) + ' ', str(round(n_nom*br_susy, 3)))
                    elif re.search('STSC', l):
                        l = l.replace('STSC' + str(k)+ ' ', str(round(unc_stat, 3)))
                    elif re.search('JESS', l):
                        l = l.replace('JESS' + str(k)+ ' ', str(round(unc_jes, 3)))
                    elif re.search('JERS', l):
                        l = l.replace('JERS' + str(k)+ ' ', str(round(unc_jer, 3)))
                    elif re.search('ESF', l):
                        l = l.replace('ESF' + str(k)+ ' ', str(round(unc_esf, 3)))
                    elif re.search('XSS', l):
                        l = l.replace('XSS' + str(k), str(round(unc_xs , 3)))
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
                        splusb_val = round(n_nom_sig * br_susy + bkg_sum, 6)
                        l = l.replace('SPLUSB' + str(k) + ' ', str(splusb_val) + ' ')
                    elif re.search('PDFS', l):
                        l = l.replace('PDFS' + str(k)+ ' ', str(round(SIG_PDF_UNC, 3)) + ' ')

                    jes_dc_min = min(jes_dc_min, unc_jes)
                    jes_dc_max = max(jes_dc_max, unc_jes)
                    if unc_jes < JEC_WARN:
                        jes_dc_min_filt = min(jes_dc_min_filt, unc_jes)
                        jes_dc_max_filt = max(jes_dc_max_filt, unc_jes)

                    jer_dc_min = min(jer_dc_min, unc_jer)
                    jer_dc_max = max(jer_dc_max, unc_jer)
                    if unc_jer < JER_WARN:
                        jer_dc_min_filt = min(jer_dc_min_filt, unc_jer)
                        jer_dc_max_filt = max(jer_dc_max_filt, unc_jer)

                    esf_dc_min = min(esf_dc_min, unc_esf)
                    esf_dc_max = max(esf_dc_max, unc_esf)

                    pu_dc_min = min(pu_dc_min, unc_PU)
                    pu_dc_max = max(pu_dc_max, unc_PU)

                    qcdscale_dc_min = min(qcdscale_dc_min, unc_qcdscale)
                    qcdscale_dc_max = max(qcdscale_dc_max, unc_qcdscale)
            file_out.write(l)

        bkg_dashes = ('{:>12s} ' * (n_processes - 1)).format(*(['-'] * (n_processes - 1)))
        sig_qcd_list = [f'QCDscale_{signal_proc}']
        for row_name in sig_qcd_list:
            file_out.write('{:30s} {:3s} {:6s}'.format(row_name, 'lnN', ''))
            for ich in range(1, n_channels + 1):
                file_out.write('{:12.3f} '.format(round(unc_qcdscale, 3)))
                file_out.write(bkg_dashes)
            file_out.write('\n') # End the row
        file_out.close()

# Integrated luminosity is applied as a flat, year-dependent uncertainty
# (not bin-wise), so its "range" collapses to a single value.
lumi_key = f'{RunYear}{whichVFP}'
flat_lumi_unc = FLAT_PU_UNC.get(lumi_key, 1.0)

print("\nSystematic uncertainties (%) on SUSY signal (T5WG):")
print("(Ranges are the min/max of the average bin-wise uncertainty across all signal regions;")
print(" luminosity and PDF are applied flat, so they show a single value.)")
print(f"  Integrated luminosity : {(flat_lumi_unc - 1)*100:.1f} %")
n_jes_excluded = len(JEC_100_cases)
n_jer_excluded = len(JER_100_cases)

def _fmt_filtered_range(lo, hi, n_excluded, warn_thresh):
    if lo > hi:
        return f"N/A - all bins flagged as low-stat (>= {(warn_thresh-1)*100:.0f}%)"
    tag = f"  ({n_excluded} low-stat bin(s) >= {(warn_thresh-1)*100:.0f}% excluded)" if n_excluded else ""
    return f"{(lo-1)*100:.1f} - {(hi-1)*100:.1f} %{tag}"

print(f"  Jet energy resolution : {_fmt_filtered_range(jer_dc_min_filt, jer_dc_max_filt, n_jer_excluded, JER_WARN)}")
print(f"  Jet energy scale      : {_fmt_filtered_range(jes_dc_min_filt, jes_dc_max_filt, n_jes_excluded, JEC_WARN)}")
print(f"  PDF                   : {(SIG_PDF_UNC - 1)*100:.1f} %")
print(f"  Pileup                : {(pu_dc_min - 1)*100:.1f} - {(pu_dc_max - 1)*100:.1f} %")
print(f"  QCD scale             : {(qcdscale_dc_min - 1)*100:.1f} - {(qcdscale_dc_max - 1)*100:.1f} %")

# ---- Diagnostics: which (year, gluino mass, neutralino mass, channel) bins hit the
#      ~100% JES/JER flag (unc >= JEC_WARN/JER_WARN), and why (instability class, n_nom).
#      These dicts were already being filled above; this just surfaces them.
from collections import Counter

def _dump_100pct_cases(label, cases_dict):
    if not cases_dict:
        print(f"\n  No {label} cases hit the >= {(JEC_WARN-1)*100:.0f}% flag.")
        return
    class_counts = Counter(v[2] for v in cases_dict.values())
    print(f"\n  {label}: {len(cases_dict)} channel/mass-point bins hit the >= {(JEC_WARN-1)*100:.0f}% flag")
    for cls, cnt in sorted(class_counts.items(), key=lambda kv: -kv[1]):
        print(f"    {cls:>16s} : {cnt}")
    worst = sorted(cases_dict.items(), key=lambda kv: -kv[1][1])[:5]
    print(f"  Worst {label} cases (year, gluino, neutralino, channel) -> n_nom, unc, class:")
    for key, (n_nom, unc, cls) in worst:
        print(f"    {key} -> n_nom={n_nom:.4g}, unc={(unc-1)*100:.1f}%, {cls}")

_dump_100pct_cases("JES", JEC_100_cases)
_dump_100pct_cases("JER", JER_100_cases)
