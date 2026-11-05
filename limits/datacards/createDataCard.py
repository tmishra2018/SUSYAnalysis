####    script to produce a format for datacard, bin wise entries for each background is stored in SignalSystematic.root
####    lnN stands for Log-normal recommended choice for multiplicative corrections (efficiencies, cross sections, ...).
####    gmN stands for Gamma, and is the recommended choice for the statistical uncertainty on a background coming from the number of events in a control region 
####
####    Usage: python createDataCard.py <n_bins> <RunYear> <preVFP> <obs_mode>
####      obs_mode:
####        unblind  -- use real data observation (reads from data log files)
####        bonly    -- set observation = total expected background  [default blinded]
####        splusb   -- set observation = signal + total expected background (signal filled in by writeT5WGcard.py via placeholder)
import re
import os
# run inside test/myResult (reads/writes its logs/: data logs, counting_exp templates)
os.chdir(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../test/myResult')))
import ROOT
import argparse
from os import system
import sys
import math

limdir = 'logs/'

n_channels = int(sys.argv[1])*2
RunYear = sys.argv[2]
preVFP = sys.argv[3]
obs_mode = sys.argv[4] if len(sys.argv) > 4 else 'bonly'   # NEW: unblind | bonly | splusb

if obs_mode not in ('unblind', 'bonly', 'splusb'):
    print(f"ERROR: obs_mode must be one of: unblind, bonly, splusb  (got '{obs_mode}')")
    sys.exit(1)

n_processes= 6
# 5 bkgrounds, 1 susy process

# rare renamed to reflect actual content: ttgamma + WWgamma + WZgamma + ZZgamma
rare_name = 'ttGam_VVGam'
pro_names = ['SUSY', 'elefakepho', 'jetfakepho', 'qcdfakelep', 'VGamma', rare_name]
# ESF split: lepton(+trigger) SF uncertainty per channel (esf_m: mu-gamma bins, esf_e: e-gamma bins)
# and the photon ID SF uncertainty (esfpho) shared by both channels
# wgnlo: LO -> NLO W+gamma ISR pT factor (V+gamma only, correlated across years and channels)
syst_names = ['jes','jer','esf_m','esf_e','esfpho','e_to_pho_syst','j_to_pho_syst','lumi','isr','wgnlo']   # V+gamma / fake-lepton normalisation: own block below
syst_hist = {'esf_m': 'esf', 'esf_e': 'esf'}   # histogram suffix when it differs from the row name

# Cross section uncertainty on rare process (ttgamma + VVgamma)
XSEC_RARE_UNC = 1.50


# Flat 5% QCD scale and PDF uncertainty on SUSY signal
SIG_QCD_PDF_UNC = 1.05

proc_to_hist = {
    'SUSY':      'SUSY',       # overridden per mass point
    'elefakepho':  'elefakepho',
    'jetfakepho':  'jetfakepho',
    'qcdfakelep':  'qcdfakelep',
    'VGamma':      'VGamma',
    rare_name:     'rare',
}

if RunYear == '2016' and preVFP == '1':
    whichVFP = 'preVFP'
elif RunYear == '2016' and preVFP == '0':
    whichVFP = 'postVFP'
else:
    whichVFP = ''

file_in = ROOT.TFile(f'/uscms_data/d3/tmishra/Output/SignalSystematic_{RunYear}{whichVFP}.root', 'read')
h_rates = {}


# Load background histograms using original ROOT names
bkg_hist_names = ['elefakepho', 'jetfakepho', 'qcdfakelep', 'VGamma', 'rare']
for hname in bkg_hist_names:
    h_rates['h_' + hname + '_norm'] = file_in.Get('h_' + hname + '_norm')
    for sys in syst_names:
        h_rates['h_' + hname + '_syserr_' + sys] = file_in.Get('h_' + hname + '_syserr_' + syst_hist.get(sys, sys))

h_rates['h_elefakepho_controlsample'] = file_in.Get('h_elefakepho_controlsample')
h_rates['h_elefakepho_transferfactor'] = file_in.Get('h_elefakepho_transferfactor')
h_rates['h_jetfakepho_controlsample'] = file_in.Get('h_jetfakepho_controlsample')
h_rates['h_jetfakepho_transferfactor'] = file_in.Get('h_jetfakepho_transferfactor')
h_rates['h_qcdfakelep_controlsample'] = file_in.Get('h_qcdfakelep_controlsample')
h_rates['h_qcdfakelep_transferfactor'] = file_in.Get('h_qcdfakelep_transferfactor')

logfile_out = open( 'logs/t5wg.log','w')

filename = limdir + 'counting_exp_XXX_YYY_{}_{}_{}{}.txt'.format(obs_mode, n_channels, RunYear, whichVFP)
file_out = open(filename, 'w')

#def compute_blinded_observation(ich):
#    """Return total expected background in channel ich."""
#    total = 0.0
#    for p in bkg_hist_names:
#        h = h_rates.get('h_' + p + '_norm')
#        print(h)
#        if h:
#            total += h.GetBinContent(ich)
#    return total


def compute_blinded_observation(ich):
    """Print each histogram bin content in a row and return total."""
    total = 0.0
    row_values = []

    for p in bkg_hist_names:
        h = h_rates.get('h_' + p + '_norm')
        if h:
            val = h.GetBinContent(ich)
        else:
            val = 0.0  # fallback if histogram missing

        row_values.append(val)
        total += val
    #print(" ".join(f"{v:8.3f}" for v in row_values) + f" | Total: {total:8.3f}")

    return total
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

# --- Observation block: three modes ---
if obs_mode == 'unblind':
    # Real data: read from the data log files (MG + EG channels)
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

elif obs_mode == 'bonly':
    # Blinded: observation = total expected background per bin
    for ich in range(1, n_channels + 1):
        file_out.write('{:6f} '.format(compute_blinded_observation(ich)))

elif obs_mode == 'splusb':
    for ich in range(1, n_channels + 1):
        file_out.write('{:>12s} '.format('SPLUSB' + str(ich)))

file_out.write('\n')
# --- end observation block ---

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
            hist_key = proc_to_hist.get(p, p)
            evtrate = h_rates['h_' + hist_key + '_norm'].GetBinContent(ich)
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
    # JER is uncorrelated across years -> year-tagged row name
    # ('isr' is year-tagged by convert_all_cards.py / convert_dict.yaml: isr -> CMS_ISR_13TeV_<year>)
    row_name = syst + '_' + RunYear + whichVFP if syst == 'jer' else syst
    file_out.write('{:15s} {:3s} {:6s}'.format(row_name,'lnN',''))
    for ich in range(1,n_channels+1):
        # bins 1..n_channels/2 are mu-gamma, the rest e-gamma
        if (syst == 'esf_m' and ich > n_channels//2) or (syst == 'esf_e' and ich <= n_channels//2):
            for p in pro_names:
                file_out.write('{:>12s} '.format('-'))
            continue
        for p in pro_names:
            if(p == 'SUSY'):
		# Just a string is saved for SUSY
                if(syst == 'jes'):
                    file_out.write('{:>12s} '.format('JESS'+ str(ich)))
                elif(syst == 'jer'):
                    file_out.write('{:>12s} '.format('JERS'+ str(ich)))
                elif(syst in ('esf_m', 'esf_e')):
                    file_out.write('{:>12s} '.format('ESF'+ str(ich)))
                elif(syst == 'esfpho'):
                    file_out.write('{:>12s} '.format('PHOSF'+ str(ich)))
                elif(syst == 'isr'):
                    file_out.write('{:>12s} '.format('ISR'+ str(ich)))
                elif(syst == 'lumi'):
                    file_out.write('{:12.3f} '.format(1.026))
                else: 
                    file_out.write('{:>12s} '.format('-'))
            else:
                syserr = 0
                hist_key = proc_to_hist.get(p, p)
                norm_val = h_rates['h_' + hist_key + '_norm'].GetBinContent(ich)
                if norm_val == 0:
		    # errors from transfer factors taken for bkgs, when central value is zero
                    if(p == 'elefakepho' and syst == 'e_to_pho_syst'):
                        syserr = h_rates['h_elefakepho_transferfactor'].GetBinError( ich) + 1
                        file_out.write('{:12.3f} '.format(syserr))
                    elif(p == 'jetfakepho' and syst == 'j_to_pho_syst'):
                        syserr = h_rates['h_jetfakepho_transferfactor'].GetBinError( ich) + 1
                        file_out.write('{:12.3f} '.format(syserr))
                    else:
                        file_out.write('{:>12s} '.format('-'))
                elif not h_rates.get('h_' + hist_key + '_syserr_' + syst):
                    # e.g. no photon-SF histogram for the data-driven backgrounds
                    file_out.write('{:>12s} '.format('-'))
                else:
                    syst_val = h_rates['h_' + hist_key + '_syserr_' + syst].GetBinContent(ich)
                    if not math.isfinite(syst_val):
                        print(f"WARNING: NaN/Inf in h_{hist_key}_syserr_{syst} "
                              f"bin {ich} norm={norm_val:.3f} writing '-'")
                        file_out.write('{:>12s} '.format('-'))
                        continue
                
                    syserr = syst_val / norm_val
                    if syserr < 0:
                        file_out.write('{:>12s} '.format('-'))
                    else:
                        syserr = syserr + 1
                        if not math.isfinite(syserr):
                            file_out.write('{:>12s} '.format('-'))
                        else:
                            file_out.write('{:12.3f} '.format(syserr))
    file_out.write('\n')

# V+gamma / fake-lepton normalisation from the dPhi template fits (VGamma/fit). The fits are independent per year,
# channel and photon-pT bin (0-100, 100-200, >200 GeV) -> one nuisance each. Both factors come from the fitted
# fake-lepton fraction f and move in opposite directions (V+gamma = (1-f)T/M, fake = fT/P), so they share the nuisance,
# both as asymmetric lnN kappa_down/kappa_up from the -1 / +1 sigma yields of pred_VGBkg.C / pred_qcdBkg.C
# (*_sfDown/sfUp_pho<k>). theta = +1: V+gamma at its +1 sigma, fake lepton at its -1 sigma.
# A fake-lepton -1 sigma yield <= 0 cannot be written as lnN: symmetric lnN with the +1 sigma error, kappa = 1/(1+r_up).
sf_name = 'CMS_SUS2411_normSF_VGamma_{}_pho{}_' + RunYear + whichVFP
nb = n_channels // 2
for chfile, ch, off in (('mg', 'mg', 0), ('egamma', 'eg', nb)):
    f_vg = ROOT.TFile(f'/uscms_data/d3/tmishra/Output/signalTree_{chfile}_VGBkg_{RunYear}{whichVFP}.root', 'read')
    f_qcd = ROOT.TFile(f'/uscms_data/d3/tmishra/Output/signalTree_{chfile}_qcd_{RunYear}{whichVFP}.root', 'read')
    vg_nom, qcd_nom = f_vg.Get(f'{ch}_VGamma_norm'), f_qcd.Get(f'{ch}_qcdfakelep_norm')
    for k in range(3):
        vg_up, vg_dn = f_vg.Get(f'{ch}_VGamma_sfUp_pho{k}'), f_vg.Get(f'{ch}_VGamma_sfDown_pho{k}')
        qcd_up, qcd_dn = f_qcd.Get(f'{ch}_qcdfakelep_sfUp_pho{k}'), f_qcd.Get(f'{ch}_qcdfakelep_sfDown_pho{k}')
        cells, used = [], False
        for ich in range(1, n_channels+1):
            b = ich - off
            inch = 1 <= b <= nb
            for p in pro_names:
                cell = '-'
                if inch and p == 'VGamma':
                    v, vu, vd = vg_nom.GetBinContent(b), vg_up.GetBinContent(b), vg_dn.GetBinContent(b)
                    if v > 0 and (vu != v or vd != v):
                        cell = '{:.3f}/{:.3f}'.format(vd/v, vu/v)
                elif inch and p == 'qcdfakelep':
                    q, qu, qd = qcd_nom.GetBinContent(b), qcd_up.GetBinContent(b), qcd_dn.GetBinContent(b)
                    if q > 0 and (qu != q or qd != q):
                        cell = '{:.3f}/{:.3f}'.format(qu/q, qd/q) if qd > 0 else '{:.3f}'.format(q/qu)
                cells.append('{:>12s} '.format(cell))
                used = used or cell != '-'
        if used:
            file_out.write('{:15s} {:3s} {:6s}'.format(sf_name.format(ch, k), 'lnN', '') + ''.join(cells) + '\n')
    f_vg.Close(); f_qcd.Close()

# Produced by pred_VGBkg.C  : h_VGamma_syserr_qcdmuR/muF/pdf
# Produced by pred_rareBkg.C: h_rare_syserr_qcdmuR/muF/pdf
# syserr_qcd* = max(|muR_up/dn - 1|, ...)  (FRACTIONAL, i.e. |kappa-1|; kappa = 1 + syserr)
# syserr_pdf  = max(|pdf_up/dn - 1|)         (FRACTIONAL, i.e. |kappa-1|; kappa = 1 + syserr)
for tag in ('VGamma', 'rare'):
    for suf in ('qcdmuR', 'qcdmuF', 'pdf'):
        h_rates[f'h_{tag}_syserr_{suf}'] = file_in.Get(f'h_{tag}_syserr_{suf}')
h_rates['h_VGamma_syserr_pu']  = file_in.Get('h_VGamma_syserr_pu')
h_rates['h_rare_syserr_pu']    = file_in.Get('h_rare_syserr_pu')
file_out.write('{:15s} {:3s} {:6s}'.format('PU','lnN',''))
for ich in range(1,n_channels+1):
    for p in pro_names:
        if p == 'SUSY':
            file_out.write('{:>12s} '.format('PUS' + str(ich)))
        elif p == 'VGamma':
            hist_key = proc_to_hist.get(p, p)
            norm_val = h_rates['h_' + hist_key + '_norm'].GetBinContent(ich)
            if norm_val > 0:
                pu_val = h_rates['h_VGamma_syserr_pu'].GetBinContent(ich)
                if pu_val and math.isfinite(pu_val) and pu_val > 0:
                    pu_kappa = pu_val / norm_val + 1
                    file_out.write('{:12.3f} '.format(pu_kappa))
                else:
                    file_out.write('{:>12s} '.format('-'))
            else:
                file_out.write('{:>12s} '.format('-'))
        elif p == rare_name:
            norm_val = h_rates['h_rare_norm'].GetBinContent(ich)
            if norm_val > 0:
                pu_val = h_rates['h_rare_syserr_pu'].GetBinContent(ich)
                if pu_val and math.isfinite(pu_val) and pu_val > 0:
                    pu_kappa = pu_val / norm_val + 1
                    file_out.write('{:12.3f} '.format(pu_kappa))
                else:
                    file_out.write('{:>12s} '.format('-'))
            else:
                file_out.write('{:>12s} '.format('-'))
        else:
            file_out.write('{:>12s} '.format('-'))
file_out.write('\n')




# QCDscale_VGamma: covers VGamma only; signal gets '-' here (signal row written per-mass-point by write*.py)
file_out.write('{:30s} {:3s} {:6s}'.format('QCDscale_VGamma', 'lnN', ''))
for ich in range(1, n_channels + 1):
    for p in pro_names:
        if p == 'VGamma':
            h_muF = h_rates.get('h_VGamma_syserr_qcdmuF')
            h_nom = h_rates.get('h_VGamma_norm')
            norm  = h_nom.GetBinContent(ich) if h_nom else 0
            if norm > 0 and h_muF:
                err_muF = max(0., h_muF.GetBinContent(ich))
                kappa   = 1.0 + err_muF            # syserr is fractional |kappa-1| (pred_VGBkg.C)
                if math.isfinite(kappa) and kappa > 1.0:
                    file_out.write('{:12.3f} '.format(kappa))
                else:
                    file_out.write('{:>12s} '.format('-'))
            else:
                file_out.write('{:>12s} '.format('-'))
        else:
            file_out.write('{:>12s} '.format('-'))
file_out.write('\n')

# QCDscale_ttGam_VVGam: covers rare only; signal gets '-' here
file_out.write('{:30s} {:3s} {:6s}'.format('QCDscale_ttGam_VVGam', 'lnN', ''))
for ich in range(1, n_channels + 1):
    for p in pro_names:
        if p == rare_name:
            # NLO sample: both muR and muF contribute; take the larger one.
            h_muR = h_rates.get('h_rare_syserr_qcdmuR')
            h_muF = h_rates.get('h_rare_syserr_qcdmuF')
            h_nom = h_rates.get('h_rare_norm')
            norm  = h_nom.GetBinContent(ich) if h_nom else 0
            if norm > 0 and h_muR and h_muF:
                err_muR = max(0., h_muR.GetBinContent(ich))
                err_muF = max(0., h_muF.GetBinContent(ich))
                kappa   = 1.0 + max(err_muR, err_muF)   # syserr is fractional |kappa-1| (pred_rareBkg.C)
                if math.isfinite(kappa) and kappa > 1.0:
                    file_out.write('{:12.3f} '.format(kappa))
                else:
                    file_out.write('{:>12s} '.format('-'))
            else:
                file_out.write('{:>12s} '.format('-'))
        else:
            file_out.write('{:>12s} '.format('-'))
file_out.write('\n')

# pdf_alphas: covers VGamma and rare; signal gets '-' here (signal row written per-mass-point by write*.py)
file_out.write('{:30s} {:3s} {:6s}'.format('pdf_alphas', 'lnN', ''))
for ich in range(1, n_channels + 1):
    for p in pro_names:
        if p == 'SUSY':
            file_out.write('{:>12s} '.format('PDFS' + str(ich)))
        elif p in ('VGamma', rare_name):
            hist_tag = 'VGamma' if p == 'VGamma' else 'rare'
            h_pdf = h_rates.get(f'h_{hist_tag}_syserr_pdf')
            h_nom = h_rates.get(f'h_{hist_tag}_norm')
            norm  = h_nom.GetBinContent(ich) if h_nom else 0
            if norm > 0 and h_pdf:
                kappa = 1.0 + max(0., h_pdf.GetBinContent(ich))   # syserr is fractional |kappa-1|
                if math.isfinite(kappa) and kappa > 1.0:
                    file_out.write('{:12.3f} '.format(kappa))
                else:
                    file_out.write('{:>12s} '.format('-'))
            else:
                file_out.write('{:>12s} '.format('-'))
        else:
            file_out.write('{:>12s} '.format('-'))
file_out.write('\n')


if RunYear in ['2016', '2017']:
    prefire_name = f'CMS_L1prefiring_13TeV_{RunYear}{whichVFP}'
    file_out.write('{:25s} {:3s} {:6s}'.format(prefire_name, 'lnN', ''))
    for ich in range(1, n_channels + 1):
        for p in pro_names:
            if p in ['SUSY', 'VGamma', rare_name]:
                file_out.write('{:12.3f} '.format(1.030))
            else:
                file_out.write('{:>12s} '.format('-'))
    file_out.write('\n')

file_out.write('{:15s} {:3s} {:6s}'.format('xsec_ttGam_VVGam', 'lnN', ''))
for ich in range(1, n_channels + 1):
    for p in pro_names:
        if p == rare_name:
            file_out.write('{:12.3f} '.format(XSEC_RARE_UNC))
        else:
            file_out.write('{:>12s} '.format('-'))
file_out.write('\n')

# V+gamma non-closure for photon pT > 200 GeV, measured in the CR (MET < 70 GeV) after the photon-pT binned SFs:
# one lnN per channel and HT range (HT < 100, HT > 100), on VGamma only, correlated across years.
# Bins: 1-18 mu-gamma, 19-36 e-gamma; within a channel 10-18 have photon pT > 200, HT order <100, 100-400, >400.
nonclosure_file = 'VGamma_nonclosure.txt'
if os.path.exists(nonclosure_file):
    nonclosure = {}
    with open(nonclosure_file) as fnc:
        for line in fnc:
            if line.strip() and not line.startswith('#'):
                ch, htrange, kappa = line.split()[:3]
                nonclosure[(ch, htrange)] = float(kappa)
    for (ch, htrange), kappa in sorted(nonclosure.items()):
        file_out.write('{:40s} {:3s} {:6s}'.format('CMS_SUS2411_VG_nonclosure_{}_HT{}'.format(ch, htrange), 'lnN', ''))
        for ich in range(1, n_channels + 1):
            ich_ch = 'mug' if ich <= n_channels // 2 else 'eg'
            local = (ich - 1) % (n_channels // 2) + 1
            ich_ht = 'lt100' if (local - 10) % 3 == 0 else 'gt100'
            for p in pro_names:
                if p == 'VGamma' and local >= 10 and ich_ch == ch and ich_ht == htrange:
                    file_out.write('{:12.3f} '.format(kappa))
                else:
                    file_out.write('{:>12s} '.format('-'))
        file_out.write('\n')

for ich in range(1,n_channels+1):
    # MC stat nuisances carry the year: MC samples of different years are independent
    file_out.write('{:15s} {:3s} {:6s}'.format('SUSY_stat'+str(ich)+'_'+RunYear+whichVFP,'lnN',''))
    for k in range(1,n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('STSC'+ str(ich),'-', '-', '-', '-', '-'))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')


for ich in range(1,n_channels+1):
    nevt = h_rates['h_elefakepho_controlsample'].GetBinContent(ich)
    fakerate = h_rates['h_elefakepho_transferfactor'].GetBinContent(ich)
    file_out.write('{:22s} {:3s} {:6d}'.format('e_to_pho_'  + RunYear + whichVFP + '_stat'+str(ich),'gmN',int(nevt)))
    for k in range(1,n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:12.8g} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', fakerate, '-', '-', '-', '-'))
            logfile_out.write('{:15} {:2d} {:10.5f}\n'.format('elefakepho', ich, nevt*fakerate))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')
for ich in range(1,n_channels+1):
    nevt = h_rates['h_jetfakepho_controlsample'].GetBinContent(ich)
    fakerate = h_rates['h_jetfakepho_transferfactor'].GetBinContent(ich)
    file_out.write('{:22s} {:3s} {:6d}'.format('j_to_pho_'  + RunYear + whichVFP + '_stat'+str(ich),'gmN',int(nevt)))
    for k in range(1,n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:>12s} {:12.8g} {:>12s} {:>12s} {:>12s} '.format('-', '-', fakerate, '-', '-', '-'))
            logfile_out.write('{:15} {:2d} {:10.5f}\n'.format('jetfakepho', ich, nevt*fakerate))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')
for ich in range(1,n_channels+1):
    nevt = h_rates['h_qcdfakelep_controlsample'].GetBinContent(ich)
    fakerate = h_rates['h_qcdfakelep_transferfactor'].GetBinContent(ich)
    file_out.write('{:22s} {:3s} {:6d}'.format('j_to_lep_'  + RunYear + whichVFP + '_stat'+str(ich),'gmN',int(nevt)))
    for k in range(1,n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:>12s} {:>12s} {:12.8g} {:>12s} {:>12s} '.format('-', '-', '-', fakerate, '-', '-'))
            logfile_out.write('{:15} {:2d} {:10.5f}\n'.format('jetfakelep', ich, nevt*fakerate))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')
# stat error for rare bkgs, factorised prediction y(bin, era) = N(era) x f(bin) (factorise_rareBkg.py):
#  - bin fraction f from all four eras: one nuisance per bin, no year tag (same f in every era)
#  - era yield N: one nuisance per era and channel on all bins of that channel (bins 1..NBIN mg, NBIN+1.. eg)
for ich in range(1,n_channels+1):
    file_out.write('{:15s} {:3s} {:6s}'.format(rare_name + '_stat'+str(ich),'lnN',''))
    if(h_rates['h_rare_norm'].GetBinContent(ich) > 0):
        staterror= 1.0 + h_rates['h_rare_norm'].GetBinError(ich)/h_rates['h_rare_norm'].GetBinContent(ich)
    else:
        sys.exit(f'rare yield in bin {ich} is {h_rates["h_rare_norm"].GetBinContent(ich)}: run factorise_rareBkg.py')
    for k in range(1,n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:12.3f} '.format('-', '-', '-', '-', '-',staterror))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')
for chfile, ch, bins in (('mg', 'mg', range(1, n_channels//2+1)), ('egamma', 'eg', range(n_channels//2+1, n_channels+1))):
    f_rare = ROOT.TFile(f'/uscms_data/d3/tmishra/Output/signalTree_{chfile}_rareBkg_{RunYear}{whichVFP}.root', 'read')
    staterror = 1.0 + f_rare.Get(f'{ch}_rare_statN').GetBinContent(1)
    f_rare.Close()
    file_out.write('{:15s} {:3s} {:6s}'.format(rare_name + '_statN_'+ch+'_'+RunYear+whichVFP,'lnN',''))
    for k in range(1,n_channels+1):
        file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-'))
        file_out.write('{:12.3f} '.format(staterror) if k in bins else '{:>12s} '.format('-'))
    file_out.write('\n')

for ich in range(1,n_channels+1):
    file_out.write('{:15s} {:3s} {:6s}'.format('VG_stat'+str(ich)+'_'+RunYear+whichVFP,'lnN',''))
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
