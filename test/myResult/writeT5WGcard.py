import re
import os
import ROOT
from os import system
import sys

#br_neu = 0.5
#br_cha = 1-br_neu
#br_susy = 2*br_neu*br_cha
#br_susy = br_susy/0.5
#print br_susy

# --- Diagnostics for large signal systematics ---
JEC_WARN = 1.99   # ~100%
JER_WARN = 1.99

JEC_100_cases = {}
JER_100_cases = {}


def classify_instability(n_nom):
    if n_nom < 1e-4:
        return "ultra-low yield"
    elif n_nom < 1e-2:
        return "low yield"
    else:
        return "shape-driven"



br_susy = 1
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
    output_file_path = f'/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_{RunYear}{whichVFP}/'
else:
    output_file_path = f'/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_{RunYear}/'

susy_in = ROOT.TFile(f'/uscms/home/tmishra/nobackup/signal_trees/signalTree_T5WG_{RunYear}{whichVFP}.root', 'read')
syst_names = ['jes','jer','esf','scale','eleshape','jetshape','qcdshape','xs','lumi']

# Global JES/JER/ESF percentage ranges for T5WG signal

jes_dc_min, jes_dc_max = 1e9, 0.0
jer_dc_min, jer_dc_max = 1e9, 0.0
esf_dc_min, esf_dc_max = 1e9, 0.0
pu_dc_min,  pu_dc_max  = 1e9, 0.0
xs_dc_min,  xs_dc_max  = 1e9, 0.0

MIN_YIELD = 0.5  # events

h_SUSYmass = susy_in.Get('SUSYMass')
h_rates = {}
for i in range(1, n_channels + 1):
    h_rates['h_chan' + str(i) + '_rate_nom'] = susy_in.Get('h_chan' + str(i) + '_rate_nom')
    for j in syst_names:
        h_rates['h_chan' + str(i) + '_syserr_' + j] = susy_in.Get('h_chan' + str(i) + '_syserr_' + j)
h_rates['t5wg_h_syserr_PU'] = susy_in.Get('t5wg_h_syserr_PU')

file_template = open('logs/counting_exp_XXX_YYY_{}_{}{}.txt'.format(n_channels,RunYear,whichVFP), 'r')
lines = [line for line in file_template.readlines()]

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

            MIN_YIELD = 0.5  # events

            if n_nom > MIN_YIELD:
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
                    e_xs  = h_rates['h_chan' + str(k) + '_syserr_xs'].GetBinContent(i, j)

                    unc_PU = h_rates['t5wg_h_syserr_PU'].GetBinContent(i, j) + 1.0
                    unc_stat = 1.0
                    unc_jes = 1.0 
                    unc_jer = 1.0
                    unc_esf = 1.0
                    unc_xs = 1.0

                    key = (RunYear, gluino_mass, neutralino_mass, k)
                    if n_nom > MIN_YIELD:
                        unc_stat = 1.0 +  e_nom/n_nom 
                        unc_jes = 1.0 +  e_jes/n_nom
                        unc_jer = 1.0 +  e_jer/n_nom 
                        unc_esf = 1.0 +  e_esf/n_nom
                        unc_xs  = 1.0 +  e_xs/n_nom

                        # --- Track large JES ---
                        if unc_jes >= JEC_WARN:
                            JEC_100_cases[key] = (n_nom, unc_jes, classify_instability(n_nom))

                        # --- Track large JER ---
                        if unc_jer >= JER_WARN:
                            JER_100_cases[key] = (n_nom, unc_jer, classify_instability(n_nom))

                    else:
                        unc_stat = 2.0
                        unc_jes = 1.0 + avg_jes 
                        unc_jer = 1.0 + avg_jer 
                        unc_esf = 1.0 + avg_esf 
                        unc_xs  = 1.0 + avg_xs

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
                        l = l.replace('XSS' + str(k), str(round(unc_xs , 3)))
                    elif re.search('PUS', l):
                        l = l.replace('PUS' + ' ', str(round(unc_PU , 3)))
                        if(n_nom > MIN_YIELD):
                          if(unc_PU < low_p):
                              low_p = unc_PU
                          if(unc_PU > high_p):
                              high_p = unc_PU
                    jes_dc_min = min(jes_dc_min, unc_jes)
                    jes_dc_max = max(jes_dc_max, unc_jes)

                    jer_dc_min = min(jer_dc_min, unc_jer)
                    jer_dc_max = max(jer_dc_max, unc_jer)

                    esf_dc_min = min(esf_dc_min, unc_esf)
                    esf_dc_max = max(esf_dc_max, unc_esf)

                    pu_dc_min = min(pu_dc_min, unc_PU)
                    pu_dc_max = max(pu_dc_max, unc_PU)

                    xs_dc_min  = min(xs_dc_min,  unc_xs)
                    xs_dc_max  = max(xs_dc_max,  unc_xs)
            file_out.write(l)
        file_out.close()
print("Total count:", count)
print("\n=== Large T5WG signal JES/JER diagnostics ===")

print(f"\nJES ≥ 100% cases: {len(JEC_100_cases)}")
for (yr, mG, mN, ch), (yld, unc, cause) in list(JEC_100_cases.items())[:10]:
    print(f"  Year {yr}, mG={mG}, mN={mN}, ch={ch}, "
          f"yield={yld:.3e}, JES={(unc-1)*100:.1f}%, cause={cause}")

print(f"\nJER ≥ 100% cases: {len(JER_100_cases)}")
for (yr, mG, mN, ch), (yld, unc, cause) in list(JER_100_cases.items())[:10]:
    print(f"  Year {yr}, mG={mG}, mN={mN}, ch={ch}, "
          f"yield={yld:.3e}, JER={(unc-1)*100:.1f}%, cause={cause}")

print("\nT5WG signal systematic uncertainty ranges (as written to datacards):")
print(f"  JES : {(jes_dc_min-1)*100:.1f} – {(jes_dc_max-1)*100:.1f} %")
print(f"  JER : {(jer_dc_min-1)*100:.1f} – {(jer_dc_max-1)*100:.1f} %")
print(f"  ESF : {(esf_dc_min-1)*100:.1f} – {(esf_dc_max-1)*100:.1f} %")
print(f"  PU  : {(pu_dc_min -1)*100:.1f} – {(pu_dc_max -1)*100:.1f} %")
print(f"  XS  : {(xs_dc_min -1)*100:.1f} – {(xs_dc_max -1)*100:.1f} %")
