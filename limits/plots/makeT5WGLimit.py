import os
import numpy as np
import re
import ROOT
import argparse

import sys
from pathlib import Path
cmsstyle_path = Path("/uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/cmsstyle/src")
if str(cmsstyle_path) not in sys.path:
    sys.path.insert(0, str(cmsstyle_path))

import cmsstyle
ROOT.gROOT.SetBatch(ROOT.kTRUE)
cmsstyle.setCMSStyle()

parser = argparse.ArgumentParser(description='Generate limit plots for T5WG.')
parser.add_argument('year', choices=["2016preVFP", "2016postVFP", "2017", "2018", "Run2"],
                    help='Year to process')
parser.add_argument('--MET1', dest='METbin1', type=int, help='First MET bin')
parser.add_argument('--MET2', dest='METbin2', type=int, help='Second MET bin')
args = parser.parse_args()

year     = args.year

# SUS-17-012 reference limits (observed only)
limdirdata_SUS17012 = '/uscms_data/d3/mengleis/CMSSW_Work/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/DataCard/T5WG_Defense/Limit/'

# Current Run-2 analysis
limdirdata = f'/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/limits/limits_{year}_unblind/'   # observed

file_in    = ROOT.TFile('/uscms/home/tmishra/nobackup/signal_trees/signalTree_T5WG_678.root', 'read')
h_SUSYmass = file_in.Get('SUSYMass')
ROOT.gStyle.SetOptStat(0)

gluino_masses    = np.array([])
neutralino_masses = np.array([])
for xbin in range(1, h_SUSYmass.GetXaxis().GetNbins() + 1):
    if int(h_SUSYmass.GetXaxis().GetBinCenter(xbin)) >= 1100:
        gluino_masses = np.append(gluino_masses, int(h_SUSYmass.GetXaxis().GetBinCenter(xbin)))
for ybin in range(1, h_SUSYmass.GetYaxis().GetNbins() + 1):
    neutralino_masses = np.append(neutralino_masses, int(h_SUSYmass.GetYaxis().GetBinCenter(ybin)))

# FIX (bug 2): the grid masses are used as bin LOW EDGES (values are filled with
# FindBin(mass)), so the highest real mass needs an edge above it. Without this
# extra edge, FindBin(top mass) returns the overflow bin and the last row/column
# of the grid is silently dropped. The last bin on each axis is a dummy
# (never filled) whose width repeats the last real grid step.
gluino_edges     = np.append(gluino_masses,
                             gluino_masses[-1] + (gluino_masses[-1] - gluino_masses[-2]))
neutralino_edges = np.append(neutralino_masses,
                             neutralino_masses[-1] + (neutralino_masses[-1] - neutralino_masses[-2]))
nbinsx = len(gluino_edges) - 1
nbinsy = len(neutralino_edges) - 1

hExpLimUp          = ROOT.TH2D('up',          'up',          nbinsx, gluino_edges, nbinsy, neutralino_edges)
hExpLimMid         = ROOT.TH2D('mid',         'mid',         nbinsx, gluino_edges, nbinsy, neutralino_edges)
hExpLimDown        = ROOT.TH2D('down',        'down',        nbinsx, gluino_edges, nbinsy, neutralino_edges)
hExp2Up            = ROOT.TH2D('up2',         'up2',         nbinsx, gluino_edges, nbinsy, neutralino_edges)
hExp2Down          = ROOT.TH2D('down2',       'down2',       nbinsx, gluino_edges, nbinsy, neutralino_edges)
hObsLim            = ROOT.TH2D('obs',         'obs',         nbinsx, gluino_edges, nbinsy, neutralino_edges)
hObsLim_SUS17012   = ROOT.TH2D('obs_SUS17012','obs_SUS17012',nbinsx, gluino_edges, nbinsy, neutralino_edges)

for i in gluino_masses:
    for j in neutralino_masses:
        if j >= i:
            continue

        # ---- expected limits (current analysis) ----------------------------
        exp_path = limdirdata + f'limit_{int(i)}_{int(j)}.txt'
        if os.path.exists(exp_path):
            with open(exp_path) as fin:
                for line in fin:
                    if re.search('16.0%', line):
                        v = float(line.split('<')[1].strip())
                        if v > 0:
                            hExpLimDown.SetBinContent(hExpLimDown.GetXaxis().FindBin(i),
                                                      hExpLimDown.GetYaxis().FindBin(j), v)
                    elif re.search('50.0%', line):
                        v = float(line.split('<')[1].strip())
                        if v > 0:
                            hExpLimMid.SetBinContent(hExpLimMid.GetXaxis().FindBin(i),
                                                     hExpLimMid.GetYaxis().FindBin(j), v)
                    elif re.search('84.0%', line):
                        v = float(line.split('<')[1].strip())
                        if v > 0:
                            hExpLimUp.SetBinContent(hExpLimUp.GetXaxis().FindBin(i),
                                                    hExpLimUp.GetYaxis().FindBin(j), v)
                    elif re.search('2.5%', line):
                        v = float(line.split('<')[1].strip())
                        if v > 0:
                            hExp2Down.SetBinContent(hExp2Down.GetXaxis().FindBin(i),
                                                    hExp2Down.GetYaxis().FindBin(j), v)
                    elif re.search('97.5%', line):
                        v = float(line.split('<')[1].strip())
                        if v > 0:
                            hExp2Up.SetBinContent(hExp2Up.GetXaxis().FindBin(i),
                                                  hExp2Up.GetYaxis().FindBin(j), v)

        # ---- observed limit (current analysis) -----------------------------
        obs_path = limdirdata + f'limit_{int(i)}_{int(j)}.txt'
        if os.path.exists(obs_path):
            with open(obs_path) as findata:
                for line in findata:
                    if re.search('Observed', line):
                        v = float(line.split('<')[1].strip())
                        if v > 0:
                            hObsLim.SetBinContent(hObsLim.GetXaxis().FindBin(i),
                                                  hObsLim.GetYaxis().FindBin(j), v)

        # ---- SUS-17-012 observed limit -------------------------------------
        sus_path = limdirdata_SUS17012 + f'limit_{int(i)}_{int(j)}.txt'
        if os.path.exists(sus_path):
            with open(sus_path) as fsus:
                for line in fsus:
                    if re.search('Observed', line):
                        v = float(line.split('<')[1].strip())
                        if v > 0:
                            hObsLim_SUS17012.SetBinContent(hObsLim_SUS17012.GetXaxis().FindBin(i),
                                                           hObsLim_SUS17012.GetYaxis().FindBin(j), v)

def read_xsec_table(filename):
    """Returns dict {mass_float -> (xsec_pb, unc_percent)}."""
    table = {}
    with open(filename) as f:
        for line in f:
            l = line.split()
            if len(l) < 3:
                continue
            table[float(l[0])] = (float(l[1]) * 1000 * 0.5, float(l[2]))
    return table


def get_xsec_entry(table, m, tol=2.5):
    """Exact match first (the common case); falls back to the nearest table
    key within `tol` GeV. The table is tabulated on a dense but fixed 5 GeV
    grid, while the signal mass grid occasionally lands slightly off of it
    (e.g. 1987 GeV vs. table entries at 1985/1990) - a real grid-alignment
    mismatch, not a genuine coverage gap. tol=2.5 (half the table's 5 GeV
    step) guarantees each mass maps to a unique nearest entry."""
    entry = table.get(m)
    if entry is not None:
        return entry
    nearest_key = min(table, key=lambda k: abs(k - m), default=None)
    if nearest_key is not None and abs(nearest_key - m) <= tol:
        return table[nearest_key]
    return None

xsec_table_t5wg = read_xsec_table('ULCrossSectionT5WG.txt')

# --- DEBUG: which gluino masses have no matching entry (exact or nearby)
# in the xsec table? ---
_missing = sorted(m for m in gluino_masses if get_xsec_entry(xsec_table_t5wg, float(m)) is None)
if _missing:
    print(f"[xsec table] {len(_missing)} gluino mass(es) with NO entry (even "
          f"within tolerance) in ULCrossSectionT5WG.txt: {_missing}")
else:
    print("[xsec table] every gluino mass in the grid has a matching entry")

def ScaleLimit(h, name, scaleByXSec, sigma):
    h_new = h.Clone(name)
    for i in range(1, h_new.GetNbinsX() + 1):
        for j in range(1, h_new.GetNbinsY() + 1):
            # gluino_masses/neutralino_masses were fed to TH2D() as bin EDGES
            # (not centers), so GetBinCenter(i) is the midpoint BETWEEN two
            # real grid masses - it essentially never matches a table key.
            # GetBinLowEdge(i) is the actual real grid value that bin i was
            # built from, and is what belongs in the xsec table lookup.
            m = h_new.GetXaxis().GetBinLowEdge(i)
            r = h_new.GetBinContent(i, j)
            entry = get_xsec_entry(xsec_table_t5wg, m)
            if entry is None:
                h_new.SetBinContent(i, j, -1.0)
                continue
            xSec, unc = entry
            if scaleByXSec:
                r *= xSec
            r /= 1.0 + sigma * unc / 100.0
            h_new.SetBinContent(i, j, r)
    return h_new


def RegularPlot(h):
    for i in range(1, h.GetNbinsX() + 1):
        for j in range(1, h.GetNbinsY() + 1):
            if i == h.GetNbinsX() or j == h.GetNbinsY():
                continue   # dummy edge bins (see bug-2 fix): never fill
            if h.GetYaxis().GetBinCenter(j) > h.GetXaxis().GetBinCenter(i) - 50:
                continue
            if h.GetBinContent(i, j) > 0:
                continue
            neighbors = []
            # search left
            for lbin in range(1, i):
                if h.GetBinContent(i - lbin, j) > 0:
                    neighbors.append(h.GetBinContent(i - lbin, j))
                    break
            # search right
            for rbin in range(1, h.GetNbinsX() - i + 1):
                if h.GetBinContent(i + rbin, j) > 0:
                    neighbors.append(h.GetBinContent(i + rbin, j))
                    break
            # search up
            for ubin in range(1, h.GetNbinsY() - j + 1):
                if h.GetBinContent(i, j + ubin) > 0:
                    neighbors.append(h.GetBinContent(i, j + ubin))
                    break
            # search down
            for dbin in range(1, j):
                if h.GetBinContent(i, j - dbin) > 0:
                    neighbors.append(h.GetBinContent(i, j - dbin))
                    break
            if neighbors:
                h.SetBinContent(i, j, sum(neighbors) / len(neighbors))

def MergePlot(h, hname='obs2'):
    h_new = ROOT.TProfile2D(hname, hname, 20, 1100, 2100, 84, 0, 2100)
    for i in range(1, h.GetNbinsX() + 1):
        for j in range(1, h.GetNbinsY() + 1):
            if h.GetYaxis().GetBinCenter(j) <= h.GetXaxis().GetBinCenter(i):
                if h.GetBinContent(i, j) > 0:
                    h_new.Fill(h.GetXaxis().GetBinCenter(i),
                               h.GetYaxis().GetBinCenter(j),
                               h.GetBinContent(i, j))
    return h_new

def MakeLimitGraph(h, name):
    """r = 1 crossing along gluino mass, for each neutralino mass.

    FIX (bug 1): values were filled at the grid mass = bin LOW EDGE, so the
    coordinates must be low edges too (GetBinCenter shifted x and y by half a
    bin, by an amount that varies on this irregular grid).
    FIX (bug 3): r falls ~exponentially with mass (sigma does), so interpolate
    log(r) linearly instead of r.
    """
    g = ROOT.TGraph()
    if name:
        g.SetName(name)
    n = 0
    max_x = max_y = None
    for i in range(1, h.GetNbinsY() + 1):
        y = h.GetYaxis().GetBinLowEdge(i)
        # stop at GetNbinsX()-1 so bin j+1 is never overflow
        for j in range(1, h.GetNbinsX()):
            xLow  = h.GetXaxis().GetBinLowEdge(j)
            xHigh = h.GetXaxis().GetBinLowEdge(j + 1)
            zLow  = h.GetBinContent(j, i)
            zHigh = h.GetBinContent(j + 1, i)
            if 0.01 < zLow < 1.0 and zHigh > 1.0:
                f = np.log(1.0 / zLow) / np.log(zHigh / zLow)
                x = xLow + (xHigh - xLow) * f
                g.SetPoint(n, x, y)
                if max_x is None or x > max_x:
                    max_x = x
                if max_y is None or y > max_y:
                    max_y = y
                n += 1
                break
    print(f"Graph '{name}': Maximum value (X, Y) = {max_x}, {max_y}")
    return g

h_xSecLimit        = ScaleLimit(hObsLim,          'obs_xsec',      True,  0.0)
h_observedLimitLow = ScaleLimit(hObsLim,          'obs_low',       False, -1.0)
h_observedLimitHigh= ScaleLimit(hObsLim,          'obs_high',      False,  1.0)

def MaskDegenerate(h, margin=0.0):
    """Zero out any bin whose center violates the physical LSP-lighter-
    than-gluino requirement (M_chargino < M_gluino), regardless of how it
    got populated upstream (bin-edge/bin-center rounding on an irregular
    mass grid, or RegularPlot's neighbor-fill creeping past its 50 GeV
    diagonal buffer). margin lets you widen the excluded band if points
    are still leaking through right at the edge."""
    for i in range(1, h.GetNbinsX() + 1):
        for j in range(1, h.GetNbinsY() + 1):
            xc = h.GetXaxis().GetBinCenter(i)
            yc = h.GetYaxis().GetBinCenter(j)
            if yc >= xc - margin:
                h.SetBinContent(i, j, 0)
 
MaskDegenerate(h_xSecLimit)
 
RegularPlot(h_xSecLimit)
 
MaskDegenerate(h_xSecLimit)
 

h_xSecMerge = MergePlot(h_xSecLimit, 'obs2')

cmsstyle.SetEnergy(13)
cmsstyle.SetLumi(138, "fb", "Run 2", 0)
cmsstyle.SetExtraText("")

ROOT.gStyle.SetNdivisions(505, 'X')

xmin, xmax = h_xSecLimit.GetXaxis().GetXmin(), h_xSecLimit.GetXaxis().GetXmax()
ymin, ymax = h_xSecLimit.GetYaxis().GetXmin(), 2400.0

c_limit = cmsstyle.cmsCanvas('c_limit', xmin, xmax, ymin, ymax,
                             'm_{#tilde{g}} (GeV)', 'm_{#tilde{#chi}} (GeV)',
                             square=True, with_z_axis=True, iPos=0,
                             yTitOffset=1.7, extraSpace=0.02)
c_limit.SetRightMargin(0.20)
c_limit.SetLogz()
cmsstyle.SetAlternative2DColor(h_xSecLimit)  

h_xSecLimit.GetZaxis().SetTitle('95% CL upper limit on #sigma (fb)')
h_xSecLimit.GetZaxis().SetTitleOffset(1.2)
h_xSecLimit.GetZaxis().SetRangeUser(0.5, 40)
h_xSecLimit.SetContour(255)
h_xSecLimit.SetTitle('')
cmsstyle.cmsObjectDraw(h_xSecLimit, 'COLZ SAME')

# Current analysis: expected
g_expMid   = MakeLimitGraph(hExpLimMid,  'hExpLimMid')
g_expDown  = MakeLimitGraph(hExpLimDown, 'hExpLimDown')
g_expUp    = MakeLimitGraph(hExpLimUp,   'hExpLimUp')
g_exp2Down = MakeLimitGraph(hExp2Down,   'hExp2Down')
g_exp2Up   = MakeLimitGraph(hExp2Up,     'hExp2Up')
# Current analysis: observed + theory unc
g_obsMid  = MakeLimitGraph(hObsLim,            'hObsLim')
g_obsLow  = MakeLimitGraph(h_observedLimitLow,  'h_observedLimitLow')
g_obsHigh = MakeLimitGraph(h_observedLimitHigh, 'h_observedLimitHigh')
# SUS-17-012: observed + theory unc
g_susMid  = MakeLimitGraph(hObsLim_SUS17012, 'hObsLim_SUS17012')

# SUS-17-012 observed (blue)
g_susMid.SetLineColor(ROOT.kBlue+1);  g_susMid.SetLineWidth(2)
g_susMid.Draw('C same')   # SUS-17-012 observed (ARC Fig. 5 comment)

# Current Run-2 observed (black)
g_obsMid.SetLineColor(ROOT.kBlack);  g_obsMid.SetLineWidth(2)
g_obsLow.SetLineColor(ROOT.kBlack);  g_obsLow.SetLineWidth(2);  g_obsLow.SetLineStyle(2)
g_obsHigh.SetLineColor(ROOT.kBlack); g_obsHigh.SetLineWidth(2); g_obsHigh.SetLineStyle(2)
g_obsMid.Draw('C same')
g_obsLow.Draw('C same')
g_obsHigh.Draw('C same')

# Current Run-2 expected (red)
g_expMid.SetLineColor(ROOT.kRed+1);  g_expMid.SetLineWidth(2)
g_expDown.SetLineColor(ROOT.kRed+1); g_expDown.SetLineWidth(2); g_expDown.SetLineStyle(2)
g_expUp.SetLineColor(ROOT.kRed+1);   g_expUp.SetLineWidth(2);   g_expUp.SetLineStyle(2)
g_expMid.Draw('C same')
g_expDown.Draw('C same')
g_expUp.Draw('C same')

# Current Run-2 expected +-2 sigma (dotted red)
#g_exp2Down.SetLineColor(ROOT.kRed+1); g_exp2Down.SetLineWidth(1); g_exp2Down.SetLineStyle(3)
#g_exp2Up.SetLineColor(ROOT.kRed+1);   g_exp2Up.SetLineWidth(1);   g_exp2Up.SetLineStyle(3)
#g_exp2Down.Draw('C same')
#g_exp2Up.Draw('C same')

le_limit = cmsstyle.cmsLeg(0.25, 0.635, 0.53, 0.83, textSize=0.025, textFont=42)

#le_limit.AddEntry(g_susLow,  'SUS-17-012 Obs. #pm1#sigma_{th}', 'l')
le_limit.AddEntry(g_obsMid,   'Obs. limit',               'l')
le_limit.AddEntry(g_obsLow,   'Obs. #pm1#sigma_{th}',     'l')
le_limit.AddEntry(g_expMid,   'Exp. limit',               'l')
le_limit.AddEntry(g_expDown,  'Exp. #pm1#sigma_{ex}',     'l')
le_limit.AddEntry(g_susMid,   'SUS-17-012 obs.',    'l')
#le_limit.AddEntry(g_exp2Down, 'Exp. #pm2#sigma_{ex}',     'l')


la_limit = ROOT.TLatex()
la_limit.SetNDC()

la_limit.SetTextFont(42); la_limit.SetTextSize(0.038)
la_limit.DrawLatex(0.26, 0.855,
    'pp #rightarrow #tilde{g}#tilde{g}, '
    '#tilde{g} #rightarrow qq#tilde{#chi}_{1}^{0/#pm}, '
    '#tilde{#chi}_{1}^{0/#pm} #rightarrow #gamma/W^{#pm}#tilde{G}')


ROOT.gPad.RedrawAxis('G')
cmsstyle.SaveCanvas(c_limit, f'/uscms_data/d3/tmishra/Output/T5WG_Run2_2016-17-2018.pdf')

file_out = ROOT.TFile(f'/uscms_data/d3/tmishra/Output/T5WG_{year}.root', 'RECREATE')
file_out.cd()
for obj in [h_xSecLimit, h_xSecMerge,
            g_expMid, g_expDown, g_expUp, g_exp2Down, g_exp2Up,
            g_obsMid, g_obsLow, g_obsHigh,
            g_susMid]:
    obj.Write()
file_out.Close()
