# for reference 
#file_in = ROOT.TFile('/uscms_data/d3/mengleis/CMSSW_Work/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/tools/signalTree_T5WG.root', 'read')

import os
import numpy as np
import re
import ROOT
import argparse
import sys
ROOT.gROOT.SetBatch(ROOT.kTRUE)

parser = argparse.ArgumentParser(description='Generate limit plots for a given year.')
parser.add_argument('year', choices=["2016preVFP", "2016postVFP", "2017", "2018", "Run2"], help='Year to process')
parser.add_argument('--MET1', dest='METbin1', type=int, help='First MET bin')
parser.add_argument('--MET2', dest='METbin2', type=int, help='Second MET bin')
args = parser.parse_args()

year = args.year
met1 = args.METbin1
met2 = args.METbin2

limdirdata = '/uscms_data/d3/mengleis/CMSSW_Work/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/DataCard/T5WG_Defense/Limit/'
limdir = f'/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/limits/limits_{year}/'

gluino_masses = np.array([])
neutralino_masses = np.array([])

file_in = ROOT.TFile(f'/uscms/home/tmishra/nobackup/signal_trees/signalTree_T5WG_{year}.root', 'read')  
h_SUSYmass = file_in.Get('SUSYMass')

ROOT.gStyle.SetOptStat(0)

for xbin in range(1, h_SUSYmass.GetXaxis().GetNbins() + 1):
    if(int(h_SUSYmass.GetXaxis().GetBinCenter(xbin)) >= 1100):
        gluino_masses = np.append(gluino_masses, int(h_SUSYmass.GetXaxis().GetBinCenter(xbin)))
for ybin in range(1, h_SUSYmass.GetYaxis().GetNbins() + 1):
    neutralino_masses = np.append(neutralino_masses, int(h_SUSYmass.GetYaxis().GetBinCenter(ybin)) )

#print(gluino_masses)
#print(neutralino_masses)
nbinsx = len(gluino_masses) - 1
nbinsy = len(neutralino_masses) - 1

hExpLimUp = ROOT.TH2D('up', 'up', nbinsx, gluino_masses, nbinsy, neutralino_masses)
hExpLimMid = ROOT.TH2D('mid', 'mid', nbinsx, gluino_masses, nbinsy, neutralino_masses)
hExpLimDown = ROOT.TH2D('down', 'down', nbinsx, gluino_masses, nbinsy, neutralino_masses)
hObsLim = ROOT.TH2D('obs', 'obs', nbinsx, gluino_masses, nbinsy, neutralino_masses)

for i in gluino_masses:
    for j in neutralino_masses:
        if j < i:
            if(not os.path.exists(limdir + 'limit_' + str(int(i)) + '_' + str(int(j)) +  '.txt')):
                continue
            with open(limdir + 'limit_' + str(int(i)) + '_' + str(int(j)) +  '.txt') as fin:
                for line in fin:

                    if re.search('16.0%', line):
                        l = line.split('<')
                        v = float(l[1].strip())
                        if v > 0:
                            hExpLimDown.SetBinContent(
                                hExpLimDown.GetXaxis().FindBin(i),
                                hExpLimDown.GetYaxis().FindBin(j), v
                                )

                    if re.search('50.0%', line):
                        l = line.split('<')
                        v = float(l[1].strip())
                        if v > 0:
                            #print(f"Setting bin ({i}, {j}) = {v}")
                            hExpLimMid.SetBinContent(
                                hExpLimMid.GetXaxis().FindBin(i),
                                hExpLimMid.GetYaxis().FindBin(j), v
                                )

                    if re.search('84.0%', line):
                        l = line.split('<')
                        v = float(l[1].strip())
                        if v > 0:
                            hExpLimUp.SetBinContent(
                                hExpLimUp.GetXaxis().FindBin(i),
                                hExpLimUp.GetYaxis().FindBin(j), v
                                )
            if(not os.path.exists(limdirdata + 'limit_' + str(int(i)) + '_' + str(int(j)) +  '.txt')):
                continue
            with open(limdirdata + 'limit_' + str(int(i)) + '_' + str(int(j)) +  '.txt') as findata:
                for line in findata:
                    if re.search('Observed', line):
                        l = line.split('<')
                        v = float(l[1].strip())
                        if v > 0:
                            #if v > 10: 
                               #print("bad point ", i, j, v) 
                            hObsLim.SetBinContent(
                                hObsLim.GetXaxis().FindBin(i),
                                hObsLim.GetYaxis().FindBin(j), v
                                )
                            # Define the original histogram


hObsLim.GetZaxis().SetRangeUser(0.0, 40.0)
hExpLimDown.GetZaxis().SetRangeUser(0.0, 40.0)
hExpLimMid.GetZaxis().SetRangeUser(0.0, 40.0)
hExpLimUp.GetZaxis().SetRangeUser(0.0, 40.0)
ROOT.gStyle.SetPadRightMargin(0.16)
#canvas = ROOT.TCanvas('can', 'can', 600,600)
#canvas.cd()
#canvas.SetLogz() 
#red2 = np.array([0.5, 0.5, 1.0, 1.0, 1.0])
#green2 = np.array([0.5, 1.0, 1.0, 0.6, 0.5])
#blue2 = np.array([1.0, 1.0, 0.5, 0.4, 0.5])
#length2 = np.array([0.0, 0.34, 0.61, 0.84, 1.0])
#table2 = ROOT.TColor.CreateGradientColorTable(5, length2, red2, green2, blue2, 255)
#limitPalette2 = [table2+i for i in range(1,255+1)]
#hObsLim.SetContour(255)
#hObsLim.Draw("colz")
#

def MakeLimitGraph(h, s):

    g = ROOT.TGraph()
    if(s != ''):
        g.SetName(s)
    n = 0
    m = 19
    max_x = None 

    for i in range(1, h.GetNbinsY()+1):

        x = 0.0
        y = h.GetYaxis().GetBinCenter(i)

        for j in range(1, h.GetNbinsX() + 1):

            xLow = h.GetXaxis().GetBinCenter(j)
            xHigh = h.GetXaxis().GetBinCenter(1+j)
            zLow = h.GetBinContent(j, i)
            zHigh = h.GetBinContent(1 + j, i)

            if zLow < 1.0 and zLow > 0.01 and zHigh > 1.0:
                x = xLow + (xHigh - xLow) * (1.0 - zLow) / (zHigh - zLow)
                g.SetPoint(n, x, y)
                
                #if h == hExpLimMid:
                #    print(f"Graph '{s}': Point {n} -> X: {x}, Y: {y}")

                if max_x is None or (x is not None and x > max_x):
                    max_x = x

                n += 1
                break
    print(f"Graph '{s}': Maximum X value = {max_x}")
    return g

def ScaleLimit(h, s, scaleByXSec, sigma):

    h_new = h.Clone(s)

    for i in range(1, h_new.GetNbinsX()+1):
        for j in range(1, h_new.GetNbinsY()+1):

            m = h_new.GetXaxis().GetBinCenter(i)
            r = h_new.GetBinContent(i, j)
            foundM = False

            if m == 1487.5:
                m = 1500.0

            with open('ULCrossSectionT5WG.txt') as f_xSec:
                for line in f_xSec:

                    l = line.split()
                    mSusy = float(l[0])
                    xSec = float(l[1])*1000*0.5
                    unc = float(l[2])

                    if mSusy == m: # careful!
                        if scaleByXSec:
                            r *= xSec
                        r /= 1.0 + sigma*unc/100.0
                        foundM = True
                        break
            if not foundM:
                r = -1.0
            #if (h_new.GetYaxis().GetBinCenter(j) > m -50):
            #    r = -1.0

            h_new.SetBinContent(i, j, r)

    return h_new

def RegularPlot(h):
    h_new = h.Clone()
    for i in range(1, h.GetNbinsX()+1):
        for j in range(1, h.GetNbinsY()+1):
            if(h.GetYaxis().GetBinCenter(j) <= h.GetXaxis().GetBinCenter(i)-50):
                #print(h.GetXaxis().GetBinCenter(i), h.GetYaxis().GetBinCenter(j), h.GetBinContent(i,j))
                if(h.GetBinContent(i,j) <= 0):
                    #print(h.GetXaxis().GetBinCenter(i), h.GetYaxis().GetBinCenter(j), h.GetBinContent(i,j))
                    r_left=0
                    r_right=0
                    r_up=0
                    r_do=0
                    for lbin in range(1, i-1):
                        if(r_left > 0):
                            continue
                        if(h.GetBinContent(i-lbin, j) > 0):
                            r_left = h.GetBinContent(i-lbin, j)
                            #print('left',lbin, h.GetBinContent(i-lbin, j))
                    for rbin in range(1, h.GetNbinsX()-i):
                        if(r_right > 0):
                            continue
                        if(h.GetBinContent(i+rbin, j) > 0):
                            r_right = h.GetBinContent(i+rbin, j)
                            #print('right',rbin, h.GetBinContent(i+rbin, j))
                    for ubin in range(1, h.GetNbinsY()-j):
                        if(r_up > 0):
                            continue
                        if(h.GetBinContent(i, j+ubin) > 0):
                            r_up = h.GetBinContent(i, j+ubin)
                            #print('up',ubin, h.GetBinContent(i, j+ubin))
                    for dbin in range(1, j-1):
                        if(r_do > 0):
                            continue
                        if(h.GetBinContent(i, j-dbin) > 0):
                            r_do = h.GetBinContent(i, j-dbin)
                            #print('down',dbin,h.GetBinContent(i, j-dbin))
                    if(r_up > 0 and r_do > 0):
                        h.SetBinContent(i,j, (r_up+r_do)/2)
                    else:
                        h.SetBinContent(i,j, (r_left+r_right)/2)
                    #print(h.GetXaxis().GetBinCenter(i), h.GetYaxis().GetBinCenter(j), h.GetBinContent(i,j))

def MergePlot(h): 
    h_new = ROOT.TProfile2D('obs2', 'obs2', 20, 1100, 2100, 84, 0, 2100);
    for i in range(1, h.GetNbinsX()+1):
        for j in range(1, h.GetNbinsY()+1):
            if(h.GetYaxis().GetBinCenter(j) <= h.GetXaxis().GetBinCenter(i)):
                if(h.GetBinContent(i,j) > 0):
                    h_new.Fill(h.GetXaxis().GetBinCenter(i), h.GetYaxis().GetBinCenter(j), h.GetBinContent(i,j))
    return h_new


ROOT.gStyle.SetPadRightMargin(0.15) 
ROOT.gStyle.SetPadLeftMargin(0.15)   
c_limit = ROOT.TCanvas('c_limit', '', 580, 620)
c_limit.SetLogz()

h_xSecLimit = ScaleLimit(hObsLim, 'obs', True, 0.0)
RegularPlot(h_xSecLimit)
h_xSecMerge = MergePlot(h_xSecLimit)
h_observedLimitLow = ScaleLimit(hObsLim, 'obs', False, -1.0)
h_observedLimitHigh = ScaleLimit(hObsLim, 'obs' , False, 1.0)


#============================================
h_dummy = h_xSecLimit.Clone("h_dummy")
h_dummy.Reset()
h_dummy.GetXaxis().SetTitle("M_{#tilde{g}} (GeV)")
h_dummy.GetYaxis().SetTitle("M_{#tilde{#chi}_{1}^{0}} (GeV)")
h_dummy.GetZaxis().SetLabelOffset(-0.003)
h_dummy.GetZaxis().SetLabelSize(0.03)
h_dummy.GetZaxis().SetTitleOffset(0.9)
h_dummy.GetZaxis().SetLabelSize(0.04)
h_dummy.GetZaxis().SetTitle("")
h_dummy.SetContour(255)
h_dummy.GetZaxis().SetRangeUser(0.5, 40)
h_dummy.SetTitle("")
h_dummy.GetXaxis().SetNdivisions(505)

h_dummy.GetYaxis().SetTitleOffset(1.8)
h_dummy.Draw("colz")
c_limit.Update()

z_title = ROOT.TLatex()
z_title.SetTextAngle(90)
z_title.SetTextSize(0.035)
z_title.SetTextAlign(31)  # Right-aligned
z_title.DrawLatexNDC(0.96, 0.6, "95% CL limit on #sigma (pb)")

c_limit.Modified()
c_limit.Update()

#============================================



red = np.array([0.5, 0.5, 1.0, 1.0, 1.0])
green = np.array([0.5, 1.0, 1.0, 0.6, 0.5])
blue = np.array([1.0, 1.0, 0.5, 0.4, 0.5])
length = np.array([0.0, 0.34, 0.61, 0.84, 1.0])

table = ROOT.TColor.CreateGradientColorTable(5, length, red, green, blue, 255)
limitPalette = [table+i for i in range(1,255+1)]


h_xSecLimit.GetXaxis().SetTitle('M_{#tilde{g}} (GeV)')
h_xSecLimit.GetYaxis().SetTitle('M_{#tilde{#chi}_{1}^{0}} (GeV)  ')
h_xSecLimit.GetZaxis().SetLabelOffset(-0.003)
h_xSecLimit.GetZaxis().SetLabelSize(0.03)
h_xSecLimit.GetZaxis().SetTitleOffset(0.9)
h_xSecLimit.GetZaxis().SetLabelSize(0.04)
h_xSecLimit.GetZaxis().SetTitle('')
h_xSecLimit.SetContour(255)
#h_xSecLimit.GetZaxis().SetRangeUser(0.001, 0.007)
h_xSecLimit.GetZaxis().SetRangeUser(0.5, 40)
h_xSecLimit.SetTitle('')

h_xSecLimit.Draw('colz same')

g_expectedLimit = MakeLimitGraph(hExpLimMid, 'hExpLimMid')
g_expectedLimitLow = MakeLimitGraph(hExpLimDown, 'hExpLimDown')
g_expectedLimitHigh = MakeLimitGraph(hExpLimUp, 'hExpLimUp')
g_observedLimit = MakeLimitGraph(hObsLim, 'hObsLim')
g_observedLimitLow = MakeLimitGraph(h_observedLimitLow, 'h_observedLimitLow')
g_observedLimitHigh = MakeLimitGraph(h_observedLimitHigh, 'h_observedLimitHigh')

g_expectedLimit.SetLineColor(ROOT.kRed+1)
g_expectedLimit.SetLineWidth(2)
g_expectedLimit.Draw('C same')

g_expectedLimitLow.SetLineColor(ROOT.kRed+1)
g_expectedLimitLow.SetLineStyle(2)
g_expectedLimitLow.SetLineWidth(2)
g_expectedLimitLow.Draw('C same')

g_expectedLimitHigh.SetLineColor(ROOT.kRed+1)
g_expectedLimitHigh.SetLineStyle(2)
g_expectedLimitHigh.SetLineWidth(2)
g_expectedLimitHigh.Draw('C same')

g_observedLimit.SetLineWidth(2)
g_observedLimit.Draw('C same')

g_observedLimitLow.SetLineStyle(2)
g_observedLimitLow.SetLineWidth(2)
#g_observedLimitLow.Draw('C same')

g_observedLimitHigh.SetLineStyle(2)
g_observedLimitHigh.SetLineWidth(2)
#g_observedLimitHigh.Draw('C same')

le_limit = ROOT.TLegend(0.2, 0.73, 0.5, 0.82)

le_limit.SetFillStyle(0)
le_limit.SetLineColor(ROOT.kWhite)
le_limit.SetTextSize(0.0325)
le_limit.SetTextFont(42)
le_limit.AddEntry(g_observedLimit, 'SUS-17-012 (Obs)', 'l')
#le_limit.AddEntry(g_observedLimitLow, 'SUS-17-012 (Obs) #pm 1 #sigma_{th}','l')
le_limit.AddEntry(g_expectedLimit, 'Exp. limit (Run 2)', 'l')
#le_limit.AddEntry(g_expectedLimitLow, 'Exp. limit #pm 1 #sigma_{ex}','l')
le_limit.Draw('same');

la_limit = ROOT.TLatex()
la_limit.SetNDC()

la_limit.SetTextAlign(11)
la_limit.SetTextSize(0.044)
la_limit.SetTextFont(42)
#la_limit.DrawLatex(0.6, 0.91, '35.8 fb^{-1} (13 TeV)')

la_limit.SetTextAlign(11)
la_limit.SetTextSize(0.044)
la_limit.SetTextFont(62)
la_limit.DrawLatex(0.21, 0.91, 'CMS (13 TeV)')

la_limit.SetTextAlign(11)
la_limit.SetTextSize(0.044)
la_limit.SetTextFont(42)
#la_limit.DrawLatex(0.2, 0.91, 'Preliminary')

la_limit.SetTextSize(0.04)
la_limit.SetTextFont(42)
la_limit.DrawLatex(0.2, 0.85, 'pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow qq#tilde{#chi}_{1}^{0/#pm}, #tilde{#chi}_{1}^{0/#pm} #rightarrow #gamma/W^{#pm}#tilde{G}')

#la_limit.DrawLatex(0.2, 0.84, 'NLO + NLL exclusion')
ROOT.gPad.RedrawAxis("G")

c_limit.Draw()
c_limit.SaveAs(f'/uscms_data/d3/tmishra/Output/T5WG_{year}_2016-17-2018.pdf')

file_out = ROOT.TFile(f'/uscms_data/d3/tmishra/Output/T5WG_{year}.root', 'RECREATE')
file_out.cd()
h_xSecLimit.Write()
g_expectedLimit.Write()
g_expectedLimitLow.Write()
g_expectedLimitHigh.Write()
g_observedLimit.Write()
g_observedLimitLow.Write()
g_observedLimitHigh.Write()
h_xSecMerge.GetXaxis().SetTitle('M_{#tilde{g}} (GeV)')
h_xSecMerge.GetYaxis().SetTitle('M_{#tilde{#chi}_{1}^{0}} (GeV)  ')
h_xSecMerge.GetZaxis().SetLabelOffset(-0.003)
h_xSecMerge.GetZaxis().SetLabelSize(0.03)
h_xSecMerge.GetZaxis().SetTitleOffset(0.9)
h_xSecMerge.GetZaxis().SetLabelSize(0.04)
h_xSecMerge.GetZaxis().SetTitle('')
h_xSecMerge.SetContour(255)
h_xSecMerge.GetZaxis().SetRangeUser(0.5, 40)
h_xSecMerge.SetTitle('')
h_xSecMerge.Write()
