#python makeT5WGLimit.py --limitFolder /eos/uscms/store/user/tmishra/CombinedLimit/Sensitive_myFile/limitsANbins/ --outputFile T5WG_ANbins

import os
import numpy as np
import re
import ROOT
import argparse
import sys
ROOT.gROOT.SetBatch(ROOT.kTRUE)


limdir = '/eos/uscms/store/user/tmishra/CombinedLimit/Sensitive_myFile/signifANbin/'
limdirdata = '/eos/uscms/store/user/tmishra/CombinedLimit/Sensitive_myFile/signifANbin/'

limdir2 = '/eos/uscms/store/user/tmishra/CombinedLimit/Sensitive_myFile/signifMT150/'
limdirdata2 = '/eos/uscms/store/user/tmishra/CombinedLimit/Sensitive_myFile/signifMT150/'


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--MET1', dest='METbin1', type=int)
parser.add_argument('--MET2', dest='METbin2', type=int)
args = parser.parse_args()

met1 = args.METbin1
met2 = args.METbin2


gluino_masses = np.array([])
neutralino_masses = np.array([])

file_in = ROOT.TFile('/uscms_data/d3/mengleis/CMSSW_Work/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/tools/signalTree_T5WG.root', 'read')
h_SUSYmass = file_in.Get('SUSYMass')

ROOT.gStyle.SetOptStat(0)

for xbin in range(1, h_SUSYmass.GetXaxis().GetNbins() + 1):
    if(int(h_SUSYmass.GetXaxis().GetBinCenter(xbin)) >= 1100):
        gluino_masses = np.append(gluino_masses, int(h_SUSYmass.GetXaxis().GetBinCenter(xbin)))
for ybin in range(1, h_SUSYmass.GetYaxis().GetNbins() + 1):
    neutralino_masses = np.append(neutralino_masses, int(h_SUSYmass.GetYaxis().GetBinCenter(ybin)) )

print(gluino_masses)
print(neutralino_masses)
nbinsx = len(gluino_masses) - 1
nbinsy = len(neutralino_masses) - 1

hExpLimUp = ROOT.TH2D('up', 'up', nbinsx, gluino_masses, nbinsy, neutralino_masses)
hExpLimMid = ROOT.TH2D('mid', 'mid', nbinsx, gluino_masses, nbinsy, neutralino_masses)
hExpLimDown = ROOT.TH2D('down', 'down', nbinsx, gluino_masses, nbinsy, neutralino_masses)
hObsLim = ROOT.TH2D('obs', 'obs', nbinsx, gluino_masses, nbinsy, neutralino_masses)

for i in gluino_masses:
    for j in neutralino_masses:
        if j < i:
            if(not os.path.exists(limdir + 'signif_' + str(int(i)) + '_' + str(int(j)) +  '.txt')):
                continue
            with open(limdir + 'signif_' + str(int(i)) + '_' + str(int(j)) +  '.txt') as fin:
                for line in fin:
                    if re.search('Significance:', line):
                        l = line.split('Significance:')
                        v = float(l[1].strip())
#                        hExpLimMid.SetBinContent(
#                                	hExpLimMid.GetXaxis().FindBin(i),
#                                	hExpLimMid.GetYaxis().FindBin(j), v
#                                	)
            if(not os.path.exists(limdir2 + 'signif_' + str(int(i)) + '_' + str(int(j)) +  '.txt')):
                continue
            with open(limdir2 + 'signif_' + str(int(i)) + '_' + str(int(j)) +  '.txt') as fin2:
                for line2 in fin2:
                    if re.search('Significance:', line2):
                        l2 = line2.split('Significance:')
                        v2 = float(l2[1].strip())
		#		print v2
                        
	    #print (v2-v) 
	    hExpLimMid.SetBinContent(
                                hExpLimMid.GetXaxis().FindBin(i),
                               	hExpLimMid.GetYaxis().FindBin(j), (v2-v)
                               	)
	    print v,v2,v2-v

	    #if((i == 1650 and j == 10) or (i == 1650 and j == 25) or (i == 1650 and j == 100) or (i == 1650 and j == 200) or (i == 1700 and j == 10) or (i == 1700 and j == 25) or (i == 1700 and j == 50) or (i == 1700 and j == 100) or (i == 1700 and j == 400) or (i == 1750 and j == 100) or (i == 1750 and j == 500) or (i == 1800 and j == 50) or (i == 1800 and j == 100) or (i == 1800 and j == 200) or (i == 1800 and j == 500) or (i == 1800 and j == 800) or (i == 1800 and j == 1000) or (i == 1800 and j == 1500) or (i == 1800 and j == 1650) or (i == 1800 and j == 1700) or(i == 1800 and j == 1750) or (i == 1800 and j == 1790) or (i == 1850 and j == 100) or (i == 1850 and j == 500) or (i == 1850 and j == 1000) or (i == 1850 and j == 1500) or (i == 1850 and j == 1750) or (i == 1850 and j == 1800)):





hObsLim.GetZaxis().SetRangeUser(0.0, 40.0)
hExpLimMid.GetZaxis().SetRangeUser(-20.0, 20.0)
ROOT.gStyle.SetPadRightMargin(0.16)
ROOT.gStyle.SetPadLeftMargin(0.1)
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

    for i in range(1, h.GetNbinsY()+1):

        x = 0.0
        y = h.GetYaxis().GetBinCenter(i)

        for j in range(1, h.GetNbinsX() + 1):

            xLow = h.GetXaxis().GetBinCenter(j)
            xHigh = h.GetXaxis().GetBinCenter(1+j)
            zLow = h.GetBinContent(j, i)
            zHigh = h.GetBinContent(1 + j, i)
	
            if zLow != 0 and zHigh == 0:
	    	print xLow,h.GetYaxis().GetBinCenter(i),zLow,zHigh
            #if zLow < 2.0 and zLow > 0.01 and zHigh > 2.0:
            #if zLow < 2.0 and zLow > 0.01 and zHigh > 3.0:
                #x = xLow + (xHigh - xLow) * (4.0 - zLow) / (zHigh - zLow)
	        #print xLow,zLow,zHigh

                #g.SetPoint(n, x, y)
                #n += 1
                break

    return g

g_expectedLimit = MakeLimitGraph(hExpLimMid, 'hExpLimMid')

def ScaleLimit(h, s, scaleByXSec, sigma):

    h_new = h.Clone(s)

    for i in range(1, h_new.GetNbinsX()+1):
        for j in range(1, h_new.GetNbinsY()+1):

            m = h_new.GetXaxis().GetBinCenter(i)
            r = h_new.GetBinContent(i, j)
            foundM = False

            if m == 1487.5:
                m = 1500.0

            #tpm: is this the right file to use?
            with open('./CrossSectionT5WG.txt') as f_xSec:
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
                print(h.GetXaxis().GetBinCenter(i), h.GetYaxis().GetBinCenter(j), h.GetBinContent(i,j))
                if(h.GetBinContent(i,j) <= 0):
                    print(h.GetXaxis().GetBinCenter(i), h.GetYaxis().GetBinCenter(j), h.GetBinContent(i,j))
                    r_left=0
                    r_right=0
                    r_up=0
                    r_do=0
                    for lbin in range(1, i-1):
                        if(r_left > 0):
                            continue
                        if(h.GetBinContent(i-lbin, j) > 0):
                            r_left = h.GetBinContent(i-lbin, j)
                            print('left',lbin, h.GetBinContent(i-lbin, j))
                    for rbin in range(1, h.GetNbinsX()-i):
                        if(r_right > 0):
                            continue
                        if(h.GetBinContent(i+rbin, j) > 0):
                            r_right = h.GetBinContent(i+rbin, j)
                            print('right',rbin, h.GetBinContent(i+rbin, j))
                    for ubin in range(1, h.GetNbinsY()-j):
                        if(r_up > 0):
                            continue
                        if(h.GetBinContent(i, j+ubin) > 0):
                            r_up = h.GetBinContent(i, j+ubin)
                            print('up',ubin, h.GetBinContent(i, j+ubin))
                    for dbin in range(1, j-1):
                        if(r_do > 0):
                            continue
                        if(h.GetBinContent(i, j-dbin) > 0):
                            r_do = h.GetBinContent(i, j-dbin)
                            print('down',dbin,h.GetBinContent(i, j-dbin))
                    if(r_up > 0 and r_do > 0):
                        h.SetBinContent(i,j, (r_up+r_do)/2)
                    else:
                        h.SetBinContent(i,j, (r_left+r_right)/2)
                    print(h.GetXaxis().GetBinCenter(i), h.GetYaxis().GetBinCenter(j), h.GetBinContent(i,j))

def MergePlot(h): 
    h_new = ROOT.TProfile2D('obs2', 'obs2', 20, 1100, 2100, 84, 0, 2100);
    for i in range(1, h.GetNbinsX()+1):
        for j in range(1, h.GetNbinsY()+1):
            if(h.GetYaxis().GetBinCenter(j) <= h.GetXaxis().GetBinCenter(i)):
                if(h.GetBinContent(i,j) > 0):
                    h_new.Fill(h.GetXaxis().GetBinCenter(i), h.GetYaxis().GetBinCenter(j), h.GetBinContent(i,j))
    return h_new



ROOT.gStyle.SetPadRightMargin(0.16)
ROOT.gStyle.SetPadLeftMargin(0.1)
c_limit = ROOT.TCanvas('c_limit', '', 580, 620)
c_limit.cd()
c_limit.SetLogz()
c_limit.SetLeftMargin(0.16)


red = np.array([0.5, 0.5, 1.0, 1.0, 1.0])
green = np.array([0.5, 1.0, 1.0, 0.6, 0.5])
blue = np.array([1.0, 1.0, 0.5, 0.4, 0.5])
length = np.array([0.0, 0.34, 0.61, 0.84, 1.0])

table = ROOT.TColor.CreateGradientColorTable(5, length, red, green, blue, 255)
limitPalette = [table+i for i in range(1,255+1)]

hExpLimMid.GetXaxis().SetTitle('M_{#tilde{g}} (GeV)')
hExpLimMid.GetYaxis().SetTitle('M_{#tilde{#chi}_{1}^{0}} (GeV)  ')
hExpLimMid.GetZaxis().SetTitle(' Difference in significance (MT>150 - MT>100)')
hExpLimMid.GetZaxis().SetLabelOffset(-0.003)
hExpLimMid.GetZaxis().SetTitleOffset(1.2)

hExpLimMid.GetZaxis().SetLabelSize(0.033)
hExpLimMid.GetXaxis().SetLabelSize(0.025)
hExpLimMid.GetYaxis().SetLabelSize(0.03)
hExpLimMid.SetContour(255)
hExpLimMid.GetZaxis().SetRangeUser(-1.0, 5.0)
#hExpLimMid.GetZaxis().SetRangeUser(0.5, 15)
hExpLimMid.SetTitle('')
hExpLimMid.Draw('colz')

g_expectedLimit.SetLineColor(ROOT.kRed+1)
g_expectedLimit.SetLineWidth(2)
g_expectedLimit.Draw('C same')

le_limit = ROOT.TLegend(0.15, 0.72, 0.4, 0.82)
le_limit.SetFillStyle(0)
le_limit.SetLineColor(ROOT.kWhite)
le_limit.SetTextSize(0.0325)
le_limit.SetTextFont(42)

la_limit = ROOT.TLatex()
la_limit.SetNDC()

la_limit.SetTextAlign(11)
la_limit.SetTextSize(0.044)
la_limit.SetTextFont(42)
#la_limit.DrawLatex(0.6, 0.91, '35.8 fb^{-1} (13 TeV)')

la_limit.SetTextAlign(11)
la_limit.SetTextSize(0.044)
la_limit.SetTextFont(62)
#la_limit.DrawLatex(0.15, 0.91, 'CMS')

la_limit.SetTextAlign(11)
la_limit.SetTextSize(0.044)
la_limit.SetTextFont(42)
#la_limit.DrawLatex(0.2, 0.91, 'Preliminary')

la_limit.SetTextSize(0.03)
la_limit.SetTextFont(42)
la_limit.DrawLatex(0.2, 0.85, 'pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrowqq#tilde{#chi}_{1}^{0/#pm}, #tilde{#chi}_{1}^{0/#pm}#rightarrow#gamma/W^{#pm}#tilde{G}')

#la_limit.DrawLatex(0.2, 0.84, 'NLO + NLL exclusion')

c_limit.Draw()
c_limit.SaveAs('T5WG_ANbin_significance_subtract.pdf')

#file_out = ROOT.TFile('T5WG_ANbins_Signif.root', 'RECREATE')
#file_out.cd()
#g_expectedLimit.Write()

