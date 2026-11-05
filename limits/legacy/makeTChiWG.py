import os
import numpy as np
import re
import ROOT

limdir = './tchiwg/'

limit_files = os.listdir(limdir)

neutralino_masses = np.array([])

file_in = ROOT.TFile('./signalTree_TChiWG_200_320.root', 'read')
h_SUSYmass = file_in.Get('p_TChiWGMASS')


for xbin in range(1, h_SUSYmass.GetXaxis().GetNbins() + 1):
    neutralino_masses = np.append(neutralino_masses, int(h_SUSYmass.GetXaxis().GetBinCenter(xbin)))

print(neutralino_masses)
nbinsx = len(neutralino_masses) - 1

hExpLimMid = ROOT.TGraph(nbinsx)
#hExpLimOneSigma = ROOT.TGraphAsymmErrors(nbinsx)
#hExpLimTwoSigma = ROOT.TGraphAsymmErrors(nbinsx)

for i in neutralino_masses:
    if(not os.path.exists(limdir + 'limit_' + str(int(i)) + '.txt')):
        continue
    with open(limdir + 'limit_' + str(int(i)) + '.txt') as fin:
        mSusy = 0
        xSec = 0
        unc = 0
        ExpLimMid = 0
        ExpLimOneUp = 0
        ExpLimOneDown = 0
        ExpLimTwoUp = 0
        ExpLimTwoDown = 0
        with open('./SusyCrossSections13TevChaNeu.txt') as f_xSec:
            for xsline in f_xSec:
                xsvalue = xsline.split()
                if(int(xsvalue[0]) == int(i)):
                    mSusy = float(xsvalue[0])
                    xSec = float(xsvalue[1])
                    unc = float(xsvalue[2])
        if(xSec==0):
            continue

        for line in fin:
            if re.search('16.0%', line):
                l = line.split('<')
                v = float(l[1].strip())
                if v > 0:
                    ExpLimOneDown = v

            if re.search('50.0%', line):
                l = line.split('<')
                v = float(l[1].strip())
                if v > 0:
                    ExpLimMid = v

            if re.search('84.0%', line):
                l = line.split('<')
                v = float(l[1].strip())
                if v > 0:
                    ExpLimOneUp = v


        hExpLimMid.SetPoint(hExpLimMid.GetXaxis().FindBin(i), i, ExpLimMid*xSec)
        #hExpLimOneSigma.SetPoint(hExpLimOneSigma.GetXaxis().FindBin(i), i, ExpLimMid*xSec)
        #hExpLimOneSigma.SetPointError(hExpLimOneSigma.GetXaxis().FindBin(i), 0,0,(ExpLimOneUp-ExpLimMid)*xSec, (ExpLimMid-ExpLimOneDown)*xSec)


#thecanvas = ROOT.TCanvas('thecanvas', 'the canvas',4,69,600,550)
#ROOT.gStyle.SetOptFit(1)
#ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetOptTitle(0)
#thecanvas.Range(5.587935e-06,-3.763583,937.5,2.11013)
#thecanvas.SetFillColor(0)
#thecanvas.SetBorderMode(0)
#thecanvas.SetBorderSize(2)
#thecanvas.SetLogy()
#thecanvas.SetGridx()
#thecanvas.SetGridy()
#thecanvas.SetTickx(1)
#thecanvas.SetTicky(1)
#thecanvas.SetLeftMargin(0.16)
#thecanvas.SetRightMargin(0.04)
#thecanvas.SetTopMargin(0.07)
#thecanvas.SetBottomMargin(0.13)
#thecanvas.SetFrameFillStyle(0)
#thecanvas.SetFrameBorderMode(0)
#thecanvas.SetFrameFillStyle(0)
#thecanvas.SetFrameBorderMode(0)
#thecanvas.cd()
#
#frame__1 = ROOT.TH1F('frame__1','',1000,300,1300)
#frame__1.SetMinimum(0.1)
#frame__1.SetMaximum(1000)
#frame__1.SetDirectory(0)
#frame__1.SetLineStyle(0)
#frame__1.SetMarkerStyle(20)
#frame__1.GetXaxis().SetTitle('Chargino mass [GeV]')
#frame__1.GetXaxis().SetNdivisions(505)
#frame__1.GetXaxis().SetLabelFont(42)
#frame__1.GetXaxis().SetLabelOffset(0.007)
#frame__1.GetXaxis().SetLabelSize(0.05)
#frame__1.GetXaxis().SetTitleSize(0.06)
#frame__1.GetXaxis().SetTitleOffset(0.9)
#frame__1.GetXaxis().SetTitleFont(42)
#frame__1.GetYaxis().SetTitle('95% CL upper limit on #sigma x BR [pb]')
#frame__1.GetYaxis().SetLabelFont(42)
#frame__1.GetYaxis().SetLabelOffset(0.007)
#frame__1.GetYaxis().SetLabelSize(0.05)
#frame__1.GetYaxis().SetTitleSize(0.06)
#frame__1.GetYaxis().SetTitleOffset(1.25)
#frame__1.GetYaxis().SetTitleFont(42)
#frame__1.GetZaxis().SetLabelFont(42)
#frame__1.GetZaxis().SetLabelOffset(0.007)
#frame__1.GetZaxis().SetLabelSize(0.05)
#frame__1.GetZaxis().SetTitleSize(0.06)
#frame__1.GetZaxis().SetTitleFont(42)
#frame__1.Draw()
#
#ci = ROOT.TColor.GetColor('#0000cc')
#hExpLimMid.SetLineColor(ci)
#hExpLimMid.SetLineStyle(9)
#hExpLimMid.SetLineWidth(3)
#hExpLimMid.SetMarkerStyle(20)
#
#Graph_Graph1 = ROOT.TH1F('Graph_Graph1','Graph',100,300,1300)
#Graph_Graph1.SetMinimum(0.1)
#Graph_Graph1.SetMaximum(1000)
#Graph_Graph1.SetDirectory(0)
#Graph_Graph1.SetStats(0)
#Graph_Graph1.SetLineStyle(0)
#Graph_Graph1.SetMarkerStyle(20)
#Graph_Graph1.GetXaxis().SetLabelFont(42)
#Graph_Graph1.GetXaxis().SetLabelOffset(0.007)
#Graph_Graph1.GetXaxis().SetLabelSize(0.05)
#Graph_Graph1.GetXaxis().SetTitleSize(0.06)
#Graph_Graph1.GetXaxis().SetTitleOffset(0.9)
#Graph_Graph1.GetXaxis().SetTitleFont(42)
#Graph_Graph1.GetYaxis().SetLabelFont(42)
#Graph_Graph1.GetYaxis().SetLabelOffset(0.007)
#Graph_Graph1.GetYaxis().SetLabelSize(0.05)
#Graph_Graph1.GetYaxis().SetTitleSize(0.06)
#Graph_Graph1.GetYaxis().SetTitleOffset(1.25)
#Graph_Graph1.GetYaxis().SetTitleFont(42)
#Graph_Graph1.GetZaxis().SetLabelFont(42)
#Graph_Graph1.GetZaxis().SetLabelOffset(0.007)
#Graph_Graph1.GetZaxis().SetLabelSize(0.05)
#Graph_Graph1.GetZaxis().SetTitleSize(0.06)
#Graph_Graph1.GetZaxis().SetTitleFont(42)
#hExpLimMid.SetHistogram(Graph_Graph1)
#
#hExpLimMid.SetLineColor(ROOT.kBlue)
#hExpLimMid.Draw('c')
#
#Graph_GraphA = ROOT.TH1F('Graph_GraphA','Graph',100,300,1300)
#Graph_GraphA.SetMinimum(0.1)
#Graph_GraphA.SetMaximum(1000)
#Graph_GraphA.SetDirectory(0)
#Graph_GraphA.SetStats(0)
#Graph_GraphA.SetLineStyle(0)
#Graph_GraphA.SetMarkerStyle(20)
#Graph_GraphA.GetXaxis().SetLabelFont(42)
#Graph_GraphA.GetXaxis().SetLabelOffset(0.007)
#Graph_GraphA.GetXaxis().SetLabelSize(0.05)
#Graph_GraphA.GetXaxis().SetTitleSize(0.06)
#Graph_GraphA.GetXaxis().SetTitleOffset(0.9)
#Graph_GraphA.GetXaxis().SetTitleFont(42)
#Graph_GraphA.GetYaxis().SetLabelFont(42)
#Graph_GraphA.GetYaxis().SetLabelOffset(0.007)
#Graph_GraphA.GetYaxis().SetLabelSize(0.05)
#Graph_GraphA.GetYaxis().SetTitleSize(0.06)
#Graph_GraphA.GetYaxis().SetTitleOffset(1.25)
#Graph_GraphA.GetYaxis().SetTitleFont(42)
#Graph_GraphA.GetZaxis().SetLabelFont(42)
#Graph_GraphA.GetZaxis().SetLabelOffset(0.007)
#Graph_GraphA.GetZaxis().SetLabelSize(0.05)
#Graph_GraphA.GetZaxis().SetTitleSize(0.06)
#Graph_GraphA.GetZaxis().SetTitleFont(42)
#hExpLimOneSigma.SetHistogram(Graph_GraphA)
#
#hExpLimOneSigma.SetName('Graph')
#hExpLimOneSigma.SetTitle('Graph')
#hExpLimOneSigma.SetFillColor(5)
#hExpLimOneSigma.SetLineWidth(3)
#hExpLimOneSigma.SetMarkerStyle(20)
#hExpLimOneSigma.Draw('3 same')
#hExpLimMid.Draw('C same')
#
#thecanvas.SaveAs('try.pdf')
