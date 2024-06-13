from ROOT import *
import numpy as np
import mmap
import time
import sys
#file_temp = TFile('plot_susyMasses_template.root')
#h_template = file_temp.Get('h_susyMasses_template')
#dictXsec={}
#dictXsecUnc={}
#
#with open('LatestXsecGluGlu.txt', 'r') as input:
#        for line in input:
#                elements = line.rstrip().split("|")
#                dictXsec[int(elements[1])]=elements[2]
#                dictXsecUnc[int(elements[1])]=elements[3]
#
#
file_temp  = TFile('../test/Result/tmp/signalTree_M1M2.root', 'read')
h_template = file_temp.Get('SUSYMass')
dictXsec={}
dictXsecUnc={}

with open('./CrossSectionM1M2.txt') as f_xSec:
    for line in f_xSec:
        l = line.split()
        if(int(l[0]) > 200 and int(l[1]) > 200):
          dictXsec[int(l[0]),int(l[1])] = float(l[8])*1000
          dictXsecUnc[int(l[0]),int(l[1])] = float(l[9])*1000

mGo=[]
mLsp=[]
limit=[]

for i in range(1, h_template.GetXaxis().GetNbins() + 1):
    for j in range(1, h_template.GetYaxis().GetNbins() + 1):
        #if h_template.GetBinContent(i, j) == 1:
        mass_gluino = int(h_template.GetXaxis().GetBinCenter(i))
        mass_neutralino = int(h_template.GetYaxis().GetBinCenter(j))
        if(mass_gluino <= 200 or mass_gluino > 1500 or mass_neutralino <= 200 or mass_neutralino > 1500):
            continue
        else:	
            mGo.append(mass_gluino)
            mLsp.append(mass_neutralino)
MuScan=TGraph2D()
MuScan.SetName("MuScan")
MuScanXsec=TGraph2D()
MuScanXsec.SetName("MuScanXsec")
MuScanExpXsec=TGraph2D()
MuScanExpXsec.SetName("MuScanExpXsec")

MuScanSup=TGraph2D()
MuScanSup.SetName("MuScanSup")
MuScanSdn=TGraph2D()
MuScanSdn.SetName("MuScanSdn")
MuScanObs=TGraph2D()
MuScanObs.SetName("MuScanObs")
MuScanObsSup=TGraph2D()
MuScanObsSup.SetName("MuScanObsSup")
MuScanObsSdn=TGraph2D()
MuScanObsSdn.SetName("MuScanObsSdn")
for m in range(len(mGo)):
  filein=TFile("../DataCard/M1M2_DataCards/Limit/higgsCombinecount_masses_%d_%d_.AsymptoticLimits.mH120.root" %(mGo[m], mLsp[m]),"READ")
  if(filein.IsOpen()):
    t = filein.Get("limit")
    t.GetEntry(2);
    ExpUL= t.limit #* float(dictXsec.get(mGo[m]))
    ExpULXSec= t.limit*float(dictXsec[mGo[m],mLsp[m]])
    t.GetEntry(5)
    ObsULXSec= t.limit*float(dictXsec[mGo[m],mLsp[m]])
    ObsUL=t.limit#*float(dictXsec.get(mGo[m]))
    t.GetEntry(1)
    ExpULSigmaUp=t.limit #*float(dictXsec.get(mGo[m]))
    t.GetEntry(3)
    ExpULSigmaDn=t.limit #*float(dictXsec.get(mGo[m]))
    shiftUp=1.0/(1-(float(dictXsecUnc[mGo[m],mLsp[m]])/dictXsec[mGo[m],mLsp[m]]));
    shiftDn=1.0/(1+(float(dictXsecUnc[mGo[m],mLsp[m]])/dictXsec[mGo[m],mLsp[m]]));
    ObsULDn=shiftDn*ObsUL
    ObsULUp=shiftUp*ObsUL
    MuScan.SetPoint(MuScan.GetN(), mGo[m], mLsp[m], ExpUL)
    MuScanSup.SetPoint(MuScanSup.GetN(), mGo[m], mLsp[m], ExpULSigmaUp)
    MuScanSdn.SetPoint(MuScanSdn.GetN(), mGo[m], mLsp[m], ExpULSigmaDn)
    MuScanObs.SetPoint( MuScanObs.GetN(), mGo[m], mLsp[m], ObsUL)
    MuScanObsSup.SetPoint( MuScanObsSup.GetN(), mGo[m], mLsp[m], ObsULUp)
    MuScanObsSdn.SetPoint( MuScanObsSdn.GetN(), mGo[m], mLsp[m], ObsULDn)
    MuScanXsec.SetPoint(MuScanXsec.GetN(), mGo[m], mLsp[m],ObsULXSec)
    MuScanExpXsec.SetPoint(MuScan.GetN(), mGo[m], mLsp[m], ExpULXSec)

MuScan.SetName("MuScan")
MuScan.SetNpx(128)
MuScan.SetNpy(160)
MuScanSup.SetNpx(128)
MuScanSup.SetNpx(160)
MuScanSdn.SetNpx(128)
MuScanSdn.SetNpx(160)
MuScanObs.SetNpx(128)
MuScanObs.SetNpy(160)
MuScanObsSup.SetNpx(128)
MuScanObsSup.SetNpy(160)
MuScanObsSdn.SetNpx(128)
MuScanObsSdn.SetNpy(160)
MuScanXsec.SetNpx(128)
MuScanXsec.SetNpy(160)
MuScanExpXsec.SetNpx(128)
MuScanExpXsec.SetNpy(160)
hExplim=MuScan.GetHistogram()
hExplimSup=MuScanSup.GetHistogram()
hExplimSdn=MuScanSdn.GetHistogram()
hObslim=MuScanObs.GetHistogram()
hObslimSup=MuScanObsSup.GetHistogram()
hObslimSdn=MuScanObsSdn.GetHistogram()
MassScan2D=MuScanXsec.GetHistogram()
MassScan2DExp=MuScanExpXsec.GetHistogram()
c=TCanvas("c","",800,800);
MuScan.Draw("colz")
MuScanSup.Draw("colz")
MuScanSdn.Draw("colz")
MuScanObs.Draw("colz")
MuScanObsSup.Draw("colz")
MuScanObsSdn.Draw("colz")
ExpLim=TGraph()
ExpLim.SetName("ExpLim")

ExpLim= MuScan.GetContourList(1.0);
ExpLimSup= MuScanSup.GetContourList(1.0);
ExpLimSdn= MuScanSdn.GetContourList(1.0);
ObsLim= MuScanObs.GetContourList(1.0);
ObsLimSup= MuScanObsSup.GetContourList(1.0);
ObsLimSdn= MuScanObsSdn.GetContourList(1.0);
fileOut=TFile("MassScanStealth.root", "RECREATE")

gr1=TGraph()
gr2=TGraph()
gr3=TGraph()
gr4=TGraph()
gr5=TGraph()
gr6=TGraph()
gr1=ObsLim.At(0)
gr2=ObsLim.At(1)
gr3=ObsLim.At(2)
gr4=ObsLim.At(3)
gr5=ObsLim.At(4)
gr6=ObsLim.At(5)


ExpLim.Write("hExpLimMid")
ExpLimSup.Write("hExpLimUp")
ExpLimSdn.Write("hExpLimDown")
ObsLim.Write("hObsLim")
ObsLimSup.Write("h_observedLimitHigh")
ObsLimSdn.Write("h_observedLimitLow")
MassScan2D.Write("obs")
MassScan2DExp.Write("MassScan2DExp")
MuScan.Write()
gr1.Write("hObsLim1")
gr2.Write("hObsLim2")
gr3.Write("hObsLim3")
gr4.Write("hObsLim4")
gr5.Write("hObsLim5")
gr6.Write("hObsLim6")
fileOut.Close()

