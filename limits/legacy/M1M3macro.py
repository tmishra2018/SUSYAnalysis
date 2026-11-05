from ROOT import *
import numpy as np
import mmap
import time
import sys
import os
flist=open("listofM1M3File.txt", 'r')
inputpath = "../DataCard/M1M3_DataCards/limit/" 

dictXsec={}
dictXsecUnc={}
with open('CrossSectionM1M3.txt', 'r') as input:
        for line in input:
                elements = line.split()
                dictXsec[(int(elements[0]),int(elements[1]))]=float(elements[8])
                dictXsecUnc[(int(elements[0]),int(elements[1]))]=float(elements[9])
limit=[]
fileOut=TFile("M1M3Scan.root", "RECREATE")

MuScan=TGraph2D()
MuScan.SetName("MuScan")
MuScanXsec=TGraph2D()
MuScanXsec.SetName("MuScanXsec")
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

for line in flist:
  fname=line.split('_')
  mGo = float(fname[3])
  mLsp = float(fname[4])
  if(not ((int(mGo),int(mLsp)) in dictXsec)):
    continue 
  if(not os.path.exists(inputpath + 'higgsCombinecount_masses_' + str(int(fname[3])) + '_' + str(int(fname[4])) + '_.AsymptoticLimits.mH120.root' )):
  	continue
  filein=TFile(inputpath + 'higgsCombinecount_masses_' + str(int(fname[3])) + '_' + str(int(fname[4])) + '_.AsymptoticLimits.mH120.root')
  t = filein.Get("limit")
  if not t:
  	continue
  t.GetEntry(1)
  limit_m1s = t.limit
  t.GetEntry(2) 
  limit_exp = t.limit
  t.GetEntry(3)
  limit_p1s = t.limit
  t.GetEntry(5)
  limit_obs = t.limit

  ExpUL= limit_exp #* float(dictXsec.get(mGo[m]))
  ExpULXSec= limit_exp* float(dictXsec[(int(mGo),int(mLsp))])
  ExpULSigmaUp = limit_p1s #*float(dictXsec.get(mGo[m]))
  ExpULSigmaDn = limit_m1s #*float(dictXsec.get(mGo[m]))
  ObsUL=limit_obs#*float(dictXsec.get(mGo[m]))
  shiftUp=1.0/(1-(float(dictXsecUnc.get((int(mGo),int(mLsp))))/100.));
  shiftDn=1.0/(1+(float(dictXsecUnc.get((int(mGo),int(mLsp))))/100.));
  ObsULDn=shiftDn*ObsUL
  ObsULUp=shiftUp*ObsUL
  if ExpUL<0.000001:continue
  MuScan.SetPoint(MuScan.GetN(), mGo, mLsp, ExpUL)
  MuScanSup.SetPoint(MuScanSup.GetN(), mGo, mLsp, ExpULSigmaUp)
  MuScanSdn.SetPoint(MuScanSdn.GetN(), mGo, mLsp, ExpULSigmaDn)
  MuScanObs.SetPoint( MuScanObs.GetN(), mGo, mLsp, ObsUL)
  MuScanObsSup.SetPoint( MuScanObsSup.GetN(), mGo, mLsp, ObsULUp)
  MuScanObsSdn.SetPoint( MuScanObsSdn.GetN(), mGo, mLsp, ObsULDn)
  MuScanXsec.SetPoint(MuScanXsec.GetN(), mGo, mLsp,ExpULXSec)
  print line
	
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
#hExplim=MuScan.GetHistogram()
#hExplimSup=MuScanSup.GetHistogram()
#hExplimSdn=MuScanSdn.GetHistogram()
hObslim=MuScanObs.GetHistogram()
#hObslimSup=MuScanObsSup.GetHistogram()
#hObslimSdn=MuScanObsSdn.GetHistogram()
#MassScan2D=MuScanXsec.GetHistogram()
print 'finish make histo'
c=TCanvas("c","",800,800);
MuScanObs.Draw("colz")
print 'finish draw'

ExpLim=TGraph()
ExpLim.SetName("ExpLim")

ExpLim= MuScan.GetContourList(1.);
ExpLimSup= MuScanSup.GetContourList(1.);
ExpLimSdn= MuScanSdn.GetContourList(1.);
ObsLim= MuScanObs.GetContourList(1.);
ObsLimSup= MuScanObsSup.GetContourList(1.);
ObsLimSdn= MuScanObsSdn.GetContourList(1.);

MuScanObs.GetZaxis().SetRangeUser(0.01,10)
c.SetLogz()
MuScanObs.Draw("colz")
#ObsLim.Draw("same")

fileOut.cd()
#	ExpLim.Write("ExpLim")
# ExpLimSup.Write("ExpLimSup")
# ExpLimSdn.Write("ExpLimSdn")
# ObsLim.Write("ObsLim")
# ObsLimSup.Write("ObsLimSup")
# ObsLimSdn.Write("ObsLimSdn")
#  MassScan2D.Write("MassScan2D")
MuScanObs.Write("MuScanObs")
MuScan.Write("MuScan")
#  
fileOut.Write()
c.Print("test.pdf")
#  fileOut.Close()
#  #if hlim is None: print "NONE"
#  #time.sleep(60)
