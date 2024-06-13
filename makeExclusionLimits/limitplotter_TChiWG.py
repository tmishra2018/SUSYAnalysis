#!/usr/bin/env python
from ROOT import *
from array import *
gROOT.SetBatch(True)
import sys, os
from os.path import isfile, join
import math

# ------------------------------
def SetHistoStyle( histo ) :
    histo.SetLineColor(2)
    histo.SetMarkerColor(1)
    histo.GetXaxis().SetTitleFont(42);
    histo.GetYaxis().SetTitleFont(42);
    histo.GetXaxis().SetTitleSize(0.05);
    histo.GetYaxis().SetTitleSize(0.05);
    histo.GetXaxis().SetTitleFont(42)
    histo.GetYaxis().SetTitleFont(42)
    histo.GetXaxis().SetLabelFont(42)
    histo.GetYaxis().SetLabelFont(42)
    histo.GetXaxis().SetTitleSize(0.055)
    histo.GetYaxis().SetTitleSize(0.055)
    histo.GetXaxis().SetTitleOffset(1.1)
    histo.GetYaxis().SetTitleOffset(1.3)
    histo.GetXaxis().SetLabelOffset(0.011)
    histo.GetYaxis().SetLabelOffset(0.011)
    histo.GetXaxis().SetLabelSize(0.045)
    histo.GetYaxis().SetLabelSize(0.045)
    histo.SetMarkerStyle(9)
    histo.SetMarkerColor(1)
    histo.SetFillColor(2)
    histo.SetFillStyle(1001)
    histo.SetTitle("")
    histo.GetYaxis().SetTitle("Count")
    histo.GetYaxis().CenterTitle()
    histo.GetXaxis().CenterTitle()
# ------------------------------

#print "\npython limitPlotter.py <path-to-dir-with-Combine-root-files>  <limitname>\n"

# Extract input arguments
# ------------------------------
inputfoldername=sys.argv[1]
limitname=sys.argv[2]


# Extract file names in target dir
# ------------------------------
filelist = [f for f in os.listdir(inputfoldername) if isfile(join(inputfoldername, f)) if ".root" in f]
#print filelist


# Signal cross sections
# ------------------------------
signalmasscrosssection = [ 
    (500 , 0.0463533),
    (525 , 0.0372636),
    (550 , 0.0301656),
    (575 , 0.0245798),
    (600 , 0.0201372),
    (625 , 0.0165706),
    (650 , 0.0137303),
    (675 , 0.0113975),
    (700 , 0.00951032),
    (725 , 0.0079595),
    (750 , 0.00669356),
    (775 , 0.00563562),
    (800 , 0.00475843),
    (825 , 0.00402646),
    (850 , 0.00342026),
    (875 , 0.00290547),
    (900 , 0.00249667),
    (925 , 0.00212907),
    (950 , 0.0018164),
    (975 , 0.00156893),
    (1000 , 0.00134352),
    (1025 , 0.00115949),
    (1050 , 0.000997903),
    (1075 , 0.00086504),
    (1100 , 0.000740372),
    (1125 , 0.000647288),
    (1150 , 0.000555594),
    (1175 , 0.000486863),
    (1200 , 0.000415851),
    (1225 , 0.000362455),
    (1250 , 0.000316975),
    (1275 , 0.000276522),
    (1300 , 0.000240739)
    ]

theorycrosssection = [ 
    (500 , 0.0463533),
    (525 , 0.0372636),
    (550 , 0.0301656),
    (575 , 0.0245798),
    (600 , 0.0201372),
    (625 , 0.0165706),
    (650 , 0.0137303),
    (675 , 0.0113975),
    (700 , 0.00951032),
    (725 , 0.0079595),
    (750 , 0.00669356),
    (775 , 0.00563562),
    (800 , 0.00475843),
    (825 , 0.00402646),
    (850 , 0.00342026),
    (875 , 0.00290547),
    (900 , 0.00249667),
    (925 , 0.00212907),
    (950 , 0.0018164),
    (975 , 0.00156893),
    (1000 , 0.00134352),
    (1025 , 0.00115949),
    (1050 , 0.000997903),
    (1075 , 0.00086504),
    (1100 , 0.000740372),
    (1125 , 0.000647288),
    (1150 , 0.000555594),
    (1175 , 0.000486863),
    (1200 , 0.000415851),
    (1225 , 0.000362455),
    (1250 , 0.000316975),
    (1275 , 0.000276522),
    (1300 , 0.000240739)
    ]
theorycrosssectionunc = [ 
    (500 , 0.00416461),
    (525 , 0.00346763),
    (550 , 0.00288065),
    (575 , 0.00240183),
    (600 , 0.00204438),
    (625 , 0.00170195),
    (650 , 0.00143519),
    (675 , 0.00121937),
    (700 , 0.00104092),
    (725 , 0.000885243),
    (750 , 0.000769988),
    (775 , 0.000654544),
    (800 , 0.000564061),
    (825 , 0.000478362),
    (850 , 0.000412315),
    (875 , 0.000366257),
    (900 , 0.000314019),
    (925 , 0.00026801),
    (950 , 0.000242682),
    (975 , 0.000217618),
    (1000 , 0.000175604),
    (1025 , 0.000155683),
    (1050 , 0.000133077),
    (1075 , 0.000117638),
    (1100 , 0.000107178),
    (1125 , 9.17503e-05),
    (1150 , 8.28113e-05),
    (1175 , 6.44736e-05),
    (1200 , 5.79252e-05),
    (1225 , 5.61888e-05),
    (1250 , 4.6491e-05),
    (1275 , 4.30197e-05),
    (1300 , 3.39561e-05)
    ]

theorymass = [ 
    (500 , 500),
    (525 , 525),
    (550 , 550),
    (575 , 575),
    (600 , 600),
    (625 , 625),
    (650 , 650),
    (675 , 675),
    (700 , 700),
    (725 , 725),
    (750 , 750),
    (775 , 775),
    (800 , 800),
    (825 , 825),
    (850 , 850),
    (875 , 875),
    (900 , 900),
    (925 , 925),
    (950 , 950),
    (975 , 975),
    (1000 , 1000),
    (1025 , 1025),
    (1050 , 1050),
    (1075 , 1075),
    (1100 , 1100),
    (1125 , 1125),
    (1150 , 1150),
    (1175 , 1175),
    (1200 , 1200),
    (1225 , 1225),
    (1250 , 1250),
    (1275 , 1275),
    (1300 , 1300)
    ]

# Extract r-values from Combine output root files
# ------------------------------
m2sigmaexparray=[]
m1sigmaexparray=[]
meanexparray=[]
p1sigmaexparray=[]
p2sigmaexparray=[]
obsarray=[]
massarray=[]
masserrorarray=[]
#print
for filename in filelist:
    #print  filename
    infile=TFile(inputfoldername+"/"+filename,"READ")
#    infile.ls()
    #print
    tr=infile.Get("limit")
#    tr.Show(0)
    for ie,event in enumerate(tr):
        #print event.quantileExpected
        #print round(event.quantileExpected,2)
        if ie==0 and round(event.quantileExpected,2) == 0.03 : m2sigmaexparray.append((round(10*event.mh,0),event.limit))
        if ie==1 and round(event.quantileExpected,2) == 0.16 : m1sigmaexparray.append((round(10*event.mh,0),event.limit))
        if ie==2 and round(event.quantileExpected,2) == 0.5  : meanexparray.append((round(10*event.mh,0),event.limit))
        if ie==3 and round(event.quantileExpected,2) == 0.84 : p1sigmaexparray.append((round(10*event.mh,0),event.limit))
        if ie==4 and round(event.quantileExpected,2) == 0.98 : p2sigmaexparray.append((round(10*event.mh,0),event.limit))
        if ie==5 and round(event.quantileExpected,2) == -1.0 : obsarray.append((round(10*event.mh,0),event.limit))
        if ie==5 and round(event.quantileExpected,2) == -1.0 : masserrorarray.append(0.0)
        #if ie==0 and round(event.quantileExpected,2) == 0.03 : m2sigmaexparray.append((round(10000*event.mh,0),event.limit))
        #if ie==1 and round(event.quantileExpected,2) == 0.16 : m1sigmaexparray.append((round(10000*event.mh,0),event.limit))
        #if ie==2 and round(event.quantileExpected,2) == 0.5  : meanexparray.append((round(10000*event.mh,0),event.limit))
        #if ie==3 and round(event.quantileExpected,2) == 0.84 : p1sigmaexparray.append((round(10000*event.mh,0),event.limit))
        #if ie==4 and round(event.quantileExpected,2) == 0.98 : p2sigmaexparray.append((round(10000*event.mh,0),event.limit))
        #if ie==5 and round(event.quantileExpected,2) == -1.0 : obsarray.append((round(10000*event.mh,0),event.limit))
        #if ie==5 and round(event.quantileExpected,2) == -1.0 : masserrorarray.append(0.0)
#print
#print
print ("got limits into array")
#for ii in m2sigmaexparray:
#    print ii

# Create r-value limit arrays
# ------------------------------
signalcrosssectionsorted = [ x[1] for x in sorted(signalmasscrosssection, key=lambda tup: tup[0]) ]
theorycrosssectionsorted = [ x[1] for x in sorted(theorycrosssection, key=lambda tup: tup[0]) ]
theorycrosssectionuncsorted = [ x[1] for x in sorted(theorycrosssectionunc, key=lambda tup: tup[0]) ]
theorymasssorted      = [ x[1] for x in sorted(theorymass,      key=lambda tup: tup[0]) ]
obsarraysorted        = [ x[1] for x in sorted(obsarray,        key=lambda tup: tup[0]) ]
m2sigmaexparraysorted = [ x[1] for x in sorted(m2sigmaexparray, key=lambda tup: tup[0]) ]
m1sigmaexparraysorted = [ x[1] for x in sorted(m1sigmaexparray, key=lambda tup: tup[0]) ]
meanexparraysorted    = [ x[1] for x in sorted(meanexparray,    key=lambda tup: tup[0]) ]
massarraysorted       = [ x[0] for x in sorted(meanexparray,    key=lambda tup: tup[0]) ]
massarraysortedV2     = [ x[0] for x in sorted(signalmasscrosssection,    key=lambda tup: tup[0]) ]
p1sigmaexparraysorted = [ x[1] for x in sorted(p1sigmaexparray, key=lambda tup: tup[0]) ]
p2sigmaexparraysorted = [ x[1] for x in sorted(p2sigmaexparray, key=lambda tup: tup[0]) ]
#print "**Before signal cross-section normalization**"
#print "signalcrosssectionsorted : ",signalcrosssectionsorted
#print "          obsarraysorted : ",obsarraysorted
#print "   m2sigmaexparraysorted : ",m2sigmaexparraysorted
#print "   m1sigmaexparraysorted : ",m1sigmaexparraysorted
#print "      meanexparraysorted : ",meanexparraysorted
#print "   p1sigmaexparraysorted : ",p1sigmaexparraysorted
#print "   p2sigmaexparraysorted : ",p2sigmaexparraysorted
#print "         massarraysorted : ",massarraysorted
#print "Th. XSec massarraysorted : ",massarraysortedV2
#print


for imass in massarraysorted:
    if imass not in massarraysortedV2 : 
        print ("\nERROR: input limit mass ",str(imass)," is not included in the signal cross-section array!")
        print ("available signal masses: ",massarraysortedV2)
        print ("..exiting!\n")
        exit(1)


# Create cross-section normalized limit arrays
# ------------------------------
obsarraysorted        = [a*b for a,b in zip(signalcrosssectionsorted,obsarraysorted)]
m2sigmaexparraysorted = [a*b for a,b in zip(signalcrosssectionsorted,m2sigmaexparraysorted)]
m1sigmaexparraysorted = [a*b for a,b in zip(signalcrosssectionsorted,m1sigmaexparraysorted)]
meanexparraysorted    = [a*b for a,b in zip(signalcrosssectionsorted,meanexparraysorted)]
p1sigmaexparraysorted = [a*b for a,b in zip(signalcrosssectionsorted,p1sigmaexparraysorted)]
p2sigmaexparraysorted = [a*b for a,b in zip(signalcrosssectionsorted,p2sigmaexparraysorted)]
#print "**After signal cross-section normalization**"
#print "signalcrosssectionsorted : ",signalcrosssectionsorted
#print "          obsarraysorted : ",obsarraysorted
#print "   m2sigmaexparraysorted : ",m2sigmaexparraysorted
#print "   m1sigmaexparraysorted : ",m1sigmaexparraysorted
#print "      meanexparraysorted : ",meanexparraysorted
#print "   p1sigmaexparraysorted : ",p1sigmaexparraysorted
#print "   p2sigmaexparraysorted : ",p2sigmaexparraysorted
#print "         massarraysorted : ",massarraysorted
#print

# Create bands for +/- 1 and 2 sigma & theory xsec uncertainty
# ------------------------------
m2bandmeanexpsorted=[]
m1bandmeanexpsorted=[]
p1bandmeanexpsorted=[]
p2bandmeanexpsorted=[]
expLimit = 0
obsLimit = 0
for im,mass in enumerate(massarraysorted) :
    m2bandmeanexpsorted.append( meanexparraysorted[im]-m2sigmaexparraysorted[im])
    m1bandmeanexpsorted.append( meanexparraysorted[im]-m1sigmaexparraysorted[im])
    p1bandmeanexpsorted.append(-meanexparraysorted[im]+p1sigmaexparraysorted[im])
    p2bandmeanexpsorted.append(-meanexparraysorted[im]+p2sigmaexparraysorted[im])
    
    # Calculate limit (obs. & exp.) automatically
    # ------------------------------
    if im+1 < len(massarraysorted) :
        # log y = a*x+b --> because of using log scale !!!
        #aExp = (TMath.Log(meanexparraysorted[im])-TMath.Log(meanexparraysorted[im+1]))/(massarraysorted[im]-massarraysorted[im+1])
        aObs = (TMath.Log(obsarraysorted[im])-TMath.Log(obsarraysorted[im+1]))/(massarraysorted[im]-massarraysorted[im+1])
        aTheo = (TMath.Log(signalcrosssectionsorted[im])-TMath.Log(signalcrosssectionsorted[im+1]))/(massarraysorted[im]-massarraysorted[im+1])
        bExp = TMath.Log(meanexparraysorted[im])-aExp*massarraysorted[im]
        bObs = TMath.Log(obsarraysorted[im])-aObs*massarraysorted[im]
        bTheo = TMath.Log(signalcrosssectionsorted[im])-aTheo*massarraysorted[im]
        
        canTF=TCanvas("canTF","canTF",0,0,700,700); canTF.cd()
        lineExp = TF1("lineExp","TMath::Exp([0]*x+[1])-TMath::Exp([2]*x+[3])",massarraysorted[im],massarraysorted[im+1])
        lineExp.SetParameter(0,aExp)
        lineExp.SetParameter(1,bExp)
        lineExp.SetParameter(2,aTheo)
        lineExp.SetParameter(3,bTheo)
        lineExp.Draw()
        lineExp.SetNpx(int(massarraysorted[len(massarraysorted)-1]-massarraysorted[0]))
        if lineExp.GetX(0,massarraysorted[0],massarraysorted[len(massarraysorted)-1]) <= massarraysorted[im+1] and lineExp.GetX(0,massarraysorted[0],massarraysorted[len(massarraysorted)-1]) >= massarraysorted[im] : expLimit = lineExp.GetX(0,massarraysorted[0],massarraysorted[len(massarraysorted)-1])
        
        lineObs = TF1("lineObs","TMath::Exp([0]*x+[1])-TMath::Exp([2]*x+[3])",massarraysorted[im],massarraysorted[im+1])
        lineObs.SetParameter(0,aObs)
        lineObs.SetParameter(1,bObs)
        lineObs.SetParameter(2,aTheo)
        lineObs.SetParameter(3,bTheo)
        lineObs.Draw()
        lineObs.SetNpx(int(massarraysorted[len(massarraysorted)-1]-massarraysorted[0]))
        if lineObs.GetX(0,massarraysorted[0],massarraysorted[len(massarraysorted)-1]) <= massarraysorted[im+1] and lineObs.GetX(0,massarraysorted[0],massarraysorted[len(massarraysorted)-1]) >= massarraysorted[im] : obsLimit = lineObs.GetX(0,massarraysorted[0],massarraysorted[len(massarraysorted)-1])
        canTF.Close()
        
print ("next: create canvas")        
    

# Create output canvas
# ------------------------------
can=TCanvas("can1","can1",0,0,700,700); can.cd()
gStyle.SetOptStat(0)
can.SetLeftMargin(0.15); can.SetRightMargin(0.10); can.SetTopMargin(0.07); can.SetBottomMargin(0.13)
can.SetTickx(1); can.SetTicky(1); can.SetLogy(0)


# Create limit curves
# ------------------------------
theorycurve             = TGraphErrors ( len(theorycrosssectionsorted), array('d',theorymasssorted), array('d',theorycrosssectionsorted) )
for im,mass in enumerate(theorymasssorted) :
    theorycurve.SetPointError(im,0,theorycrosssectionuncsorted[im])
    
observedlimit           = TGraph( len(massarraysorted), array('d',massarraysorted), array('d',obsarraysorted) )
m2sigmaexpectedlimit    = TGraph( len(massarraysorted), array('d',massarraysorted), array('d',m2sigmaexparraysorted) )
m1sigmaexpectedlimit    = TGraph( len(massarraysorted), array('d',massarraysorted), array('d',m1sigmaexparraysorted) )
expectedlimit           = TGraph( len(massarraysorted), array('d',massarraysorted), array('d',meanexparraysorted)    )
expected1sigmalimitarea = TGraphAsymmErrors( len(massarraysorted), array('d',massarraysorted), array('d',meanexparraysorted), array('d',masserrorarray), array('d',masserrorarray), array('d',m1bandmeanexpsorted), array('d',p1bandmeanexpsorted)   )
expected2sigmalimitarea = TGraphAsymmErrors( len(massarraysorted), array('d',massarraysorted), array('d',meanexparraysorted), array('d',masserrorarray), array('d',masserrorarray), array('d',m2bandmeanexpsorted), array('d',p2bandmeanexpsorted)   )
p1sigmaexpectedlimit    = TGraph( len(massarraysorted), array('d',massarraysorted), array('d',p1sigmaexparraysorted) )
p2sigmaexpectedlimit    = TGraph( len(massarraysorted), array('d',massarraysorted), array('d',p2sigmaexparraysorted) )

# Create canvas
# ------------------------------
can.cd(); can.SetLogy(1)
can.SetLeftMargin(0.16)
can.SetRightMargin(0.08)
can.SetTopMargin(0.07)
can.SetBottomMargin(0.13)

# Create legends
# ------------------------------
legend2 = TLegend(0.52,0.84,0.83,0.9)
legend2.SetFillColor(0); legend2.SetLineColor(0); legend2.SetFillColor(0); legend2.SetFillStyle(0)
legend2.SetTextSize(0.035); legend2.SetTextFont(42); legend2.SetTextAlign(12)
legend2.AddEntry(theorycurve,"#sigma_{EWKino} with unc.","FL")
legend = TLegend(0.52,0.57,0.83,0.82)
legend.SetFillColor(0); legend.SetLineColor(0); legend.SetFillColor(0); legend.SetFillStyle(0)
legend.SetTextSize(0.035); legend.SetTextFont(42); legend.SetTextAlign(12)
legend.SetHeader("95% CL upper limits")
legend.AddEntry(observedlimit,"Observed","PL")
legend.AddEntry(expectedlimit,"Expected","PL")
legend.AddEntry(expected1sigmalimitarea,"Expected #pm 1 s.d.","F")
legend.AddEntry(expected2sigmalimitarea,"Expected #pm 2 s.d.","F")


# Create labels
# ------------------------------
cmslabel      = TLatex(); cmslabel.SetTextFont(61);      cmslabel.SetTextAlign(11);      cmslabel.SetNDC();      cmslabel.SetTextSize(0.042)
cmsextralabel = TLatex(); cmsextralabel.SetTextFont(52); cmsextralabel.SetTextAlign(11); cmsextralabel.SetNDC(); cmsextralabel.SetTextSize(0.042)
lumilabel     = TLatex(); lumilabel.SetTextFont(42);     lumilabel.SetTextAlign(31);     lumilabel.SetNDC();     lumilabel.SetTextSize(0.042)


# Draw & save output 
# ------------------------------
hFrame = TH1D("","",1,500,1260)
hFrame.Draw("HIST")
SetHistoStyle(hFrame)
hFrame.GetYaxis().SetRangeUser(0.0001,2.)
hFrame.GetYaxis().SetTitle("#sigma (pb)")
hFrame.GetXaxis().SetTitle("EWKino Mass (GeV)")
expectedlimit.SetLineStyle(3)
expectedlimit.SetLineWidth(2)
observedlimit.SetMarkerStyle(9)
observedlimit.SetLineStyle(7)
observedlimit.SetLineWidth(3)
expected1sigmalimitarea.SetFillColor(kGreen+1)
expected1sigmalimitarea.SetLineColor(kGreen+1)
expected1sigmalimitarea.SetFillStyle(1001)
expected2sigmalimitarea.SetFillColor(kOrange)
expected2sigmalimitarea.SetLineColor(kOrange)
expected2sigmalimitarea.SetFillStyle(1001)
theorycurve.SetFillColor(kRed-9)
theorycurve.SetFillStyle(1001)
theorycurve.SetLineColor(kRed+1)
theorycurve.SetLineWidth(2)
expected2sigmalimitarea.Draw("3same")
expected1sigmalimitarea.Draw("3same")
#m1sigmaexpectedlimit.Draw("Lsame")
#p1sigmaexpectedlimit.Draw("Lsame")
#m2sigmaexpectedlimit.Draw("Lsame")
#p2sigmaexpectedlimit.Draw("Lsame")
theorycurve.Draw("3same")
theorycurve.Draw("LXsame")
expectedlimit.Draw("PLsame")
observedlimit.Draw("PLsame")
legend2.Draw()
legend.Draw()
cmslabel.DrawLatex(0.16,0.942,"CMS")
cmsextralabel.DrawLatex(0.26,0.942,"Preliminary")
lumilabel.DrawLatex(0.92,0.942,"35.9 fb^{-1} (13 TeV)")
gPad.RedrawAxis()
outFile = TFile("limitPlotter_TChiWG.root","recreate")
theorycurve.Write()
expectedlimit.Write()
observedlimit.Write()
can.SaveAs("limitPlotter_"+limitname+".pdf")

#print 
if expLimit <= massarraysorted[0] or expLimit >= massarraysorted[len(massarraysorted)-1] : print ("No expected limit calculation possible!")
elif expLimit != 0 : 
    print ("Expected Limit: %.0f GeV" % expLimit)
print ("Expected Limit: %.0f GeV" % expLimit)
if obsLimit <= massarraysorted[0] or obsLimit >= massarraysorted[len(massarraysorted)-1] : print ("No observed limit calculation possible!")
elif obsLimit != 0 : 
    print ("Observed Limit: %.0f GeV" % obsLimit)
print ("Observed Limit: %.0f GeV" % obsLimit)
