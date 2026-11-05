import ROOT, sys
ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetOptTitle(0)
O=sys.argv[1]
R="/uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/SUSYAnalysis/test/ScaleFactor/tnp/results/%s/SF_%s_%s_hybrid.root"
OLD="/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/R9_SF.root"
years=["2016preVFP","2016postVFP","2017","2018"]; cols=[ROOT.kBlue,ROOT.kCyan+1,ROOT.kRed,ROOT.kGreen+2]
etab=[(-2.5,-2.0),(-2.0,-1.566),(-1.4442,-0.8),(-0.8,0.0),(0.0,0.8),(0.8,1.4442),(1.566,2.0),(2.0,2.5)]
keep=[]
def get(fn):
    f=ROOT.TFile.Open(fn); h=f.Get("EGamma_SF2D"); h.SetDirectory(0); keep.append(f); return h
def val(h,eta,pt):
    b=h.FindBin(eta,pt); return h.GetBinContent(b),h.GetBinError(b)
old=get(OLD)
txt=open(O+"/SF_comparison_tables.txt","w")
for sf in ("R9","miniIso"):
    new={y:get(R%(y,sf,y)) for y in years}
    ptedges=[new["2017"].GetYaxis().GetBinLowEdge(i) for i in range(1,new["2017"].GetNbinsY()+2)]
    c=ROOT.TCanvas("c","",1600,800); c.Divide(4,2)
    txt.write(f"\n===== {sf} SF  (old = {'mengleis 2016 R9_SF.root' if sf=='R9' else 'none applied (=1)'}) =====\n")
    txt.write(f"{'eta bin':18s} {'pT bin':10s} {'old':>14s} "+" ".join(f"{y:>16s}" for y in years)+"\n")
    for ip,(ea,eb) in enumerate(etab):
        pad=c.cd(ip+1); pad.SetLogx(); pad.SetGridy(); ec=0.5*(ea+eb)
        fr=pad.DrawFrame(25,0.85 if sf=="R9" else 0.92,500,1.08); fr.GetXaxis().SetTitle("p_{T} [GeV]"); fr.GetYaxis().SetTitle(sf+" SF")
        fr.GetXaxis().SetMoreLogLabels()
        leg=ROOT.TLegend(0.15,0.15,0.6,0.42); leg.SetBorderSize(0); leg.SetHeader(f"{ea} < #eta_{{SC}} < {eb}")
        go=ROOT.TGraphErrors()
        if sf=="R9":
            for i in range(old.GetNbinsY()):
                lo,hi=old.GetYaxis().GetBinLowEdge(i+1),old.GetYaxis().GetBinUpEdge(i+1)
                if hi<=25: continue
                v,e=val(old,ec,0.5*(max(lo,25)+hi)); n=go.GetN(); go.SetPoint(n,0.5*(max(lo,25)+hi),v); go.SetPointError(n,0.5*(hi-max(lo,25)),e)
        else:
            go.SetPoint(0,262.5,1.0); go.SetPointError(0,237.5,0)
        go.SetMarkerStyle(20); go.SetLineWidth(2); go.Draw("PE same"); leg.AddEntry(go,"old (2016 file)" if sf=="R9" else "old: not applied","lp"); keep.append(go)
        for iy,y in enumerate(years):
            g=ROOT.TGraphErrors(); h=new[y]
            for i in range(len(ptedges)-1):
                lo,hi=ptedges[i],ptedges[i+1]; pc=0.5*(lo+hi)*(1+0.03*(iy-1.5))
                v,e=val(h,ec,0.5*(lo+hi)); n=g.GetN(); g.SetPoint(n,pc,v); g.SetPointError(n,0,e)
            g.SetLineColor(cols[iy]); g.SetMarkerColor(cols[iy]); g.SetMarkerStyle(24+iy); g.SetLineWidth(2); g.Draw("PE same"); leg.AddEntry(g,"new "+y,"lp"); keep.append(g)
        leg.Draw(); keep.append(leg)
        for i in range(len(ptedges)-1):
            lo,hi=ptedges[i],ptedges[i+1]; pm=0.5*(lo+hi)
            if sf=="R9":
                ov=val(old,ec,min(pm,199))
                os=f"{ov[0]:.3f}±{ov[1]:.3f}"
            else: os="1 (none)"
            txt.write(f"[{ea:7.3f},{eb:7.3f}] [{lo:3.0f},{hi:3.0f}] {os:>14s} "+" ".join(f"{val(new[y],ec,pm)[0]:.3f}±{val(new[y],ec,pm)[1]:.3f}".rjust(16) for y in years)+"\n")
    c.SaveAs(f"{O}/SF_{sf}_old_vs_new.png"); c.SaveAs(f"{O}/SF_{sf}_old_vs_new.pdf")
txt.close()
