#!/usr/bin/env python3
# Generator-level ISR pT of LO vs NLO V+gamma samples (UL18 NanoAODv9), one NanoAOD file per call.
# ISR pT = |vector sum| of GenJets (pT > 30, |eta| < 2.5) with dR > 0.4 to the prompt leptons and the photon.
# Selection: >= 1 prompt e/mu (pT > 25, |eta| < 2.5) and a prompt photon (pT > 35, |eta| < 1.4442),
# dR(photon, lepton) > 0.8. Histograms are filled with weight 1 (unsigned) and sign(genWeight) (signed).
#   python3 compare_VG_LOvsNLO_nano.py <input root:// file> <output.root>
import sys
import ROOT
ROOT.gROOT.SetBatch(True)

ROOT.gInterpreter.Declare("""
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"
using ROOT::VecOps::RVec;
RVec<bool> promptLep(const RVec<int>& id, const RVec<int>& st, const RVec<int>& fl, const RVec<float>& pt, const RVec<float>& eta){
  RVec<bool> r(id.size());
  for(size_t i=0;i<id.size();i++) r[i] = st[i]==1 && (abs(id[i])==11 || abs(id[i])==13) && (fl[i]&1) && pt[i]>25 && fabs(eta[i])<2.5;
  return r;
}
RVec<bool> promptPho(const RVec<int>& id, const RVec<int>& st, const RVec<int>& fl, const RVec<float>& pt, const RVec<float>& eta){
  RVec<bool> r(id.size());
  for(size_t i=0;i<id.size();i++) r[i] = st[i]==1 && id[i]==22 && (fl[i]&1) && pt[i]>35 && fabs(eta[i])<1.4442;
  return r;
}
float dR(float e1, float p1, float e2, float p2){ float dp = TVector2::Phi_mpi_pi(p1-p2); return sqrt((e1-e2)*(e1-e2)+dp*dp); }
// photon: leading prompt photon with dR > 0.8 to every selected lepton; returns index or -1
int pickPho(const RVec<bool>& isPho, const RVec<bool>& isLep, const RVec<float>& pt, const RVec<float>& eta, const RVec<float>& phi){
  int best=-1;
  for(size_t i=0;i<isPho.size();i++){ if(!isPho[i])continue;
    bool ok=true; for(size_t j=0;j<isLep.size();j++) if(isLep[j] && dR(eta[i],phi[i],eta[j],phi[j])<0.8) ok=false;
    if(ok && (best<0 || pt[i]>pt[best])) best=i; }
  return best;
}
float isrPt(const RVec<float>& jpt, const RVec<float>& jeta, const RVec<float>& jphi, const RVec<float>& jm,
            const RVec<bool>& isLep, int ipho, const RVec<float>& eta, const RVec<float>& phi){
  ROOT::Math::PtEtaPhiMVector s(0,0,0,0);
  for(size_t k=0;k<jpt.size();k++){
    if(jpt[k]<30 || fabs(jeta[k])>2.5)continue;
    bool clean = dR(jeta[k],jphi[k],eta[ipho],phi[ipho])>0.4;
    for(size_t j=0;j<isLep.size();j++) if(isLep[j] && dR(jeta[k],jphi[k],eta[j],phi[j])<0.4) clean=false;
    if(clean) s += ROOT::Math::PtEtaPhiMVector(jpt[k],jeta[k],jphi[k],jm[k]);
  }
  return s.Pt();
}
""")

inp, out = sys.argv[1], sys.argv[2]
df = ROOT.RDataFrame("Events", inp)
df = (df.Define("w", "genWeight < 0 ? -1.f : 1.f")
        .Define("isLep", "promptLep(GenPart_pdgId, GenPart_status, GenPart_statusFlags, GenPart_pt, GenPart_eta)")
        .Define("isPho", "promptPho(GenPart_pdgId, GenPart_status, GenPart_statusFlags, GenPart_pt, GenPart_eta)"))
cnt_all = df.Sum("w")
n_all = df.Count()
sel = (df.Filter("Sum(isLep) >= 1")
         .Define("ipho", "pickPho(isPho, isLep, GenPart_pt, GenPart_eta, GenPart_phi)")
         .Filter("ipho >= 0")
         .Define("isr", "isrPt(GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, isLep, ipho, GenPart_eta, GenPart_phi)")
         .Define("phopt", "GenPart_pt[ipho]")
         .Define("ht", "Sum(GenJet_pt[GenJet_pt > 30 && abs(GenJet_eta) < 2.5])")
         .Define("one", "1.f"))
edges = [0, 50, 100, 150, 200, 250, 300, 800]
import array
e = array.array('d', edges)
hs = []
for tag, wv in (("signed", "w"), ("unsigned", "one")):
    hs.append(sel.Histo1D(ROOT.RDF.TH1DModel("isr_" + tag, ";ISR p_{T} (GeV)", len(edges) - 1, e), "isr", wv))
    hs.append(sel.Histo1D(ROOT.RDF.TH1DModel("isrfine_" + tag, ";ISR p_{T} (GeV)", 40, 0, 800), "isr", wv))
    hs.append(sel.Histo1D(ROOT.RDF.TH1DModel("ht_" + tag, ";H_{T} (GeV)", 30, 0, 1500), "ht", wv))
    hs.append(sel.Histo1D(ROOT.RDF.TH1DModel("phopt_" + tag, ";photon p_{T} (GeV)", 30, 0, 900), "phopt", wv))
f = ROOT.TFile(out, "RECREATE")
for h in hs:
    h.GetValue().Write()
hc = ROOT.TH1D("h_count", "bin1: events, bin2: sum sign(genWeight)", 2, 0, 2)
hc.SetBinContent(1, n_all.GetValue()); hc.SetBinContent(2, cnt_all.GetValue()); hc.Write()
f.Close()
print("[DONE]", inp, n_all.GetValue())
