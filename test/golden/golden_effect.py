import ROOT, json, sys, glob, os
ROOT.gROOT.SetBatch(True)
G="/uscms/homes/t/tmishra/work/CMSSW_14_1_0_pre4/src/SUSYAnalysis/test/golden/golden_%s.json"
ROOT.gInterpreter.Declare("""
#include <unordered_set>
std::unordered_set<long long> gmask;
bool inGolden(int run, int lumi){ return gmask.count((long long)run*1000000LL+lumi)>0; }
""")
loaded=None
def load(year):
    global loaded
    if loaded==year: return
    ROOT.gmask.clear()
    for r,rngs in json.load(open(G%year)).items():
        for a,b in rngs:
            for l in range(a,b+1): ROOT.gmask.insert(int(r)*1000000+l)
    loaded=year
def year_of(fn):
    for y in ("2016","2017","2018"):
        if y in os.path.basename(fn): return y
files=sys.argv[1:]
print(f"{'file':55s} {'tree':12s} {'total':>10s} {'removed':>9s} {'frac%':>6s}")
for f in files:
    y=year_of(f); load(y)
    tf=ROOT.TFile.Open(f)
    if not tf or tf.IsZombie(): print(f,"cannot open"); continue
    for k in tf.GetListOfKeys():
        if "TTree" not in k.GetClassName(): continue
        t=tf.Get(k.GetName())
        br=[b.GetName() for b in t.GetListOfBranches()]
        if "run" not in br or "lumis" not in br:
            print(f"{os.path.basename(f):55s} {k.GetName():12s} {t.GetEntries():>10d}   (no run/lumis)"); continue
        df=ROOT.RDataFrame(t)
        n=df.Count(); nb=df.Filter("!inGolden(run,lumis)").Count()
        n=n.GetValue(); nb=nb.GetValue()
        print(f"{os.path.basename(f):55s} {k.GetName():12s} {n:>10d} {nb:>9d} {100.*nb/max(n,1):6.2f}")
