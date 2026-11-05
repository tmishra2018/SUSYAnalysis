// Runs toys first..last of one template fit in a single ROOT session, instead of one ROOT process per toy
// (which spends ~4 s per toy on ROOT start-up). Same macro and arguments as the per-toy call
//   Fitfraction<ch>.C+(i,metlow,methigh,leplow,lephigh,4,RunYear,preVFP,htlow,hthigh,1,pholow,phohigh)
// so the output (one line per toy appended to VGamma_scalefactor_<ch>_<era>[_<bin>].txt, the toy-0 MINOS file)
// is the same; checked 2026-10-31 against one process per toy for all 32 fits (24 photon-pT + 8 inclusive, 1000 toys):
// identical nominal fits, 1 sigma SF errors within 1.05% of the SF. The macro leaks ~5 MB per call, so keep
// last-first <= ~100 per session.
//   root -b -q 'run_toys.C("eg",0,99,2018,0,0,100)'                    photon pT 0-100 GeV
//   root -b -q 'run_toys.C("eg",0,99,2018,0,-1,-1)'                    inclusive
//   root -b -q 'run_toys.C("eg",0,99,2018,0,-1,-1,40,70,0,1000,0,100)' HT 0-100 GeV
#include "TROOT.h"
#include "TString.h"

void run_toys(const char *ch, int first, int last, int RunYear, int preVFP, int pholow, int phohigh,
              int metlow = 40, int methigh = 70, int leplow = 0, int lephigh = 1000, int htlow = -1, int hthigh = -1){
	gROOT->ProcessLine(Form(".L Fitfraction%s.C+", ch));
	for(int i = first; i <= last; i++)
		gROOT->ProcessLine(Form("Fitfraction%s(%d,%d,%d,%d,%d,4,%d,%d,%d,%d,1,%d,%d)", ch, i, metlow, methigh, leplow, lephigh,
		                        RunYear, preVFP, htlow, hthigh, pholow, phohigh));
}
