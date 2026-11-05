#!/bin/bash
# =============================================================================
set -e

CARDDIR="/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards"
OUTDIR="postfit_unblind_T5WG"
MASS="1700_1100"
MASS_TAG="1700/1100 GeV"
TOPOLOGY="T5WG"
NBINS=36

declare -A DATACARDS=(
  ["2016preVFP"]="${CARDDIR}/cards_2016preVFP_unblind/counting_t5wg_${MASS}.txt"
  ["2016postVFP"]="${CARDDIR}/cards_2016postVFP_unblind/counting_t5wg_${MASS}.txt"
  ["2017"]="${CARDDIR}/cards_2017_unblind/counting_t5wg_${MASS}.txt"
  ["2018"]="${CARDDIR}/cards_2018_unblind/counting_t5wg_${MASS}.txt"
  ["Run2"]="${CARDDIR}/cards_Run2_unblind/counting_t5wg_${MASS}.txt"
)

# Luminosities per year (for plot label)
declare -A LUMIS=(
  ["2016preVFP"]="19.5"
  ["2016postVFP"]="16.8"
  ["2017"]="41.5"
  ["2018"]="59.8"
  ["Run2"]="137.6"
)

YEAR_ORDER="2016preVFP 2016postVFP 2017 2018 Run2"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"
echo "Working in: $(pwd)"

# ---------------------------------------------------------------------------
# STEP 1: FitDiagnostics
# Per-year: observed (--expectSignal 0) + Asimov (-t -1 --expectSignal 0)
# Run2:     observed only (combined card, jointly constrained nuisances)
#
# NOTE on postfit_exp:
#   Asimov with --expectSignal 0 generates data from the prefit model,
#   so shapes_fit_b ~ shapes_prefit (nuisances at nominal, no pulls).
#   It is still useful to show alongside postfit_obs because it demonstrates
#   the reduction in uncertainty that real data provides.
#   Relabel it clearly as "B-only expected (prefit model)" not "postfit exp".
# ---------------------------------------------------------------------------
for YEAR in ${YEAR_ORDER}; do

  DATACARD=${DATACARDS[$YEAR]}
  TAG="${TOPOLOGY}_${YEAR}"

  echo ""
  echo "=============================================="
  echo " FitDiagnostics: ${YEAR}"
  echo "=============================================="

  echo "  [obs] Running combine (B-only, observed) ..."
  combine -M FitDiagnostics "${DATACARD}" \
    --saveShapes \
    --saveWithUncertainties \
    --saveNorm \
    --robustFit 1 \
    -n _${TAG}_obs \
    2>&1 | grep -E "Best fit|ERROR|FitDiag"

  OBSFILE="fitDiagnostics_${TAG}_obs.root"
  [ ! -f "${OBSFILE}" ] && { echo "ERROR: ${OBSFILE} not found!"; continue; }
  echo "  [obs] Done -> ${OBSFILE}"

  # Asimov only for per-year (not Run2 combined)
  if [ "${YEAR}" != "Run2" ]; then
    echo "  [exp] Running combine (B-only, Asimov) ..."
    combine -M FitDiagnostics "${DATACARD}" \
      --expectSignal 0 \
      -t -1 \
      --saveShapes \
      --saveWithUncertainties \
      --saveNormalizations \
      --robustFit 1 \
      -n _${TAG}_exp \
      2>&1 | grep -E "Best fit|ERROR|FitDiag"

    EXPFILE="fitDiagnostics_${TAG}_exp.root"
    [ ! -f "${EXPFILE}" ] && { echo "ERROR: ${EXPFILE} not found!"; continue; }
    echo "  [exp] Done -> ${EXPFILE}"
  fi

done

# ---------------------------------------------------------------------------
# STEP 2: Extract NormPostFit txt files
# ---------------------------------------------------------------------------
cat << 'ROOTEOF' > extract_normpostfit.C
#include "TFile.h"
#include "TH1.h"
#include "TKey.h"
#include "TDirectory.h"
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

std::vector<TString> getBinDirs(TDirectory* d){
  std::vector<TString> v;
  if(!d) return v;
  TIter next(d->GetListOfKeys());
  TKey* key;
  while((key=(TKey*)next())){
    TObject* obj = d->Get(key->GetName());
    if(obj && obj->InheritsFrom("TDirectory"))
      v.push_back(TString(key->GetName()));
  }
  return v;
}

int binIndex(const TString& name){
  TString s = name;
  if(s.Contains("_")) s = s(s.Last('_')+1, s.Length());
  s.ReplaceAll("bin","");
  return s.IsDigit() ? s.Atoi() : 0;
}
bool cmpBinDir(const TString& a, const TString& b){
  return binIndex(a) < binIndex(b);
}

void extract_normpostfit(){
  const int NBINS = 36;
  std::vector<std::string> years_peryear = {
    "2016preVFP","2016postVFP","2017","2018"
  };
  std::vector<std::string> processes = {
    "total_background","elefakepho","jetfakepho","qcdfakelep","VGamma","ttGam_VVGam"
  };

  for(auto& year : years_peryear){
    TString obsfile = Form("fitDiagnostics_T5WG_%s_obs.root", year.c_str());
    TFile* fobs = TFile::Open(obsfile);
    if(!fobs || fobs->IsZombie()){
      printf("WARN: %s not found\n", obsfile.Data()); continue;
    }
    TDirectory* d_pre  = (TDirectory*)fobs->Get("shapes_prefit");
    TDirectory* d_post = (TDirectory*)fobs->Get("shapes_fit_b");
    if(!d_pre || !d_post){
      printf("ERROR: shapes dirs missing in %s\n", obsfile.Data());
      fobs->Close(); continue;
    }
    std::vector<TString> bindirs = getBinDirs(d_pre);
    std::sort(bindirs.begin(), bindirs.end(), cmpBinDir);
    printf("\n=== %s : %d bins ===\n", year.c_str(), (int)bindirs.size());

    for(auto& proc : processes){
      std::ofstream outf(Form("NormPostFit_%s_%s.txt", proc.c_str(), year.c_str()));
      outf << std::fixed << std::setprecision(4);
      for(int i=0; i<(int)bindirs.size() && i<NBINS; i++){
        TH1* hp = (TH1*)d_pre ->Get(bindirs[i]+"/"+proc.c_str());
        TH1* hf = (TH1*)d_post->Get(bindirs[i]+"/"+proc.c_str());
        double pre_y(0),pre_e(0),post_y(0),post_e(0);
        if(hp){ pre_y=hp->GetBinContent(1); pre_e=hp->GetBinError(1); }
        if(hf){ post_y=hf->GetBinContent(1); post_e=hf->GetBinError(1); }
        outf << std::setw(10) << pre_y  << "  " << std::setw(8) << pre_e  << "  "
             << std::setw(10) << post_y << "  " << std::setw(8) << post_e << "\n";
      }
      outf.close();
      printf("  Written: NormPostFit_%s_%s.txt\n", proc.c_str(), year.c_str());
    }

    printf("  Pull summary for %s:\n", year.c_str());
    for(int i=0; i<(int)bindirs.size() && i<NBINS; i++){
      TH1* hp=(TH1*)d_pre ->Get(bindirs[i]+"/total_background");
      TH1* hf=(TH1*)d_post->Get(bindirs[i]+"/total_background");
      if(!hp||!hf) continue;
      double pre=hp->GetBinContent(1), post=hf->GetBinContent(1);
      double ratio=(pre>0)?post/pre:-1;
      const char* flag=(ratio>0&&fabs(ratio-1.0)>0.3)?"  <-- LARGE PULL":"";
      printf("    bin%d  pre=%.2f post=%.2f ratio=%.3f%s\n",
             i+1, pre, post, ratio, flag);
    }
    fobs->Close();
  }

  // Run2: sum ch1..ch4 per bin from combined fitDiagnostics
  printf("\n=== Run2 (combined, 144 ch -> 36 bins) ===\n");
  TFile* frun2 = TFile::Open("fitDiagnostics_T5WG_Run2_obs.root");
  if(!frun2 || frun2->IsZombie()){
    printf("WARN: Run2 file not found\n"); return;
  }
  TDirectory* d_pre_r2  = (TDirectory*)frun2->Get("shapes_prefit");
  TDirectory* d_post_r2 = (TDirectory*)frun2->Get("shapes_fit_b");
  if(!d_pre_r2 || !d_post_r2){
    printf("ERROR: shapes dirs missing in Run2\n"); frun2->Close(); return;
  }
  std::vector<TString> all_dirs = getBinDirs(d_pre_r2);
  printf("  Total channels: %d\n", (int)all_dirs.size());

  for(auto& proc : processes){
    double pre_y[NBINS]={}, pre_e2[NBINS]={}, post_y[NBINS]={}, post_e2[NBINS]={};
    for(auto& d : all_dirs){
      int bidx = binIndex(d) - 1;
      if(bidx<0||bidx>=NBINS) continue;
      TH1* hp = (TH1*)d_pre_r2 ->Get(d+"/"+proc.c_str());
      TH1* hf = (TH1*)d_post_r2->Get(d+"/"+proc.c_str());
      if(hp){ pre_y[bidx]+=hp->GetBinContent(1); pre_e2[bidx]+=pow(hp->GetBinError(1),2); }
      if(hf){ post_y[bidx]+=hf->GetBinContent(1); post_e2[bidx]+=pow(hf->GetBinError(1),2); }
    }
    std::ofstream outf(Form("NormPostFit_%s_Run2.txt", proc.c_str()));
    outf << std::fixed << std::setprecision(4);
    for(int i=0; i<NBINS; i++)
      outf << std::setw(10) << pre_y[i]       << "  " << std::setw(8) << sqrt(pre_e2[i])  << "  "
           << std::setw(10) << post_y[i]       << "  " << std::setw(8) << sqrt(post_e2[i]) << "\n";
    outf.close();
    printf("  Written: NormPostFit_%s_Run2.txt\n", proc.c_str());
  }
  frun2->Close();
}
ROOTEOF

echo ""
echo "=== [Step 2] Extracting NormPostFit_*.txt ==="
root -l -b -q extract_normpostfit.C
rm -f extract_normpostfit.C

echo ""
echo "=== PART 1 DONE ==="
echo ""
echo "Outputs written to: $(pwd)"
echo "                                                  proc in {total_background, elefakepho,"
echo "                                                  jetfakepho, qcdfakelep, VGamma, ttGam_VVGam}"
echo "                                                  year in {2016preVFP, 2016postVFP, 2017, 2018, Run2}"
echo ""
echo "Next: run run_plot_eventct.sh (Part 2), pointing POSTFITDIR at this directory,"
echo "      to feed these into plot_eventct.C."
