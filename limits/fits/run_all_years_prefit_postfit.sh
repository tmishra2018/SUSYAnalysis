#!/bin/bash

# Years and datacards
declare -A DATACARDS=(
  ["2016preVFP"]="/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_2016preVFP_unblind/counting_t5wg_1700_1100.txt"
  ["2016postVFP"]="/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_2016postVFP_unblind/counting_t5wg_1700_1100.txt"
  ["2017"]="/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_2017_unblind/counting_t5wg_1700_1100.txt"
  ["2018"]="/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_2018_unblind/counting_t5wg_1700_1100.txt"
  ["Run2"]="/eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_Run2_unblind/counting_t5wg_1700_1100.txt"
)

NBINS=36

for YEAR in "${!DATACARDS[@]}"; do

  DATACARD=${DATACARDS[$YEAR]}
  TAG="T5WG_${YEAR}"

  echo "=============================================="
  echo "Running FitDiagnostics for $YEAR"
  echo "=============================================="

  combine -M FitDiagnostics $DATACARD \
    --saveShapes \
    --saveWithUncertainties \
    -n _$TAG

  ROOTFILE="fitDiagnostics_$TAG.root"

  if [ ! -f "$ROOTFILE" ]; then
    echo "ERROR: $ROOTFILE not found!"
    continue
  fi

  # ROOT macro
  cat << EOF > extract_${YEAR}.C
#include <iomanip>
void extract_${YEAR}() {

  TFile *f = TFile::Open("$ROOTFILE");

  TDirectory *prefit  = (TDirectory*)f->Get("shapes_prefit");
  TDirectory *postfit = (TDirectory*)f->Get("shapes_fit_s");

  double pre_err[$NBINS];
  double post_err[$NBINS];

  std::cout << "\n=== $YEAR ===\n";
  std::cout << "bin  prefit_yield  prefit_err  postfit_yield  postfit_err\n";
  std::cout << "---------------------------------------------------------\n";

  for (int i = 1; i <= $NBINS; ++i) {

    TString bin = Form("bin%d", i);

    TH1 *h_pre  = (TH1*)prefit ->Get(bin + "/total_background");
    TH1 *h_post = (TH1*)postfit->Get(bin + "/total_background");

    if (!h_pre || !h_post) {
      pre_err[i-1]  = 0;
      post_err[i-1] = 0;
      continue;
    }

    double pre_y  = h_pre ->GetBinContent(1);
    double pre_e  = h_pre ->GetBinError(1);
    double post_y = h_post->GetBinContent(1);
    double post_e = h_post->GetBinError(1);

    pre_err[i-1]  = pre_e;
    post_err[i-1] = post_e;

    std::cout << std::fixed << std::setprecision(2)
              << "bin" << i << "  "
              << pre_y  << "  " << pre_e << "  "
              << post_y << "  " << post_e << std::endl;
  }

  // Print arrays
  std::cout << "\n// $YEAR\n";
  std::cout << "double totalerror_prefit_$YEAR[$NBINS] = {";
  for (int i = 0; i < $NBINS; ++i) {
    std::cout << pre_err[i];
    if (i != $NBINS-1) std::cout << ", ";
  }
  std::cout << "};\n";

  std::cout << "double totalerror_postfit_$YEAR[$NBINS] = {";
  for (int i = 0; i < $NBINS; ++i) {
    std::cout << post_err[i];
    if (i != $NBINS-1) std::cout << ", ";
  }
  std::cout << "};\n";
}
EOF

  root -l -b -q extract_${YEAR}.C
  rm extract_${YEAR}.C

done
