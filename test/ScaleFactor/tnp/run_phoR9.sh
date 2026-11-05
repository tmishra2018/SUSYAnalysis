#!/bin/bash
# Derive the photon R9 SF (R9 > 0.5, EB and EE) for one year from the photon-probe TnP trees
# (convert_year.sh <tag> pho). Denominator: probe passing the analysis photon selection (Et > 35, Fall17v2 loose).
# usage: ./run_phoR9.sh <tag>      e.g. ./run_phoR9.sh 2017   (needs fileset_<tag>.json and pileup_<tag>.json here)
# Same tools and layout as run_sf.sh; outputs go to results/<tag>/phoR9 and results/<tag>/SF_phoR9_<tag>*.root
set -e
TAG=$1
HERE=$(cd $(dirname $0) && pwd)
EGAMMA_TNP_SRC=${EGAMMA_TNP_SRC:-/uscms_data/d3/tmishra/tnp_tools/egamma-tnp}
PYLIB=${PYLIB:-/uscms_data/d3/tmishra/tnp_tools/pylib}
TNP_WORKDIR=${TNP_WORKDIR:-/uscms_data/d3/tmishra/tnp_work}
IMAGE=/cvmfs/unpacked.cern.ch/registry.cern.ch/cms-egamma/egamma-tnp:lxplus-el9-latest
OUT=$HERE/results/$TAG
W=$TNP_WORKDIR/$TAG/phoR9
mkdir -p $OUT/phoR9 $W

# electron fileset -> photon-probe trees
sed -e 's#/tnp_#/tnpPho_#' -e 's#tnpEleIDs/fitter_tree#tnpPhoIDs/fitter_tree#' $HERE/fileset_$TAG.json > $W/fileset_pho_$TAG.json

in_container() {
  apptainer exec -B /uscms -B /uscms/homes -B /uscms_data -B /tmp -B /eos/uscms -B /cvmfs $IMAGE bash -c \
    "export PYTHONPATH=$EGAMMA_TNP_SRC/src:$PYLIB MPLCONFIGDIR=$W/.mpl; cd $HERE; $1"
}

if ! ls $W/hists/data_*/get_nd_tnp_histograms_1/passingR9_histos.pkl >/dev/null 2>&1; then
  in_container "python3 -m egamma_tnp.scripts.run_analysis --config config_phoR9.json --settings settings_phoR9.json \
      --fileset $W/fileset_pho_$TAG.json --binning binning_pho.json --output $W/hists --executor threads"
fi
DATA=$(ls -d $W/hists/data_*/get_nd_tnp_histograms_1)
MC=$(ls -d $W/hists/DY_*/get_nd_tnp_histograms_1)
in_container "python3 -m egamma_tnp.scripts.prepare_histograms --source $DATA/passingR9_histos.pkl --target $W/fitter_data.root \
    --binning $W/fitter_binning.pkl --axes ph_eta ph_et && \
  python3 -m egamma_tnp.scripts.prepare_histograms --source $MC/passingR9_histos.pkl --target $W/fitter_mc.root \
    --binning $W/fitter_binning.pkl --axes ph_eta ph_et && \
  python3 fit_all.py --data $W/fitter_data.root --mc $W/fitter_mc.root --binning $W/fitter_binning.pkl --outdir $W/fits \
    --workers 4 --axes ph_eta ph_et"
cp $W/fits/sf_table.txt $W/fits/sf_summary.json $OUT/phoR9/

# TH2F maps (PyROOT from the analysis CMSSW release)
source /cvmfs/cms.cern.ch/cmsset_default.sh
(cd $HERE/../../.. && eval `scramv1 runtime -sh` 2>/dev/null
 python3 $HERE/make_sf_hist.py $OUT/phoR9/sf_summary.json $OUT/SF_phoR9_$TAG.root
 python3 $HERE/make_sf_hist.py $OUT/phoR9/sf_summary.json $OUT/SF_phoR9_${TAG}_cutAndCount.root --cnc
 python3 $HERE/make_sf_hist.py $OUT/phoR9/sf_summary.json $OUT/SF_phoR9_${TAG}_hybrid.root --hybrid)
echo "done: $OUT (work files: $W)"
