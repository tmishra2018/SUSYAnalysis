#!/bin/bash
# Derive electron miniIso and R9 SFs for one year from the converted TnP trees (makeTnPTree.exe output).
# usage: ./run_sf.sh <tag>        e.g. ./run_sf.sh 2017   (needs fileset_<tag>.json and pileup_<tag>.json here)
# EGAMMA_TNP_SRC: checkout of https://gitlab.cern.ch/mamaxout/egamma-tnp ; PYLIB: dir with numba_stats (pip install --target)
# Histograms, fits and plots go to TNP_WORKDIR (nobackup, large); the SF maps and tables go to results/<tag> here.
set -e
TAG=$1
HERE=$(cd $(dirname $0) && pwd)
EGAMMA_TNP_SRC=${EGAMMA_TNP_SRC:-/uscms_data/d3/tmishra/tnp_tools/egamma-tnp}
PYLIB=${PYLIB:-/uscms_data/d3/tmishra/tnp_tools/pylib}
TNP_WORKDIR=${TNP_WORKDIR:-/uscms_data/d3/tmishra/tnp_work}
IMAGE=/cvmfs/unpacked.cern.ch/registry.cern.ch/cms-egamma/egamma-tnp:lxplus-el9-latest
OUT=$HERE/results/$TAG
WORK=$TNP_WORKDIR/$TAG
mkdir -p $OUT $WORK

in_container() {
  apptainer exec -B /uscms -B /uscms/homes -B /uscms_data -B /tmp -B /eos/uscms -B /cvmfs $IMAGE bash -c \
    "export PYTHONPATH=$EGAMMA_TNP_SRC/src:$PYLIB MPLCONFIGDIR=$WORK/.mpl; cd $HERE; $1"
}

for M in miniIso R9; do
  F=passing${M^}
  W=$WORK/$M
  if ! ls $W/hists/data_*/get_nd_tnp_histograms_1/${F}_histos.pkl >/dev/null 2>&1; then
    in_container "python3 -m egamma_tnp.scripts.run_analysis --config config_$M.json --settings settings_$M.json \
        --fileset fileset_$TAG.json --binning binning.json --output $W/hists --executor threads"
  fi
  DATA=$(ls -d $W/hists/data_*/get_nd_tnp_histograms_1)
  MC=$(ls -d $W/hists/DY_*/get_nd_tnp_histograms_1)
  in_container "python3 -m egamma_tnp.scripts.prepare_histograms --source $DATA/${F}_histos.pkl --target $W/fitter_data.root \
      --binning $W/fitter_binning.pkl --axes el_eta el_pt && \
    python3 -m egamma_tnp.scripts.prepare_histograms --source $MC/${F}_histos.pkl --target $W/fitter_mc.root \
      --binning $W/fitter_binning.pkl --axes el_eta el_pt && \
    python3 fit_all.py --data $W/fitter_data.root --mc $W/fitter_mc.root --binning $W/fitter_binning.pkl --outdir $W/fits --workers 4"
  mkdir -p $OUT/$M
  cp $W/fits/sf_table.txt $W/fits/sf_summary.json $OUT/$M/
done

# TH2F maps (PyROOT from the analysis CMSSW release)
source /cvmfs/cms.cern.ch/cmsset_default.sh
(cd $HERE/../../.. && eval `scramv1 runtime -sh` 2>/dev/null
 for M in miniIso R9; do
   python3 $HERE/make_sf_hist.py $OUT/$M/sf_summary.json $OUT/SF_${M}_$TAG.root
   python3 $HERE/make_sf_hist.py $OUT/$M/sf_summary.json $OUT/SF_${M}_${TAG}_cutAndCount.root --cnc
   python3 $HERE/make_sf_hist.py $OUT/$M/sf_summary.json $OUT/SF_${M}_${TAG}_hybrid.root --hybrid
 done)
echo "done: $OUT (work files: $WORK)"
