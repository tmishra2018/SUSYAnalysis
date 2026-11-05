#!/usr/bin/env python
"""
diagnose_postfit_vs_limit.py
SUS-24-011  T5WG  --  pin down the postfit-pull vs. limit-gap "inconsistency"

WHY THIS SCRIPT EXISTS
-----------------------
run_postfit.sh produces ONE Run2-combined B-only postfit picture
(SRbins_T5WG_Run2_postfit_obs.pdf) for a SINGLE benchmark datacard
(MASS="1700_1100"). runLimit.sh + makeT5WGLimit.py instead scan MANY
(gluino, chargino) mass points and report, for each one, an
observed-vs-expected limit gap.

Those are not the same statistical object:
  - the postfit plot is a background-ONLY fit (r fixed to 0): any data
    excess can only be explained via nuisance-parameter pulls in
    specific SR bins.
  - the limit gap comes from a PROFILED fit (r floating): an excess
    only pulls the observed limit away from expected at mass points
    whose signal acceptance actually lives in the same SR bins where
    the excess sits.

So a real, single background mismodeling feature can show up as
"mostly absorbed" in the aggregate postfit plot while still producing
a large observed>expected gap at some mass points and none at others.
This script tests that directly: it builds a per-SR-bin data/postfit
pull vector from the SAME fitDiagnostics ROOT file run_postfit.sh
already produces, builds each mass point's per-bin signal-acceptance
shape from its combine datacard (the same ones runLimit.sh feeds to
AsymptoticLimits), and correlates a signal-acceptance-weighted "pull
overlap" against that mass point's actual observed/expected limit
ratio.

A strong positive correlation across the mass plane confirms: the
apparent inconsistency is explained by different mass points being
sensitive to different SR bins (hypothesis #2). A weak/no correlation
instead points at the other candidates discussed alongside this
script (asymptotic-formula breakdown at low stats, signal
contamination in the fake-rate control regions, or differing
combine fit configurations between FitDiagnostics and
AsymptoticLimits).

INPUTS (all already produced by your existing pipeline)
---------------------------------------------------------
  --postfit-root   fitDiagnostics_T5WG_Run2_obs.root
                    (written by run_postfit.sh, Step 1, Run2 entry)
  --card-dir        cards_Run2_unblind directory
                    (same CARDDIR/cards_Run2_unblind used by both
                     run_postfit.sh and runLimit.sh)
  --limit-dir        limits_Run2_unblind directory
                    (where runLimit.sh xrdcp's limit_<g>_<c>.txt)
  --model            datacard/file prefix, e.g. t5wg
  --nbins            number of true SR bins (36 for T5WG, per
                     run_postfit.sh's NBINS)

OUTPUTS
--------
  <out-csv>   one row per mass point: gluino, chargino, pull_overlap,
              observed limit, expected(50%) limit, obs/exp ratio
  <out-pdf>   scatter of pull_overlap vs. obs/exp ratio (ROOT canvas)
  stdout      per-bin pull table for the Run2 combined postfit,
              Pearson correlation coefficient + plain-English verdict

USAGE
------
  python diagnose_postfit_vs_limit.py \\
      --postfit-root postfit_Bonly_T5WG/fitDiagnostics_T5WG_Run2_obs.root \\
      --card-dir     /eos/uscms/store/user/tmishra/CombinedLimit/T5WG/cards/cards_Run2_unblind \\
      --limit-dir    /eos/uscms/store/user/tmishra/CombinedLimit/T5WG/limits/limits_Run2_unblind \\
      --model        t5wg \\
      --out-csv      pull_vs_limit_gap.csv \\
      --out-pdf      pull_vs_limit_gap.pdf

NOTE: this has not been run against real fitDiagnostics/datacard files
(none were available in the environment this was written in) - the
ROOT-shapes-directory parsing and datacard parsing below follow the
exact conventions already used by run_postfit.sh's extract_normpostfit.C
and by makeT5WGLimit.py's limit-file parsing, but please sanity-check
the first few printed per-bin rows against a NormPostFit_*_Run2.txt
file you already trust before relying on the correlation number.
"""

import argparse
import glob
import math
import os
import re
import sys

import ROOT

ROOT.gROOT.SetBatch(ROOT.kTRUE)


# ---------------------------------------------------------------------------
# Shared bin-naming convention, mirrored from run_postfit.sh's
# extract_normpostfit.C (getBinDirs / binIndex) so this script buckets
# channels into the same 36 true SR bins as the existing pipeline.
# ---------------------------------------------------------------------------
def bin_index_from_name(name):
    s = name
    if "_" in s:
        s = s.rsplit("_", 1)[-1]
    s = s.replace("bin", "")
    return int(s) if s.isdigit() else None


def get_bin_dirs(tdirectory):
    out = []
    if not tdirectory:
        return out
    for key in tdirectory.GetListOfKeys():
        obj = tdirectory.Get(key.GetName())
        if obj and obj.InheritsFrom("TDirectory"):
            out.append(key.GetName())
    return out


def get_data_obs_yield(chdir):
    """The data graph in combine's saved per-channel shape directories
    is keyed "data" (not "data_obs" - that name is only used at the
    top-level Shape_ObserveData style objects in some other combine
    outputs). It's usually a TGraphAsymmErrors (Poisson points),
    occasionally a TH1. Handle both; fall back to sqrt(N) for the
    statistical error if the graph's own errors are degenerate (can
    happen for empty bins).

    NOTE: TDirectory::Get() returns a null-but-not-None ROOT proxy
    when the key is missing -- "is None" never catches that in
    PyROOT, only "not obj" does."""
    obj = chdir.Get("data")
    if not obj:
        return None, None
    if obj.InheritsFrom("TH1"):
        y = obj.GetBinContent(1)
        e = obj.GetBinError(1)
        return y, (e if e > 0 else math.sqrt(max(y, 0.0)))
    if obj.InheritsFrom("TGraph"):
        n = obj.GetN()
        if n == 0:
            return 0.0, 0.0
        y = float(obj.GetPointY(0))
        try:
            eh = obj.GetErrorYhigh(0)
            el = obj.GetErrorYlow(0)
            e = 0.5 * (eh + el) if (eh > 0 or el > 0) else math.sqrt(max(y, 0.0))
        except Exception:
            e = math.sqrt(max(y, 0.0))
        return y, e
    return None, None


def extract_bin_pulls(fitdiag_path, nbins, postfit_dirname="shapes_fit_b"):
    """Returns dict: bin_index -> dict(data, data_err, pre, pre_err,
    post, post_err, pull). Sums all channels (across years) sharing
    the same trailing bin index, exactly like run_postfit.sh's Run2
    aggregation loop. postfit_dirname lets this be reused for either
    the b-only postfit (shapes_fit_b) or the s+b postfit (shapes_fit_s,
    which combine -M FitDiagnostics already computes in the same run)."""
    f = ROOT.TFile.Open(fitdiag_path)
    if not f or f.IsZombie():
        sys.exit(f"ERROR: cannot open {fitdiag_path}")

    d_pre = f.Get("shapes_prefit")
    d_post = f.Get(postfit_dirname)
    if not d_pre or not d_post:
        sys.exit(f"ERROR: shapes_prefit/{postfit_dirname} missing in {fitdiag_path} "
                  f"(shapes_fit_s only exists if the datacard used for this "
                  f"FitDiagnostics run contained a signal process)")

    all_dirs = get_bin_dirs(d_pre)

    data_y = [0.0] * (nbins + 1)
    data_e2 = [0.0] * (nbins + 1)
    pre_y = [0.0] * (nbins + 1)
    pre_e2 = [0.0] * (nbins + 1)
    post_y = [0.0] * (nbins + 1)
    post_e2 = [0.0] * (nbins + 1)
    seen = [False] * (nbins + 1)

    for chname in all_dirs:
        bidx = bin_index_from_name(chname)
        if bidx is None or bidx < 1 or bidx > nbins:
            continue
        seen[bidx] = True

        chdir_pre = d_pre.Get(chname)
        if chdir_pre:
            dy, de = get_data_obs_yield(chdir_pre)
            if dy is not None:
                data_y[bidx] += dy
                data_e2[bidx] += de * de

        hp = chdir_pre.Get("total_background") if chdir_pre else None
        if hp:
            pre_y[bidx] += hp.GetBinContent(1)
            pre_e2[bidx] += hp.GetBinError(1) ** 2

        chdir_post = d_post.Get(chname)
        hf = chdir_post.Get("total_background") if chdir_post else None
        if hf:
            post_y[bidx] += hf.GetBinContent(1)
            post_e2[bidx] += hf.GetBinError(1) ** 2

    result = {}
    for b in range(1, nbins + 1):
        if not seen[b]:
            continue
        dsig = math.sqrt(data_e2[b])
        psig = math.sqrt(post_e2[b])
        denom = math.sqrt(dsig * dsig + psig * psig)
        pull = (data_y[b] - post_y[b]) / denom if denom > 0 else 0.0
        result[b] = dict(
            data=data_y[b], data_err=dsig,

            pre=pre_y[b], pre_err=math.sqrt(pre_e2[b]),
            post=post_y[b], post_err=psig,
            pull=pull,
        )
    f.Close()
    return result


# ---------------------------------------------------------------------------
# Combine text-datacard parsing: pull the signal process's per-bin rate.
# Standard counting-datacard block:
#   bin      ch1   ch1   ch1   ch2  ...
#   process  SIG   BKG1  BKG2  SIG  ...
#   process   0     1     2    0   ...
#   rate     r1    r2    r3   r4   ...
# Signal = the process column(s) with process-number <= 0.
# ---------------------------------------------------------------------------
def _find_rate_block(lines):
    """Returns (header_line_index, bin_line, proc_name_line, proc_num_line,
    rate_line) for the bin/process/process/rate block, or None if not found.
    Shared by parse_signal_yield_per_bin and parse_gmN_entries so both use
    the exact same header parsing."""
    for i, ln in enumerate(lines):
        toks = ln.split()
        if not toks:
            continue
        if toks[0] == "bin" and i > 0:
            if i + 3 < len(lines) and lines[i + 1].split()[0:1] == ["process"] \
               and lines[i + 2].split()[0:1] == ["process"] \
               and lines[i + 3].split()[0:1] == ["rate"]:
                return (i, toks[1:], lines[i + 1].split()[1:],
                        lines[i + 2].split()[1:], lines[i + 3].split()[1:])
    return None


def parse_gmN_entries(datacard_path, nbins, proc_names, target_bins=None):
    """Scans the systematics block (everything after the rate line) for
    gmN entries affecting any of proc_names, in any of target_bins (or all
    bins if target_bins is None). Returns list of dicts: bin, process,
    syst_name, N (control-region count), kappa (scale factor entry)."""
    with open(datacard_path) as fh:
        lines = fh.readlines()

    block = _find_rate_block(lines)
    if block is None:
        return []
    header_idx, bin_line, proc_name_line, proc_num_line, rate_line = block

    # map column index -> (bin_index, process_name), for columns we care about
    cols_of_interest = {}
    for col, (chname, pname) in enumerate(zip(bin_line, proc_name_line)):
        bidx = bin_index_from_name(chname)
        if bidx is None:
            continue
        if target_bins is not None and bidx not in target_bins:
            continue
        if pname not in proc_names:
            continue
        cols_of_interest[col] = (bidx, pname)

    results = []
    for ln in lines[header_idx + 4:]:
        toks = ln.split()
        if len(toks) < 3 or toks[1] != "gmN":
            continue
        syst_name, N = toks[0], toks[2]
        values = toks[3:]
        for col, (bidx, pname) in cols_of_interest.items():
            if col < len(values) and values[col] != "-":
                results.append(dict(bin=bidx, process=pname, syst_name=syst_name,
                                     N=N, kappa=values[col]))
    return results


def parse_signal_yield_per_bin(datacard_path, nbins):
    with open(datacard_path) as fh:
        lines = fh.readlines()

    block = _find_rate_block(lines)
    if block is None:
        sys.exit(f"ERROR: could not locate bin/process/process/rate block in {datacard_path}")
    _, bin_line, proc_name_line, proc_num_line, rate_line = block

    yields = {}
    for chname, pnum, rate in zip(bin_line, proc_num_line, rate_line):
        try:
            pnum_i = int(float(pnum))
            rate_f = float(rate)
        except ValueError:
            continue
        if pnum_i > 0:
            continue  # background process, skip
        bidx = bin_index_from_name(chname)
        if bidx is None or bidx < 1 or bidx > nbins:
            continue
        yields[bidx] = yields.get(bidx, 0.0) + rate_f

    return yields


def get_fitted_r(fitdiag_path):
    """Reads the best-fit signal strength r and its asymmetric errors
    from the fit_s RooFitResult that combine -M FitDiagnostics stores
    in the same file (this is the s+b best fit already performed --
    no extra combine run needed)."""
    f = ROOT.TFile.Open(fitdiag_path)
    if not f or f.IsZombie():
        return None
    fit_s = f.Get("fit_s")
    if not fit_s:
        f.Close()
        return None
    pars = fit_s.floatParsFinal()
    r_par = pars.find("r")
    if not r_par:
        f.Close()
        return None
    val, elo, ehi = r_par.getVal(), r_par.getErrorLo(), r_par.getErrorHi()
    f.Close()
    return val, elo, ehi


NORMPOSTFIT_PROCESSES = ["total_background", "elefakepho", "jetfakepho",
                          "qcdfakelep", "VGamma", "ttGam_VVGam"]


def read_normpostfit_bin(normpostfit_dir, proc, year, bin_index):
    """Reads one line (1 bin) from NormPostFit_<proc>_<year>.txt, written
    by run_postfit.sh's extract_normpostfit.C: each line is
    'pre_y  pre_e  post_y  post_e', one line per bin, in bin order."""
    path = os.path.join(normpostfit_dir, f"NormPostFit_{proc}_{year}.txt")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        lines = fh.readlines()
    if bin_index < 1 or bin_index > len(lines):
        return None
    pre_y, pre_e, post_y, post_e = (float(x) for x in lines[bin_index - 1].split())
    return dict(pre=pre_y, pre_err=pre_e, post=post_y, post_err=post_e)


def report_bin_process_breakdown(normpostfit_dir, bins, year="Run2", sample_datacard=None, nbins=36):
    print(f"\n=== Per-process postfit breakdown, {year}, bins {bins} ===")
    print(f"{'bin':>4} {'process':>16} {'pre':>10} {'post':>10} {'ratio':>7} {'delta':>9}")
    for b in bins:
        for proc in NORMPOSTFIT_PROCESSES:
            info = read_normpostfit_bin(normpostfit_dir, proc, year, b)
            if info is None:
                continue
            ratio = info["post"] / info["pre"] if info["pre"] > 0 else float("nan")
            delta = info["post"] - info["pre"]
            flag = "  <-- largest mover" if abs(ratio - 1.0) > 0.15 and proc != "total_background" else ""
            print(f"{b:4d} {proc:>16} {info['pre']:10.3f} {info['post']:10.3f} "
                  f"{ratio:7.3f} {delta:9.3f}{flag}")
    print("\nThe non-total_background process with the largest |ratio-1| (or "
          "largest |delta| relative to its own size) in these bins is the most "
          "likely source of the systematic mismodeling -- that process's "
          "fake-rate/estimation method should be checked specifically in this "
          "kinematic region.")

    if sample_datacard:
        procs = [p for p in NORMPOSTFIT_PROCESSES if p != "total_background"]
        gmn = parse_gmN_entries(sample_datacard, nbins, procs, target_bins=set(bins))
        print(f"\n=== gmN control-region counts (from {sample_datacard}) ===")
        if not gmn:
            print("  No gmN systematics found for these processes/bins -- "
                  "check whether these processes use a different systematic "
                  "type (lnN, shape, or a rateParam) instead.")
        else:
            print(f"{'bin':>4} {'process':>16} {'syst_name':>20} {'N':>5} {'kappa':>8}")
            for e in sorted(gmn, key=lambda d: (d["bin"], d["process"])):
                flag = "  <-- LOW STATS control region" if e["N"].isdigit() and int(e["N"]) <= 2 else ""
                print(f"{e['bin']:4d} {e['process']:>16} {e['syst_name']:>20} "
                      f"{e['N']:>5} {e['kappa']:>8}{flag}")
            print("\n  N is the observed control-region event count backing that "
                  "process's estimate in that bin. N<=2 means a very weakly "
                  "constrained Gamma nuisance -- exactly the condition that lets "
                  "a correlated pull from other bins drive this bin's postfit "
                  "yield to (or near) zero, even though prefit was non-zero.")


def scan_boundary_collapses(normpostfit_dir, nbins, year="Run2", collapse_thresh=0.1):
    """Scans every bin for a background process whose postfit yield collapsed
    to near-zero relative to prefit (ratio < collapse_thresh) -- the signature
    of an under-constrained rateParam/systematic hitting its boundary, rather
    than a normal, bounded nuisance-parameter pull."""
    procs = [p for p in NORMPOSTFIT_PROCESSES if p != "total_background"]
    hits = {p: [] for p in procs}
    for b in range(1, nbins + 1):
        for proc in procs:
            info = read_normpostfit_bin(normpostfit_dir, proc, year, b)
            if info is None or info["pre"] <= 0:
                continue
            ratio = info["post"] / info["pre"]
            if ratio < collapse_thresh:
                hits[proc].append((b, info["pre"], info["post"], ratio))

    print(f"\n=== Scan for postfit boundary-collapse (post/pre < {collapse_thresh}), "
          f"{year}, bins 1-{nbins} ===")
    any_hits = False
    for proc, bin_hits in hits.items():
        if not bin_hits:
            continue
        any_hits = True
        bins_str = ", ".join(str(b) for b, *_ in bin_hits)
        print(f"  {proc}: collapsed in {len(bin_hits)}/{nbins} bins -> [{bins_str}]")
    if not any_hits:
        print(f"  No process collapsed below post/pre={collapse_thresh} in any bin -- "
              "the pattern seen in your checked bins is not widespread.")
    else:
        print("\n  If a process collapses in many bins concentrated in the high-bin-"
              "number range (roughly the high-MET/high-pT tail), that's the region "
              "driving the gluino>1400 GeV exclusion gap, and this process's "
              "rateParam/systematic is the most likely root cause: check the "
              "datacard for how its normalization is implemented (a freely-"
              "floating rateParam with no Gaussian constraint term will let the "
              "fit zero it out with no real penalty; a properly constrained lnN "
              "or shape systematic tied to the actual control-region precision "
              "would not).")
    return hits


def parse_limit_file(path):
    obs = exp50 = None
    with open(path) as fh:
        for line in fh:
            if "Observed" in line:
                obs = float(line.split("<")[1].strip())
            elif "50.0%" in line:
                exp50 = float(line.split("<")[1].strip())
    return obs, exp50


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0 or syy <= 0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--postfit-root", required=True,
                    help="fitDiagnostics_T5WG_Run2_obs.root from run_postfit.sh")
    ap.add_argument("--card-dir", required=True,
                    help="cards_Run2_unblind directory (same one runLimit.sh reads from)")
    ap.add_argument("--limit-dir", required=True,
                    help="directory of limit_<gluino>_<chargino>.txt from runLimit.sh")
    ap.add_argument("--model", default="t5wg",
                    help="lowercase model tag used in file names (default: t5wg)")
    ap.add_argument("--nbins", type=int, default=36)
    ap.add_argument("--out-csv", default="pull_vs_limit_gap.csv")
    ap.add_argument("--out-pdf", default="pull_vs_limit_gap.pdf")
    ap.add_argument("--compare-sb", action="store_true",
                    help="Also extract shapes_fit_s (s+b postfit, already computed "
                         "by the same FitDiagnostics run) and the fitted signal "
                         "strength r for the benchmark used in --postfit-root, and "
                         "print a bin-by-bin b-only vs. s+b pull comparison. Only "
                         "meaningful for that one benchmark's mass point.")
    ap.add_argument("--gluino-min", type=float, default=None,
                    help="If set, also report an aggregate diagnostic for all mass "
                         "points with gluino mass >= this value: the average "
                         "per-bin acceptance shape across that subset, weighted "
                         "against the b-only pulls, to find which bins are "
                         "collectively responsible for a systematic (not "
                         "single-point) observed<expected gap across that region "
                         "of the mass plane -- e.g. --gluino-min 1400 to match a "
                         "visible gap for m_gluino > 1400 GeV in the exclusion plot.")
    ap.add_argument("--exp-near-one", type=float, nargs=2, default=None,
                    metavar=("LOW", "HIGH"),
                    help="Restrict the --gluino-min aggregate comparison to mass "
                         "points whose EXPECTED limit falls in [LOW, HIGH] (e.g. "
                         "0.3 1.5). Points far from expected r=1 sit far from the "
                         "actual exclusion contour and their obs/exp ratio is not "
                         "representative of the visible gap between the plotted "
                         "curves -- use this to restrict the comparison to "
                         "points that actually sit near the boundary.")
    ap.add_argument("--normpostfit-dir", default=None,
                    help="Directory containing NormPostFit_<proc>_<year>.txt "
                         "(the OUTDIR from run_postfit.sh). If set together with "
                         "--check-bins, prints a per-process pre/post breakdown "
                         "for those bins to identify which background component "
                         "is responsible for the pull there.")
    ap.add_argument("--check-bins", default=None,
                    help="Comma-separated bin numbers to break down by process "
                         "(requires --normpostfit-dir), e.g. --check-bins 18,36")
    ap.add_argument("--sample-datacard", default=None,
                    help="A single representative datacard to inspect for gmN "
                         "systematics on the --check-bins/processes (background "
                         "structure is the same across mass points, so any one "
                         "datacard in --card-dir works). If omitted but "
                         "--check-bins is given, the first datacard found in "
                         "--card-dir is used automatically.")
    ap.add_argument("--scan-collapses", action="store_true",
                    help="Scan ALL bins (not just --check-bins) for any "
                         "background process whose postfit yield collapsed to "
                         "near-zero relative to prefit -- the boundary-collapse "
                         "signature discussed for jetfakepho/qcdfakelep. "
                         "Requires --normpostfit-dir.")
    args = ap.parse_args()

    # -----------------------------------------------------------------
    # 1) Per-bin data/postfit-background pulls, Run2 combined
    # -----------------------------------------------------------------
    pulls = extract_bin_pulls(args.postfit_root, args.nbins)

    print(f"\n=== Per-bin Run2 postfit pulls  ({args.postfit_root}) ===")
    print(f"{'bin':>4} {'data':>8} {'post-bkg':>10} {'post-err':>9} {'pull':>7}")
    for b in sorted(pulls):
        p = pulls[b]
        flag = "  <-- |pull|>2" if abs(p["pull"]) > 2 else ""
        print(f"{b:4d} {p['data']:8.1f} {p['post']:10.2f} {p['post_err']:9.2f} "
              f"{p['pull']:7.2f}{flag}")

    if args.compare_sb:
        r_info = get_fitted_r(args.postfit_root)
        sb_pulls = None
        try:
            sb_pulls = extract_bin_pulls(args.postfit_root, args.nbins,
                                          postfit_dirname="shapes_fit_s")
        except SystemExit as e:
            print(f"\n[--compare-sb] {e}")

        if r_info:
            rval, elo, ehi = r_info
            print(f"\nBest-fit signal strength for this benchmark: "
                  f"r = {rval:.3f}  (+{ehi:.3f} / {elo:.3f})")
            if abs(rval) < 1e-6:
                print("r is consistent with zero: this benchmark's signal shape "
                      "does not meaningfully improve the fit over background-only.")

        if sb_pulls:
            print(f"\n=== b-only vs. s+b postfit pulls, same benchmark ===")
            print(f"{'bin':>4} {'pull(b)':>9} {'pull(s+b)':>10} {'reduced?':>9}")
            for b in sorted(pulls):
                pb = pulls[b]["pull"]
                ps = sb_pulls.get(b, {}).get("pull", float("nan"))
                improved = "yes" if (ps == ps and abs(ps) < abs(pb) - 0.1) else \
                           ("no" if ps == ps else "?")
                print(f"{b:4d} {pb:9.2f} {ps:10.2f} {improved:>9}")
            print("\nIf pull(s+b) is much smaller than pull(b) specifically in the "
                  "bins that were large under b-only, the excess there is well "
                  "described by this benchmark's signal shape. If pulls barely "
                  "move (or move in bins unrelated to the original large-pull "
                  "bins), the excess is not well explained by this particular "
                  "signal hypothesis -- look at background mismodeling instead.")

    # -----------------------------------------------------------------
    # 2) For every mass point with both a datacard and a limit file:
    #    build the signal-acceptance-weighted pull overlap, and the
    #    observed/expected limit ratio.
    # -----------------------------------------------------------------
    limit_files = glob.glob(os.path.join(args.limit_dir, "limit_*.txt"))
    rows = []
    for lf in sorted(limit_files):
        m = re.search(r"limit_(\d+)_(\d+)\.txt$", os.path.basename(lf))
        if not m:
            continue
        gluino, chargino = m.group(1), m.group(2)

        card_path = os.path.join(args.card_dir, f"counting_{args.model}_{gluino}_{chargino}.txt")
        if not os.path.exists(card_path):
            continue

        obs, exp50 = parse_limit_file(lf)
        if obs is None or exp50 is None or exp50 <= 0:
            continue

        sig_yield = parse_signal_yield_per_bin(card_path, args.nbins)
        total_sig = sum(sig_yield.values())
        if total_sig <= 0:
            continue

        overlap = sum(
            (sig_yield.get(b, 0.0) / total_sig) * pulls[b]["pull"]
            for b in pulls
        )
        dominant_bin = max(sig_yield, key=sig_yield.get) if sig_yield else None
        dominant_frac = sig_yield.get(dominant_bin, 0.0) / total_sig if dominant_bin else 0.0
        sig_frac = {b: y / total_sig for b, y in sig_yield.items()}

        rows.append(dict(
            gluino=int(gluino), chargino=int(chargino),
            overlap=overlap, obs=obs, exp=exp50, ratio=obs / exp50,
            dominant_bin=dominant_bin, dominant_frac=dominant_frac,
            sig_frac=sig_frac,
        ))

    if not rows:
        sys.exit("ERROR: no mass points had both a datacard and a limit file -- "
                 "check --card-dir/--limit-dir/--model.")

    # -----------------------------------------------------------------
    # 3) Write CSV, make the diagnostic scatter, print the verdict
    # -----------------------------------------------------------------
    with open(args.out_csv, "w") as fh:
        fh.write("gluino,chargino,pull_overlap,observed_limit,expected_limit,"
                 "obs_over_exp,dominant_bin,dominant_bin_frac\n")
        for r in rows:
            fh.write(f"{r['gluino']},{r['chargino']},{r['overlap']:.4f},"
                     f"{r['obs']:.4f},{r['exp']:.4f},{r['ratio']:.4f},"
                     f"{r['dominant_bin']},{r['dominant_frac']:.3f}\n")
    print(f"\nWrote {len(rows)} mass points -> {args.out_csv}")

    top_gap = sorted(rows, key=lambda r: r["ratio"], reverse=True)[:10]
    print(f"\n=== Top 10 mass points by observed/expected gap ===")
    print(f"{'gluino':>7} {'chargino':>9} {'obs/exp':>8} {'dom.bin':>8} "
          f"{'dom.frac':>9} {'pull(dom.bin)':>14}")
    for r in top_gap:
        db = r["dominant_bin"]
        db_pull = pulls.get(db, {}).get("pull", float("nan"))
        print(f"{r['gluino']:7d} {r['chargino']:9d} {r['ratio']:8.2f} "
              f"{db:8d} {r['dominant_frac']:9.2f} {db_pull:14.2f}")
    print("\nIf the same dominant_bin (and a large positive pull there) keeps "
          "showing up across these top-gap mass points, that single bin -- not "
          "a real signal anywhere on this grid -- is what's driving their "
          "excess. Worth checking that bin's background composition "
          "(which of elefakepho/jetfakepho/qcdfakelep/VGamma/ttGam_VVGam "
          "dominates it) and whether it sits in a low-stat corner of phase "
          "space where the fake-rate estimate is less reliable.")

    if args.gluino_min is not None:
        pool = rows
        if args.exp_near_one is not None:
            lo, hi = args.exp_near_one
            pool = [r for r in rows if lo <= r["exp"] <= hi]
            print(f"\n[--exp-near-one] restricting to {len(pool)}/{len(rows)} "
                  f"points with expected limit in [{lo}, {hi}]")
            if not pool:
                print("  No points survived this window -- widen it, or check "
                      "that --exp-near-one is in the same r units as your "
                      "limit_*.txt files.")
        group = [r for r in pool if r["gluino"] >= args.gluino_min]
        rest = [r for r in pool if r["gluino"] < args.gluino_min]
        if not group:
            print(f"\n[--gluino-min] no mass points with gluino >= {args.gluino_min}")
        else:
            def _median(vals):
                s = sorted(vals)
                n = len(s)
                return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])

            print(f"\n=== Aggregate check: gluino >= {args.gluino_min} GeV "
                  f"({len(group)} points) vs. gluino < {args.gluino_min} GeV "
                  f"({len(rest)} points) ===")
            print(f"  median obs/exp, gluino >= {args.gluino_min}: "
                  f"{_median([r['ratio'] for r in group]):.2f}")
            if rest:
                print(f"  median obs/exp, gluino <  {args.gluino_min}: "
                      f"{_median([r['ratio'] for r in rest]):.2f}")

            # Average (unweighted-by-point) acceptance shape across the group,
            # then weight against the b-only pulls to see which bins are
            # collectively responsible for a *systematic* gap across this
            # region of the plane, as opposed to any single point's own
            # dominant bin.
            avg_frac = {b: 0.0 for b in pulls}
            for r in group:
                for b, f in r["sig_frac"].items():
                    if b in avg_frac:
                        avg_frac[b] += f / len(group)

            contrib = {b: avg_frac[b] * pulls[b]["pull"] for b in avg_frac}
            ranked = sorted(contrib, key=lambda b: abs(contrib[b]), reverse=True)

            print(f"\n  Bins ranked by contribution to the average pull-overlap "
                  f"across this group (avg.acceptance_fraction x pull):")
            print(f"  {'bin':>4} {'avg.frac':>9} {'pull(b)':>8} {'contribution':>13}")
            for b in ranked[:10]:
                print(f"  {b:4d} {avg_frac[b]:9.3f} {pulls[b]['pull']:8.2f} "
                      f"{contrib[b]:13.4f}")
            print(f"\n  Sum of contributions (= average pull-overlap for this "
                  f"group): {sum(contrib.values()):.3f}")

            from collections import Counter
            dom_counts = Counter(r["dominant_bin"] for r in group)
            print(f"\n  Individual dominant-bin frequency within this group "
                  f"({len(group)} points) -- how many points actually have "
                  f"THIS bin as their single largest-acceptance bin, not just "
                  f"on average:")
            for b, cnt in dom_counts.most_common(10):
                print(f"    bin {b:2d}: dominant for {cnt}/{len(group)} points "
                      f"({100.0*cnt/len(group):.0f}%)  pull={pulls[b]['pull']:.2f}")
            print("\n  A high avg.frac driven by this bin being dominant for "
                  "MOST points confirms a broad structural effect. A high "
                  "avg.frac coming from a few points with very large fractions "
                  "and many with small ones is a different, more localized "
                  "story -- check which case this is before concluding the "
                  "whole gluino>=1400 region is uniformly affected.")
            print("\n  A gap that's spread across SEVERAL of the top bins here "
                  "(rather than one single bin dominating, as in the top-10 "
                  "table above) points to a broader, systematic effect -- most "
                  "likely a background component that is mismodeled across the "
                  "whole high-MET/high-pT tail these bins correspond to, rather "
                  "than a single fluctuating bin. Check the per-process "
                  "NormPostFit_<proc>_Run2.txt files (elefakepho, jetfakepho, "
                  "qcdfakelep, VGamma, ttGam_VVGam) for exactly these bins to "
                  "see which process needed the largest postfit correction "
                  "there -- that process's fake-rate/estimation method in this "
                  "kinematic region is the most likely root cause.")

    if args.normpostfit_dir and args.check_bins:
        try:
            bins_to_check = [int(x) for x in args.check_bins.split(",")]
        except ValueError:
            sys.exit(f"ERROR: --check-bins must be comma-separated integers, "
                     f"got '{args.check_bins}'")

        sample_datacard = args.sample_datacard
        if sample_datacard is None:
            found = sorted(glob.glob(os.path.join(args.card_dir, f"counting_{args.model}_*.txt")))
            sample_datacard = found[0] if found else None
            if sample_datacard:
                print(f"\n[gmN check] auto-picked sample datacard: {sample_datacard}")
            else:
                print(f"\n[gmN check] no datacard found in {args.card_dir} to inspect "
                      f"for gmN entries -- pass --sample-datacard explicitly.")

        report_bin_process_breakdown(args.normpostfit_dir, bins_to_check,
                                      sample_datacard=sample_datacard, nbins=args.nbins)

    if args.scan_collapses:
        if not args.normpostfit_dir:
            print("\n[--scan-collapses] requires --normpostfit-dir")
        else:
            scan_boundary_collapses(args.normpostfit_dir, args.nbins)

    xs = [r["overlap"] for r in rows]
    ys = [math.log(r["ratio"]) for r in rows]
    corr = pearson(xs, ys)

    g = ROOT.TGraph(len(rows))
    for i, r in enumerate(rows):
        g.SetPoint(i, r["overlap"], math.log(r["ratio"]))
    g.SetMarkerStyle(20)
    g.SetTitle(";signal-acceptance-weighted pull overlap;log(observed/expected limit)")

    c = ROOT.TCanvas("c_diag", "c_diag", 700, 600)
    g.Draw("AP")
    line = ROOT.TLine(g.GetXaxis().GetXmin(), 0, g.GetXaxis().GetXmax(), 0)
    line.SetLineStyle(2)
    line.Draw("same")
    c.SaveAs(args.out_pdf)
    print(f"Wrote scatter plot -> {args.out_pdf}")

    print(f"\nPearson correlation (pull overlap vs. log obs/exp): {corr:.3f}")
    if corr != corr:  # NaN
        print("Could not compute a correlation (degenerate overlap or ratio values).")
    elif corr > 0.5:
        print("STRONG positive correlation: mass points most exposed to the excess "
              "bins are exactly the ones with the largest observed>expected gap. "
              "This confirms the postfit-vs-limit 'inconsistency' is simply "
              "different SR-bin weighting across the mass plane (hypothesis #2) "
              "-- not a separate bug.")
    elif corr > 0.2:
        print("Weak-to-moderate positive correlation: SR-bin weighting explains "
              "PART of the effect, but something else (asymptotic-limit breakdown "
              "at low stats, signal contamination in the fake-rate control "
              "regions, or differing FitDiagnostics/AsymptoticLimits fit "
              "configurations) is likely contributing too.")
    else:
        print("No meaningful correlation: SR-bin weighting does NOT explain the "
              "mismatch. Look instead at: (a) AsymptoticLimits' asymptotic "
              "formula reliability in your lowest-stat bins -- cross-check a "
              "few mass points with HybridNew toys; (b) whether the fake-rate "
              "control regions (elefakepho/jetfakepho/qcdfakelep) are "
              "signal-depleted for the mass points showing the largest gaps; "
              "(c) --robustFit/tolerance differences between the FitDiagnostics "
              "call in run_postfit.sh and the AsymptoticLimits call in "
              "runLimit.sh.")


if __name__ == "__main__":
    main()
