#!/usr/bin/env python3
"""Bin-by-bin prefit vs postfit comparison of all backgrounds and data.

Prefit  : ../../../test/myResult/logs/SignalRegionTable_<tag>.tex   (written by test/myResult/plot_eventct.C)
Postfit : <postfitDir>/NormPostFit_<proc>_<year>.txt      (pre_y pre_e post_y post_e per bin)
Data    : logs/SignalRegionTable_<tag>.txt                (written by plot_eventct.C here)

Usage: python3 compare_prefit_postfit.py [year ...]   (default: all eras)
"""
import os
import re
import sys

os.chdir(os.path.dirname(os.path.abspath(__file__)))   # all paths below are relative to this folder

POSTFITDIR = "../../fits/postfit_unblind_T5WG"   # written by limits/fits/run_postfit.sh
PREFITDIR = "../../../test/myResult/logs"
NBINS = 36

# (short name, combine process name, column index in the prefit tex table)
PROCS = [
    ("e->g",    "elefakepho",       5),
    ("j->g",    "jetfakepho",       6),
    ("misid-l", "qcdfakelep",       7),
    ("Vg",      "VGamma",           8),
    ("ttg/VVg", "ttGam_VVGam",      9),
    ("Total",   "total_background", 10),
]
YEARS = ["2016preVFP", "2016postVFP", "2017", "2018", "Run2"]


def file_tag(year):
    # plot_eventct.C names Run2 outputs with RunYear=678
    return "678" if year == "Run2" else year


def read_prefit_tex(path):
    """Return {bin: [(val, err) per tex column 5..10]}"""
    rows = {}
    pm = re.compile(r"([-\d.]+)\s*\$\\pm\$\s*([-\d.]+)")
    with open(path) as f:
        for line in f:
            if not re.match(r"^\d+ &", line):
                continue
            cols = [c.strip() for c in line.split("&")]
            vals = []
            for c in cols[5:11]:
                m = pm.search(c)
                vals.append((float(m.group(1)), float(m.group(2))) if m else (0.0, 0.0))
            rows[int(cols[0])] = vals
    return rows


def read_normpostfit(proc, year):
    path = os.path.join(POSTFITDIR, "NormPostFit_%s_%s.txt" % (proc, year))
    out = []
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) == 4:
                out.append(tuple(float(x) for x in p))
    return out[:NBINS]


def read_data(tag):
    """Obs per bin from the postfit text table."""
    obs = {}
    with open("logs/SignalRegionTable_%s.txt" % tag) as f:
        for line in f:
            p = line.split()
            if p and p[0].isdigit():
                obs[int(p[0])] = float(p[-2])
    return obs


def ratio(a, b):
    return a / b if b > 0 else float("nan")


def compare(year):
    tag = file_tag(year)
    pre = read_prefit_tex(os.path.join(PREFITDIR, "SignalRegionTable_%s.tex" % tag))
    post = {name: read_normpostfit(proc, year) for name, proc, _ in PROCS}
    obs = read_data(tag)

    out = open("logs/PrefitPostfitComparison_%s.txt" % tag, "w")

    def w(s=""):
        out.write(s + "\n")

    w("# %s: prefit (myResult/plot_eventct.C) vs postfit (combine) per search bin" % year)
    w("# columns per process: prefit +- unc | postfit +- unc | post/pre")
    w("# pulls: (post - pre) / prefit unc, using the prefit unc from the macro table")
    w()

    sums_pre = {n: [0.0, 0.0] for n, _, _ in PROCS}   # [mug, eg]
    sums_post = {n: [0.0, 0.0] for n, _, _ in PROCS}
    sums_obs = [0.0, 0.0]
    mismatches = []

    for name, proc, col in PROCS:
        w("=" * 96)
        w("%s  (%s)" % (name, proc))
        w("%-4s %-4s %20s %20s %9s %8s" % ("Bin", "Ch", "Prefit", "Postfit", "Post/Pre", "Pull"))
        for ibin in range(1, NBINS + 1):
            ch = 0 if ibin <= NBINS // 2 else 1
            pv, pe = pre[ibin][col - 5]
            cpre_y, cpre_e, qv, qe = post[name][ibin - 1]
            sums_pre[name][ch] += pv
            sums_post[name][ch] += qv
            pull = (qv - pv) / pe if pe > 0 else float("nan")
            w("%-4d %-4s %10.2f +- %6.2f %10.2f +- %6.2f %9.3f %8.2f"
              % (ibin, "mug" if ch == 0 else "eg", pv, pe, qv, qe, ratio(qv, pv), pull))
            # macro prefit should equal combine's prefit central value
            if abs(pv - cpre_y) > max(0.01, 0.005 * abs(cpre_y)):
                mismatches.append((name, ibin, pv, cpre_y))
        w()

    for ibin in range(1, NBINS + 1):
        sums_obs[0 if ibin <= NBINS // 2 else 1] += obs.get(ibin, 0.0)

    # Compact per-bin summary: total bkg prefit/postfit vs data
    w("=" * 96)
    w("Total background vs data")
    w("%-4s %-4s %20s %20s %6s %9s %9s" % ("Bin", "Ch", "Prefit", "Postfit", "Obs", "Obs/Pre", "Obs/Post"))
    for ibin in range(1, NBINS + 1):
        pv, pe = pre[ibin][5]
        qv, qe = post["Total"][ibin - 1][2:]
        o = obs.get(ibin, 0.0)
        w("%-4d %-4s %10.2f +- %6.2f %10.2f +- %6.2f %6.0f %9.3f %9.3f"
          % (ibin, "mug" if ibin <= NBINS // 2 else "eg", pv, pe, qv, qe, o, ratio(o, pv), ratio(o, qv)))
    w()

    # Channel totals per process (central values only)
    w("=" * 96)
    w("Channel totals (central values)")
    w("%-9s %12s %12s %9s %12s %12s %9s" % ("Process", "mug pre", "mug post", "post/pre",
                                              "eg pre", "eg post", "post/pre"))
    for name, _, _ in PROCS:
        a, b = sums_pre[name], sums_post[name]
        w("%-9s %12.2f %12.2f %9.3f %12.2f %12.2f %9.3f"
          % (name, a[0], b[0], ratio(b[0], a[0]), a[1], b[1], ratio(b[1], a[1])))
    w("%-9s %12.0f %12s %9s %12.0f" % ("Data", sums_obs[0], "", "", sums_obs[1]))
    w()

    if mismatches:
        w("WARNING: macro prefit central value differs from combine prefit (proc, bin, macro, combine):")
        for m in mismatches:
            w("  %-8s bin %2d  %10.3f  %10.3f" % m)
    else:
        w("Macro prefit central values agree with combine prefit in all bins.")
    out.close()

    print("%-12s written logs/PrefitPostfitComparison_%s.txt  (%d prefit mismatches vs combine)"
          % (year, tag, len(mismatches)))


if __name__ == "__main__":
    for y in (sys.argv[1:] or YEARS):
        compare(y)
