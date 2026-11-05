# Check that the FitKer condor jobs of one era/type wrote their outputs to EOS.
# A bin counts as done only if its txt file is newer than the batch submit time
# (mtime of condorSubFile/condor_sub_<tag>_1.txt), so leftovers from older runs are ignored.
# Missing jobs are written to RunnerScript_missing_<tag>.txt, to resubmit with
#   python3 command_HTcondor.py year preVFP isData doEB RunnerScript_missing_<tag>.txt
# Exits 1 if any DY-ker pt/eta/vtx bin is missing. ht is reported but not used by fitFakeFunc_pT_nVtx.C;
# makeTxt.py fills a missing Bw-ker bin from DY-ker or Bw-expo, and a missing Bw-expo bin from Bw-ker.
#
# Examples:
# python3 checkFits.py --Year 2018 --Type DATA
# python3 checkFits.py --Year 2016preVFP --Type DY

import os
import sys
import argparse
import datetime

EOS = "/eos/uscms/store/user/tmishra/elefakepho"
FITFUNC = {1: "Bw-ker", 2: "DY-ker", 3: "Bw-expo"}
BINTYPE = {0: "pt", 1: "eta", 2: "vtx", 3: "ht"}
FITTYPE = {0: "den", 1: "num"}
ERAS = {"2016preVFP": (2016, 1), "2016postVFP": (2016, 0), "2017": (2017, 0), "2018": (2018, 0)}

parser = argparse.ArgumentParser(description='Check FitKer outputs on EOS for one era and type')
parser.add_argument('--Year', type=str, required=True, choices=list(ERAS), help="The run era, e.g. 2016postVFP")
parser.add_argument('--Type', type=str, required=True, choices=["DATA", "DY"], help="DATA or DY")
parser.add_argument('--doEB', type=int, default=1, help="1 = barrel, 0 = endcap")
parser.add_argument('--Runner', type=str, default="RunnerScript.txt", help="Job list that was submitted")
args = parser.parse_args()

year, preVFP = ERAS[args.Year]
isData = 1 if args.Type == "DATA" else 0
process = "Data" if isData else "DY"
tag = "{}_{}_{}_{}".format(year, preVFP, isData, args.doEB)

subfile = "condorSubFile/condor_sub_{}_1.txt".format(tag)
if not os.path.exists(subfile):
    sys.exit("No submit file {}: batch {} was never submitted".format(subfile, tag))
since = os.path.getmtime(subfile)


def fmt(x):
    # Same as ostream's default float output in FitKer.C (0.0 -> 0, 0.05 -> 0.05)
    return "%g" % float(x)


def folder_suffix(bintype, low, high):
    # Same as the CROSS_TAG logic of command_HTcondor.sh
    var = {0: "nVtx", 2: "pT"}[bintype]
    if float(high) >= 10000:
        return "_{}-gt{}".format(var, fmt(low))
    return "_{}{}-{}".format(var, fmt(low), fmt(high))


expected = []  # (folder, filename, runner line, bintype)
with open(args.Runner) as f:
    for line in f:
        cols = line.split()
        if not cols:
            continue
        ch, bintype, fittype = int(cols[0]), int(cols[1]), int(cols[2])
        low, high, fitlow, fithigh = cols[3], cols[4], cols[5], cols[6]
        # FitKer cuts |eta| <= 1.4442 for EB, so higher eta bins have no entries and write nothing
        if bintype == 1 and args.doEB and float(low) >= 1.4442:
            continue
        name = "EleFakeRate-{}-{}-{}-{}-{}-{}-{}-{}".format(
            process, FITFUNC[ch], BINTYPE[bintype], FITTYPE[fittype], fmt(low), fmt(high), fitlow, fithigh)
        folder = "{}Result{}".format(args.Type, args.Year)
        if len(cols) >= 9:
            name += "-{}{}-{}".format("nvtx" if bintype == 0 else "pt", fmt(cols[7]), fmt(cols[8]))
            folder += folder_suffix(bintype, cols[7], cols[8])
        expected.append((folder, name + ".txt", line.strip(), BINTYPE[bintype]))

print("{} {}: batch submitted {}".format(args.Year, args.Type,
      datetime.datetime.fromtimestamp(since).strftime("%Y-%m-%d %H:%M")))

missing = []
summary = {}
for folder, name, line, var in expected:
    path = os.path.join(EOS, folder, name)
    done = os.path.exists(path) and os.path.getmtime(path) >= since
    key = (folder, var)
    n_done, n_all = summary.get(key, (0, 0))
    summary[key] = (n_done + done, n_all + 1)
    if not done:
        missing.append((folder, name, line, var))

for (folder, var), (n_done, n_all) in sorted(summary.items()):
    flag = "" if n_done == n_all else "   <-- missing {}".format(n_all - n_done)
    print("  {:40s} {:4s} {:4d}/{:<4d}{}".format(folder, var, n_done, n_all, flag))

missfile = "RunnerScript_missing_{}.txt".format(tag)
if missing:
    print("Missing outputs:")
    for folder, name, line, var in missing:
        print("  {}/{}".format(folder, name))
    with open(missfile, "w") as f:
        for _, _, line, _ in missing:
            f.write(line + "\n")
    print("Wrote {} jobs to {}".format(len(missing), missfile))
elif os.path.exists(missfile):
    os.remove(missfile)

blocking = [m for m in missing if m[3] != "ht" and "-DY-ker-" in m[1]]
for fit, fill in [("Bw-ker", "DY-ker or Bw-expo"), ("Bw-expo", "Bw-ker")]:
    n = len([m for m in missing if m[3] != "ht" and "-{}-".format(fit) in m[1]])
    if n:
        print("{} {} pt/eta/vtx bins missing, makeTxt.py fills them from {}".format(n, fit, fill))
if blocking:
    print("NOT READY: {} DY-ker pt/eta/vtx bins missing".format(len(blocking)))
    sys.exit(1)
print("READY: all DY-ker pt/eta/vtx bins done")
