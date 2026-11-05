#!/usr/bin/env python3
# Split V+gamma skim on LPC condor (DAGMan). For each sample x era:
#   norm parts (NORM_CHUNK events each) -> merge (hadd = full-sample ntotalevent and signed MC PU) -> range chunks (CHUNK events)
# Outputs: /eos/uscms/store/user/tmishra/VGamma/<tag>/<sample>_<era>/chunks/resTree_VGamma_<sample>_<era>_<first>.root,
# merged afterwards with merge_skim.sh. Run from this directory, in the CMSSW environment, with a valid grid proxy:
#   python3 submit_skim.py [--tag zgisr] [--samples ZGToLLG,DYJetsToLL] [--eras 2018] [--chunk N] [--norm-chunk N] [--dry-run]
import argparse
import os
import subprocess
import ROOT

SAMPLES = ["WGToLNuG", "WGJet40", "WGJet130", "TTJets", "DYJetsToLL", "TTGJets", "WWG", "WZG", "WW", "WZ", "ZZ", "ZGToLLG"]
ERAS = {"2016postVFP": (2016, 0), "2016preVFP": (2016, 1), "2017": (2017, 0), "2018": (2018, 0)}
GROUP = ("DYJetsToLL", "GJets_DoubleEM", "TTJets", "ZGToLLG", "QCDMuEnriched", "QCD_DoubleEM")  # as analysis_VGamma_Run2.C
CHUNK = 10_000_000
NORM_CHUNK = 10_000_000
HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.abspath(os.path.join(HERE, "../../.."))  # SUSYAnalysis


def input_file(sample, era):
    base = "/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni" if sample in GROUP else "/eos/uscms/store/user/tmishra/InputFilesMC"
    return f"{base}/{sample}/{sample}_{era}.root"


def n_events(path):
    f = ROOT.TFile.Open(path)
    n = f.Get("ggNtuplizer/EventTree").GetEntries()
    f.Close()
    return n


def ranges(n, size):
    return [(a, min(a + size, n)) for a in range(0, n, size)]


ap = argparse.ArgumentParser()
ap.add_argument("--tag", default="zgisr")
ap.add_argument("--samples", default=",".join(SAMPLES))
ap.add_argument("--eras", default=",".join(ERAS))
ap.add_argument("--chunk", type=int, default=CHUNK, help="events per range job")
ap.add_argument("--norm-chunk", type=int, default=NORM_CHUNK, help="events per norm job")
ap.add_argument("--dry-run", action="store_true")
a = ap.parse_args()

work = os.path.join("/uscms_data/d3/tmishra/vg_condor", a.tag)
os.makedirs(os.path.join(work, "logs"), exist_ok=True)
exe = os.path.join(HERE, "../src/analysis_VGamma_Run2_condor.exe")
inputs = ",".join([exe, f"{TOP}/lib/libAnaClasses.so", f"{TOP}/lib/analysis_photonDict_rdict.pcm", f"{TOP}/test/PUfiles"])
proxy = subprocess.run(["voms-proxy-info", "-path"], capture_output=True, text=True).stdout.strip()

with open(os.path.join(work, "skim.jdl"), "w") as j:
    j.write(f"""universe = vanilla
executable = {HERE}/skim_job.sh
arguments = $(args)
transfer_input_files = {inputs}
initialdir = {work}
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
output = {work}/logs/$(name).out
error = {work}/logs/$(name).err
log = {work}/logs/skim.log
request_memory = 2500
+DesiredOS = "EL9"
x509userproxy = {proxy}
use_x509userproxy = true
queue
""")

lines, njobs, total = [], 0, 0
for era in a.eras.split(","):
    year, vfp = ERAS[era]
    for s in a.samples.split(","):
        n = n_events(input_file(s, era))
        total += n
        d = f"/store/user/tmishra/VGamma/{a.tag}/{s}_{era}"
        os.makedirs(f"/eos/uscms{d}/norm", exist_ok=True)
        os.makedirs(f"/eos/uscms{d}/chunks", exist_ok=True)
        key = f"{s}_{era}"
        norms = [f"norm_{key}_{f}" for f, _ in ranges(n, a.norm_chunk)]
        for (f, l), name in zip(ranges(n, a.norm_chunk), norms):
            lines.append(f'JOB {name} {work}/skim.jdl\nVARS {name} name="{name}" args="norm {year} {vfp} {s} {f} {l} {d}"')
        merge = f"merge_{key}"
        lines.append(f'JOB {merge} {work}/skim.jdl\nVARS {merge} name="{merge}" args="merge {d}"')
        lines.append(f"PARENT {' '.join(norms)} CHILD {merge}")
        rng = []
        for f, l in ranges(n, a.chunk):
            name = f"range_{key}_{f}"
            rng.append(name)
            lines.append(f'JOB {name} {work}/skim.jdl\nVARS {name} name="{name}" args="range {year} {vfp} {s} {f} {l} {d}"')
        lines.append(f"PARENT {merge} CHILD {' '.join(rng)}")
        njobs += len(norms) + 1 + len(rng)
        print(f"{key:24s} {n:>11d} events: {len(norms)} norm, {len(rng)} range jobs")
lines.append("RETRY ALL_NODES 2")
dag = os.path.join(work, "skim.dag")
open(dag, "w").write("\n".join(lines) + "\n")
print(f"{njobs} jobs, {total / 1e6:.0f}M events; DAG {dag}")
if not a.dry_run:
    # submitted from the work area on /uscms_data, so nothing a job leaves behind lands in the 3 GB home area
    subprocess.run(["condor_submit_dag", "-f", dag], check=True, cwd=work)
