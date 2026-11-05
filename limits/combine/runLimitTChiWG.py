# Submit TChiWG Run 2 observed limits to condor, one job per mass point.
# To run (from a fresh run directory holding runLimitTChiWG.py and runLimitTChiWG.sh):
#     python3 runLimitTChiWG.py            (unblind cards)
#     python3 runLimitTChiWG.py unblind dry   (write the submit file only)
# Output: /eos/uscms/store/user/tmishra/CombinedLimit/TChiWG/limits/limits_Run2/limit_observed_<mass>.txt

import glob
import os
import re
import subprocess
import sys

mode = sys.argv[1].lower() if len(sys.argv) > 1 else "unblind"
dry = len(sys.argv) > 2 and sys.argv[2].lower() == "dry"
if mode not in ("splusb", "unblind"):
    sys.exit("Error: mode must be 'splusb' or 'unblind'")

cardDir = f"/eos/uscms/store/user/tmishra/CombinedLimit/TChiWG/cards/cards_Run2_{mode}"
limitDir = "/eos/uscms/store/user/tmishra/CombinedLimit/TChiWG/limits/limits_Run2"

masses = sorted(int(re.search(r"counting_tchiwg_(\d+)\.txt$", f).group(1))
                for f in glob.glob(f"{cardDir}/counting_tchiwg_*.txt"))
if not masses:
    sys.exit(f"Error: no cards in {cardDir}")

if os.path.exists(limitDir) and not os.path.isdir(limitDir):
    sys.exit(f"Error: {limitDir} is a file, not a directory")
os.makedirs(limitDir, exist_ok=True)
old = glob.glob(f"{limitDir}/limit_observed_*.txt")
if old:
    sys.exit(f"Error: {limitDir} already holds {len(old)} limit files; move them aside first")

logDir = f"logs_{mode}_tchiwg"
os.makedirs(logDir, exist_ok=True)
with open("masses.txt", "w") as f:
    f.write("\n".join(str(m) for m in masses) + "\n")

sub = f"""universe              = vanilla
Executable            = runLimitTChiWG.sh
Arguments             = $(mass) {mode}
request_memory        = 4096
request_cpus          = 4
Should_Transfer_Files = YES
WhenToTransferOutput  = ON_EXIT_OR_EVICT
x509userproxy         = $ENV(X509_USER_PROXY)
Output                = {logDir}/job_$(mass).stdout
Error                 = {logDir}/job_$(mass).stderr
Log                   = {logDir}/jobs.log
queue mass from masses.txt
"""
with open("condor_sub_tchiwg.txt", "w") as f:
    f.write(sub)

print(f"{len(masses)} mass points ({masses[0]}-{masses[-1]}), limits -> {limitDir}")
if dry:
    print("dry run: condor_sub_tchiwg.txt written, nothing submitted")
else:
    subprocess.run("condor_submit condor_sub_tchiwg.txt", shell=True, check=True)  # LPC condor_submit is a wrapper: needs a shell
