# To run: python3 runLimit.py T5WG splusb
#     or: python3 runLimit.py T5WG unblind
#     or: python3 runLimit.py T6WG splusb
#     or: python3 runLimit.py T6WG unblind
# Add "mask" as a 3rd argument to mask bins 13,25,30,34,35,36 (test run), e.g.:
#         python3 runLimit.py T5WG unblind mask

import sys
import os
from subprocess import check_output

if len(sys.argv) < 3:
    print("Usage: python3 runLimit.py [T5WG|T6WG] [splusb|unblind] [mask]")
    sys.exit(1)

model = sys.argv[1]        # T5WG or T6WG
mode = sys.argv[2].lower() # splusb or unblind

if mode not in ("splusb", "unblind"):
    print("Error: mode must be 'splusb' or 'unblind'")
    sys.exit(1)

maskBins = len(sys.argv) > 3 and sys.argv[3].lower() == "mask"  # mask bins boolean
maskTag = "_masked_MET_200_400" if maskBins else ""

model_lower = model.lower()  # t5wg or t6wg

logDir = f"logs_{mode}_{model_lower}{maskTag}"

submitFileTT = f"""
universe              = vanilla
Executable            = runLimit.sh
request_memory        = 4096
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT_OR_EVICT
Transfer_Input_Files  = cardsSplitted
x509userproxy = $ENV(X509_USER_PROXY)
+JobFlavour            = "nextweek"
"""

# Output directory of the limit files; a plain file with this name (left by an earlier run) would block it
limitDir = f"/eos/uscms/store/user/tmishra/CombinedLimit/{model}/limits/limits_Run2_{mode}{maskTag}"
if os.path.exists(limitDir) and not os.path.isdir(limitDir):
    print(f"Error: {limitDir} is a file, not a directory; remove it before submitting")
    sys.exit(1)
os.makedirs(limitDir, exist_ok=True)
print(f"Limits will be written to: {limitDir}")

os.system(f"rm -rf {logDir}")
os.system(f"mkdir -p {logDir}")

os.system(f"rm -rf /eos/uscms/store/user/tmishra/susy/{model}/{mode}{maskTag}")
os.system(f"mkdir -p /eos/uscms/store/user/tmishra/susy/{model}/{mode}{maskTag}")

os.system("rm -rf cardsSplitted DatacardList.txt")
os.system("mkdir -p cardsSplitted")

path = f"/eos/uscms/store/user/tmishra/CombinedLimit/{model}/cards/cards_Run2_{mode}/"
os.system(f"ls -d -1 {path}*  >> DatacardList.txt")
os.system("cd cardsSplitted && split -l 2 -d -a 3 ../DatacardList.txt")

for i in range(0, 10):
    os.system(f"cd cardsSplitted && mv x00{i} x{i}")
for i in range(10, 100):
    os.system(f"cd cardsSplitted && mv x0{i} x{i}")

njobs = int(check_output("ls -l cardsSplitted | wc -l", shell=True).strip())
print(f"Number of jobs: {njobs}")

jobNumber = 0

while jobNumber < njobs - 1:
    fileParts = [submitFileTT]
    fileParts.append(f"Output    = {logDir}/job{jobNumber}.stdout\n")
    fileParts.append(f"error     = {logDir}/job{jobNumber}.stderr\n")
    fileParts.append(f"Log       = {logDir}/job{jobNumber}.log\n")
    fileParts.append(f"Arguments = {jobNumber} {model} {mode} {str(maskBins).lower()}\n")
    fileParts.append("Queue \n \n")

    with open(f"/eos/uscms/store/user/tmishra/susy/{model}/{mode}{maskTag}/condor_sub_{jobNumber}.txt", "w") as fout:
        fout.write(''.join(fileParts))

    command = f"condor_submit /eos/uscms/store/user/tmishra/susy/{model}/{mode}{maskTag}/condor_sub_{jobNumber}.txt"
    os.system(command)

    jobNumber += 1
