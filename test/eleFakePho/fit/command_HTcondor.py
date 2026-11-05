#  python3 command_HTcondor.py
submitFileTT="""
universe = vanilla
executable            = command_HTcondor.sh
transfer_output_files   = ""
+JobFlavour            = "tomorrow"
Transfer_Input_Files  = files.tar.gz
#Request_Disk = 1024
Request_Memory = 4096
#Request_cpus = 1
#Requirements = OpSys == "LINUX" && (Arch != "DUMMY" )
Should_Transfer_Files = YES
x509userproxy = $ENV(X509_USER_PROXY)
"""
import os
import sys
from os import system, environ

# usage: python3 command_HTcondor.py year preVFP isData doEB [runnerFile]
#   runnerFile defaults to RunnerScript.txt; any other file (e.g. RunnerScript_missing_<tag>.txt
#   from checkFits.py) is a resubmission and gets its own submit-file tag, so the original
#   condor_sub_<tag>_1.txt that checkFits.py uses as the batch submit time is not touched.
if len(sys.argv) not in (5, 6):
    sys.exit("usage: python3 command_HTcondor.py year preVFP isData doEB [runnerFile]")
year, preVFP, isData, doEB = [int(a) for a in sys.argv[1:5]]
runner = sys.argv[5] if len(sys.argv) == 6 else "RunnerScript.txt"
tag = "{}_{}_{}_{}".format(year, preVFP, isData, doEB)
if runner != "RunnerScript.txt":
    tag += "_resub"

files = open(runner,"r") 
count = 0
system("mkdir -p logs condorSubFile")
for ij in files:
    cols = ij.strip().split(' ')
    fileParts = [submitFileTT]
    count += 1
    if len(cols) >= 9:
        # useCrossBinning line: 7 default args + crossbinlow + crossbinhigh
        fileParts.append("Arguments = {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(
            year, preVFP, isData, doEB,
            int(cols[0]), int(cols[1]), int(cols[2]), float(cols[3]), float(cols[4]),
            int(cols[5]), int(cols[6]), float(cols[7]), float(cols[8])))
    else:
        # default 7-argument line (useCrossBinning off)
        fileParts.append("Arguments = {} {} {} {} {} {} {} {} {} {} {}\n".format(
            year, preVFP, isData, doEB,
            int(cols[0]), int(cols[1]), int(cols[2]), float(cols[3]), float(cols[4]),
            int(cols[5]), int(cols[6])))
    fileParts.append("Queue \n \n")
    fout = open("condorSubFile/condor_sub_{}_{}.txt".format(tag, count),"w")
    fout.write(''.join(fileParts))
    fout.close()
    fileParts.pop(-1)
    fileParts.pop(-1)
    system('condor_submit condorSubFile/condor_sub_{}_{}.txt'.format(tag, count))
    del fileParts
files.close()
