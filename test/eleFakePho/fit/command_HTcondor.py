#  python3 command_HTcondor.py
submitFileTT="""
universe = vanilla
executable            = command_HTcondor.sh
transfer_output_files   = ""
+JobFlavour            = "tomorrow"
Transfer_Input_Files  = files.tar.gz
#Request_Disk = 1024
#Request_Memory = 1024
#Request_cpus = 1
#Requirements = OpSys == "LINUX" && (Arch != "DUMMY" )
Should_Transfer_Files = YES
x509userproxy = $ENV(X509_USER_PROXY)
"""
import os
from os import system, environ

files = open("RunnerScript.txt","r") 
count = 0
system("mkdir -p logs condorSubFile")
for ij in files:
    cols = ij.strip().split(' ')
    fileParts = [submitFileTT]
    count += 1
#    fileParts.append("Output    = logs/job{}_$(Cluster)_$(Process).stdout\n".format(count))
#    fileParts.append("error     = logs/job{}_$(Cluster)_$(Process).stderr\n".format(count))
#    fileParts.append("Log       = logs/job{}_$(Cluster)_$(Process).log\n".format(count))
    fileParts.append("Arguments = {} {} {} {} {} {} {}\n".format(int(cols[0]),int(cols[1]),int(cols[2]),float(cols[3]),float(cols[4]),int(cols[5]),int(cols[6])))
    fileParts.append("Queue \n \n")
    fout = open("condorSubFile/condor_sub_{}.txt".format(count),"w")
    fout.write(''.join(fileParts))
    fout.close()
    fileParts.pop(-1)
    fileParts.pop(-1)
    system('condor_submit condorSubFile/condor_sub_{}.txt'.format(count))
    del fileParts
files.close()
