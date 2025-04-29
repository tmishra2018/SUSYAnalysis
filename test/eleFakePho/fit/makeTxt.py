# usage : python makeTxt.py --Year 2016preVFP --Type Data
#!/usr/bin/python
import os
from os import system, environ
import argparse
parser = argparse.ArgumentParser(description='pass RunYear')
parser.add_argument('--Year',type=str,dest="Year",default=1, help="the year")
parser.add_argument('--Type',type=str,dest="Type",default=1, help="data or mc")

args = parser.parse_args()
Year = args.Year
Type = args.Type

path = "/eos/uscms/store/user/tmishra/elefakepho/{}Result{}/".format(Type,Year)

system("mkdir -p {}Result{}".format(Type,Year))
outDir = "{}Result{}".format(Type,Year)


system("cat {}*Bw-expo-pt*.txt >> {}/EleFakeRate-{}-Bw-expo-pt-60-120.txt".format(path,outDir,Type))
system("cat {}*Bw-ker-pt*.txt >> {}/EleFakeRate-{}-Bw-ker-pt-60-120.txt".format(path,outDir,Type))
system("cat {}*DY-ker-pt*.txt >> {}/EleFakeRate-{}-DY-ker-pt-60-120.txt".format(path,outDir,Type))

system("cat {}*Bw-expo-eta*.txt >> {}/EleFakeRate-{}-Bw-expo-eta-60-120.txt".format(path,outDir,Type))
system("cat {}*Bw-ker-eta*.txt >> {}/EleFakeRate-{}-Bw-ker-eta-60-120.txt".format(path,outDir,Type))
system("cat {}*DY-ker-eta*.txt >> {}/EleFakeRate-{}-DY-ker-eta-60-120.txt".format(path,outDir,Type))

system("cat {}*Bw-expo-vtx*.txt >> {}/EleFakeRate-{}-Bw-expo-vtx-60-120.txt".format(path,outDir,Type))
system("cat {}*Bw-ker-vtx*.txt >> {}/EleFakeRate-{}-Bw-ker-vtx-60-120.txt".format(path,outDir,Type))
system("cat {}*DY-ker-vtx*.txt >> {}/EleFakeRate-{}-DY-ker-vtx-60-120.txt".format(path,outDir,Type))
