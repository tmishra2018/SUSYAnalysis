#!/usr/bin/python
import argparse, os, tempfile, shutil, sys,math,pickle,itertools
from subprocess import call, PIPE, STDOUT, Popen
import argparse
import datetime
import os
from os import system, environ

parser = argparse.ArgumentParser(description='Submit job to batch')
parser.add_argument('--Arg0',type=int,dest="Arg0",default=1, help="Arg0")
parser.add_argument('--Arg1',type=int,dest="Arg1",default=0, help="Arg1")
parser.add_argument('--Arg2',type=int,dest="Arg2",default=0, help="Arg2")
parser.add_argument('--Arg3',type=int,dest="Arg3",default=30, help="Arg3")
parser.add_argument('--Arg4',type=int,dest="Arg4",default=35, help="Arg4")
parser.add_argument('--Arg5',type=int,dest="Arg5",default=60, help="Arg5")
parser.add_argument('--Arg6',type=int,dest="Arg6",default=120, help="Arg6")
args = parser.parse_args()

system("./FitKer.exe {} {} {} {} {} {} {}".format(args.Arg0,args.Arg1,args.Arg2,args.Arg3,args.Arg4,args.Arg5,args.Arg6))   
system("xrdcp EleFakeRate*.txt root://cmseos.fnal.gov//store/user/tmishra/Fitting/")
system("xrdcp *.png root://cmseos.fnal.gov//store/user/tmishra/Fitting/")
