#!/bin/bash
# fitJetFunc(int ichannel,int  RunYear,bool ISpreVFP)

#root -l -q "fitJetFunc.C+(1,2016,1)"
root -l -q "fitJetFunc.C+(2,2016,1)"

#root -l -q "fitJetFunc.C+(1,2016,0)"
root -l -q "fitJetFunc.C+(2,2016,0)"

#root -l -q "fitJetFunc.C+(1,2017,1)"
root -l -q "fitJetFunc.C+(2,2017,1)"

#root -l -q "fitJetFunc.C+(1,2018,1)"
root -l -q "fitJetFunc.C+(2,2018,1)"
