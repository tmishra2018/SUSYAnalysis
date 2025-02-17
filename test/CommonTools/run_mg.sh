# int RunYear, const char *Era /*for Data*/, bool preVFP, bool isMC, const char *Sample /*for MC*/

g++ `root-config --cflags` ../../lib/libAnaClasses.so sigtree_mg.C -o sigtree_mg.exe `root-config --libs`

# MC
# 2016 preVFP
./sigtree_mg.exe 2016 "" 1 1 WGJet40
./sigtree_mg.exe 2016 "" 1 1 WGJet130
./sigtree_mg.exe 2016 "" 1 1 WGToLNuG
./sigtree_mg.exe 2016 "" 1 1 TTJets
./sigtree_mg.exe 2016 "" 1 1 TTGJets
./sigtree_mg.exe 2016 "" 1 1 WWG
./sigtree_mg.exe 2016 "" 1 1 WZG
./sigtree_mg.exe 2016 "" 1 1 WW
./sigtree_mg.exe 2016 "" 1 1 WZ
./sigtree_mg.exe 2016 "" 1 1 DYJetsToLL
./sigtree_mg.exe 2016 "" 1 1 ZGToLLG

# 2016 postVFP
./sigtree_mg.exe 2016 "" 0 1 WGJet40
./sigtree_mg.exe 2016 "" 0 1 WGJet130
./sigtree_mg.exe 2016 "" 0 1 WGToLNuG
./sigtree_mg.exe 2016 "" 0 1 TTJets
./sigtree_mg.exe 2016 "" 0 1 TTGJets 
./sigtree_mg.exe 2016 "" 0 1 WWG
./sigtree_mg.exe 2016 "" 0 1 WZG
./sigtree_mg.exe 2016 "" 0 1 WW
./sigtree_mg.exe 2016 "" 0 1 WZ 
./sigtree_mg.exe 2016 "" 0 1 DYJetsToLL
./sigtree_mg.exe 2016 "" 0 1 ZGToLLG

# 2017
./sigtree_mg.exe 2017 "" 0 1 WGJet40
./sigtree_mg.exe 2017 "" 0 1 WGJet130
./sigtree_mg.exe 2017 "" 0 1 WGToLNuG
./sigtree_mg.exe 2017 "" 0 1 TTJets
./sigtree_mg.exe 2017 "" 0 1 TTGJets
./sigtree_mg.exe 2017 "" 0 1 WWG
./sigtree_mg.exe 2017 "" 0 1 WZG
./sigtree_mg.exe 2017 "" 0 1 WW
./sigtree_mg.exe 2017 "" 0 1 WZ
./sigtree_mg.exe 2017 "" 0 1 DYJetsToLL
./sigtree_mg.exe 2017 "" 0 1 ZGToLLG

# 2018
./sigtree_mg.exe 2018 "" 0 1 WGJet40
./sigtree_mg.exe 2018 "" 0 1 WGJet130
./sigtree_mg.exe 2018 "" 0 1 WGToLNuG
./sigtree_mg.exe 2018 "" 0 1 TTJets
./sigtree_mg.exe 2018 "" 0 1 TTGJets
./sigtree_mg.exe 2018 "" 0 1 WWG
./sigtree_mg.exe 2018 "" 0 1 WZG
./sigtree_mg.exe 2018 "" 0 1 WW
./sigtree_mg.exe 2018 "" 0 1 WZ
./sigtree_mg.exe 2018 "" 0 1 DYJetsToLL
./sigtree_mg.exe 2018 "" 0 1 ZGToLLG
