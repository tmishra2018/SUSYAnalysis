cd "$(dirname "$0")"
g++ `root-config --cflags` MuonTrigger_ExtraPlots.C -o MuonTrigger_ExtraPlots.exe `root-config --libs`
./MuonTrigger_ExtraPlots.exe 2016 1
./MuonTrigger_ExtraPlots.exe 2016 0
./MuonTrigger_ExtraPlots.exe 2017 1
./MuonTrigger_ExtraPlots.exe 2018 1
