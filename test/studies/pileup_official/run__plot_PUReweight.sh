g++ `root-config --cflags` plot_PUReweight.C -o plot_PUReweight.exe `root-config --libs`
./plot_PUReweight.exe 2016 preVFP
./plot_PUReweight.exe 2016 postVFP
./plot_PUReweight.exe 2017 ""
./plot_PUReweight.exe 2018 ""
