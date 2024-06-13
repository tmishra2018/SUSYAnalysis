#(int RunYear, int preVFP, int isQCD, int ichannel)
rm scaleFactors.txt
g++ `root-config --cflags` plotGauss.C -o plotGauss.exe `root-config --libs`

#./plotGauss.exe 2016 1 1 1
#./plotGauss.exe 2016 1 1 2
#./plotGauss.exe 2016 1 0 1
#./plotGauss.exe 2016 1 0 2
#
#./plotGauss.exe 2016 0 1 1
#./plotGauss.exe 2016 0 1 2
./plotGauss.exe 2016 0 0 1
#./plotGauss.exe 2016 0 0 2
#
#./plotGauss.exe 2017 1 1 1
#./plotGauss.exe 2017 1 1 2
./plotGauss.exe 2017 1 0 1
#./plotGauss.exe 2017 1 0 2
#
./plotGauss.exe 2018 1 1 1
#./plotGauss.exe 2018 1 1 2
./plotGauss.exe 2018 1 0 1
#./plotGauss.exe 2018 1 0 2
