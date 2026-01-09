#(int RunYear, int preVFP, int ichannel)

# Compile
g++ `root-config --cflags` plotFit.C -o plotFit.exe `root-config --libs`


./plotFit.exe 2016 1  2
./plotFit.exe 2016 0  2
./plotFit.exe 2017 1  2
./plotFit.exe 2018 1  2
  
./plotFit.exe 2016 1  1
./plotFit.exe 2016 0  1
./plotFit.exe 2017 1  1
./plotFit.exe 2018 1  1
