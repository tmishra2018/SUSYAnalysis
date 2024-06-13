rm temp.txt
root -l -q "plot_bkg.C(1,2016,1)" >> temp.txt
root -l -q "plot_bkg.C(2,2016,1)" >> temp.txt

root -l -q "plot_bkg.C(1,2016,0)" >> temp.txt
root -l -q "plot_bkg.C(2,2016,0)" >> temp.txt

root -l -q "plot_bkg.C(1,2017,1)" >> temp.txt
root -l -q "plot_bkg.C(2,2017,1)" >> temp.txt
 
root -l -q "plot_bkg.C(1,2018,1)" >> temp.txt
root -l -q "plot_bkg.C(2,2018,1)" >> temp.txt

grep -F "2016  preVFP" temp.txt
grep -F "2016  postVFP" temp.txt
grep -F "2017 " temp.txt
grep -F "2018 " temp.txt

#scp -r /eos/uscms/store/user/tmishra/Background/plots/* trmishra@lxplus.cern.ch:/eos/home-t/trmishra/www/Plots/SUSYAnalysis/ControlRegion
