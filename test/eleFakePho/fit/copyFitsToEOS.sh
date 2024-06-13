#!/bin/bash

cd /eos/uscms/store/user/tmishra/elefakepho/DataFitting16preVFP
mkdir -p DY_ker Bw_ker Bw_expo
mkdir -p DY_ker/pt DY_ker/eta DY_ker/vtx
mkdir -p Bw_ker/pt Bw_ker/eta Bw_ker/vtx
mkdir -p Bw_expo/pt Bw_expo/eta Bw_expo/vtx
mkdir -p Toy 

mv Data_DY_ker_eta_*.png DY_ker/eta && cp /eos/uscms/store/user/tmishra/index.php DY_ker/eta
mv Data_Bw_ker_eta_*.png Bw_ker/eta && cp /eos/uscms/store/user/tmishra/index.php Bw_ker/eta 
mv Data_Bw_expo_eta_*.png Bw_expo/eta && cp /eos/uscms/store/user/tmishra/index.php Bw_expo/eta

mv Data_DY_ker_pt_*.png DY_ker/pt && cp /eos/uscms/store/user/tmishra/index.php DY_ker/pt
mv Data_Bw_ker_pt_*.png Bw_ker/pt && cp /eos/uscms/store/user/tmishra/index.php Bw_ker/pt
mv Data_Bw_expo_pt_*.png Bw_expo/pt && cp /eos/uscms/store/user/tmishra/index.php Bw_expo/pt

mv Data_DY_ker_vtx_*.png DY_ker/vtx && cp /eos/uscms/store/user/tmishra/index.php DY_ker/vtx
mv Data_Bw_ker_vtx_*.png Bw_ker/vtx && cp /eos/uscms/store/user/tmishra/index.php Bw_ker/vtx
mv Data_Bw_expo_vtx_*.png Bw_expo/vtx && cp /eos/uscms/store/user/tmishra/index.php Bw_expo/vtx

mv *ToyMC*.png Toy 

ls DY_ker/pt | wc -l
ls DY_ker/eta | wc -l
ls DY_ker/vtx | wc -l
ls Bw_ker/pt | wc -l
ls Bw_ker/eta | wc -l
ls Bw_ker/vtx | wc -l
ls Bw_expo/pt | wc -l
ls Bw_expo/eta | wc -l
ls Bw_expo/vtx | wc -l

cd ../
scp -r DataFitting16preVFP trmishra@lxplus.cern.ch:/eos/home-t/trmishra/www/Plots/
