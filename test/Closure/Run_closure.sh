#!/bin/bash
# closure_efakepho(int ichannel, int RunYear, bool preVFP)
#root -l -q "closure_efakepho.C(1, 2016, 1)"
#root -l -q "closure_efakepho.C(2, 2016, 1)"

#root -l -q "closure_efakepho.C(1, 2016, 0)"
#root -l -q "closure_efakepho.C(2, 2016, 0)"

#root -l -q "closure_efakepho.C(1, 2017, 1)"
#root -l -q "closure_efakepho.C(2, 2017, 1)"

#root -l -q "closure_efakepho.C(1, 2018, 1)"
#root -l -q "closure_efakepho.C(2, 2018, 1)"
# closure_jetfakepho(int ichannel, int RunYear, bool preVFP)
#root -l -q "closure_jetfakepho.C(1, 2016, 1)"
#root -l -q "closure_jetfakepho.C(2, 2016, 1)"
#
#root -l -q "closure_jetfakepho.C(1, 2016, 0)"
#root -l -q "closure_jetfakepho.C(2, 2016, 0)"

#root -l -q "closure_jetfakepho.C(1, 2017, 1)"
#root -l -q "closure_jetfakepho.C(2, 2017, 1)"

#root -l -q "closure_jetfakepho.C(1, 2018, 1)"
#root -l -q "closure_jetfakepho.C(2, 2018, 1)"
# closure_jetfakelep(int ichannel, int RunYear, bool preVFP)
root -l -q "closure_jetfakelep.C(1, 2016, 1)"
root -l -q "closure_jetfakelep.C(2, 2016, 1)"

root -l -q "closure_jetfakelep.C(1, 2016, 0)"
root -l -q "closure_jetfakelep.C(2, 2016, 0)"

root -l -q "closure_jetfakelep.C(1, 2017, 1)"
root -l -q "closure_jetfakelep.C(2, 2017, 1)"

root -l -q "closure_jetfakelep.C(1, 2018, 1)"
root -l -q "closure_jetfakelep.C(2, 2018, 1)"
