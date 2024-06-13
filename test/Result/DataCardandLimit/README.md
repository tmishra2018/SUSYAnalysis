# LepPhoLimit
limit setting

Scripts to create datacard.

analysis_SUSY.C:
  1. Find the mass of the gluino of each event and store this information in a "SUSYTree".
  2. Select e+gamma, mu+gamma events and store their kinamatic information in the "egTree", "mgTree".
  
analysis_T5WG.C:
  1. read the trees created by analysis_SUSY.C
  2. bin the events by their MET, HT and phoEt, and fill the number of events into TH2D histograms.
  3. scale the number of events by Lumi and cross section.
  
counting_t5Wg_XXX_YYY.txt:
  datacard template
 
testT5WGcard.py:
   create a datacard for each Gluino-Chargino mass point.
