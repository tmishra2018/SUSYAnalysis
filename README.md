# SUSYAnalysis: SUS-24-011 / AN-21-007

CMS Run 2 search for SUSY in events with one photon, one lepton (e or μ) and missing transverse momentum
(eγ and μγ channels; T5WG, T6WG and TChiWG interpretations).

The code is plain ROOT/C++ macros plus bash and python scripts. It is not a CMSSW package and nothing is built
with `scram b`. Big inputs and outputs live on EOS and in `/uscms_data/d3/tmishra`, never in the repository.

## Layout

```
CMSSW_14_1_0_pre4/src/
  HiggsAnalysis/CombinedLimit/   Combine v10.3.0 (used by limits/)
  CombineHarvester/              combineTool.py, plotImpacts.py
  cmsstyle/                      CMS plot style (used by the limit and post-fit plots)
  SUSYAnalysis/                  this repository
    include/ src/ lib/           shared classes (recoPhoton, recoEle, recoMuon, rawData, esfScaleFactor, RooFit shapes)
    SFs/ sf/                     ID/iso scale-factor maps read by esfScaleFactor
    macro/                       old (2021) plotting macros, kept for reference
    test/                        analysis chain, one folder per step (see test/README.md)
    limits/                      everything after the signal-region prediction:
      datacards/                   datacard writing and conversion to CMS nuisance names
      combine/                     Run 2 card combination and limit jobs (condor)
      fits/                        post-fit, impacts, goodness of fit, likelihood scan, masked-bin limits
      plots/                       exclusion contours; postfit/ = post-fit signal-region plot
      TChiWG/                      self-contained TChiWG card -> limit -> plot chain
      notes/                       replies to the reviewers
```

## Setup

```bash
# once: CMSSW release with Combine v10
cmsrel CMSSW_14_1_0_pre4 && cd CMSSW_14_1_0_pre4/src && cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
(cd HiggsAnalysis/CombinedLimit && git checkout v10.3.0)
git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
git clone https://github.com/cms-cat/cmsstyle.git cmsstyle
scram b -j 8
git clone -b SUS24011-Repo git@github.com:tmishra2018/SUSYAnalysis.git

# every session (the login shell is tcsh; the run scripts are bash)
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd CMSSW_14_1_0_pre4/src && eval "$(scramv1 runtime -sh)"

# shared library (rebuild after editing include/analysis_{photon,ele,muon,rawData,scalefactor}.h or src/*.cc)
cd SUSYAnalysis/include && make ana       # -> lib/libAnaClasses.so, lib/libRooFitClasses.so
```

Condor steps (e→γ fits, Vγ skims, limits) need a grid proxy: `voms-proxy-init -voms cms -valid 192:00`.

## Rules that apply to every step

- **Run each script from its own folder.** Most scripts use relative paths (`../../lib`, `logs/`).
- **One era at a time per folder.** `test/Background/runBkg.sh`, `test/myResult/runLimit.sh` and
  `test/ISRReweight/run_all_ISR.sh` rewrite the folder's config file (`BkgPredConfig.txt`, `SigConfig.txt`,
  `ISRConfig.txt`) before every job, and the macros read the year from it. Never run two of these loops in the same
  folder at the same time.
- **How the config file is found.** `esfScaleFactor` (`include/analysis_scalefactor.h`) picks the config from the end
  of the working directory (`/test/myResult`, `/test/Background`, `/test/ISRReweight`), so it works in any copy.
  `SetRunConfig()` and `SetSignalConfig()` in `include/analysis_commoncode.h` always read
  `test/Background/BkgPredConfig.txt`, `test/myResult/SigConfig.txt` and `test/myResult/binConfig.txt` of the area
  that contains the working directory (`anaBaseDir()` walks up to the folder holding `include/analysis_commoncode.h`),
  so macros in `test/studies/...`, `test/myResult/SRplots` or `limits/plots/postfit` also see those configs.
  So the config files of a copy (for example on lxplus) are its own. Some other files still contain absolute paths
  to this area (`include/analysis_tools.h`, `include/analysis_fakes.h`, `test/ISRReweight/plot_ISRweight.C`, TnP
  filesets, `limits/fits/phase2_2026-10-31/*.sh`): `git grep -l "$PWD"` from the top folder lists them.
- **Many steps overwrite production files** on EOS or in `/uscms_data/d3/tmishra/Output` (named below). Back them up
  first if you only want to test.
- Eras: `2016preVFP`, `2016postVFP`, `2017`, `2018` (scripts usually take `<year> <preVFP 1/0>`); `678` means Run 2.
- Inside a folder, scripts that are not named here are old or one-off studies; check before using them.

## Running the analysis

The steps are in the order they depend on each other.

### 1. Object scale factors, triggers, pileup

**Tag-and-probe SFs** (electron R9 and mini-isolation, photon R9), `test/ScaleFactor/`:
```bash
cd test/ScaleFactor
g++ $(root-config --cflags) ../../lib/libAnaClasses.so makeTnPTree.C    -o makeTnPTree.exe    $(root-config --libs)
g++ $(root-config --cflags) ../../lib/libAnaClasses.so makeTnPPhoTree.C -o makeTnPPhoTree.exe $(root-config --libs)
bash tnp/convert_year.sh <era> [pho]     # TnP trees -> /eos/uscms/store/user/tmishra/TnP_ggNtuple/<era>/
bash tnp/run_sf.sh <era>                 # electron R9 + miniIso SFs (egamma-tnp container)
bash tnp/run_phoR9.sh <era>              # photon R9 SF (after convert_year.sh <era> pho)
```
Output: `test/ScaleFactor/tnp/results/<era>/SF_*_<era>_hybrid.root`, read directly by `esfScaleFactor`.
The SF scripts reuse cached histograms in `/uscms_data/d3/tmishra/tnp_work/<era>`; remove that folder after making
new trees. Pileup JSON for TnP (inside the container): `python3 tnp/make_pu_json.py <era> <path to test/> pileup_<era>.json`.

**Trigger efficiencies and SFs**, `test/Trigger/`:
```bash
cd test/Trigger
setsid nohup bash rerun_trigger_data.sh > /eos/uscms/store/user/tmishra/Trigger/logs_rerun/runner.log 2>&1 &   # data ntuples, all eras
bash Run-mgTrig-mc.sh            # μγ DY MC (and the MET dataset)
bash runEGTrigger.sh             # eγ SFs -> /eos/.../Trigger/files/diphoton_{pholeg,eleg}_<era>.root
bash runMuonTrigger.sh           # μγ SFs -> /eos/.../Trigger/files/muonphoton_trigger_<era>.root
```
`esfScaleFactor` reads these files directly, so a rerun replaces the production SFs. `Run-egTrig.sh` and
`Run-mgTrig.sh` are superseded; do not use them.

**Data pileup profiles**, run from `test/`: `bash PUfiles/makeDataPUDistributions.sh`
→ `test/PUfiles/dataPU_<era>{,_PUUp,_PUDown}.root` (read as `../../PUfiles` by the skims).

### 2. Skims

| What | Command (folder) | Output |
|---|---|---|
| Data signal/control trees | `MAX_JOBS=8 bash run_egmg_parallel.sh`, then `bash merge_egmg_years.sh` (`test/dataSkim`) | `/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_{egsignal_DoubleEG,mgsignal_MuonEG}_<era>.root` |
| Vγ MC (condor) | build `src/analysis_VGamma_Run2_condor.exe`; `python3 submit_skim.py --tag <tag>`, then `bash merge_skim.sh <tag>` (`test/VGamma/condor`) | `/eos/.../VGamma/<tag>/`; copy by hand to `/eos/.../VGamma/`, which Background and myResult read |
| Signal trees | `bash run_Signal.sh` (`test/myResult/Signal-samples`) | `/uscms_data/d3/tmishra/signal_trees/resTree_<model>_<era>.root` |
| MC for the closure tests | `bash Run_MC_all.sh` (`test/Closure`) | `/eos/uscms/store/user/tmishra/{egMC,mgMC}/` |

`merge_egmg_years.sh` overwrites the per-era files that all later steps read (`DRYRUN=1` and `BACKUP=1` are available).

### 3. e→γ misidentification rate, `test/eleFakePho/`
```bash
cd test/eleFakePho/src && bash run_elefakepho.sh && bash run_bgTemplate.sh   # TnP trees and background templates
#   then hadd the data eras per period by hand: plot_elefakepho_DataTnP_dR05_<era>.root
cd ../fit
bash RunFits.sh "2016preVFP 2016postVFP 2017 2018" "DY Data" 1    # about 660 template fits per era on condor
bash processEleFake.sh <era>                                        # checkFits.py, makeTxt.py, runFakeFunc.sh
root -l -b -q 'fitFakeFunc_pT_nVtx.C++(<year>,<preVFP>,true)'      # DY function, needed by the closure test
```
Output: `/eos/uscms/store/user/tmishra/elefakepho/{DATAResult,DrellYanResult}<era>/`. `checkFits.py` uses the time
stamps in `fit/condorSubFile/`, so keep that folder.

### 4. jet→γ misidentification rate, `test/jetFakePho/`
```bash
cd test/jetFakePho/src && bash RunHad.sh && bash Run_mgGJet.sh     # hadron templates (Run_mgGJet is for the ISR step)
cd ../fit
rm -f /eos/uscms/store/user/tmishra/jetfakepho/txt*/JetFakeRate-{DoubleEG,MuonEG}-EB_pT_nVtx.txt   # the fit APPENDS to these
bash run_FitJetFake_pT_nVtx.sh && bash run_fitJetFunc.sh
```
Output: `/eos/uscms/store/user/tmishra/jetfakepho/txt<era>/JetFakeRate-transferfactor-*-EB_pT_nVtx.txt`.

### 5. ISR weights for Vγ, `test/ISRReweight/` (details in `test/ISRReweight/README.md`)
```bash
nohup bash test/ISRReweight/run_all_ISR.sh >& /eos/uscms/store/user/tmishra/ISRweighting/logs_ISR/run_all_ISR.log &
```
Four steps: skims, the jet fake rate in the ISR region (about 7 h), weights per era, `format_isrweights.py`.
Paste the printed blocks into `test/Background/analysis_VGBkg.C` and `test/myResult/pred_VGBkg.C`, then redo steps 6–8.

### 6. Vγ and misidentified-lepton normalisation
```bash
cd test/Background        # set ToDeriveScale=true (the other flags false) at the top of runBkg.sh
bash runBkg.sh && bash merge_controlTree_Run2.sh     # control trees -> /uscms_data/d3/tmishra/Background/
cd ../VGamma/fit
bash run_fit_ALL_HT.sh    # photon-pT binned SFs: set doPhoScan=true and RunYears at the top
bash run_fit_ALL.sh       # inclusive SFs (1000 toys) -> /eos/.../VGamma/VGamma_scalefactor_{eg,mg}_<era>.txt
bash RunDependence.sh     # -> /eos/.../VGamma/<era>/VGammaScaleHist_{eg,mg}_<era>.root, used by steps 7 and 8
```
The inclusive fallback factors in `include/analysis_commoncode.h` are updated by hand.

### 7. Control- and validation-region predictions, `test/Background/`
Set `makeCRPlots=true` or `makeVRPlots=true` in `runBkg.sh` and run `bash runBkg.sh` (4 eras × 2 channels).
The jobs run in the background, so check `logs/` for failures. Output: `/uscms_data/d3/tmishra/Background/`.

### 8. Signal-region prediction, `test/myResult/`
```bash
cd test/myResult
bash runLimit.sh     # despite the name: SR predictions for all eras, signal trees, factorise_rareBkg.py, event counts
bash makePlot.sh     # SR plots and the Run 2 (_678) merges needed by the Run 2 cards
```
Output: `/uscms_data/d3/tmishra/Output/signalTree_{egamma,mg}_<bkg>_<era>.root`,
`/uscms_data/d3/tmishra/signal_trees/signalTree_<model>_<era>.root`, `logs/SignalRegionTable_<era>.tex` and
`SR_counts_SUS-24-011_<era>.txt`.

### 9. Closure tests, `test/Closure/`
`bash Run_closure_temp.sh` (e→γ closure; needs the DY function from step 3). Do not use `Run_closure.sh`: it writes
the same output names with the older method.

### 10. Datacards, `limits/datacards/`
```bash
bash limits/datacards/makeCards.sh             # set years/bins/modes at the top -> /eos/.../CombinedLimit/<model>/cards/cards_<era>_<mode>/
python3 limits/datacards/convert_all_cards.py --dry-run   # then without --dry-run: CMS nuisance names, luminosity scheme
cd limits/combine && ./combine_cards_Run2.sh T5WG   # and T6WG -> cards_Run2_unblind/
```
`createDataCard.py` and `write*card.py` work inside `test/myResult` (they read its `logs/`) wherever they are started
from. TChiWG: `limits/TChiWG/combine_cards_Run2.sh`.

### 11. Limits, `limits/combine/` (condor)
```bash
cd limits/combine
mkdir -p T5WG_run_<date> && cp runLimit.py runLimit.sh T5WG_run_<date>/ && cd T5WG_run_<date>
python3 runLimit.py T5WG unblind     # T6WG the same; TChiWG with runLimitTChiWG.py / runLimitTChiWG.sh
```
Output: `/eos/uscms/store/user/tmishra/CombinedLimit/<model>/limits/limits_Run2_<mode>/`.

### 12. Fits and checks, `limits/fits/`
```bash
cd limits/fits
bash run_postfit.sh                          # FitDiagnostics per era and for Run 2 -> postfit_unblind_T5WG/
bash runImpacts_Run2_unblind.sh              # impacts and goodness of fit (choose tasks and mass points at the top)
bash phase2_2026-10-31/scan_r_1700_200.sh    # 1D likelihood scan
(cd masked_rerun && python3 submit.py T5WG)  # masked-bin limits (condor; run from masked_rerun/)
```

### 13. Plots, `limits/plots/`
```bash
cd limits/plots
python3 makeT5WGLimit.py Run2                # and makeT6WGLimit.py -> /uscms_data/d3/tmishra/Output/T5WG_Run2_2016-17-2018.pdf
bash postfit/run_plot_eventct_postFit.sh     # post-fit SR plots and SR_counts (reads ../fits/postfit_unblind_T5WG)
python3 postfit/compare_prefit_postfit.py    # pre-fit vs post-fit per bin
cd ../TChiWG && bash makelimit.sh && root -l -b -q plot_limit.C
```
These overwrite the figures used in the paper and the AN.

### 14. Figures into the AN
`bash test/notes/CopyingToAN/copy_to_AN.sh --dry-run`, then without `--dry-run`: copies the figures listed in
`an_figures_manifest.txt` to the AN area on EOS and on lxplus.

## Where things are stored

| What | Where |
|---|---|
| Input ntuples | `/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/`, `/eos/uscms/store/user/tmishra/InputFilesMC/` |
| Intermediate results (fake rates, Vγ, trigger, ISR, data trees) | `/eos/uscms/store/user/tmishra/{elefakepho,jetfakepho,VGamma,Trigger,ISRweighting,eg_mg_treesData}/` |
| Predictions, signal trees, plots | `/uscms_data/d3/tmishra/{Output,Background,signal_trees}/` |
| Datacards and limits | `/eos/uscms/store/user/tmishra/CombinedLimit/<model>/{cards,limits}/` |
| Archived old output and code | `/eos/uscms/store/user/tmishra/SUSYAnalysis_out/` (one `.tgz` per folder; restore with `tar xzf`) |

## References

- e/γ HLT: https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary, https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgHLTPathDetails
- HLT paths: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HLTPathsRunIIList, https://fwyzard.web.cern.ch/fwyzard/hlt/2016/summary.html
- Muon HLT: https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2016, https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2017,
  https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2018
- MET filters and uncertainties: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2,
  https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription
- SUSY recommendations: https://twiki.cern.ch/twiki/bin/view/CMS/SUSRecommendationsRun2Legacy,
  https://twiki.cern.ch/twiki/bin/view/CMS/SUSRecommendations18
- Photon ID and e/γ POG: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2,
  https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPOG, https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations#Fall17v2
- L1 prefiring: https://twiki.cern.ch/twiki/bin/view/CMS/L1PrefiringWeightRecipe
- Systematics naming: https://cms-analysis.docs.cern.ch/guidelines/systematics/
