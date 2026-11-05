# SUSYAnalysis: SUS-24-011 / AN-21-007

CMS Run 2 search for supersymmetry in events with one photon, one lepton (e or μ) and missing transverse momentum,
in the eγ and μγ channels, interpreted in the T5WG, T6WG and TChiWG simplified models.

The code is plain ROOT/C++ macros with bash and python drivers. It is not a CMSSW package: nothing here is built with
`scram b`. Large inputs and outputs live on EOS and in `/uscms_data/d3/tmishra`, never in the repository.

**Contents:** [Workflow](#workflow) · [Repository layout](#repository-layout) · [Setup](#setup) ·
[General rules](#general-rules) · [Steps 1–14](#running-the-analysis) · [Where outputs are stored](#where-outputs-are-stored) ·
[References](#references)

## Workflow

Each step uses the outputs of the steps before it.

```
 INPUTS (ntuples on EOS)
   │
   ├─ 1  Object SFs, trigger SFs, pileup profiles ─────────────┐
   │                                                            │
   ├─ 2  Skims: data trees, Vγ MC, signal trees, closure MC     │
   │                                                            │
   ├─ 3  e→γ misidentification rate ─┐                          │
   ├─ 4  jet→γ misidentification rate ┼─ 5  ISR weights for Vγ  │
   │                                  │                         │
   ├──────────────────────────────────┴─ 6  Vγ / misid.-lepton normalisation
   │                                        │
   │                                        ├─ 7  Control- and validation-region predictions
   │                                        ├─ 8  Signal-region prediction ── 9  Closure tests
   │                                        │
   └─────────────────────────────────────── 10 Datacards ── 11 Limits ── 12 Fits and checks
                                                                  │
                                                       13 Plots ── 14 Figures into the AN
```

| Step | What | Folder | Main command |
|---|---|---|---|
| 1 | Object, trigger and pileup corrections | `test/ScaleFactor`, `test/Trigger`, `test/PUfiles` | `tnp/run_sf.sh`, `runEGTrigger.sh`, `makeDataPUDistributions.sh` |
| 2 | Skims | `test/dataSkim`, `test/VGamma/condor`, `test/myResult/Signal-samples`, `test/Closure` | `run_egmg_parallel.sh` (`analysis_{eg,mg}.C`), `submit_skim.py`, `run_Signal.sh` |
| 3 | e→γ rate | `test/eleFakePho` | `RunFits.sh`, `processEleFake.sh` |
| 4 | jet→γ rate | `test/jetFakePho` | `run_FitJetFake_pT_nVtx.sh`, `run_fitJetFunc.sh` |
| 5 | ISR weights | `test/ISRReweight` | `run_all_ISR.sh` |
| 6 | Vγ normalisation | `test/Background`, `test/VGamma/fit` | `runBkg.sh`, `run_fit_ALL_HT.sh` |
| 7 | CR/VR predictions | `test/Background` | `runBkg.sh` |
| 8 | SR prediction | `test/myResult` | `runLimit.sh`, `makePlot.sh` |
| 9 | Closure tests | `test/Closure` | `run_closure_efakepho.sh` |
| 10 | Datacards | `limits/datacards`, `limits/combine` | `makeCards.sh`, `convert_all_cards.py`, `combine_cards_Run2.sh` |
| 11 | Limits | `limits/combine` | `runLimit.py` |
| 12 | Fits and checks | `limits/fits` | `run_postfit.sh`, `runImpacts_Run2_unblind.sh` |
| 13 | Plots | `limits/plots`, `limits/TChiWG` | `makeT5WGLimit.py`, `run_plot_eventct_postFit.sh` |
| 14 | AN figures | `test/notes/CopyingToAN` | `copy_to_AN.sh` |

## Repository layout

```
CMSSW_14_1_0_pre4/src/
  HiggsAnalysis/CombinedLimit/   Combine (interactive fits)
  CombineHarvester/              combineTool.py, plotImpacts.py
  cmsstyle/                      CMS plot style (limit and post-fit plots)
  SUSYAnalysis/                  this repository
    include/ src/ lib/           shared classes: recoPhoton, recoEle, recoMuon, rawData, esfScaleFactor, RooFit shapes
    SFs/ sf/                     ID/isolation scale-factor maps read by esfScaleFactor
    macro/                       2021 plotting macros, kept for reference
    test/                        analysis chain, one folder per step (overview in test/README.md)
      studies/ tools/ notes/ legacy/   side studies, helper scripts, AN figure copying, retired code
    limits/                      everything after the signal-region prediction
      datacards/                   datacard writing and conversion to CMS nuisance names
      combine/                     Run 2 card combination and limit jobs (condor)
      fits/                        post-fit, impacts, goodness of fit, likelihood scan, masked-bin limits
      plots/                       exclusion contours; postfit/ = post-fit signal-region plot
      TChiWG/                      self-contained TChiWG card → limit → plot chain
```

## Setup

**Once:**
```bash
cmsrel CMSSW_14_1_0_pre4 && cd CMSSW_14_1_0_pre4/src && cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
(cd HiggsAnalysis/CombinedLimit && git checkout v10.3.0)
git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
git clone https://github.com/cms-cat/cmsstyle.git cmsstyle
scram b -j 8
git clone -b SUS24011-Repo git@github.com:tmishra2018/SUSYAnalysis.git
cd SUSYAnalysis/include && make ana        # -> lib/libAnaClasses.so, lib/libRooFitClasses.so
```

**Every session** (the login shell is tcsh; the run scripts are bash):
```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd CMSSW_14_1_0_pre4/src && eval "$(scramv1 runtime -sh)"
```

- Rebuild the shared library (`make ana`) after editing `include/analysis_{photon,ele,muon,rawData,scalefactor}.h`
  or `src/*.cc`.
- Condor steps (e→γ fits, Vγ skims, limits) need a grid proxy: `voms-proxy-init -voms cms -valid 192:00`.
- Combine versions: the local install (v10.3.0) is used for the interactive fits in step 12. The limit jobs of step 11
  build Combine from the GitHub master branch on the worker node; the production limits (31 Oct 2026) used v11.1.0.

## General rules

- **Run each script from its own folder.** Most scripts use relative paths (`../../lib`, `logs/`).
- **One era at a time per folder.** `test/Background/runBkg.sh`, `test/myResult/runLimit.sh` and
  `test/ISRReweight/run_all_ISR.sh` rewrite their folder's config file (`BkgPredConfig.txt`, `SigConfig.txt`,
  `ISRConfig.txt`) before every job, and the macros read the year from it. Never run two of these loops in the same
  folder at the same time.
- **Production files are overwritten.** Many steps write to EOS or to `/uscms_data/d3/tmishra/Output` (named in each
  step). Back these up before running a test.
- **Eras:** `2016preVFP`, `2016postVFP`, `2017`, `2018`; scripts usually take `<year> <preVFP 1/0>`. `678` means Run 2.
- **Only the scripts named below are part of the chain.** Other scripts in a folder are studies or helpers.

<details>
<summary>How the macros find their config files</summary>

- `esfScaleFactor` (`include/analysis_scalefactor.h`) picks the config from the end of the working directory
  (`/test/myResult`, `/test/Background`, `/test/ISRReweight`), so it works in any copy of the area.
- `SetRunConfig()` and `SetSignalConfig()` (`include/analysis_commoncode.h`) read `test/Background/BkgPredConfig.txt`,
  `test/myResult/SigConfig.txt` and `test/myResult/binConfig.txt` of the area that contains the working directory
  (`anaBaseDir()` walks up to the folder holding `include/analysis_commoncode.h`). Macros in `test/studies/...` or
  `limits/plots/postfit` therefore see the same configs, and a copy of the area (for example on lxplus) uses its own.
- Some files still contain absolute paths to this area (`include/analysis_tools.h`, `include/analysis_fakes.h`,
  `test/ISRReweight/plot_ISRweight.C`, the TnP filesets, `limits/fits/phase2_2026-10-31/*.sh`).
  `git grep -l "$PWD"` from the top folder lists them.
</details>

## Running the analysis

### Step 1. Object scale factors, triggers, pileup

**1a. Tag-and-probe SFs** (electron R9 and mini-isolation, photon R9), folder `test/ScaleFactor/`
```bash
g++ $(root-config --cflags) ../../lib/libAnaClasses.so makeTnPTree.C    -o makeTnPTree.exe    $(root-config --libs)
g++ $(root-config --cflags) ../../lib/libAnaClasses.so makeTnPPhoTree.C -o makeTnPPhoTree.exe $(root-config --libs)
bash tnp/convert_year.sh <era> [pho]     # TnP trees -> /eos/uscms/store/user/tmishra/TnP_ggNtuple/<era>/
bash tnp/run_sf.sh <era>                 # electron R9 + mini-isolation SFs (egamma-tnp container)
bash tnp/run_phoR9.sh <era>              # photon R9 SF (after convert_year.sh <era> pho)
```
- Output: `tnp/results/<era>/SF_*_<era>_hybrid.root`, read directly by `esfScaleFactor`.
- The SF scripts reuse cached histograms in `/uscms_data/d3/tmishra/tnp_work/<era>`; remove that folder after making
  new trees.
- Pileup JSON for TnP (inside the container): `python3 tnp/make_pu_json.py <era> <path to test/> pileup_<era>.json`.

**1b. Trigger efficiencies and SFs**, folder `test/Trigger/`
```bash
setsid nohup bash rerun_trigger_data.sh > /eos/uscms/store/user/tmishra/Trigger/logs_rerun/runner.log 2>&1 &   # data, all eras
bash Run-mgTrig-mc.sh            # μγ DY MC (and the MET dataset)
bash runEGTrigger.sh             # eγ SFs
bash runMuonTrigger.sh           # μγ SFs
```
- Output: `/eos/uscms/store/user/tmishra/Trigger/files/{diphoton_pholeg,diphoton_eleg,muonphoton_trigger}_<era>.root`,
  read directly by `esfScaleFactor`, so a rerun replaces the production SFs.

**1c. Data pileup profiles**, from `test/`
```bash
bash PUfiles/makeDataPUDistributions.sh
```
- Output: `test/PUfiles/dataPU_<era>{,_PUUp,_PUDown}.root`, read as `../../PUfiles` by the skims.

### Step 2. Skims

**2a. Data signal/control trees**, folder `test/dataSkim/`
- `analysis_eg.C` (eγ, DoubleEG) and `analysis_mg.C` (μγ, MuonEG) select events from the ntuples of one run period
  and write `resTree_{egsignal_DoubleEG,mgsignal_MuonEG}_<year><run>.root` to `/eos/uscms/store/user/tmishra/eg_mg_treesData/`.
```bash
MAX_JOBS=8 bash run_egmg_parallel.sh     # compiles both macros, runs all years and run periods
bash merge_egmg_years.sh                 # run periods -> per-era files read by all later steps (DRYRUN=1, BACKUP=1)
# one run period by hand: ./analysis_eg.exe <year> <run>, e.g. ./analysis_eg.exe 2018 A
```
- Output: `/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_{egsignal_DoubleEG,mgsignal_MuonEG}_<era>.root`.
- `merge_egmg_years.sh` overwrites these per-era files; use `DRYRUN=1` to check first and `BACKUP=1` to keep the old ones.

**2b. Simulation skims**

| Sample | Folder | Commands | Output |
|---|---|---|---|
| Vγ MC (condor) | `test/VGamma/condor` | build `../src/analysis_VGamma_Run2_condor.exe`; `python3 submit_skim.py --tag <tag>`, then `bash merge_skim.sh <tag>` | `/eos/.../VGamma/<tag>/`, then copied by hand to `/eos/.../VGamma/` |
| Signal trees | `test/myResult/Signal-samples` | `bash run_Signal.sh` | `/uscms_data/d3/tmishra/signal_trees/resTree_<model>_<era>.root` |
| MC for the closure tests | `test/Closure` | `bash Run_MC_all.sh` | `/eos/uscms/store/user/tmishra/{egMC,mgMC}/` |

### Step 3. e→γ misidentification rate, folder `test/eleFakePho/`
```bash
cd src && bash run_elefakepho.sh && bash run_bgTemplate.sh   # TnP trees and background templates
#   then hadd the data eras per period: plot_elefakepho_DataTnP_dR05_<era>.root
cd ../fit
bash RunFits.sh "2016preVFP 2016postVFP 2017 2018" "DY Data" 1    # about 660 template fits per era (condor)
bash processEleFake.sh <era>                                        # checkFits.py, makeTxt.py, runFakeFunc.sh
root -l -b -q 'fitFakeFunc_pT_nVtx.C++(<year>,<preVFP>,true)'      # DY function, needed by the closure test (step 9)
```
- Output: `/eos/uscms/store/user/tmishra/elefakepho/{DATAResult,DrellYanResult}<era>/`.
- `checkFits.py` uses the time stamps in `fit/condorSubFile/`, so keep that folder.

### Step 4. jet→γ misidentification rate, folder `test/jetFakePho/`
```bash
cd src && bash RunHad.sh && bash Run_mgGJet.sh     # hadron templates (Run_mgGJet.sh is for step 5)
cd ../fit
rm -f /eos/uscms/store/user/tmishra/jetfakepho/txt*/JetFakeRate-{DoubleEG,MuonEG}-EB_pT_nVtx.txt   # the fit appends to these
bash run_FitJetFake_pT_nVtx.sh && bash run_fitJetFunc.sh
```
- Output: `/eos/uscms/store/user/tmishra/jetfakepho/txt<era>/JetFakeRate-transferfactor-*-EB_pT_nVtx.txt`.

### Step 5. ISR weights for Vγ, folder `test/ISRReweight/`
```bash
nohup bash run_all_ISR.sh >& /eos/uscms/store/user/tmishra/ISRweighting/logs_ISR/run_all_ISR.log &
```
- Four sub-steps: skims, jet fake rate in the ISR region (about 7 h), weights per era, `format_isrweights.py`.
- Paste the printed blocks into `test/Background/analysis_VGBkg.C` and `test/myResult/pred_VGBkg.C`, then redo
  steps 6–8.
- Details: `test/ISRReweight/README.md`.

### Step 6. Vγ and misidentified-lepton normalisation
```bash
cd test/Background        # set ToDeriveScale=true (other flags false) at the top of runBkg.sh
bash runBkg.sh && bash merge_controlTree_Run2.sh     # control trees -> /uscms_data/d3/tmishra/Background/
cd ../VGamma/fit
bash run_fit_ALL_HT.sh    # photon-pT binned SFs: set doPhoScan=true and RunYears at the top
bash run_fit_ALL.sh       # inclusive SFs (1000 toys) -> /eos/.../VGamma/VGamma_scalefactor_{eg,mg}_<era>.txt
bash RunDependence.sh     # -> /eos/.../VGamma/<era>/VGammaScaleHist_{eg,mg}_<era>.root, used by steps 7 and 8
```
- The inclusive fallback factors in `include/analysis_commoncode.h` are updated by hand.

### Step 7. Control- and validation-region predictions, folder `test/Background/`
```bash
bash runBkg.sh            # with makeCRPlots=true or makeVRPlots=true at the top (4 eras × 2 channels)
```
- Output: `/uscms_data/d3/tmishra/Background/`.
- The jobs run in the background; check `logs/` for failures.

### Step 8. Signal-region prediction, folder `test/myResult/`
```bash
bash runLimit.sh          # SR predictions for all eras, signal trees, factorise_rareBkg.py, event counts
bash makePlot.sh          # SR plots and the Run 2 (_678) merges needed by the Run 2 cards
```
- Output: `/uscms_data/d3/tmishra/Output/signalTree_{egamma,mg}_<bkg>_<era>.root`,
  `/uscms_data/d3/tmishra/signal_trees/signalTree_<model>_<era>.root`, `logs/SignalRegionTable_<era>.tex`,
  `SR_counts_SUS-24-011_<era>.txt`.
- `runLimit.sh` makes the predictions only; the limits are step 11.

### Step 9. Closure tests, folder `test/Closure/`
```bash
bash run_closure_efakepho.sh     # e→γ closure for all eras (needs the DY function from step 3)
```

### Step 10. Datacards, folders `limits/datacards/` and `limits/combine/`
```bash
bash limits/datacards/makeCards.sh             # set years/bins/modes at the top -> /eos/.../CombinedLimit/<model>/cards/cards_<era>_<mode>/
python3 limits/datacards/convert_all_cards.py --dry-run   # then without --dry-run: CMS nuisance names, luminosity scheme
cd limits/combine && ./combine_cards_Run2.sh T5WG          # and T6WG -> cards_Run2_unblind/
```
- `createDataCard.py` and `write*card.py` work inside `test/myResult` (they read its `logs/`) from wherever they are
  started.
- TChiWG: `limits/TChiWG/combine_cards_Run2.sh`.

### Step 11. Limits (condor), folder `limits/combine/`
```bash
mkdir -p T5WG_run_<date> && cp runLimit.py runLimit.sh T5WG_run_<date>/ && cd T5WG_run_<date>
python3 runLimit.py T5WG unblind     # T6WG the same; TChiWG with runLimitTChiWG.py / runLimitTChiWG.sh
```
- Output: `/eos/uscms/store/user/tmishra/CombinedLimit/<model>/limits/limits_Run2_<mode>/`.

### Step 12. Fits and checks, folder `limits/fits/`
```bash
bash run_postfit.sh                          # FitDiagnostics per era and for Run 2 -> postfit_unblind_T5WG/
bash runImpacts_Run2_unblind.sh              # impacts and goodness of fit (tasks and mass points at the top)
bash phase2_2026-10-31/scan_r_1700_200.sh    # 1D likelihood scan
(cd masked_rerun && python3 submit.py T5WG)  # masked-bin limits (condor)
```

### Step 13. Plots, folders `limits/plots/` and `limits/TChiWG/`
```bash
cd limits/plots
python3 makeT5WGLimit.py Run2                # and makeT6WGLimit.py -> /uscms_data/d3/tmishra/Output/T5WG_Run2_2016-17-2018.pdf
bash postfit/run_plot_eventct_postFit.sh     # post-fit SR plots and SR_counts (reads ../fits/postfit_unblind_T5WG)
python3 postfit/compare_prefit_postfit.py    # pre-fit vs post-fit per bin
cd ../TChiWG && bash makelimit.sh && root -l -b -q plot_limit.C
```
- These overwrite the figures used in the paper and the AN.

### Step 14. Figures into the AN, folder `test/notes/CopyingToAN/`
```bash
bash copy_to_AN.sh --dry-run     # then without --dry-run
```
- Copies the figures listed in `an_figures_manifest.txt` to the AN area on EOS and on lxplus.

## Where outputs are stored

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
