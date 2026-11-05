# test/: analysis chain for AN-21-007 (eγ / μγ)

How to run every step is in `../README.md`. Set up the environment and build `lib/` first, and run every script from
its own folder.

## Order of the chain

```
ScaleFactor/  Trigger/  PUfiles/          object SFs (TnP), trigger SFs, data pileup profiles
        │
dataSkim/  VGamma/condor  myResult/Signal-samples  Closure/(MC)      skims: data trees, Vγ MC, signal, closure MC
        │
eleFakePho/  →  jetFakePho/               e→γ and jet→γ misidentification rates (src/ skims, fit/ fits)
        │
ISRReweight/                              ISR weights for Vγ (paste into Background/ and myResult/ VGBkg macros)
        │
Background/ (SCALE)  →  VGamma/fit        Vγ and misidentified-lepton normalisation (photon-pT binned SFs)
        │
Background/ (CR, VR)                      control- and validation-region predictions
        │
myResult/                                 signal-region predictions, signal trees, SR tables and plots
        │
../limits/                                datacards, Run 2 combination, limits, fits, exclusion plots
```

## Folders

| Folder | What it holds |
|---|---|
| `ScaleFactor/` | TnP trees and SFs. `tnp/results/` is production: `esfScaleFactor` reads the SF files from there |
| `Trigger/` | trigger ntuples (`rerun_trigger_data.sh`, `Run-mgTrig-mc.sh`) and trigger SFs (`runEGTrigger.sh`, `runMuonTrigger.sh`) |
| `PUfiles/` | data pileup profiles in use (`dataPU_*.root`), read as `../../PUfiles` by the skims |
| `golden/` | golden JSONs (read by `include/analysis_tools.h`) and the golden-filter effect check |
| `cross/` | signal cross sections (`ULsusyCrossSection.root` is read by the signal macros) |
| `dataSkim/` | data signal/control trees: `analysis_{eg,mg}.C`, `run_egmg_parallel.sh`, `merge_egmg_years.sh` |
| `eleFakePho/` | e→γ rate. `fit/condorSubFile/` is still used by `checkFits.py` |
| `jetFakePho/` | jet→γ rate |
| `VGamma/` | Vγ skims (`condor/`, `src/`) and normalisation fits (`fit/`) |
| `ISRReweight/` | ISR weights (`run_all_ISR.sh`; see `ISRReweight/README.md`) |
| `Background/` | CR/VR predictions and the SCALE trees for the Vγ fit; `BkgPredConfig.txt` selects channel/era/cuts |
| `myResult/` | SR predictions (`runLimit.sh`), SR plots (`makePlot.sh`), `SigConfig.txt`, `binConfig.txt`; signal trees in `Signal-samples/` |
| `Closure/` | closure tests of the misidentification methods |
| `studies/` | side studies: `HEM/`, `dR/`, `pileup_official/` (PU POG recipe; its `PUfiles/PUWeights_*` feed the TnP pileup JSON), `MCStudy/`, `QCD/`, `DoubleFake/`, `EndCap/`, `Supplementary/` |
| `tools/` | helpers: `CommonTools/` (tables, comparisons), `get_nevents.sh`, `compare_root_shapes.py` |
| `notes/` | review notes, `CHANGES_*.txt`, method-change slides, `CopyingToAN/` (copies the AN figures) |
| `legacy/` | old code kept for reference: `Result/` (2016-style chain), `PostFit_2016/`, `macro_old/`, Sep-2026 top-level macros. Some of it reads inputs that no longer exist |

## Archived output

Old plots, logs and backups are archived as one `.tgz` per folder in
`/eos/uscms/store/user/tmishra/SUSYAnalysis_out/` (`test/<same path>/`, `limits/`, `toplevel_unused_dirs_2026-11-04.tgz`).
To restore one: `xrdcp root://cmseos.fnal.gov//store/user/tmishra/SUSYAnalysis_out/<path>.tgz . && tar xzf <file>.tgz`.
