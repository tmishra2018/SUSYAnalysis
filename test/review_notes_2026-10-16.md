# Analysis review notes (2026-10-16)

Nothing below has been changed in your code unless marked DONE.

## A. Data ntuples (/eos/uscms/store/group/lpcsusyphotons/SoftPhoton/Tribeni/)
1. SingleElectron_2016B-ver2.root and SingleElectron_2016C.root contain SingleMuon events
   (same first run/event as SingleMuon files, ~99% events with a muon). Real SingleElectron
   2016B-ver2 (246.4M) / 2016C (97.3M) never ntuplized -> remake on lxplus.
   Used by: Trigger/analysis_egTrigger.C, eleFakePho/src/analysis_elefakepho.C
2. SingleMuon_2016D only 56% complete (54.8M of 98.0M DAS events) -> remake.
3. Partially processed (95-99.9%): SingleElectron 2016D/F_pre/G/H; SingleMuon 2016B-ver2/C/E/F_pre/G/H,
   2017C, 2018D; MET 2016postVFP/2017/2018. Lumi for these must come from processed lumis AND golden JSON.
   DoubleEG and MuonEG: 100% complete (fine).

## B. Analysis code
4. DONE: golden JSON filter. loadGolden/passGolden in include/analysis_tools.h,
   applied in analysis_eg.C and analysis_mg.C (JSONs in test/golden/). Other data macros not yet.
   Lumi values in analysis_tools.h are correct (official normtag_BRIL golden lumis).
5. include/analysis_scalefactor.h:94-99 loads pre-UL 2016 SFs from /uscms/home/mengleis/... for all years;
   src/analysis_scalefactor.cc:34 multiplies electron SF by 2016 R9 SF.
   New per-year maps: test/ScaleFactor/tnp/results/<year>/SF_{miniIso,R9}_<year>_hybrid.root (EGamma_SF2D).
6. PU weights: getPUESF16/17/18(nVtx) tables (analysis_tools.h:517) used in VGamma/src/analysis_DY.C:332-335,
   Reweight/analysis_ISRMC.C:226-229. Official truePU weights exist in
   PileUp_with_official_recipe/PUfiles/PUWeights_<year>.root but are unused. Compare before switching.
7. R9 cuts inconsistent:
   - analysis_eg.C:408 photon cut R9 < R9EBCut (0.5) for both EB and EE; electrons use 0.5/0.8 (line 472).
   - analysis_mg.C:410, :489 hard-code 0.5/0.8 instead of R9EBCut/R9EECut.
   - Electrons use eleR9 (not full5x5), photons phoR9Full5x5. New electron SFs match eleR9;
     switching to full5x5 would need SF re-derivation.
8. analysis_eg.C:367 event-level trigger `//if(!raw.passHLT())continue;` commented out (output named _noHLT);
   only per-object passHLTSelection applied. Intentional?
9. isMC check (analysis_eg.C:18-21) only knows 2016 MC RunTypes; harmless in data macros, wrong if copied to MC macros.
10. Trigger/analysis_R9.C: hard-coded 2016 trigger bits, first 1M events only, DY path in store/user.

## C. Electron SFs (egamma-tnp), all years done
- miniIso SF ~1 (within ~1-2%) everywhere.
- R9 SF ~1 in barrel; endcap 1.566<|eta|<2.0 below 1 (2017: 0.965-0.982, 2018 ~0.99).
- High-pT bins (>100 GeV) use cut-and-count; 200-500 GeV error = 0 -> assign flat 1-2%.
- Fit plots not yet inspected (2016 endcap R9 0.950 / 1.04 bins especially).
- 2016preVFP SFs use SingleElectron B-ver1, D, E, F_pre only (because of item 1).
- Tools: /uscms_data/d3/tmishra/tnp_tools; work: /uscms_data/d3/tmishra/tnp_work; script test/ScaleFactor/tnp/run_sf.sh

## Decisions needed
(a) remake ntuples on lxplus (items 1, 2)
(b) switch analysis_scalefactor.h to new per-year SFs (item 5)
(c) fix R9 cuts / trigger check (items 7, 8)
(d) PU method (item 6)
