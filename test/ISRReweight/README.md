# ISR reweighting (`test/ISRReweight`)

`bash run_all_ISR.sh [first] [last]` runs four steps and stops at the first failure:
1. `run_Reweight_step1.sh`: builds and runs `analysis_ISR` (data), `analysis_ISRMC` (MC) and
   `analysis_Hadron` (jet→γ proxy), then merges the outputs with hadd. Uses `MAX_JOBS` parallel jobs.
2. `runFit_step2.sh`: `FitJetFake.C` template fits per channel × EB/EE × era × photon-pT bin, written to
   `JetFakeRate-ISR-{EB,EE}_{mg,egamma}_<era>.txt` on EOS (about 6–7 h).
3. `run_ISRweight_step3.sh`: `plot_ISRweight.exe` derives the weights and then applies them as a check;
   it writes `ISRConfig.txt` for each era.
4. `format_isrweights.py` writes `isr_blocks.txt`, to paste into `Background/analysis_VGBkg.C` and
   `myResult/pred_VGBkg.C`. After that, re-derive the VGamma scale (`ToDeriveScale=true` in `runBkg.sh`) before `pred_VGBkg.C`.

Logs go to `/eos/uscms/store/user/tmishra/ISRweighting/logs_ISR/`.
The method follows the 2016 analysis (`/uscms/homes/m/mengleis/work/SUSY2016/SUSYAnalysis/test/Reweight/`,
another user's area, read-only). It uses EB+EE photons, a fixed 0.03 fake rate for EE photons with
150 ≤ pT < 250, and no fake subtraction above 250 GeV. A fake rate that is ≤ 0, ≥ 1 or NaN falls back to 0.03.
As in 2016, the data/Zγ ISR weights are applied to all Vγ MC (Wγ, Zγ, DY); the LO→NLO factor
(`LOtoNLOFactor`) is applied to Wγ (mcType 1–3) only.
