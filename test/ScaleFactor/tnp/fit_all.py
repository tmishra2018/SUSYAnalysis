# Fit every (eta, pt) bin of the prepared fitter inputs with egamma-tnp's run_fitter and collect data/MC SFs.
# Run inside the egamma-tnp container with PYTHONPATH pointing at the mamaxout fork (needs numba_stats).
# python3 fit_all.py --data data.root --mc mc.root --binning binning.pkl --outdir <dir> [--workers 4]
import argparse
import json
import math
import pickle
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import uproot

FIT_TYPES = ("dcb_cms", "dcb_exp")  # nominal, alternative (fit-model systematic)


def write_config(outdir, bin_name, fit_type, data, mc):
    cfg = {
        "info_level": "INFO",
        "mass": "Z",
        "input": {"root_files_DATA": {"DATA": str(data)}, "root_files_MC": {"MC": str(mc)}},
        "fit": {
            "bin_ranges": [[0, 1]],
            "bin": ["bin0"],
            "fit_type": fit_type,
            "use_cdf": False,
            "sigmoid_eff": False,
            "interactive": False,
            "x_min": 60,
            "x_max": 120,
            "hist_pass_name": f"{bin_name}_Pass",
            "hist_fail_name": f"{bin_name}_Fail",
        },
        "output": {
            "plot_dir": str(outdir / "plots" / fit_type / bin_name),
            "results_file": str(outdir / "results" / fit_type / f"{bin_name}.json"),
        },
        "scale_factors": {"data_mc_pair": {"SF": ["DATA", "MC"]}},
    }
    path = outdir / "configs" / fit_type / f"{bin_name}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    Path(cfg["output"]["results_file"]).parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cfg, indent=1))
    return path


def run_fit(cfg_path):
    log = cfg_path.with_suffix(".log")
    with open(log, "w") as f:
        ret = subprocess.run(["python3", "-m", "egamma_tnp.scripts.run_fitter", "--config", str(cfg_path)], stdout=f, stderr=subprocess.STDOUT)
    return cfg_path, ret.returncode


def read_eff(results_file):
    eff = {}
    try:
        for entry in json.load(open(results_file)):
            eff[entry["sample_type"]] = (entry["epsilon"], entry["epsilon_err"], entry.get("fit_quality_status", "UNKNOWN"))
    except (OSError, ValueError, KeyError):
        pass
    return eff


def ratio(d, de, m, me):
    if not m:
        return float("nan"), float("nan")
    return d / m, math.hypot(de / m, d * me / m**2)


def cut_and_count(root_file, bin_name, lo=81.0, hi=101.0):
    f = uproot.open(root_file)
    vals = []
    for cat in ("Pass", "Fail"):
        h = f[f"{bin_name}_{cat}"]
        centers = h.axis().centers()
        mask = (centers > lo) & (centers < hi)
        vals.append((h.values()[mask].sum(), (h.variances()[mask]).sum()))
    (p, pv), (fl, fv) = vals
    tot = p + fl
    if tot <= 0:
        return float("nan"), float("nan")
    eff = p / tot
    err = math.sqrt(fl**2 * pv + p**2 * fv) / tot**2
    return eff, err


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--mc", required=True)
    ap.add_argument("--binning", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()

    outdir = Path(args.outdir).resolve()
    bins = pickle.load(open(args.binning, "rb"))["bins"]
    jobs = [write_config(outdir, b["name"], ft, Path(args.data).resolve(), Path(args.mc).resolve()) for b in bins for ft in FIT_TYPES]
    with ThreadPoolExecutor(args.workers) as pool:
        for cfg, rc in pool.map(run_fit, jobs):
            if rc != 0:
                print(f"FIT FAILED (rc={rc}): {cfg}")

    rows = []
    for b in bins:
        name = b["name"]
        row = {
            "bin": name,
            "eta": [float(b["vars"]["el_eta"]["min"]), float(b["vars"]["el_eta"]["max"])],
            "pt": [float(b["vars"]["el_pt"]["min"]), float(b["vars"]["el_pt"]["max"])],
        }
        for ft in FIT_TYPES:
            eff = read_eff(outdir / "results" / ft / f"{name}.json")
            if "DATA" in eff and "MC" in eff:
                sf, sferr = ratio(eff["DATA"][0], eff["DATA"][1], eff["MC"][0], eff["MC"][1])
                row[ft] = {"eff_data": eff["DATA"][0], "eff_mc": eff["MC"][0], "sf": sf, "sf_stat": sferr,
                           "quality_data": eff["DATA"][2], "quality_mc": eff["MC"][2]}
        cd, cde = cut_and_count(args.data, name)
        cm, cme = cut_and_count(args.mc, name)
        csf, csfe = ratio(cd, cde, cm, cme)
        row["cut_and_count"] = {"eff_data": cd, "eff_mc": cm, "sf": csf, "sf_stat": csfe}
        nom, alt = row.get(FIT_TYPES[0]), row.get(FIT_TYPES[1])
        if nom:
            syst = abs(alt["sf"] - nom["sf"]) if alt else 0.0
            row["sf"] = nom["sf"]
            row["sf_err"] = math.hypot(nom["sf_stat"], syst)
            row["sf_syst_fitmodel"] = syst
        rows.append(row)

    (outdir / "sf_summary.json").write_text(json.dumps(rows, indent=1))
    with open(outdir / "sf_table.txt", "w") as f:
        f.write(f"{'eta':>18} {'pt':>12} | {'effD':>6} {'effMC':>6} {'SF':>6} {'err':>6} {'qD/qMC':>22} | {'CnC SF':>6}\n")
        for r in rows:
            nom = r.get(FIT_TYPES[0], {})
            f.write(f"[{r['eta'][0]:7.3f},{r['eta'][1]:7.3f}] [{r['pt'][0]:4.0f},{r['pt'][1]:4.0f}] | "
                    f"{nom.get('eff_data', float('nan')):6.3f} {nom.get('eff_mc', float('nan')):6.3f} "
                    f"{r.get('sf', float('nan')):6.3f} {r.get('sf_err', float('nan')):6.3f} "
                    f"{nom.get('quality_data', '-') + '/' + nom.get('quality_mc', '-'):>22} | "
                    f"{r['cut_and_count']['sf']:6.3f}\n")
    print((outdir / "sf_table.txt").read_text())


if __name__ == "__main__":
    main()
