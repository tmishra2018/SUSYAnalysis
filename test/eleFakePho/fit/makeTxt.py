# Examples:
# python3 makeTxt.py --Year 2017 --Type DY --Range pT-gt80
# python3 makeTxt.py --Year 2018 --Type DY --Range pT30-50
# python3 makeTxt.py --Year 2018 --Type DY --Range pT50-80
# python3 makeTxt.py --Year 2018 --Type DATA --Range nVtx-gt30
# python3 makeTxt.py --Year 2018 --Type DATA --Range nVtx20-30
# python3 makeTxt.py --Year 2018 --Type DATA --Range nVtx8-20
# python3 makeTxt.py --Year 2016postVFP --Type DY --Range pT-gt80
# python3 makeTxt.py --Year 2018 --Type DATA                      (main folder: pt, eta, vtx)
# python3 makeTxt.py --Year 2018 --Type DATA --Since "2026-10-19 12:07"

import os
import sys
import argparse
import datetime
import re
from glob import glob

# Valid range suffixes -> which folder "family" they belong to.
PT_RANGES = ["pT30-50", "pT50-80", "pT-gt80"]
NVTX_RANGES = ["nVtx-gt30", "nVtx20-30", "nVtx8-20"]

# For folders binned in pT, the fits performed inside them are vs vtx.
# For folders binned in nVtx, the fits performed inside them are vs pt.
PATTERNS_FOR_PT_FOLDER = ["Bw-expo-vtx", "Bw-ker-vtx", "DY-ker-vtx"]
PATTERNS_FOR_NVTX_FOLDER = ["Bw-expo-pt", "Bw-ker-pt", "DY-ker-pt"]
# Without --Range the main folder is merged, with the nominal pt, eta and vtx fits.
PATTERNS_FOR_MAIN_FOLDER = [f"{fit}-{axis}" for axis in ["pt", "eta", "vtx"]
                            for fit in ["Bw-expo", "Bw-ker", "DY-ker"]]

parser = argparse.ArgumentParser(description='Merge electron fake rate txt files by year, type, and pT/nVtx range')
parser.add_argument('--Year', type=str, required=True, help="The run year, e.g. 2016postVFP")
parser.add_argument('--Type', type=str, required=True, help="Data type: DATA, DY, etc.")
parser.add_argument('--Range', type=str, default=None, choices=PT_RANGES + NVTX_RANGES,
                     help=f"Binning folder suffix. One of: {', '.join(PT_RANGES + NVTX_RANGES)}. "
                          "Omit it to merge the main folder.")
parser.add_argument('--Since', type=str, default=None,
                     help="Only merge files written after this time, \"YYYY-mm-dd HH:MM\" "
                          "(the batch submit time), to skip leftovers from older runs")

args = parser.parse_args()
Year = args.Year
Type = args.Type
Range = args.Range
since = datetime.datetime.strptime(args.Since, "%Y-%m-%d %H:%M").timestamp() if args.Since else None

if Range is None:
    patterns = PATTERNS_FOR_MAIN_FOLDER
elif Range in PT_RANGES:
    patterns = PATTERNS_FOR_PT_FOLDER
elif Range in NVTX_RANGES:
    patterns = PATTERNS_FOR_NVTX_FOLDER
else:
    # argparse's choices already guards this, but keep a defensive fallback.
    raise ValueError(f"Unrecognized --Range value: {Range}")

folder = f"{Type}Result{Year}" + (f"_{Range}" if Range else "")
input_dir = f"/eos/uscms/store/user/tmishra/elefakepho/{folder}/"
out_dir = f"{folder}_tmp"
if not os.path.isdir(input_dir):
    sys.exit(f"No input folder {input_dir}")
os.makedirs(out_dir, exist_ok=True)


def extract_axis_bin_and_type(line):
    """Extracts ('den' or 'num'), bin value (float), and axis type ('pt', 'eta', 'vtx') from the line."""
    match = re.search(r'(pt|eta|vtx|ht) (den|num) (\d+(?:\.\d+)?)', line)
    if match:
        axis = match.group(1)
        which = match.group(2)
        bin_val = float(match.group(3))
        return axis, which, bin_val, line.strip()
    return None, None, float('inf'), line.strip()


def fit_failed(line):
    """A fit failed if its toy mean is far from the fitted signal, its spread is huge, or a value is nan/inf.
    Columns: axis den/num bin signal signalErr mean rms"""
    try:
        sig, sigerr, mean, rms = [float(x) for x in line.split()[3:7]]
    except (ValueError, IndexError):
        return True
    if any(x != x or abs(x) == float('inf') for x in (sig, sigerr, mean, rms)):
        return True
    return sig <= 0 or not 0.9 < mean / sig < 1.1 or rms > 0.2 * abs(mean)


# merged[pattern] = {"den": {bin: line}, "num": {bin: line}}
merged = {}
for pattern in patterns:
    # Only the per-bin FitKer outputs (...-<pattern>-den-... / -num-...), never a merged file
    files = (glob(os.path.join(input_dir, f"*-{pattern}-den-*.txt")) +
             glob(os.path.join(input_dir, f"*-{pattern}-num-*.txt")))
    if since is not None:
        files = [f for f in files if os.path.getmtime(f) >= since]

    bins = {"den": {}, "num": {}}
    for f_path in sorted(files, key=os.path.getmtime):  # newest file wins for a repeated bin
        try:
            with open(f_path, 'r') as f:
                raw_line = f.readline().strip()
                axis, which, bin_val, formatted_line = extract_axis_bin_and_type(raw_line)
                if which in bins:
                    bins[which][bin_val] = formatted_line
        except Exception as e:
            print(f"Warning: Failed to read {f_path}: {e}")
    merged[pattern] = bins

def rel_spread(line):
    cols = line.split()
    return float(cols[6]) / abs(float(cols[5]))


# A failed or missing Bw-ker fit takes the DY-ker or Bw-expo value of the same bin: whichever of
# the two did not fail, and if both are fine the one with the smaller relative spread (rms/mean).
# If both alternatives failed too, the Bw-ker value is kept and copied to both of them.
for pattern in patterns:
    if not pattern.startswith("Bw-ker-"):
        continue
    alternatives = [alt for alt in (pattern.replace("Bw-ker-", "DY-ker-"), pattern.replace("Bw-ker-", "Bw-expo-"))
                    if alt in merged]
    for which in ("den", "num"):
        all_bins = set(merged[pattern][which])
        for alt in alternatives:
            all_bins |= set(merged[alt][which])
        for bin_val in sorted(all_bins):
            ker_line = merged[pattern][which].get(bin_val)
            if ker_line is not None and not fit_failed(ker_line):
                continue
            good = [(rel_spread(merged[alt][which][bin_val]), i, alt) for i, alt in enumerate(alternatives)
                    if bin_val in merged[alt][which] and not fit_failed(merged[alt][which][bin_val])]
            state = 'missing' if ker_line is None else 'failed (' + ker_line + ')'
            if not good:
                if ker_line is None:
                    print(f"WARNING: {folder} {pattern} {which} {bin_val:g}: missing, and no good DY-ker/Bw-expo fit to use")
                    continue
                # Keep the Bw-ker value and copy it to both alternatives, so this bin gets no fit-function systematic
                print(f"{folder} {pattern} {which} {bin_val:g}: {state}, no good DY-ker/Bw-expo fit "
                      f"-> Bw-ker value kept and copied to {', '.join(alternatives)}")
                for alt in alternatives:
                    merged[alt][which][bin_val] = ker_line
                continue
            best = min(good)[2]
            print(f"{folder} {pattern} {which} {bin_val:g}: {state} -> {best.split('-')[0]}-{best.split('-')[1]} value")
            merged[pattern][which][bin_val] = merged[best][which][bin_val]

# A failed or missing Bw-expo fit takes the Bw-ker value of the same bin, so it adds no
# fit-function systematic and the three files the macro reads in parallel stay aligned.
for pattern in patterns:
    if not pattern.startswith("Bw-expo-"):
        continue
    ker = merged.get(pattern.replace("Bw-expo-", "Bw-ker-"))
    if ker is None:
        continue
    for which in ("den", "num"):
        for bin_val, ker_line in ker[which].items():
            expo_line = merged[pattern][which].get(bin_val)
            if expo_line is None or fit_failed(expo_line):
                print(f"{folder} {pattern} {which} {bin_val:g}: "
                      f"{'missing' if expo_line is None else 'failed (' + expo_line + ')'} -> Bw-ker value")
                merged[pattern][which][bin_val] = ker_line

for pattern in patterns:
    sorted_den_lines = sorted(merged[pattern]["den"].items())
    sorted_num_lines = sorted(merged[pattern]["num"].items())

    # Write to output file
    if len(sorted_den_lines) != len(sorted_num_lines) or not sorted_den_lines:
        print(f"WARNING: {folder} {pattern}: {len(sorted_den_lines)} den bins, {len(sorted_num_lines)} num bins")
    else:
        print(f"{folder} {pattern}: {len(sorted_den_lines)} bins")

    if Range:
        output_file = os.path.join(out_dir, f"EleFakeRate-{Type}-{Range}-{pattern}-60-120.txt")
    else:
        output_file = os.path.join(out_dir, f"EleFakeRate-{Type}-{pattern}-60-120.txt")
    if os.path.exists(output_file):
        os.remove(output_file)
    with open(output_file, "w") as out_f:
        for _, line in sorted_den_lines:
            out_f.write(line + "\n")
        for _, line in sorted_num_lines:
            out_f.write(line + "\n")

# Move merged files back into the matching source folder
# (e.g. DATAResult2018_pT30-50/, DYResult2018_nVtx-gt30/, ...)
os.system(f"mv {out_dir}/* {input_dir}")
os.rmdir(out_dir)
