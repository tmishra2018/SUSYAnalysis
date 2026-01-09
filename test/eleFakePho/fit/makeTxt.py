# python3 makeTxt.py --Year 2017 --Type DY 
# python3 makeTxt.py --Year 2018 --Type DY 
# python3 makeTxt.py --Year 2016postVFP --Type DY

import os
import argparse
import re
from glob import glob

parser = argparse.ArgumentParser(description='Merge electron fake rate txt files by year and type')
parser.add_argument('--Year', type=str, required=True, help="The run year, e.g. 2016postVFP")
parser.add_argument('--Type', type=str, required=True, help="Data type: DATA, DY, etc.")

args = parser.parse_args()
Year = args.Year
Type = args.Type

input_dir = f"/eos/uscms/store/user/tmishra/elefakepho/{Type}Result{Year}/"
out_dir = f"{Type}Result{Year}"
os.makedirs(out_dir, exist_ok=True)

patterns = [
    "Bw-expo-pt", "Bw-ker-pt", "DY-ker-pt",
    "Bw-expo-eta", "Bw-ker-eta", "DY-ker-eta",
    "Bw-expo-vtx", "Bw-ker-vtx", "DY-ker-vtx"
]

def extract_axis_bin_and_type(line):
    """Extracts ('den' or 'num'), bin value (float), and axis type ('pt', 'eta', 'vtx') from the line."""
    match = re.search(r'(pt|eta|vtx) (den|num) (\d+(?:\.\d+)?)', line)
    if match:
        axis = match.group(1)
        which = match.group(2)
        bin_val = float(match.group(3))
        return axis, which, bin_val, line.strip()
    return None, None, float('inf'), line.strip()

for pattern in patterns:
    files = glob(os.path.join(input_dir, f"*{pattern}*.txt"))

    den_lines_set = set()
    num_lines_set = set()

    for f_path in files:
        try:
            with open(f_path, 'r') as f:
                raw_line = f.readline().strip()
                axis, which, bin_val, formatted_line = extract_axis_bin_and_type(raw_line)
                if which == "den":
                    den_lines_set.add((bin_val, formatted_line))
                elif which == "num":
                    num_lines_set.add((bin_val, formatted_line))
        except Exception as e:
            print(f"Warning: Failed to read {f_path}: {e}")

    # Sort lines by bin value
    sorted_den_lines = sorted(den_lines_set, key=lambda x: x[0])
    sorted_num_lines = sorted(num_lines_set, key=lambda x: x[0])

    # Write to output file
    output_file = os.path.join(out_dir, f"EleFakeRate-{Type}-{pattern}-60-120.txt")
    if os.path.exists(output_file):
        os.remove(output_file)
    with open(output_file, "w") as out_f:
        for _, line in sorted_den_lines:
            out_f.write(line + "\n")
        for _, line in sorted_num_lines:
            out_f.write(line + "\n")

# Move to EOS directory
if Type.lower() == "data":
    os.system(f"mv {out_dir}/* /eos/uscms/store/user/tmishra/elefakepho/DATAResult{Year}/")
elif Type.lower() in ["drell", "dy"]:
    os.system(f"mv {out_dir}/* /eos/uscms/store/user/tmishra/elefakepho/DrellYanResult{Year}/")

os.rmdir(out_dir)
