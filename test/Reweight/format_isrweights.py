#!/usr/bin/env python3
"""Format ISRweights_*.txt files into C++ if/else blocks.

Expected input file format:
  # normalization <value>
  # ptLow ptHigh weight
  <ptLow> <ptHigh> <weight>
  ...
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Tuple


def parse_isr_file(path: Path) -> Tuple[float, List[Tuple[float, float, float]]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")

    normalization = None
    bins: List[Tuple[float, float, float]] = []

    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            if line.startswith("# normalization"):
                parts = line.split()
                if len(parts) >= 3:
                    normalization = float(parts[2])
            continue

        parts = line.split()
        if len(parts) < 3:
            continue
        low, high, weight = float(parts[0]), float(parts[1]), float(parts[2])
        bins.append((low, high, weight))

    if normalization is None:
        raise ValueError(f"No '# normalization <value>' found in {path}")
    if len(bins) == 0:
        raise ValueError(f"No ISR bins found in {path}")

    return normalization, bins


def fmt_num(x: float) -> str:
    return f"{x:.6g}"


def format_block(run_year: int, pre_vfp: int | None, normalization: float, bins: List[Tuple[float, float, float]], is_first: bool) -> str:
    if run_year == 2016:
        head = f"{'if' if is_first else 'else if'}(RunYear==2016 && preVFP=={pre_vfp}){{"
    else:
        head = f"{'if' if is_first else 'else if'} (RunYear=={run_year}){{"

    lines = [head]
    for i, (low, high, w) in enumerate(bins):
        if i == 0:
            lines.append(f"\tif(ISRJetPt < {fmt_num(high)})reweightF = {fmt_num(w)};")
        elif i < len(bins) - 1:
            lines.append(
                f"\telse if(ISRJetPt >= {fmt_num(low)} && ISRJetPt < {fmt_num(high)})reweightF  = {fmt_num(w)};"
            )
        else:
            lines.append(f"\telse if(ISRJetPt >= {fmt_num(low)})reweightF = {fmt_num(w)};")

    lines.append(f"\tNormalization = {fmt_num(normalization)};    }}")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="Format ISRweights text files into C++ if/else blocks")
    parser.add_argument("--dir", default=".", help="Directory containing ISRweights_*.txt files")
    args = parser.parse_args()

    base = Path(args.dir)
    mapping = [
        (2016, 1, base / "ISRweights_2016preVFP.txt"),
        (2016, 0, base / "ISRweights_2016postVFP.txt"),
        (2017, None, base / "ISRweights_2017.txt"),
        (2018, None, base / "ISRweights_2018.txt"),
    ]

    blocks = []
    for i, (year, prevfp, p) in enumerate(mapping):
        norm, bins = parse_isr_file(p)
        blocks.append(format_block(year, prevfp, norm, bins, is_first=(i == 0)))

    print("\n\n".join(blocks))


if __name__ == "__main__":
    main()
