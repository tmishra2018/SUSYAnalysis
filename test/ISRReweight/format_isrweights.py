#!/usr/bin/env python3
"""Format channel-specific ISRweights_*.txt files into C++ if/else blocks.

Expected input file format:
  # normalization <value>
  # ptLow ptHigh weight
  <ptLow> <ptHigh> <weight>
  ...
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional, Sequence, Tuple


Bin = Tuple[float, float, float]
EraMapping = Tuple[int, Optional[int], Path]
ChannelMapping = Tuple[int, str, Sequence[EraMapping]]


def parse_isr_file(path: Path) -> Tuple[float, List[Bin]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")

    normalization = None
    bins: List[Bin] = []

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
    if not bins:
        raise ValueError(f"No ISR bins found in {path}")

    return normalization, bins


def fmt_num(x: float) -> str:
    return f"{x:.6g}"


def format_year_block(
    run_year: int,
    pre_vfp: int | None,
    normalization: float,
    bins: List[Bin],
    is_first: bool,
    indent: str = "\t",
) -> str:
    if run_year == 2016:
        head = f"{indent}{'if' if is_first else 'else if'}(RunYear==2016 && preVFP=={pre_vfp}){{"
    else:
        head = f"{indent}{'if' if is_first else 'else if'}(RunYear=={run_year}){{"

    lines = [head]
    for i, (low, high, weight) in enumerate(bins):
        if i == 0:
            lines.append(f"{indent}\tif(ISRJetPt < {fmt_num(high)})reweightF = {fmt_num(weight)};")
        elif i < len(bins) - 1:
            lines.append(
                f"{indent}\telse if(ISRJetPt >= {fmt_num(low)} && ISRJetPt < {fmt_num(high)})reweightF = {fmt_num(weight)};"
            )
        else:
            lines.append(f"{indent}\telse if(ISRJetPt >= {fmt_num(low)})reweightF = {fmt_num(weight)};")

    lines.append(f"{indent}\tNormalization = {fmt_num(normalization)};")
    lines.append(f"{indent}}}")
    return "\n".join(lines)


def format_channel_block(channel_type: int, channel_name: str, mappings: Sequence[EraMapping], is_first: bool) -> str:
    lines = [f"{'if' if is_first else 'else if'}(channelType == {channel_type}){{"]

    for i, (year, pre_vfp, path) in enumerate(mappings):
        normalization, bins = parse_isr_file(path)
        lines.append(
            format_year_block(
                year,
                pre_vfp,
                normalization,
                bins,
                is_first=(i == 0),
            )
        )

    lines.append("}")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Format channel-specific ISRweights text files into C++ if/else blocks"
    )
    parser.add_argument(
        "--dir",
        default="/eos/uscms/store/user/tmishra/ISRweighting/",
        help="Directory containing ISRweights_*.txt files",
    )
    args = parser.parse_args()
    base = Path(args.dir)
    channel_mappings: List[ChannelMapping] = [
        (
            1,
            "egamma",
            [
                (2016, 1, base / "ISRweights_2016preVFP_egamma.txt"),
                (2016, 0, base / "ISRweights_2016postVFP_egamma.txt"),
                (2017, None, base / "ISRweights_2017_egamma.txt"),
                (2018, None, base / "ISRweights_2018_egamma.txt"),
            ],
        ),
        (
            2,
            "mgamma",
            [
                (2016, 1, base / "ISRweights_2016preVFP_mgamma.txt"),
                (2016, 0, base / "ISRweights_2016postVFP_mgamma.txt"),
                (2017, None, base / "ISRweights_2017_mgamma.txt"),
                (2018, None, base / "ISRweights_2018_mgamma.txt"),
            ],
        ),
    ]

    blocks = [
        format_channel_block(channel_type, channel_name, mappings, is_first=(i == 0))
        for i, (channel_type, channel_name, mappings) in enumerate(channel_mappings)
    ]
    print("\n\n".join(blocks))


if __name__ == "__main__":
    main()
