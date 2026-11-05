#!/usr/bin/env python3
"""
check_names.py  --  Validate and rename systematics in CMS Combine datacards
                   to follow the CMS recommended nomenclature.

Usage
-----
# Single year card:
python3 systematics/check_names.py \\
    --input   counting_t6wg_1500_1350_2018.txt \\
    --year    2018 \\
    --systematics-dict  syst_SUSxxxx.yaml \\
    --rename-dict       convert_dict.yaml \\
    --output  datacard_Signal_2018_new.txt

# Loop over all years and mass points:
for year in 2016preVFP 2016postVFP 2017 2018; do
  for card in counting_t6wg_*_${year}.txt; do
    base=$(basename $card .txt)
    python3 systematics/check_names.py \\
        --input  ${card} \\
        --year   ${year} \\
        --systematics-dict syst_SUSxxxx.yaml \\
        --rename-dict      convert_dict.yaml \\
        --output ${base}_new.txt
  done
done

Notes on lumi treatment
-----------------------
The single "lumi" line in the original card is expanded into 2-3 lines
depending on the year, each with the correct correlated/uncorrelated value:

  2016preVFP / 2016postVFP:
      CMS_lumi_13TeV_correlated   lnN  1.006  ...
      CMS_lumi_13TeV_2016         lnN  1.010  ...

  2017:
      CMS_lumi_13TeV_correlated   lnN  1.009  ...
      CMS_lumi_13TeV_1718         lnN  1.006  ...
      CMS_lumi_13TeV_2017         lnN  1.020  ...

  2018:
      CMS_lumi_13TeV_correlated   lnN  1.020  ...
      CMS_lumi_13TeV_1718         lnN  1.002  ...
      CMS_lumi_13TeV_2018         lnN  1.015  ...

When combineCards.py merges per-year cards, nuisances with the same
name are automatically 100% correlated. Different names = uncorrelated.

Author   : auto-generated scaffold for SUSxxxx SUSY analysis
Reference: https://cms-analysis.docs.cern.ch/guidelines/systematics/
"""

import argparse
import re
import sys
import yaml
from pathlib import Path


VALID_YEARS = {"2016preVFP", "2016postVFP", "2017", "2018"}

# Combine nuisance types that mark a line as a systematic
COMBINE_TYPES = {
    "lnN", "lnU", "gmN", "shape", "rateParam",
    "flatParam", "param", "discrete", "extArg", "autoMCStats",
}


def load_yaml(path: str) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def infer_year(input_path: str) -> str | None:
    """Try to guess the era from the filename if --year is not supplied."""
    name = Path(input_path).stem
    for era in ("2016preVFP", "2016postVFP", "2017", "2018"):
        if era in name:
            return era
    return None


def parse_syst_lines(card_lines: list[str]) -> list[tuple[int, str, str, list[str]]]:
    """
    Return (line_index, syst_name, syst_type, tokens) for every nuisance line.
    """
    result = []
    for i, line in enumerate(card_lines):
        tokens = line.rstrip("\n").split()
        if len(tokens) >= 2 and tokens[1] in COMBINE_TYPES:
            result.append((i, tokens[0], tokens[1], tokens))
    return result


def build_rename_rules(rename_cfg: dict, year: str) -> list[dict]:
    """
    Convert the rename section of convert_dict.yaml into a flat list of rules.

    Each rule is a dict with keys:
      old      : exact old name to match (str) or None for regex rules
      regex    : compiled regex or None
      new      : new name string, or None for regex/lumi rules
      template : replacement template string (for regex rules)
      lumi     : list of {name, value} dicts (for the lumi expansion)
    """
    rules = []
    for key, val in rename_cfg.get("rename", {}).items():

        if key == "lumi":
            if str(year) not in val:
                print(f"  [WARN] No lumi entry for year '{year}' in convert_dict.yaml")
                continue
            rules.append({
                "old":      "lumi",
                "regex":    None,
                "new":      None,
                "template": None,
                "lumi":     val[str(year)],   # list of {name, value}
            })

        elif isinstance(val, dict) and "regex" in val:
            rules.append({
                "old":      None,
                "regex":    re.compile(val["regex"]),
                "new":      None,
                "template": val["replace"],
                "lumi":     None,
            })

        else:
            rules.append({
                "old":      key,
                "regex":    None,
                "new":      val,
                "template": None,
                "lumi":     None,
            })

    return rules


def apply_rename(name: str, rules: list[dict]) -> tuple:
    """
    Apply the first matching rule to `name`.

    Returns one of:
      (name,           False, None)        -- no match, unchanged
      (new_name_str,   True,  None)        -- simple rename
      (new_name_str,   True,  None)        -- regex rename
      (None,           True,  lumi_list)   -- lumi expansion
    """
    for rule in rules:
        # regex rule
        if rule["regex"] is not None:
            m = rule["regex"].match(name)
            if m:
                new_name = rule["template"]
                for idx, grp in enumerate(m.groups(), start=1):
                    new_name = new_name.replace(f"{{{idx}}}", grp)
                return new_name, True, None

        # lumi expansion
        elif rule["lumi"] is not None:
            if rule["old"] == name:
                return None, True, rule["lumi"]

        # simple rename
        elif rule["old"] == name:
            return rule["new"], True, None

    return name, False, None


def build_approved_set(syst_cfg: dict) -> tuple[set, list]:
    """Return (approved_exact_names, approved_compiled_regexes)."""
    approved = set()
    for entry in syst_cfg.get("systematics", []):
        approved.add(entry["name"])
        for alias in entry.get("old_names", []):
            approved.add(alias)
    regexes = [
        re.compile(p["regex"])
        for p in syst_cfg.get("approved_patterns", [])
    ]
    return approved, regexes


def validate(names: list[str], approved: set, regexes: list) -> dict:
    ok, unknown = [], []
    for n in names:
        if n in approved or any(rx.match(n) for rx in regexes):
            ok.append(n)
        else:
            unknown.append(n)
    return {"ok": ok, "unknown": unknown}


def expand_lumi_line(original_line: str, old_name: str,
                     lumi_components: list[dict]) -> list[str]:
    """
    Given the original "lumi lnN 1.026 - 1.026 ..." line,
    produce one replacement line per lumi component, each with
    the correct lnN value substituted everywhere the original value appeared.

    Strategy:
      - Find the original lnN value(s) in the line (the non-"-" numeric tokens
        after the type token).
      - Replace each such token with the component's canonical value.
      - Replace the name token with the component's name.
    """
    tokens = original_line.rstrip("\n").split()
    # tokens[0] = old_name, tokens[1] = "lnN", tokens[2:] = values

    # Collect the set of numeric values present (to replace them)
    orig_values = set()
    for tok in tokens[2:]:
        if tok != "-":
            try:
                float(tok)
                orig_values.add(tok)
            except ValueError:
                pass

    new_lines = []
    for comp in lumi_components:
        new_tok = list(tokens)
        # Replace the name
        new_tok[0] = comp["name"]
        # Replace every numeric value with this component's value
        for i in range(2, len(new_tok)):
            if new_tok[i] in orig_values:
                new_tok[i] = comp["value"]
        new_line = "  ".join(new_tok) + "\n"  # preserve rough spacing
        new_lines.append(new_line)

    return new_lines


def convert_card(input_path: str,
                 syst_cfg: dict,
                 rename_cfg: dict,
                 year: str,
                 output_path: str,
                 dry_run: bool = False) -> int:
    """
    Returns exit code: 0 = all OK after renaming, 1 = unknowns remain.
    """
    card_lines = Path(input_path).read_text().splitlines(keepends=True)
    rules      = build_rename_rules(rename_cfg, year)
    approved, regexes = build_approved_set(syst_cfg)

    nuisances      = parse_syst_lines(card_lines)
    names_in_card  = [n for (_, n, _, _) in nuisances]

    val_before = validate(names_in_card, approved, regexes)
    print(f"\n{'='*60}")
    print(f"  Year   : {year}")
    print(f"  Input  : {input_path}")
    print(f"{'='*60}")
    print(f"\n[BEFORE] Recognised: {len(val_before['ok'])}  |  "
          f"Unknown: {len(val_before['unknown'])}")
    for u in val_before["unknown"]:
        print(f"    [UNKNOWN] {u}")

    renamed_log = []
    # Work on a mutable list; build a dict of lines to insert after each idx
    new_lines   = list(card_lines)
    inserts     = {}   # {original_line_idx: [list_of_new_line_strings]}

    for (idx, old_name, stype, tokens) in nuisances:
        new_name, was_renamed, lumi_comps = apply_rename(old_name, rules)

        if not was_renamed:
            continue

        renamed_log.append((old_name, lumi_comps or new_name))

        if lumi_comps:
            # Expand lumi -> multiple lines, each with correct value
            expanded = expand_lumi_line(card_lines[idx], old_name, lumi_comps)
            # Replace the original line with the first component
            new_lines[idx] = expanded[0]
            # Queue remaining components to insert after
            if len(expanded) > 1:
                inserts[idx] = expanded[1:]
        else:
            # Simple rename: replace first occurrence of the name token
            new_lines[idx] = new_lines[idx].replace(old_name, new_name, 1)

    # Insert extra lines in reverse order so indices stay valid
    for insert_idx in sorted(inserts.keys(), reverse=True):
        pos = insert_idx + 1
        for extra_line in reversed(inserts[insert_idx]):
            new_lines.insert(pos, extra_line)

    nuisances_after = parse_syst_lines(new_lines)
    names_after     = [n for (_, n, _, _) in nuisances_after]
    val_after       = validate(names_after, approved, regexes)

    print(f"\n[AFTER]  Recognised: {len(val_after['ok'])}  |  "
          f"Unknown: {len(val_after['unknown'])}")
    for u in val_after["unknown"]:
        print(f"    [STILL UNKNOWN] {u}  <- needs manual review")

    print("  Rename summary:")
    if renamed_log:
        for old, new in renamed_log:
            if isinstance(new, list):
                names = " + ".join(c["name"] for c in new)
                print(f"    {old:38s} ->  {names}")
            else:
                print(f"    {old:38s} ->  {new}")
    else:
        print("    (no renames applied)")

    if dry_run:
        print(f"\n  [DRY RUN] Would write: {output_path}")
    else:
        Path(output_path).write_text("".join(new_lines))
        print(f"\n  Output written: {output_path}")

    return 0 if not val_after["unknown"] else 1


def parse_args():
    p = argparse.ArgumentParser(
        description="Check and rename CMS datacard systematics to CMS-recommended "
                    "nomenclature, with year-aware Run2 lumi decomposition.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single card:
  python3 check_names.py -i counting_t6wg_1500_1350_2018.txt \\
      -y 2018 -s syst_SUSxxxx.yaml -r convert_dict.yaml \\
      -o datacard_Signal_2018_new.txt

  # Year not in filename -- supply explicitly:
  python3 check_names.py -i card.txt -y 2016preVFP \\
      -s syst_SUSxxxx.yaml -r convert_dict.yaml -o card_new.txt

  # Validate only (no output written):
  python3 check_names.py -i card.txt -y 2017 \\
      -s syst_SUSxxxx.yaml -r convert_dict.yaml -o /dev/null --check-only
""")
    p.add_argument("--input",  "-i", required=True,
                   help="Input datacard (.txt)")
    p.add_argument("--year",   "-y",
                   choices=list(VALID_YEARS),
                   help="Era: 2016preVFP | 2016postVFP | 2017 | 2018 "
                        "(auto-detected from filename if omitted)")
    p.add_argument("--systematics-dict", "-s", required=True,
                   help="syst_SUSxxxx.yaml -- approved systematics definitions")
    p.add_argument("--rename-dict", "-r", required=True,
                   help="convert_dict.yaml -- old->new name mapping")
    p.add_argument("--output", "-o", required=True,
                   help="Output datacard path")
    p.add_argument("--dry-run", action="store_true",
                   help="Print report but do not write output file")
    p.add_argument("--check-only", action="store_true",
                   help="Validate current names only; no renaming, no output")
    return p.parse_args()


def main():
    args = parse_args()

    # Resolve year
    year = args.year or infer_year(args.input)
    if year is None:
        print("ERROR: Cannot determine year from filename. "
              "Please supply --year explicitly.", file=sys.stderr)
        sys.exit(2)
    if year not in VALID_YEARS:
        print(f"ERROR: Unrecognised year '{year}'. "
              f"Choose from {sorted(VALID_YEARS)}.", file=sys.stderr)
        sys.exit(2)

    syst_cfg   = load_yaml(args.systematics_dict)
    rename_cfg = load_yaml(args.rename_dict)

    if args.check_only:
        card_lines = Path(args.input).read_text().splitlines(keepends=True)
        nuisances  = parse_syst_lines(card_lines)
        names      = [n for (_, n, _, _) in nuisances]
        approved, regexes = build_approved_set(syst_cfg)
        val = validate(names, approved, regexes)
        print(f"\nCheck-only: {args.input}  (year={year})")
        print(f"  Recognised : {len(val['ok'])}")
        print(f"  Unknown    : {len(val['unknown'])}")
        for u in val["unknown"]:
            print(f"    [UNKNOWN] {u}")
        sys.exit(0 if not val["unknown"] else 1)

    rc = convert_card(
        input_path  = args.input,
        syst_cfg    = syst_cfg,
        rename_cfg  = rename_cfg,
        year        = year,
        output_path = args.output,
        dry_run     = args.dry_run,
    )
    sys.exit(rc)


if __name__ == "__main__":
    main()
