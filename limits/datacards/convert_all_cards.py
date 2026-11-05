#!/usr/bin/env python3
"""
convert_all_cards.py
====================
Batch-convert CMS Combine datacards for all SUSY signal topologies and
Run2 years, applying CMS-recommended systematic nomenclature.

Directory structure
-------------------
  <base_dir>/<TOPOLOGY>/cards/cards_<YEAR>/counting_*.txt   (input, replaced in-place)

Conversion strategy
-------------------
  1. Write converted cards to a temporary folder  cards_<YEAR>_tmp_<pid>/
  2. Rename the original folder to               cards_<YEAR>_backup_<pid>/
  3. Move the temp folder into place as          cards_<YEAR>/
  4. Remove the backup (unless --keep-backup is set)

  If any step fails the original folder is left untouched.

Topologies : T5WG  T6WG  TChiWG
Years       : 2016preVFP  2016postVFP  2017  2018

Modes
-----
  (default)     Convert all cards in-place (no *_new folders)
  --dry-run     Simulate conversion, print what would be done, write nothing
  --check-only  Simulate conversion in memory, report post-rename validation
                status for every card WITHOUT writing any files.
  --keep-backup Keep the cards_<YEAR>_backup_<pid>/ folder after conversion.

Usage
-----
  python3 convert_all_cards.py                          # full run
  python3 convert_all_cards.py --dry-run                # preview
  python3 convert_all_cards.py --check-only             # validate only
  python3 convert_all_cards.py --keep-backup            # keep originals
  python3 convert_all_cards.py --topologies T6WG        # one topology
  python3 convert_all_cards.py --years 2017 2018        # two years
  python3 convert_all_cards.py --verbose                # show per-card details
  python3 convert_all_cards.py \\
      --base-dir /eos/uscms/store/user/tmishra/CombinedLimit \\
      --systematics-dict /path/to/syst_SUS2411.yaml \\
      --rename-dict      /path/to/convert_dict.yaml

Author   : scaffold for SUS2411 SUSY analysis
Reference: https://cms-analysis.docs.cern.ch/guidelines/systematics/
"""

import argparse
import os
import re
import shutil
import sys
import yaml
from pathlib import Path

DEFAULT_BASE_DIR   = "/eos/uscms/store/user/tmishra/CombinedLimit"
DEFAULT_TOPOLOGIES = ["T5WG", "T6WG", "TChiWG"]
DEFAULT_YEARS      = ["2016preVFP", "2016postVFP", "2017", "2018"]

INPUT_FOLDER_PATTERN = "cards_{year}_unblind"

DEFAULT_SYST_YAML   = Path(__file__).parent / "syst_SUS2411.yaml"
DEFAULT_RENAME_YAML = Path(__file__).parent / "convert_dict.yaml"

COMBINE_TYPES = {
    "lnN", "lnU", "gmN", "shape", "rateParam",
    "flatParam", "param", "discrete", "extArg", "autoMCStats",
}
VALID_YEARS = {"2016preVFP", "2016postVFP", "2017", "2018"}

PID = os.getpid()


def load_yaml(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def parse_syst_lines(card_lines):
    result = []
    for i, line in enumerate(card_lines):
        tokens = line.rstrip("\n").split()
        if len(tokens) >= 2 and tokens[1] in COMBINE_TYPES:
            result.append((i, tokens[0], tokens[1], tokens))
    return result


def build_rename_rules(rename_cfg, year):
    """Build flat list of rename rules for a given year."""
    rules = []
    lumi_found = False
    for key, val in rename_cfg.get("rename", {}).items():
        if key == "lumi":
            if str(year) in val:
                rules.append({
                    "old": "lumi", "regex": None, "new": None,
                    "template": None, "lumi": val[str(year)],
                })
                lumi_found = True
        elif key == "isr":
            yr = str(year)
            if yr in val:
                rules.append({
                    "old": "isr", "regex": None,
                    "new": val[yr], "template": None, "lumi": None,
                })
        elif isinstance(val, dict) and "regex" in val:
            rules.append({
                "old": None, "regex": re.compile(val["regex"]),
                "new": None, "template": val["replace"], "lumi": None,
            })
        else:
            rules.append({
                "old": key, "regex": None, "new": val,
                "template": None, "lumi": None,
            })
    return rules, lumi_found


def apply_rename(name, rules):
    """
    Returns (new_name_or_None, was_renamed, lumi_components_or_None).
    """
    for rule in rules:
        if rule["regex"] is not None:
            m = rule["regex"].match(name)
            if m:
                t = rule["template"]
                for idx, grp in enumerate(m.groups(), start=1):
                    t = t.replace(f"{{{idx}}}", grp)
                return t, True, None
        elif rule["lumi"] is not None:
            if rule["old"] == name:
                return None, True, rule["lumi"]
        elif rule["old"] == name:
            return rule["new"], True, None
    return name, False, None


def build_approved_set(syst_cfg):
    """Return (exact_names_set, compiled_regex_list).

    Supports two YAML formats:

    Format A — old list format (syst_SUS2411.yaml before CI migration):
        systematics:
          - name: CMS_scale_j
            old_names: [jes]
          ...
        approved_patterns:
          - regex: "^CMS_SUS2411_VG_stat_bin\\d+$"

    Format B — new flat-dict format required by check_names.py CI:
        CMS_scale_j:
          class: custom
          description: "..."
        CMS_SUS2411_VG_stat_bin[0-9]+:   # key contains [ → treated as regex
          class: custom
          description: "..."
    """
    approved = set()
    regexes  = []

    # ── Format A: list under "systematics" key ────────────────────────
    if "systematics" in syst_cfg:
        for entry in syst_cfg.get("systematics", []):
            approved.add(entry["name"])
            for alias in entry.get("old_names", []):
                approved.add(alias)
        for p in syst_cfg.get("approved_patterns", []):
            regexes.append(re.compile(p["regex"]))
        return approved, regexes

    # ── Format B: flat dict — keys are exact names or regex patterns ──
    # A key is treated as a regex if it contains any of:  *  [  (
    # This matches the logic in check_names.py (find_np).
    for key in syst_cfg:
        if key in ("rename",):          # skip non-NP top-level keys
            continue
        val = syst_cfg[key]
        if not isinstance(val, dict):   # skip non-entry values
            continue
        if any(c in key for c in ('*', '[', '(')):
            # Regex key — anchor it for fullmatch
            pattern = key if key.startswith('^') else '^' + key + '$'
            try:
                regexes.append(re.compile(pattern))
            except re.error:
                pass                    # skip malformed patterns
        else:
            approved.add(key)

    return approved, regexes


def validate(names, approved, regexes):
    ok, unknown = [], []
    for n in names:
        if n in approved or any(rx.match(n) for rx in regexes):
            ok.append(n)
        else:
            unknown.append(n)
    return {"ok": ok, "unknown": unknown}


def expand_lumi_line(original_line, old_name, lumi_components):
    """Replace a single lumi lnN line with one line per component."""
    tokens = original_line.rstrip("\n").split()
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
        new_tok[0] = comp["name"]
        for i in range(2, len(new_tok)):
            if new_tok[i] in orig_values:
                new_tok[i] = comp["value"]
        new_lines.append("  ".join(new_tok) + "\n")
    return new_lines


def process_card(input_path, year, approved, regexes, rename_rules):
    """
    Apply all rename rules to a card in memory.

    Returns dict:
      before_ok   : int   -- recognised names BEFORE renaming
      before_unk  : list  -- unknown names BEFORE renaming
      after_ok    : int   -- recognised names AFTER renaming
      after_unk   : list  -- unknown names AFTER renaming (should be [])
      renamed     : list of (old, new_or_lumi_list)
      new_lines   : list[str]  -- transformed card content
    """
    card_lines = Path(input_path).read_text().splitlines(keepends=True)
    nuisances  = parse_syst_lines(card_lines)
    names_in   = [n for (_, n, _, _) in nuisances]
    val_before = validate(names_in, approved, regexes)

    renamed_log = []
    new_lines   = list(card_lines)
    inserts     = {}   # {line_idx: [extra_line_strings]}

    for (idx, old_name, stype, tokens) in nuisances:
        new_name, was_renamed, lumi_comps = apply_rename(old_name, rename_rules)
        if not was_renamed:
            continue
        renamed_log.append((old_name, lumi_comps or new_name))
        if lumi_comps:
            expanded = expand_lumi_line(card_lines[idx], old_name, lumi_comps)
            new_lines[idx] = expanded[0]
            if len(expanded) > 1:
                inserts[idx] = expanded[1:]
        else:
            new_lines[idx] = new_lines[idx].replace(old_name, new_name, 1)

    # Insert extra lumi lines in reverse order so indices stay valid
    for insert_idx in sorted(inserts.keys(), reverse=True):
        pos = insert_idx + 1
        for extra_line in reversed(inserts[insert_idx]):
            new_lines.insert(pos, extra_line)

    nuisances_after = parse_syst_lines(new_lines)
    names_after     = [n for (_, n, _, _) in nuisances_after]
    val_after       = validate(names_after, approved, regexes)

    return {
        "before_ok":  len(val_before["ok"]),
        "before_unk": val_before["unknown"],
        "after_ok":   len(val_after["ok"]),
        "after_unk":  val_after["unknown"],
        "renamed":    renamed_log,
        "new_lines":  new_lines,
    }


# Status flag characters
FLAG = {"OK": "✓", "WARN": "⚠", "ERROR": "✗"}


def atomic_replace(in_folder, tmp_folder, keep_backup=False):
    """
    Atomically replace in_folder with tmp_folder:
      1. Rename in_folder  ->  in_folder_backup_<pid>
      2. Rename tmp_folder ->  in_folder
      3. Remove backup (unless keep_backup)

    Returns the backup path (str) for logging.
    Raises on any failure; in_folder is left untouched if step 1 fails.
    """
    backup = in_folder.parent / f"{in_folder.name}_backup_{PID}"
    in_folder.rename(backup)
    try:
        tmp_folder.rename(in_folder)
    except Exception:
        # Roll back: restore original
        backup.rename(in_folder)
        raise
    if not keep_backup:
        shutil.rmtree(backup)
        return None
    return str(backup)


def run_batch(base_dir, topologies, years, syst_cfg, rename_cfg,
              dry_run=False, check_only=False, verbose=False,
              keep_backup=False):
    """
    Walk all topology × year combinations, process every card found.

    Conversion flow (write mode only):
      - Write converted cards into  cards_<YEAR>_tmp_<pid>/
      - Atomically swap it into place as  cards_<YEAR>/
      - No *_new folders are ever left on disk.

    Returns (n_ok, n_warn, n_error, n_total).
    """
    approved, regexes = build_approved_set(syst_cfg)

    n_ok = n_warn = n_error = n_total = 0
    summary_rows = []

    write_files = not dry_run and not check_only

    for topology in topologies:
        for year in years:
            in_folder  = (Path(base_dir) / topology / "cards"
                          / INPUT_FOLDER_PATTERN.format(year=year))
            tmp_folder = (Path(base_dir) / topology / "cards"
                          / f"cards_{year}_tmp_{PID}")

            if not in_folder.exists():
                print(f"  [SKIP] Not found: {in_folder}")
                continue

            cards = sorted(in_folder.glob("counting_*.txt"))
            if not cards:
                print(f"  [SKIP] No counting_*.txt in {in_folder}")
                continue

            # Build rename rules once per year
            rename_rules, lumi_found = build_rename_rules(rename_cfg, year)
            lumi_warn = (None if lumi_found else
                         f"[WARN] No lumi entry for year '{year}' in rename dict")

            print(f"\n{'━'*70}")
            print(f"  Topology : {topology}   Year : {year}")
            print(f"  Input    : {in_folder}")
            if write_files:
                print(f"  Temp dir : {tmp_folder}")
                print(f"  Final    : {in_folder}  [replaced in-place]")
            elif dry_run:
                print(f"  Mode     : DRY RUN — no files written")
            else:
                print(f"  Mode     : CHECK ONLY — simulating rename, no files written")
            print(f"  Cards    : {len(cards)}")
            print(f"{'━'*70}")

            if write_files:
                tmp_folder.mkdir(parents=True, exist_ok=True)

            card_results = []   # collect before atomic swap

            for card in cards:
                n_total += 1
                tmp_path = tmp_folder / card.name

                result = process_card(card, year, approved, regexes, rename_rules)

                has_unk  = bool(result["after_unk"])
                has_warn = bool(lumi_warn)

                if has_unk:
                    status = "ERROR"; n_error += 1
                elif has_warn:
                    status = "WARN";  n_warn  += 1
                else:
                    status = "OK";    n_ok    += 1

                flag_char = FLAG[status]

                print(f"  {flag_char} {card.name:50s}  "
                      f"ok={result['after_ok']:3d}  "
                      f"unk={len(result['after_unk']):3d}")

                if result["after_unk"]:
                    for u in result["after_unk"]:
                        print(f"      [?] {u}")

                if verbose and result["renamed"]:
                    for old, new in result["renamed"]:
                        tag = (" + ".join(c["name"] for c in new)
                               if isinstance(new, list) else new)
                        print(f"      {old:40s} -> {tag}")

                if lumi_warn:
                    print(f"    {lumi_warn}")

                summary_rows.append((topology, year, card.name, status,
                                     result["after_ok"],
                                     len(result["after_unk"])))

                # Stage converted card into temp folder
                if write_files:
                    tmp_path.write_text("".join(result["new_lines"]))
                    card_results.append(card.name)

            # ── Atomic swap: tmp -> in_folder ─────────────────────────────
            if write_files:
                if n_error > 0:
                    # Don't swap if there are errors; clean up temp folder
                    shutil.rmtree(tmp_folder)
                    print(f"\n  ✗ Skipping in-place swap for {in_folder.name} "
                          f"due to {n_error} error(s). Temp folder removed.")
                else:
                    backup_path = atomic_replace(in_folder, tmp_folder,
                                                 keep_backup=keep_backup)
                    if keep_backup and backup_path:
                        print(f"\n  ✓ Converted {len(card_results)} card(s) in-place.")
                        print(f"    Original backed up to: {backup_path}")
                    else:
                        print(f"\n  ✓ Converted {len(card_results)} card(s) in-place.")
                        print(f"    Original folder replaced: {in_folder}")

    # ── Final summary ─────────────────────────────────────────────────────
    print(f"\n{'═'*70}")
    print("  FINAL SUMMARY")
    print(f"{'═'*70}")

    if len(summary_rows) <= 100:
        col_w = max((len(r[2]) for r in summary_rows), default=20)
        col_w = max(col_w, 30)
        print(f"  {'Topo':<8} {'Year':<14} {'Card':<{col_w}}  St   OK  UNK")
        print(f"  {'-'*7} {'-'*13} {'-'*col_w}  --  ---  ---")
        for (topo, yr, cname, st, ok, unk) in summary_rows:
            print(f"  {FLAG[st]} {topo:<7} {yr:<14} {cname:<{col_w}}  "
                  f"{ok:3d}  {unk:3d}")
        print()

    print(f"  Total cards  : {n_total}")
    print(f"  ✓ OK         : {n_ok}")
    print(f"  ⚠ Warnings   : {n_warn}")
    print(f"  ✗ Errors     : {n_error}")

    if check_only:
        print("\n  Mode: CHECK ONLY — no files were written.")
        print("  All counts above reflect the state AFTER renaming (simulated).")
        if n_error == 0:
            print("  → All cards are ready for conversion. Run without --check-only to proceed.")
        else:
            print("  → Fix the unknown systematics above, then re-run --check-only.")
    elif dry_run:
        print("\n  [DRY RUN] No files were written.")
        if n_error == 0:
            print("  → All cards look good. Remove --dry-run to write output files.")
    else:
        if n_error == 0:
            print(f"\n  → All cards replaced in-place under cards_<YEAR>/.")
            print(f"     No *_new or *_tmp folders remain on disk.")
        else:
            print(f"\n  → {n_error} card(s) have residual unknown systematics.")
            print(f"     Check the [?] lines above and update convert_dict.yaml.")
            print(f"     Affected folders were NOT modified.")

    print(f"{'═'*70}")

    return n_ok, n_warn, n_error, n_total


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Batch-convert CMS Combine datacards for all topologies and years "
            "to CMS-recommended systematic nomenclature. "
            "Converted cards replace the originals in-place; no *_new folders."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Modes:
  (default)     Convert all cards and replace cards_<YEAR>/ in-place
  --dry-run     Simulate conversion, print results, write nothing
  --check-only  Simulate conversion, validate post-rename names, write nothing
                (Use this FIRST to confirm all systematics are handled)
  --keep-backup Retain cards_<YEAR>_backup_<pid>/ after successful conversion

Examples:
  python3 convert_all_cards.py --check-only
  python3 convert_all_cards.py --dry-run
  python3 convert_all_cards.py
  python3 convert_all_cards.py --keep-backup
  python3 convert_all_cards.py --topologies T6WG --years 2017 2018
  python3 convert_all_cards.py --verbose
  python3 convert_all_cards.py \\
      --base-dir /eos/uscms/store/user/tmishra/CombinedLimit \\
      --systematics-dict /path/to/syst_SUS2411.yaml \\
      --rename-dict      /path/to/convert_dict.yaml
""")
    p.add_argument("--base-dir", "-b", default=DEFAULT_BASE_DIR,
                   help=f"Base directory (default: {DEFAULT_BASE_DIR})")
    p.add_argument("--topologies", "-t", nargs="+",
                   default=DEFAULT_TOPOLOGIES,
                   choices=["T5WG", "T6WG", "TChiWG"],
                   help="Signal topologies to process (default: all three)")
    p.add_argument("--years", "-y", nargs="+",
                   default=DEFAULT_YEARS,
                   choices=list(VALID_YEARS),
                   help="Years to process (default: all four)")
    p.add_argument("--systematics-dict", "-s",
                   default=str(DEFAULT_SYST_YAML),
                   help="Path to syst_SUS2411.yaml")
    p.add_argument("--rename-dict", "-r",
                   default=str(DEFAULT_RENAME_YAML),
                   help="Path to convert_dict.yaml")
    p.add_argument("--dry-run", action="store_true",
                   help="Simulate but do not write output files")
    p.add_argument("--check-only", action="store_true",
                   help="Validate post-rename names only; no files written")
    p.add_argument("--keep-backup", action="store_true",
                   help="Keep cards_<YEAR>_backup_<pid>/ after conversion")
    p.add_argument("--verbose", "-v", action="store_true",
                   help="Print per-card rename details")
    return p.parse_args()


def main():
    args = parse_args()

    for label, path in [("--systematics-dict", args.systematics_dict),
                        ("--rename-dict",       args.rename_dict)]:
        if not Path(path).exists():
            print(f"ERROR: {label} not found: {path}", file=sys.stderr)
            sys.exit(2)

    syst_cfg   = load_yaml(args.systematics_dict)
    rename_cfg = load_yaml(args.rename_dict)

    mode = ("CHECK ONLY" if args.check_only
            else ("DRY RUN" if args.dry_run else "CONVERT (in-place)"))
    print(f"Base dir          : {args.base_dir}")
    print(f"Topologies        : {args.topologies}")
    print(f"Years             : {args.years}")
    print(f"Systematics dict  : {args.systematics_dict}")
    print(f"Rename dict       : {args.rename_dict}")
    print(f"Mode              : {mode}")
    if args.keep_backup:
        print(f"Keep backup       : yes  (cards_<YEAR>_backup_{PID}/)")

    n_ok, n_warn, n_error, n_total = run_batch(
        base_dir    = args.base_dir,
        topologies  = args.topologies,
        years       = args.years,
        syst_cfg    = syst_cfg,
        rename_cfg  = rename_cfg,
        dry_run     = args.dry_run,
        check_only  = args.check_only,
        verbose     = args.verbose,
        keep_backup = args.keep_backup,
    )

    sys.exit(0 if n_error == 0 else 1)


if __name__ == "__main__":
    main()
