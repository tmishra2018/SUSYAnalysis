#!/usr/bin/env python3
"""Derive proxy/sample composition ratios from SUSYAnalysis ROOT trees.

Computes, for each input file:
  proxyTree/signalTree
  jetTree/signalTree
  fakeLepTree/signalTree
  proxyTree/p_eventcount(bin1)
  jetTree/p_eventcount(bin1)
  fakeLepTree/p_eventcount(bin1)

Supports the per-year mg/eg files and the legacy reference files.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Tuple


MG_YEAR_FILES = {
    "2016preVFP": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016preVFP_doubleFake.root",
    "2016postVFP": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2016postVFP_doubleFake.root",
    "2017": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2017_doubleFake.root",
    "2018": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_mgsignal_MuonEG_2018_doubleFake.root",
}

EG_YEAR_FILES = {
    "2016preVFP": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2016preVFP_doubleFake.root",
    "2016postVFP": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2016postVFP_doubleFake.root",
    "2017": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2017_doubleFake.root",
    "2018": "/eos/uscms/store/user/tmishra/eg_mg_treesData/resTree_egsignal_DoubleEG_2018_doubleFake.root",
}

REF_FILES = {
    "mg_reference": "/uscms_data/d3/mengleis/FullStatusOct/resTree_mgsignal_MuonEG_FullEcal.root",
    "eg_reference": "/uscms_data/d3/mengleis/FullStatusOct/resTree_egsignal_DoubleEG_ReMiniAOD_FullEcal.root",
}


@dataclass
class EntrySummary:
    signal: float
    proxy: float
    jet: float
    fakelep: float
    eventcount_bin1: float


class ReaderError(RuntimeError):
    pass


def _load_backend(name: str):
    """Return (backend_name, module)."""
    if name in ("auto", "uproot"):
        try:
            import uproot  # type: ignore

            return ("uproot", uproot)
        except Exception:
            if name == "uproot":
                raise ReaderError("Requested backend 'uproot' but uproot is not available.")
    if name in ("auto", "pyroot"):
        try:
            import ROOT  # type: ignore

            return ("pyroot", ROOT)
        except Exception:
            if name == "pyroot":
                raise ReaderError("Requested backend 'pyroot' but ROOT is not available.")
    raise ReaderError(
        "No supported ROOT-reading backend found. Install 'uproot' (preferred) "
        "or use a Python with PyROOT available."
    )


def _read_with_uproot(uproot_mod, path: str) -> EntrySummary:
    try:
        with uproot_mod.open(path) as f:
            signal = float(f["signalTree"].num_entries)
            proxy = float(f["proxyTree"].num_entries)
            jet = float(f["jetTree"].num_entries)
            fakelep = float(f["fakeLepTree"].num_entries)
            h = f["p_eventcount"]
            vals, _edges = h.to_numpy()
            eventcount_bin1 = float(vals[0])
            return EntrySummary(signal, proxy, jet, fakelep, eventcount_bin1)
    except Exception as exc:
        raise ReaderError(f"Failed reading '{path}' with uproot: {exc}") from exc


def _read_with_pyroot(root_mod, path: str) -> EntrySummary:
    f = root_mod.TFile.Open(path)
    if not f or f.IsZombie():
        raise ReaderError(f"Failed to open ROOT file '{path}' with PyROOT")
    try:
        def _entries(obj_name: str) -> float:
            obj = f.Get(obj_name)
            if not obj:
                raise ReaderError(f"Object '{obj_name}' missing in '{path}'")
            return float(obj.GetEntries())

        signal = _entries("signalTree")
        proxy = _entries("proxyTree")
        jet = _entries("jetTree")
        fakelep = _entries("fakeLepTree")

        h = f.Get("p_eventcount")
        if not h:
            raise ReaderError(f"Histogram 'p_eventcount' missing in '{path}'")
        eventcount_bin1 = float(h.GetBinContent(1))
        return EntrySummary(signal, proxy, jet, fakelep, eventcount_bin1)
    finally:
        f.Close()


def read_summary(path: str, backend: str, module) -> EntrySummary:
    if backend == "uproot":
        return _read_with_uproot(module, path)
    if backend == "pyroot":
        return _read_with_pyroot(module, path)
    raise ReaderError(f"Unknown backend '{backend}'")


def _safe_ratio(num: float, den: float) -> float:
    if den == 0:
        return math.nan
    return num / den


def _fmt(v: float) -> str:
    if math.isnan(v):
        return "nan"
    return f"{v:.6g}"


def _print_one(label: str, path: str, s: EntrySummary) -> None:
    r_proxy_sig = _safe_ratio(s.proxy, s.signal)
    r_jet_sig = _safe_ratio(s.jet, s.signal)
    r_fake_sig = _safe_ratio(s.fakelep, s.signal)

    r_proxy_evt = _safe_ratio(s.proxy, s.eventcount_bin1)
    r_jet_evt = _safe_ratio(s.jet, s.eventcount_bin1)
    r_fake_evt = _safe_ratio(s.fakelep, s.eventcount_bin1)

    print(f"\n[{label}] {path}")
    print(f"  signalTree entries      : {_fmt(s.signal)}")
    print(f"  proxyTree entries       : {_fmt(s.proxy)}")
    print(f"  jetTree entries         : {_fmt(s.jet)}")
    print(f"  fakeLepTree entries     : {_fmt(s.fakelep)}")
    print(f"  p_eventcount bin1       : {_fmt(s.eventcount_bin1)}")
    print("  -- requested ratios --")
    print(f"  proxyTree/signalTree    : {_fmt(r_proxy_sig)}")
    print(f"  jetTree/signalTree      : {_fmt(r_jet_sig)}")
    print(f"  fakeLepTree/signalTree  : {_fmt(r_fake_sig)}")
    print(f"  proxyTree/eventcount1   : {_fmt(r_proxy_evt)}")
    print(f"  jetTree/eventcount1     : {_fmt(r_jet_evt)}")
    print(f"  fakeLepTree/eventcount1 : {_fmt(r_fake_evt)}")


def _iter_targets(channel: str) -> Iterable[Tuple[str, str]]:
    if channel in ("mg", "both"):
        for year, path in MG_YEAR_FILES.items():
            yield (f"mg_{year}", path)
        yield ("mg_reference", REF_FILES["mg_reference"])
    if channel in ("eg", "both"):
        for year, path in EG_YEAR_FILES.items():
            yield (f"eg_{year}", path)
        yield ("eg_reference", REF_FILES["eg_reference"])


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--channel",
        choices=["mg", "eg", "both"],
        default="both",
        help="Which channel files to process (default: both)",
    )
    parser.add_argument(
        "--backend",
        choices=["auto", "uproot", "pyroot"],
        default="auto",
        help="ROOT reader backend (default: auto)",
    )
    parser.add_argument(
        "--strict-missing",
        action="store_true",
        help="Fail immediately if an input file is missing (default: skip missing with warning).",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    try:
        backend, module = _load_backend(args.backend)
    except ReaderError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    print(f"Using backend: {backend}")

    missing_any = False
    failed_any = False
    for label, path in _iter_targets(args.channel):
        if not os.path.exists(path):
            missing_any = True
            msg = f"WARNING: missing file for {label}: {path}"
            if args.strict_missing:
                print(msg, file=sys.stderr)
                return 3
            print(msg)
            continue

        try:
            summary = read_summary(path, backend, module)
            _print_one(label, path, summary)
        except ReaderError as exc:
            failed_any = True
            print(f"ERROR [{label}]: {exc}", file=sys.stderr)

    if failed_any:
        return 4
    if missing_any:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
