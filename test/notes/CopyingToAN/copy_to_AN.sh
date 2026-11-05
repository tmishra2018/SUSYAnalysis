#!/bin/bash
# Copy the AN figures listed in an_figures_manifest.txt (source <TAB> destination relative to the AN directory) to
#   1. the AN area on cmslpc EOS:  /eos/uscms/store/user/tmishra/AN-21-007/<destination>
#   2. the AN area on lxplus:      trmishra@lxplus.cern.ch:/afs/cern.ch/work/t/trmishra/AnalysisNote/AN-21-007/<destination>
# The manifest holds only figures whose source is current (reviewed 2026-10-31); figures that still need a rerun are
# not in it. The previous wildcard-based version (listing.sh) is kept as copy_to_AN.sh.before_manifest_2026-10-31.
#   bash copy_to_AN.sh [--dry-run] [--eos-only] [--lxplus-only] [manifest]

set -euo pipefail

AN_EOS="/eos/uscms/store/user/tmishra/AN-21-007"
DEST_BASE="/afs/cern.ch/work/t/trmishra/AnalysisNote/AN-21-007"
USER="trmishra@lxplus.cern.ch"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${HERE}/an_figures_manifest.txt"
DRY=0; DO_EOS=1; DO_LXPLUS=1
for arg in "$@"; do
    case "$arg" in
        --dry-run)     DRY=1 ;;
        --eos-only)    DO_LXPLUS=0 ;;
        --lxplus-only) DO_EOS=0 ;;
        *)             MANIFEST="$arg" ;;
    esac
done

function green() { echo -e "\033[0;32m$1\033[0m"; }
function yellow() { echo -e "\033[1;33m$1\033[0m"; }
function red() { echo -e "\033[0;31m$1\033[0m"; }

[[ -f "$MANIFEST" ]] || { red "Manifest $MANIFEST not found"; exit 1; }
LIST=$(mktemp)
trap 'rm -f "$LIST"' EXIT

missing=0; n=0
while IFS=$'\t' read -r src dest || [[ -n "${src:-}" ]]; do
    [[ -z "${src:-}" || "$src" == \#* ]] && continue
    if [[ ! -f "$src" ]]; then
        yellow "  Warning: source $src not found, skipping"
        missing=$((missing + 1)); continue
    fi
    echo "$dest" >> "$LIST"; n=$((n + 1))
    if [[ $DO_EOS -eq 1 ]]; then
        if [[ $DRY -eq 1 ]]; then
            echo "  cp $src -> ${AN_EOS}/${dest}"
        elif [[ "$(readlink -f "$src")" != "$(readlink -f "${AN_EOS}/${dest}")" ]]; then
            mkdir -p "$(dirname "${AN_EOS}/${dest}")"
            cp -p "$src" "${AN_EOS}/${dest}"
        fi
    fi
done < "$MANIFEST"
green "${n} files in the manifest, ${missing} missing sources"
[[ $DO_EOS -eq 1 && $DRY -eq 0 ]] && green "✅ Copied to ${AN_EOS}"

if [[ $DO_LXPLUS -eq 1 ]]; then
    green "Syncing the same files to ${USER}:${DEST_BASE}/ ..."
    opts=(-av --files-from="$LIST")
    [[ $DRY -eq 1 ]] && opts+=(--dry-run)
    rsync "${opts[@]}" "${AN_EOS}/" "${USER}:${DEST_BASE}/"
    green "✅ Synced to lxplus"
fi
