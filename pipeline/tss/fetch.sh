#!/usr/bin/env bash
# Download the two FANTOM5 files the TSS-weight pipeline needs (CC BY 4.0, see README.md).
# ~116 MB + ~5 MB. They are NOT committed; re-run this to rebuild from scratch.
set -euo pipefail
DEST="${1:-${TSS_CACHE:-$HOME/.cache/utrrace/fantom5}}"
BASE=https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra
mkdir -p "$DEST"
curl -fL --retry 3 -o "$DEST/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz" \
  "$BASE/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz"
curl -fL --retry 3 -o "$DEST/hg38_fair+new_CAGE_peaks_phase1and2_counts.osc.txt.gz" \
  "$BASE/CAGE_peaks_expression/hg38_fair+new_CAGE_peaks_phase1and2_counts.osc.txt.gz"
echo "downloaded into $DEST"
