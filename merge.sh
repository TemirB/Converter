#!/usr/bin/env bash
set -euo pipefail

# Preserve the original defaults; optional arguments make the script usable elsewhere.
SRC="${1:-/mnt/LHEP_2TB/Temir/3gev_FermiMotionFrozen_lso}"
OUT="${2:-/home/ubuntu/Temir/smash/build/data/Converter/lso-result.mcDst.root}"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

if (( $# > 2 )); then
    echo "Usage: $0 [INPUT_DIRECTORY [OUTPUT.root]]" >&2
    exit 1
fi
if [[ ! -x "$SCRIPT_DIR/mcdst_merge" ]]; then
    echo "Build the merger first: make -j4 mcdst_merge" >&2
    exit 1
fi

exec "$SCRIPT_DIR/mcdst_merge" "$OUT" "$SRC"
