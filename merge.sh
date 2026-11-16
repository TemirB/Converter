
#!/usr/bin/env bash
set -euo pipefail

SRC="/mnt/LHEP_2TB/Temir/3gev_FermiMotionFrozen_lso"
OUT="/home/ubuntu/Temir/smash/build/data/Converter/lso-result.mcDst.root"
TMP_BASE="/home/ubuntu/Temir/smash/build/data/Converter/tmp"

BATCH=50
JOBS=4

# Find all McDst files, regardless of case
mapfile -d '' -t files < <(
    find "$SRC" -maxdepth 1 -type f \
        -iname '*mcdst.root' -print0 | sort -z
)

total=${#files[@]}

echo "Found $total ROOT files"

if (( total == 0 )); then
    echo "ERROR: No McDst files found"
    exit 1
fi

if [[ -e "$OUT" ]]; then
    echo "ERROR: Output already exists: $OUT"
    exit 1
fi

# Intermediate files on fast storage
mkdir -p "$TMP_BASE"
TMP=$(mktemp -d "$TMP_BASE/hadd_merge.XXXXXXXX")

echo "Temporary directory: $TMP"

stage=0

while (( ${#files[@]} > 1 )); do
    stage=$((stage + 1))
    echo "===== Stage $stage: ${#files[@]} files ====="

    next_files=()
    batch_id=0

    for ((i=0; i<${#files[@]}; i+=BATCH)); do
        chunk=( "${files[@]:i:BATCH}" )

        part="$TMP/stage_${stage}_$(printf '%05d' "$batch_id").root"

        echo "Merging batch $batch_id (${#chunk[@]} files)"

        if (( ${#chunk[@]} == 1 )); then
            cp -- "${chunk[0]}" "$part"
        else
            hadd -fk -j "$JOBS" "$part" "${chunk[@]}"
        fi

        next_files+=( "$part" )
        batch_id=$((batch_id + 1))
    done

    files=( "${next_files[@]}" )
done

echo "===== Saving final result ====="

mv -- "${files[0]}" "$OUT"

echo "SUCCESS: $OUT"
du -h "$OUT"

# Only clean up after a successful merge
rm -rf -- "$TMP"
