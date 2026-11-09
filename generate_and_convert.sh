#!/usr/bin/env bash
# Example: ./generate_and_convert.sh 16 125 job
set -u

if [[ $# -lt 2 || $# -gt 3 || ! $1 =~ ^[1-9][0-9]*$ || ! $2 =~ ^[1-9][0-9]*$ ]]; then
    echo "Usage: $0 <tasks> <runs> [job_name]" >&2
    exit 1
fi
N_TASKS=$1
N_RUNS=$2
JOB_NAME=${3:-job}
if [[ ! $JOB_NAME =~ ^[A-Za-z0-9_.-]+$ ]]; then
    echo "Job name may contain only letters, numbers, underscores, dots and dashes." >&2
    exit 1
fi

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
CONVERTER_DIR=${CONVERTER_DIR:-$SCRIPT_DIR}
SMASH_BUILD_DIR=${SMASH_BUILD_DIR:-/home/tahea/smash/build}
TMP_DIR=${TMP_DIR:-${CONVERTER_DIR}/tmp}
OSCAR_DST_DIR=${OSCAR_DST_DIR:-${CONVERTER_DIR}/oscarDst}
LOG_DIR=${LOG_DIR:-${CONVERTER_DIR}/data/generator/logs}
export ROOTSYS=${ROOTSYS:-/home/tahea/root_install}
if [[ -f $ROOTSYS/bin/thisroot.sh ]]; then
    source "$ROOTSYS/bin/thisroot.sh"
fi
if [[ ! -x $SMASH_BUILD_DIR/smash || ! -x $CONVERTER_DIR/converter ]]; then
    echo "Set SMASH_BUILD_DIR and CONVERTER_DIR to directories with executable smash/converter." >&2
    exit 1
fi
mkdir -p "$TMP_DIR" "$OSCAR_DST_DIR" "$LOG_DIR" || exit 1
SMASH_LOG=${LOG_DIR}/${JOB_NAME}_smash.log
CONVERT_LOG=${LOG_DIR}/${JOB_NAME}_convert.log

run_task() {
    local task=$1 run=$2
    local output_dir=${TMP_DIR}/${JOB_NAME}_${task}_${run}
    local root_file=${OSCAR_DST_DIR}/${JOB_NAME}_${task}_${run}.McDst.root
    if [[ -e $output_dir ]]; then
        echo "Preserving existing directory: $output_dir" >> "$SMASH_LOG"
        return 1
    fi
    mkdir -p "$output_dir" || return 1
    if ! (cd "$SMASH_BUILD_DIR" && ./smash -o "$output_dir" >> "$SMASH_LOG" 2>&1); then
        echo "Generation failed: task=$task run=$run; preserved $output_dir" >> "$SMASH_LOG"
        return 1
    fi
    if ! "$CONVERTER_DIR/converter" "$output_dir/full_event_history.oscar" "$root_file" >> "$CONVERT_LOG" 2>&1; then
        echo "Conversion failed: task=$task run=$run; preserved $output_dir" >> "$CONVERT_LOG"
        return 1
    fi
    # This directory was created by this task; delete only after successful conversion.
    rm -rf -- "$output_dir"
}

for ((run = 1; run <= N_RUNS; ++run)); do
    echo "Starting $JOB_NAME run $run/$N_RUNS"
    pids=()
    for ((task = 1; task <= N_TASKS; ++task)); do
        run_task "$task" "$run" &
        pids+=("$!")
    done
    failed=0
    for pid in "${pids[@]}"; do
        if ! wait "$pid"; then
            failed=1
        fi
    done
    if ((failed)); then
        echo "Run failed; original OSCAR files for failed tasks were preserved. See $LOG_DIR." >&2
        exit 1
    fi
done
