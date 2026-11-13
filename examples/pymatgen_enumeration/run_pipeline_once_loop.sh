#!/bin/bash
# automag pipeline launcher - one-shot loop (manual-restart simulator).
#
# The pipeline is advanced one pass at a time with `run_automag --once`: every
# pass starts whatever can be started, records the jobs in
# `automag_state.sqlite` and exits.  This is exactly what you would do by hand,
#
#     ./run_pipeline_once.sh          # one pass
#     ... wait for the calculations, run it again ...
#
# only automated: the script repeats the pass, sleeps between passes and stops
# when the pipeline reports that every step is complete.  Because all progress
# lives in the state file, the loop is resumable and several invocations (or a
# cron job) never do the same work twice.
#
# Usage (from this directory):
#
#     export AUTOMAG_HOME=/path/to/automag      # source checkout
#     # or: export AUTOMAG_EXE=/shared/bin/run_automag   # compiled binary
#     ./run_pipeline_once_loop.sh
#
# Optional environment variables:
#
#     AUTOMAG_ONCE_INTERVAL     seconds between passes (default 60)
#     AUTOMAG_ONCE_MAX_ROUNDS   stop after this many passes, 0 = no limit
#                               (default 0; useful for cron, where one pass per
#                               invocation is enough: `*/20 * * * * ./run_pipeline_once.sh`)
#     AUTOMAG_HOME / AUTOMAG_EXE / AUTOMAG_CONDA_ENV / AUTOMAG_CONDA_SH
#         as in run_pipeline.sh
#
# Extra arguments are forwarded to every pass, e.g.
#
#     ./run_pipeline_once_loop.sh --dry-run
#     ./run_pipeline_once_loop.sh --stages coll --steps submit
#
# For the "run everything and wait here" behaviour use run_pipeline.sh (or
# run_pipeline_slurm.sh to put the supervisor into a SLURM job).
set -euo pipefail

cd "$(dirname "$0")"

# --- Python environment (only needed for a source checkout) ----------------
if [ -z "${AUTOMAG_EXE:-}" ]; then
    AUTOMAG_CONDA_ENV="${AUTOMAG_CONDA_ENV-automag}"
    if [ -n "${AUTOMAG_CONDA_ENV}" ]; then
        # conda's own activation hooks are not `set -u` safe (deactivate.d
        # scripts reference unset variables), so relax nounset while activating.
        set +u
        if [ -n "${AUTOMAG_CONDA_SH:-}" ] && [ -f "${AUTOMAG_CONDA_SH}" ]; then
            # shellcheck disable=SC1090
            source "${AUTOMAG_CONDA_SH}"
        elif command -v conda >/dev/null 2>&1; then
            eval "$(conda shell.bash hook)"
        elif [ -f /path/to/miniconda3/etc/profile.d/conda.sh ]; then
            # EDIT THIS PATH if conda is not on PATH.
            source /path/to/miniconda3/etc/profile.d/conda.sh
        fi
        conda activate "${AUTOMAG_CONDA_ENV}" || true
        set -u
    fi
fi

# --- locate the pipeline: binary, checkout, PATH, or ../../ ---------------
COMMAND=()
if [ -n "${AUTOMAG_EXE:-}" ]; then
    if [ ! -x "${AUTOMAG_EXE}" ]; then
        echo "ERROR: '${AUTOMAG_EXE}' is not an executable file." >&2
        exit 1
    fi
    COMMAND=("${AUTOMAG_EXE}")
elif command -v run_automag >/dev/null 2>&1; then
    COMMAND=("$(command -v run_automag)")
elif [ -n "${AUTOMAG_HOME:-}" ]; then
    COMMAND=(python "${AUTOMAG_HOME}/run_automag.py")
elif [ -f ../../run_automag.py ]; then
    COMMAND=(python "$(cd ../.. && pwd)/run_automag.py")
else
    echo "ERROR: set AUTOMAG_HOME (source checkout) or AUTOMAG_EXE (compiled binary)." >&2
    exit 1
fi

# --- the controlling file: input_automag.py wins over .yaml, which wins .json
CONFIG=""
for candidate in input_automag.py input_automag.yaml input_automag.yml input_automag.json; do
    if [ -f "$candidate" ]; then
        CONFIG="$candidate"
        break
    fi
done
if [ -z "$CONFIG" ]; then
    echo "ERROR: no input_automag.py/.yaml/.json in $PWD" >&2
    exit 1
fi

INTERVAL="${AUTOMAG_ONCE_INTERVAL:-60}"
MAX_ROUNDS="${AUTOMAG_ONCE_MAX_ROUNDS:-0}"

echo "automag:  ${COMMAND[*]}"
echo "config:   ${PWD}/${CONFIG}"
echo "mode:     one pass every ${INTERVAL}s (max rounds: ${MAX_ROUNDS:-unlimited})"

round=0
while true; do
    round=$((round + 1))
    echo
    echo "=== pass ${round}: $(date -Is) ==="
    set +e
    output="$("${COMMAND[@]}" "$CONFIG" --once "$@" 2>&1)"
    status=$?
    set -e
    printf '%s\n' "$output"
    if [ "$status" -ne 0 ]; then
        echo "ERROR: automag exited with ${status}; stopping the loop." >&2
        exit "$status"
    fi
    if printf '%s' "$output" | grep -q "every step is complete"; then
        echo "Pipeline finished after ${round} pass(es): $(date -Is)"
        break
    fi
    if [ "${MAX_ROUNDS}" -gt 0 ] && [ "$round" -ge "${MAX_ROUNDS}" ]; then
        echo "Stopped after ${round} pass(es) (AUTOMAG_ONCE_MAX_ROUNDS)."
        echo "Run the same command again to continue."
        break
    fi
    sleep "${INTERVAL}"
done