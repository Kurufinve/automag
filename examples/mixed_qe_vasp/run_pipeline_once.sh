#!/bin/bash
# automag pipeline launcher - one pass (manual restart).
#
# Advances the pipeline by exactly one pass and exits: whatever can be started
# is started and recorded in `automag_state.sqlite`, everything that still has
# to wait for external programs is left for the next invocation.  Nothing is
# recomputed, so running this command again (after the calculations finished,
# or tomorrow) simply continues where the previous pass stopped:
#
#     ./run_pipeline_once.sh          # pass 1: submits the first batch
#     ... the cluster runs the jobs ...
#     ./run_pipeline_once.sh          # pass 2: analysis + the next batch
#     ./run_pipeline_once.sh          # ... until it reports "every step is complete"
#
# `run_pipeline_once_loop.sh` does the same automatically (sleeps between the
# passes), and `run_pipeline.sh` combines both into one blocking command.
#
# Usage (from this directory):
#
#     export AUTOMAG_HOME=/path/to/automag      # source checkout
#     # or: export AUTOMAG_EXE=/shared/bin/run_automag   # compiled binary
#     ./run_pipeline_once.sh
#
# Suitable for cron, where a new pass every few minutes is all that is needed:
#
#     */20 * * * * cd /work/project && /path/to/automag/examples/vasp/run_pipeline_once.sh >> automag.log 2>&1
#
# Extra arguments are forwarded to the pipeline, e.g.
#
#     ./run_pipeline_once.sh --dry-run
#     ./run_pipeline_once.sh --stages coll --steps submit
#
# Environment variables: AUTOMAG_HOME / AUTOMAG_EXE, AUTOMAG_CONDA_ENV,
# AUTOMAG_CONDA_SH - as in run_pipeline.sh.
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

echo "automag: ${COMMAND[*]}"
echo "config:  ${PWD}/${CONFIG}"
echo "pass:    $(date -Is)"

set +e
"${COMMAND[@]}" "$CONFIG" --once "$@"
status=$?
set -e
if [ "$status" -ne 0 ]; then
    exit "$status"
fi
echo
echo "One pass finished.  Run ./run_pipeline_once.sh again to continue"
echo "(or use ./run_pipeline.sh to wait for the calculations in one command)."