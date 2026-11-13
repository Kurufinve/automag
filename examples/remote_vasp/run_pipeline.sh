#!/bin/bash
# automag pipeline launcher.
#
# Copy this whole directory anywhere on the cluster and run:
#
#     ./run_pipeline.sh
#
# The script changes into its own directory first, so every result
# (0_conv/, 2_coll/, 2_relax/, 3_monte_carlo/, 4_mae/, CalcFold/, geometries/)
# is written next to input_automag.py and the automag checkout stays clean.
#
# Tell the launcher where automag is, either:
#
#     export AUTOMAG_HOME=/path/to/automag      # python source checkout
#     export AUTOMAG_EXE=/path/to/run_automag   # compiled single-file build
#
# and optionally
#
#     export AUTOMAG_CONDA_ENV=automag          # conda env with the dependencies
#
# Extra arguments are forwarded to the pipeline, e.g.
#
#     ./run_pipeline.sh --dry-run            # show the plan, run nothing
#     ./run_pipeline.sh --instances coll_vasp
set -euo pipefail

cd "$(dirname "$0")"

AUTOMAG_CONDA_ENV="${AUTOMAG_CONDA_ENV:-automag}"

# Locate the pipeline: an explicit binary, an explicit checkout, a binary on
# PATH, or the checkout this example was copied from.
if [ -z "${AUTOMAG_EXE:-}" ] && [ -z "${AUTOMAG_HOME:-}" ]; then
    if command -v run_automag >/dev/null 2>&1; then
        AUTOMAG_EXE="$(command -v run_automag)"
    elif [ -f "../../run_automag.py" ]; then
        AUTOMAG_HOME="$(cd ../.. && pwd)"
    fi
fi

AUTOMAG_CONDA_ENV="${AUTOMAG_CONDA_ENV-automag}"
if [ -n "${AUTOMAG_CONDA_ENV}" ]; then
    # conda's own activation hooks are not `set -u` safe (deactivate.d scripts
    # reference unset variables), so relax nounset while activating.
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

# Pick the controlling file: input_automag.py wins over .yaml, which wins
# over .json (same priority as the pipeline itself).
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

if [ -n "${AUTOMAG_EXE:-}" ]; then
    exec "${AUTOMAG_EXE}" "$CONFIG" --wait "$@"
fi

if [ -z "${AUTOMAG_HOME:-}" ]; then
    echo "ERROR: set AUTOMAG_HOME (source checkout) or AUTOMAG_EXE (compiled binary)." >&2
    exit 1
fi

exec python "${AUTOMAG_HOME}/run_automag.py" "$CONFIG" --wait "$@"
