#!/bin/bash
# automag pipeline launcher - compiled (single-file) build.
#
# Identical to run_pipeline.sh, except that the pipeline is started from the
# compiled automag binary (dist/run_automag) instead of the Python sources.
# The binary bundles Python, the dependencies, the pipeline and the stage
# scripts, so no conda environment is required on the machine that runs it.
#
# Where is the binary?
#
#     export AUTOMAG_EXE=/shared/bin/run_automag      # this is what you need
#     ./run_pipeline_compiled.sh
#
# A `run_automag` found on PATH is used automatically.  Build the binary in the
# automag checkout with `./build.sh` (see docs/INSTALL.md, section 5); the
# result is dist/run_automag.
#
# Without AUTOMAG_EXE the launcher stops with an error instead of falling back
# to the sources - use run_pipeline.sh for a source checkout.
#
# Optional environment variables:
#
#     AUTOMAG_EXTERNAL_CODES=/path/to/external/codes.json
#         Registry of the external codes (enumlib/VAMPIRE/ESpinS).  A
#         codes.json next to the binary is picked up automatically.
#     AUTOMAG_CONDA_ENV=""            # deleted/ignored: the binary needs no env
#
# Extra arguments are forwarded to the pipeline, e.g.
#
#     ./run_pipeline_compiled.sh --dry-run          # show the plan, run nothing
#     ./run_pipeline_compiled.sh --stages coll --steps submit
#     ./run_pipeline_compiled.sh --once             # one pass, then exit
set -euo pipefail

cd "$(dirname "$0")"

# --- locate the compiled pipeline ------------------------------------------
if [ -z "${AUTOMAG_EXE:-}" ]; then
    if command -v run_automag >/dev/null 2>&1; then
        AUTOMAG_EXE="$(command -v run_automag)"
    else
        echo "ERROR: set AUTOMAG_EXE to the compiled automag binary." >&2
        echo "       export AUTOMAG_EXE=/path/to/automag/dist/run_automag" >&2
        echo "       (build it with './build.sh' in the automag checkout," >&2
        echo "        or use ./run_pipeline.sh for a source checkout)" >&2
        exit 1
    fi
fi
if [ ! -x "${AUTOMAG_EXE}" ]; then
    echo "ERROR: '${AUTOMAG_EXE}' is not an executable file." >&2
    exit 1
fi
# a relative AUTOMAG_EXE would break after `cd` (and the SLURM job needs an
# absolute path), so resolve it once, here
case "${AUTOMAG_EXE}" in
    /*) ;;
    *) AUTOMAG_EXE="$(cd "$(dirname "${AUTOMAG_EXE}")" && pwd)/$(basename "${AUTOMAG_EXE}")" ;;
esac
export AUTOMAG_EXE

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

echo "automag: ${AUTOMAG_EXE}"
echo "config:  ${PWD}/${CONFIG}"

exec "$AUTOMAG_EXE" "$CONFIG" --wait "$@"
