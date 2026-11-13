#!/bin/bash
#SBATCH --job-name=automag-driver
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --output=automag-driver-%j.out
#SBATCH --error=automag-driver-%j.err

# automag pipeline launcher - SLURM driver job (compiled build).
#
# Same as run_pipeline_slurm.sh, but the supervisor is the compiled automag
# binary (dist/run_automag) instead of the Python sources: it bundles Python,
# the dependencies, the pipeline and the stage scripts, so the compute node
# needs no conda environment.  The calculations it launches are ordinary SLURM
# jobs, started with sbatch because the step headers in input_automag.py
# contain #SBATCH directives (set `"scheduler": "slurm"` in PIPELINE["hosts"]
# if your headers do not).
#
# Usage (from this directory, on a machine with sbatch):
#
#     export AUTOMAG_EXE=/shared/bin/run_automag   # compiled build
#     ./run_pipeline_slurm_compiled.sh             # submits itself, prints the job id
#     squeue -u $USER                               # find it
#     tail -f automag-driver-*.out                  # follow it
#
# `run_automag` on PATH is used automatically.  Build the binary in the automag
# checkout with `./build.sh` (docs/INSTALL.md, section 5).
#
# Extra arguments are forwarded to the pipeline and therefore to the job:
#
#     ./run_pipeline_slurm_compiled.sh --dry-run
#     ./run_pipeline_slurm_compiled.sh --stages coll --steps submit
#
# Resources: the defaults are in the #SBATCH lines above.  Override them at
# submit time without editing the file:
#
#     AUTOMAG_SLURM_OPTS="--partition=compute --account=myproj --time=2-00:00:00" \
#         ./run_pipeline_slurm_compiled.sh
#
# Environment (all optional):
#
#     AUTOMAG_EXE              the compiled binary (required)
#     AUTOMAG_EXTERNAL_CODES   registry of enumlib/VAMPIRE/ESpinS; a
#                              codes.json next to the binary is picked up
#                              automatically
#     AUTOMAG_PRE_COMMAND      shell line run at the start of the job, e.g.
#                              "module load intel" (executed with eval)
#     AUTOMAG_SLURM_OPTS       extra sbatch options (see above)
#
# The driver job must be allowed to call sbatch itself.  If your site forbids
# that, run ./run_pipeline_compiled.sh on the login node (nohup) instead.
set -euo pipefail

# ---------------------------------------------------------------------------
# Phase 1 - started by hand: submit this very file as a SLURM job.
# ---------------------------------------------------------------------------
if [ -z "${AUTOMAG_SLURM_JOB:-}" ]; then
    cd "$(dirname "$0")"
    if ! command -v sbatch >/dev/null 2>&1; then
        echo "ERROR: sbatch not found: there is no SLURM here." >&2
        echo "       Use ./run_pipeline_compiled.sh to run the pipeline directly." >&2
        exit 1
    fi
    if [ -z "${AUTOMAG_EXE:-}" ]; then
        if command -v run_automag >/dev/null 2>&1; then
            AUTOMAG_EXE="$(command -v run_automag)"
        else
            echo "ERROR: set AUTOMAG_EXE to the compiled automag binary." >&2
            echo "       export AUTOMAG_EXE=/path/to/automag/dist/run_automag" >&2
            echo "       (build it with './build.sh', or use run_pipeline_slurm.sh" >&2
            echo "        for a source checkout)" >&2
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
    # shellcheck disable=SC2086  # AUTOMAG_SLURM_OPTS is a list of options
    exec sbatch ${AUTOMAG_SLURM_OPTS:-} \
         --export=ALL,AUTOMAG_SLURM_JOB=1,AUTOMAG_SLURM_DIR="$PWD",AUTOMAG_EXE="$AUTOMAG_EXE" \
         --chdir="$PWD" \
         "$0" "$@"
fi

# ---------------------------------------------------------------------------
# Phase 2 - running as a SLURM job: start the supervisor.
# ---------------------------------------------------------------------------
cd "${AUTOMAG_SLURM_DIR:-$PWD}"
echo "automag driver (compiled): job ${SLURM_JOB_ID:-?} on ${SLURMD_NODENAME:-$(hostname)}"
echo "working directory: $PWD"
echo "started: $(date -Is)"

# Site-specific setup (modules, licences, ...).
if [ -n "${AUTOMAG_PRE_COMMAND:-}" ]; then
    eval "${AUTOMAG_PRE_COMMAND}"
fi

if [ -z "${AUTOMAG_EXE:-}" ] || [ ! -x "${AUTOMAG_EXE}" ]; then
    echo "ERROR: AUTOMAG_EXE is not set to an executable binary." >&2
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

echo "automag: ${AUTOMAG_EXE}"
echo "config:  ${PWD}/${CONFIG}"

set +e
"${AUTOMAG_EXE}" "$CONFIG" --wait "$@"
status=$?
set -e
echo "finished: $(date -Is) (exit ${status})"
exit "${status}"