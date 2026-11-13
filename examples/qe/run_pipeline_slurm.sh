#!/bin/bash
#SBATCH --job-name=automag-driver
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --output=automag-driver-%j.out
#SBATCH --error=automag-driver-%j.err

# automag pipeline launcher - SLURM driver job (Python sources).
#
# The pipeline is a long-running supervisor: it submits the calculations of a
# step, waits for them, runs the analysis and continues with the next step.
# This script makes *the supervisor itself* a SLURM job, so it survives the
# logout and does not occupy the login node.  The calculations it launches are
# ordinary SLURM jobs, started with sbatch because the step headers in
# input_automag.py contain #SBATCH directives (set `"scheduler": "slurm"` in
# PIPELINE["hosts"] if your headers do not).
#
# Usage (from this directory, on a machine with sbatch):
#
#     export AUTOMAG_HOME=/path/to/automag      # source checkout
#     ./run_pipeline_slurm.sh                   # submits itself, prints the job id
#     squeue -u $USER                            # find it
#     tail -f automag-driver-*.out               # follow it
#
# Extra arguments are forwarded to the pipeline and therefore to the job:
#
#     ./run_pipeline_slurm.sh --dry-run          # the job only prints the plan
#     ./run_pipeline_slurm.sh --stages coll --steps submit
#
# Resources: the defaults are in the #SBATCH lines above.  Override them at
# submit time without editing the file:
#
#     AUTOMAG_SLURM_OPTS="--partition=compute --account=myproj --time=2-00:00:00" \
#         ./run_pipeline_slurm.sh
#
# Environment (all optional):
#
#     AUTOMAG_HOME        automag checkout (default: this example's ../../)
#     AUTOMAG_CONDA_ENV   conda env with the dependencies (default: automag,
#                         empty string = do not activate anything)
#     AUTOMAG_CONDA_SH    path to conda.sh when `conda` is not on PATH
#     AUTOMAG_PRE_COMMAND shell line run before activating the env, e.g.
#                         "module load python/3.12" (executed with eval)
#     AUTOMAG_SLURM_OPTS  extra sbatch options (see above)
#
# The driver job must be allowed to call sbatch itself.  If your site forbids
# that, run ./run_pipeline.sh on the login node (nohup) instead.
set -euo pipefail

# ---------------------------------------------------------------------------
# Phase 1 - started by hand: submit this very file as a SLURM job.
# ---------------------------------------------------------------------------
if [ -z "${AUTOMAG_SLURM_JOB:-}" ]; then
    cd "$(dirname "$0")"
    if ! command -v sbatch >/dev/null 2>&1; then
        echo "ERROR: sbatch not found: there is no SLURM here." >&2
        echo "       Use ./run_pipeline.sh to run the pipeline directly." >&2
        exit 1
    fi
    # shellcheck disable=SC2086  # AUTOMAG_SLURM_OPTS is a list of options
    exec sbatch ${AUTOMAG_SLURM_OPTS:-} \
         --export=ALL,AUTOMAG_SLURM_JOB=1,AUTOMAG_SLURM_DIR="$PWD" \
         --chdir="$PWD" \
         "$0" "$@"
fi

# ---------------------------------------------------------------------------
# Phase 2 - running as a SLURM job: start the supervisor.
# ---------------------------------------------------------------------------
cd "${AUTOMAG_SLURM_DIR:-$PWD}"
echo "automag driver: job ${SLURM_JOB_ID:-?} on ${SLURMD_NODENAME:-$(hostname)}"
echo "working directory: $PWD"
echo "started: $(date -Is)"

# Site-specific setup (modules, spack, ...) before the Python environment.
if [ -n "${AUTOMAG_PRE_COMMAND:-}" ]; then
    eval "${AUTOMAG_PRE_COMMAND}"
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

# --- locate the pipeline: an explicit checkout or this example's ../../ ----
if [ -z "${AUTOMAG_HOME:-}" ] && [ -f ../../run_automag.py ]; then
    AUTOMAG_HOME="$(cd ../.. && pwd)"
fi
if [ -z "${AUTOMAG_HOME:-}" ]; then
    echo "ERROR: set AUTOMAG_HOME to the automag checkout." >&2
    exit 1
fi
# the job may start in another directory, so store an absolute path
AUTOMAG_HOME="$(cd "${AUTOMAG_HOME}" && pwd)"
export AUTOMAG_HOME

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

echo "automag: ${AUTOMAG_HOME}"
echo "config:  ${PWD}/${CONFIG}"

set +e
python "${AUTOMAG_HOME}/run_automag.py" "$CONFIG" --wait "$@"
status=$?
set -e
echo "finished: $(date -Is) (exit ${status})"
exit "${status}"