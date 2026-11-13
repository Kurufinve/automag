#!/usr/bin/env python
"""
automag.tools.execution_test
============================

End-to-end test of the execution engine: launch modes, parallel calculations,
the SQLite state and the resumable ``--once`` mode.

No DFT code and no real cluster are needed:

* the "calculation" is a **fake calculator command** that only writes the files
  a finished run would leave behind, so the generated job scripts complete
  instantly;
* the real stage script (``4_mae/1_submit.py``) generates the real job scripts,
  so the task manifest, the ``z``-before-grid ordering and the parallel limit
  are exercised for real;
* **fake ``sbatch``/``squeue``** scripts on ``PATH`` really start the job
  scripts in the background and report their state, so the scheduler code path
  is exercised too.

Scenarios
---------
1. ``--once`` starts the 4_mae reference calculation first, never more than
   ``max_parallel`` calculations at a time, records host + job id in the state
   file and exits.
2. Continuing with ``--once`` finishes the grid; every calculation is started
   exactly once and the step ends up ``done``.
3. Changing a setting resets only that step and recomputes it (the stale
   results are not reused).
4. A ``header`` with ``#SBATCH`` directives sends the calculations through the
   task manager (fake ``sbatch``) instead of starting plain processes.
5. A step can run on another "machine" (``transport: shell`` with a wrapper).

    source /path/to/miniconda3/bin/activate automag
    python tools/execution_test.py           # exit code 0 on success
    python tools/execution_test.py --keep    # keep the temporary projects

.. codeauthor:: Automag contributors
"""

import argparse
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FORMULA = "Fe4O6"

#: Fake "calculator": the job script already changed into ``singlepoint``.
FAKE_CALC = (
    'sh -c \'echo "General timing and accounting informations for this job:" '
    "> OUTCAR; "
    'echo "   1 F= -.1E+02 E0= -.1E+02  d E =0.0  mag= 4.0 4.0 -4.0 -4.0" '
    "> OSZICAR'"
)

FAKE_SBATCH = """#!/bin/sh
# fake sbatch: start the job script in the background and remember its state
job="$1"
dir=$(pwd)
id=$(( $(cat "$AUTOMAG_FAKE_QUEUE.count" 2>/dev/null || echo 1000) + 1 ))
echo "$id" > "$AUTOMAG_FAKE_QUEUE.count"
printf '%s %s\\n' "$id" "$dir/$job" >> "$AUTOMAG_FAKE_QUEUE"
( cd "$dir" && sh "$job" > /dev/null 2>&1 ) &
printf 'Submitted batch job %s\\n' "$id"
"""

FAKE_SQUEUE = """#!/bin/sh
# fake squeue: -j <id> -> state, empty output when the job is unknown
id=""
while [ $# -gt 0 ]; do
    case "$1" in
        -j) shift; id="$1" ;;
    esac
    shift
done
entry=$(grep "^$id " "$AUTOMAG_FAKE_QUEUE" 2>/dev/null | tail -1)
[ -z "$entry" ] && exit 0
script=$(echo "$entry" | cut -d' ' -f2-)
marker="$(dirname "$script")/.automag_done"
if [ -f "$marker" ]; then
    echo "COMPLETED"
elif [ -f "$AUTOMAG_FAKE_QUEUE.done.$id" ]; then
    echo "COMPLETED"
else
    echo "RUNNING"
fi
"""

CONFIG = '''"""Temporary pipeline written by tools/execution_test.py."""

PARAMS = {
    "xc": "PBE", "prec": "Accurate", "ncore": 1, "encut": %(encut)s,
    "ediff": 1e-6, "ismear": 1, "sigma": 0.1, "nelm": 20, "kpts": 2,
    "lcharg": False, "lwave": False, "isym": 0,
    "ldau": True, "ldautype": 2, "ldaul": [2, -1], "ldauu": [4.0, 0],
    "ldauj": [0, 0], "ldauprint": 2, "lreal": "Auto",
}

PIPELINE = {
    "workdir": ".",
    "structure": "structure.vasp",
    "state_file": "state.sqlite",
    "poll_interval": 1,
    "wait_timeout": 600,
    %(hosts)s
    "run": [
        {"id": "mae_vasp", "stage": "mae", "calculator": "vasp",
         "coll_calculator": "prophet_spin", "configuration": "afm1",
         "Nph": 1, "Nth": 1, "N_MAE": 1,
         "standardize_cell": True, "use_primitive_cell": True,
         "spin_values": {"Fe": [4]}, "params": PARAMS,
         "max_parallel": %(max_parallel)s,
         "execution": %(execution)s%(host)s},
    ],
    "submit": {
        "use_fireworks": False,
        "auto_submit": True,
        "jobheader": %(header)s,
        "calculator_command": %(command)s,
        "environment_activate": "",
        "environment_deactivate": "",
        "parallel_over_configurations": True,
    },
}
'''


def build_project(project, *, max_parallel=1, execution='"process"',
                  header='""', command=None, encut=300, host=None,
                  hosts=""):
    structure = ROOT / "examples" / "vasp" / "structure.vasp"
    (project / "geometries").mkdir(parents=True)
    (project / "CalcFold").mkdir()
    shutil.copy2(structure, project / "structure.vasp")
    shutil.copy2(structure, project / "geometries" / "structure.vasp")

    # a fabricated 2_coll bundle, exactly like the smoke test does
    bundle = project / "2_coll" / "coll_prophet" / f"{FORMULA}_prophet_spin"
    trials = bundle / f"trials_{FORMULA}"
    trials.mkdir(parents=True)
    shutil.copy2(structure, trials / "setting001.vasp")
    (trials / "configurations001.txt").write_text("afm1 1 -1 1 -1\n")
    (bundle / f"{FORMULA}_prophet_spin_most_stable_configuration.txt").write_text(
        "afm1")
    (project / "CalcFold" / f"{FORMULA}_singlepoint_prophet_spin.txt").write_text(
        "afm1       1  singlepoint=convergedconverged  Fe4O6       P1  "
        "energy=-0.1\n"
        "          initial_magmoms=[ 4. -4.  4. -4.]  "
        "final_magmoms=[ 3.9 -3.8  3.7 -3.6]\n")

    (project / "input_automag.py").write_text(CONFIG % {
        "max_parallel": max_parallel,
        "execution": execution,
        "header": header,
        "command": command or repr(FAKE_CALC),
        "encut": encut,
        "hosts": hosts,
        "host": f',\n         "host": "{host}"' if host else "",
    })


def fake_scheduler(bin_dir):
    bin_dir.mkdir(parents=True, exist_ok=True)
    for name, text in (("sbatch", FAKE_SBATCH), ("squeue", FAKE_SQUEUE)):
        path = bin_dir / name
        path.write_text(text)
        path.chmod(0o755)


def run_pipeline(project, bin_dir=None, *args, timeout=900):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT)
    env["AUTOMAG_FAKE_QUEUE"] = str(bin_dir / "queue.txt") if bin_dir else ""
    if bin_dir:
        env["PATH"] = str(bin_dir) + os.pathsep + env["PATH"]
    return subprocess.run(
        [sys.executable, str(ROOT / "run_automag.py"), *args],
        cwd=str(project), env=env, capture_output=True, text=True,
        timeout=timeout)


def _rows(path, table):
    if not Path(path).is_file():
        return []
    try:
        with sqlite3.connect(str(path)) as connection:
            connection.row_factory = sqlite3.Row
            return [dict(row) for row in
                    connection.execute(f"SELECT * FROM {table}")]
    except sqlite3.Error:
        return []


def items(project):
    return _rows(project / "state.sqlite", "items")


def tasks(project):
    return _rows(project / "state.sqlite", "tasks")


def steps(project):
    return {row["step_key"]: row
            for row in _rows(project / "state.sqlite", "steps")}


def markers(project):
    return sorted(project.glob("CalcFold/**/.automag_done"))


# ---------------------------------------------------------------------------
# Scenario 1: --once barrier, parallel limit, resume without recomputation
# ---------------------------------------------------------------------------
def scenario_once_and_resume(project, failures, max_parallel=2):
    bin_dir = project / "bin"
    fake_scheduler(bin_dir)

    first = run_pipeline(project, bin_dir, "--once")
    output = first.stdout + first.stderr
    if first.returncode != 0:
        failures.append("--once failed: " + output[-1200:])
        return

    registered = items(project)
    if len(registered) != 5:
        failures.append(f"expected 5 calculations (z + 2x2 grid), "
                        f"got {len(registered)}")
    started = markers(project)
    if not started:
        failures.append("--once started no calculation")
    if len(started) > max_parallel:
        failures.append(f"--once started {len(started)} calculations, "
                        f"max_parallel is {max_parallel}")
    recorded = tasks(project)
    if not recorded or not any(task["external_id"] for task in recorded):
        failures.append("the started job has no process id in the state file")
    elif not json.loads(recorded[0]["items"])[0].endswith("/z"):
        failures.append("the first calculation started was not the reference "
                        f"one: {json.loads(recorded[0]['items'])}")
    if "has to wait for external programs" not in output:
        failures.append("--once did not report that it waits for the jobs")

    # continue; the second call must start the grid but not repeat anything
    time.sleep(1.0)
    second = run_pipeline(project, bin_dir, "--once")
    if second.returncode != 0:
        failures.append("the second --once failed: "
                        + (second.stderr or second.stdout)[-1200:])
    if len(markers(project)) > 2 * max_parallel:
        failures.append(f"the parallel limit was exceeded: "
                        f"{len(markers(project))} calculations started")

    # finish the step
    result = run_pipeline(project, bin_dir, "--steps", "submit", "--wait",
                         timeout=900)
    if result.returncode != 0:
        failures.append("--wait did not finish the step: "
                        + (result.stderr or result.stdout)[-1500:])
        return
    if len(markers(project)) != 5:
        failures.append(f"expected 5 finished calculations, "
                        f"got {len(markers(project))}")
    if steps(project)["mae_vasp.submit"]["status"] != "done":
        failures.append("mae_vasp.submit is not 'done' after all jobs finished")
    attempts = {item["path"]: item["attempts"] for item in items(project)}
    if any(value > 1 for value in attempts.values()):
        failures.append(f"some calculation was started more than once: "
                        f"{attempts}")


# ---------------------------------------------------------------------------
# Scenario 2: a changed setting resets the step and recomputes it
# ---------------------------------------------------------------------------
def scenario_config_change(project, failures):
    """A changed setting resets the step; the same calculations are recomputed.

    ``nelm`` does not appear in the MAE directory name, so the calculation
    directories stay the same and the test proves that the stale results are
    *not* reused (there is no new directory to hide behind).
    """
    bin_dir = project / "bin"
    before = {path: path.stat().st_mtime for path in markers(project)}
    if len(before) != 5:
        failures.append("the project did not finish before the settings changed")
        return
    if steps(project)["mae_vasp.submit"]["status"] != "done":
        failures.append("the step was not 'done' before the settings changed")

    time.sleep(0.1)
    changed_at = time.time()
    (project / "input_automag.py").write_text(
        (project / "input_automag.py").read_text().replace(
            '"nelm": 20', '"nelm": 25'))
    result = run_pipeline(project, bin_dir, "--once", "--steps", "submit")
    output = result.stdout + result.stderr
    if "resetting step mae_vasp.submit" not in output:
        failures.append("a changed setting did not reset the step: "
                        + output[-900:])
        return
    rows = items(project)
    if len(rows) != 5:
        failures.append(f"the reset step has {len(rows)} calculations, "
                        "expected 5")
    if any(row["status"] == "done" for row in rows):
        failures.append("stale results were marked done after the reset "
                        "(they are older than the reset)")
    result = run_pipeline(project, bin_dir, "--steps", "submit", "--wait",
                          timeout=900)
    if result.returncode != 0:
        failures.append("--wait after the settings change failed: "
                        + (result.stderr or result.stdout)[-1200:])
        return
    stale = [path for path in markers(project)
             if path.stat().st_mtime < changed_at]
    if stale:
        failures.append(f"{len(stale)} of 5 calculation(s) were not recomputed "
                        "after the settings change")
    if steps(project)["mae_vasp.submit"]["status"] != "done":
        failures.append("the reset step did not finish")


# ---------------------------------------------------------------------------
# Scenario 3: batch directives -> the task manager starts the calculations
# ---------------------------------------------------------------------------
def scenario_scheduler(project, failures):
    bin_dir = project / "bin"
    fake_scheduler(bin_dir)
    result = run_pipeline(project, bin_dir, "--steps", "submit", "--wait",
                         timeout=900)
    output = result.stdout + result.stderr
    if result.returncode != 0:
        failures.append("--wait with a SLURM header failed: "
                        + output[-1500:])
        return
    kinds = {task["kind"] for task in tasks(project)}
    if kinds != {"slurm"}:
        failures.append(f"the calculations did not go through sbatch "
                        f"(task kinds: {kinds})")
    if not (bin_dir / "queue.txt").is_file():
        failures.append("sbatch was never called")
    if steps(project)["mae_vasp.submit"]["status"] != "done":
        failures.append("the step is not 'done' after --wait")


# ---------------------------------------------------------------------------
# Scenario 4: a step on another "machine" (shell transport)
# ---------------------------------------------------------------------------
def scenario_remote(project, failures):
    result = run_pipeline(project, None, "--once")
    output = result.stdout + result.stderr
    if result.returncode != 0:
        failures.append("--once on a shell transport failed: "
                        + output[-1200:])
        return
    entry = steps(project).get("mae_vasp.submit") or {}
    if entry.get("host") != "box":
        failures.append(f"the step did not run on the configured host "
                        f"(recorded '{entry.get('host')}')")
    if not items(project):
        failures.append("the step on the remote host produced no calculations")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--keep", action="store_true",
                        help="keep the temporary projects for inspection")
    args = parser.parse_args(argv)

    base = Path(tempfile.mkdtemp(prefix="automag_exec_"))
    failures = []
    try:
        once = base / "once"
        build_project(once, max_parallel=2)
        scenario_once_and_resume(once, failures, max_parallel=2)
        scenario_config_change(once, failures)

        scheduler = base / "scheduler"
        build_project(
            scheduler, max_parallel=2,
            header='"""#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --ntasks=1\n"""')
        scenario_scheduler(scheduler, failures)

        remote = base / "remote"
        build_project(remote, host="box",
                      hosts='"hosts": {"box": {"transport": "shell", '
                            '"wrapper": "nice -n 0"}},')
        scenario_remote(remote, failures)

        if failures:
            print("FAILED:")
            for item in failures:
                print(" -", item)
            return 1
        print("execution test OK: --once barrier, resume, config change, "
              "scheduler launcher, remote host")
        return 0
    finally:
        if args.keep:
            print(f"projects kept in {base}")
        else:
            shutil.rmtree(base, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
