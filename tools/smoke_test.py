#!/usr/bin/env python
"""
automag.tools.smoke_test
========================

Fast regression test that needs **no DFT code and no enumlib**.

It builds a throw-away project in a temporary directory, fabricates a
``2_coll`` bundle that looks like the output of a ``prophet_spin`` run, and
drives the real pipeline over it.  This exercises the parts that break most
easily and that a ``--dry-run`` cannot see:

* stage scripts really execute (``from input import *`` resolution);
* a ``relax``/``monte_carlo``/``mae`` instance reads **another** calculator's
  ``2_coll`` results via ``coll_calculator`` (``common.results``);
* ``reuse`` skips a finished stage;
* ``on_missing = "skip"`` skips a stage whose calculator produced no data and
  prints the advisory note listing the calculators that did;
* ``auto_submit = False`` writes inputs without launching jobs;
* the Monte Carlo backend is launched through the external-code registry
  (``common.external``), including its pre/post commands, and an unavailable
  code fails with an actionable message.

    source /path/to/miniconda3/bin/activate automag
    python tools/smoke_test.py            # exit code 0 on success
    python tools/smoke_test.py --keep     # keep the temporary project

.. codeauthor:: Automag contributors
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FORMULA = "Fe4O6"

CONFIG = '''"""Temporary pipeline written by tools/smoke_test.py."""

VASP_PARAMS = {
    "xc": "PBE", "prec": "Accurate", "ncore": 4, "encut": 400,
    "ediff": 1e-6, "ismear": 1, "sigma": 0.1, "nelm": 50, "kpts": 2,
    "lmaxmix": 4, "lcharg": False, "lwave": False, "isym": 0,
    "ldau": True, "ldautype": 2, "ldaul": [2, -1], "ldauu": [4.0, 0],
    "ldauj": [0, 0], "ldauprint": 2, "lreal": "Auto",
}

MC = {
    "configuration": "lowest",
    "cutoff_radius": 4.1,
    "control_group_size": 0.4,
    "append_coupling_constants": True,
    "magnetic_atoms": ["Fe"],
    "mc_backend": "vampire",
    "vampire_run_command": "vampire",
}

PIPELINE = {
    "workdir": ".",
    "structure": "structure.vasp",
    "submit": {
        "use_fireworks": False,
        "auto_submit": False,
        "jobheader": "#!/bin/bash",
        "calculator_command": "mpirun vasp_std",
        "environment_activate": "",
        "environment_deactivate": "",
        "parallel_over_configurations": True,
    },
    "run": [
        {"id": "coll_prophet", "stage": "coll", "calculator": "prophet_spin",
         "reuse": True, "params": {}, "supercell_size": 1,
         "spin_values": {"Fe": [4]}, "magnetic_atoms": ["Fe"]},
        {"id": "mae_vasp", "stage": "mae", "calculator": "vasp",
         "coll_calculator": "prophet_spin", "on_missing": "skip",
         "configuration": "lowest", "Nph": 2, "Nth": 2, "N_MAE": 2,
         "standardize_cell": True, "use_primitive_cell": True,
         "spin_values": {"Fe": [4]}, "params": VASP_PARAMS},
        {"id": "mc_vasp", "stage": "monte_carlo", "calculator": "vasp",
         "coll_calculator": "vasp", "on_missing": "skip", **MC},
    ],
}
'''


def build_project(project):
    structure = ROOT / "examples" / "vasp" / "structure.vasp"
    (project / "geometries").mkdir(parents=True)
    (project / "CalcFold").mkdir()
    shutil.copy2(structure, project / "structure.vasp")
    shutil.copy2(structure, project / "geometries" / "structure.vasp")

    bundle = project / "2_coll" / "coll_prophet" / f"{FORMULA}_prophet_spin"
    trials = bundle / f"trials_{FORMULA}"
    trials.mkdir(parents=True)
    shutil.copy2(structure, trials / "setting001.vasp")
    (trials / "configurations001.txt").write_text(
        "afm1 1 -1 1 -1\nfm1 1 1 1 1\n")
    (bundle / f"{FORMULA}_prophet_spin_setting001.vasp").write_text(
        structure.read_text())
    (bundle / f"{FORMULA}_prophet_spin_states001.txt").write_text(
        json.dumps([[1, -1, 1, -1], [1, 1, 1, 1]]))
    (bundle / f"{FORMULA}_prophet_spin_energies001.txt").write_text(
        json.dumps([0.0, 0.35]))
    (bundle / f"{FORMULA}_prophet_spin_most_stable_configuration.txt").write_text(
        "afm1")

    (project / "CalcFold" / f"{FORMULA}_singlepoint_prophet_spin.txt").write_text(
        "afm1       1  singlepoint=convergedconverged  Fe4O6       P1  "
        "energy=-0.1\n"
        "          initial_magmoms=[ 4. -4.  4. -4.]  "
        "final_magmoms=[ 3.9 -3.8  3.7 -3.6]\n")
    (project / "input_automag.py").write_text(CONFIG)


FAKE_VAMPIRE = """#!/bin/sh
echo "fake-vampire ran in $PWD" >> "$PWD/ran.log"
printf '# T M\\n1.0 1.0\\n2.0 0.5\\n' > output
"""

MC_INPUT = """mc_backend = 'vampire'
mc_timeout = 60
vampire_run_command = None
vampire_input_file = None
vampire_material_file = None
vampire_output_file = 'output'
external_codes = {
    'vampire': {
        'command': '%(fake)s',
        'pre_command': 'echo pre >> pre.log',
        'post_command': 'echo post >> post.log',
    }
}
"""


def check_external_codes(project):
    """Run 3_monte_carlo/5_run_mc.py against a fake VAMPIRE."""
    failures = []
    work = project / "external_codes"
    (work / "vampire").mkdir(parents=True)
    (work / "vampire" / "vamp.ucf").write_text("# fake ucf\n")
    (work / "CalcFold").mkdir()

    fake = work / "fake_vampire.sh"
    fake.write_text(FAKE_VAMPIRE)
    fake.chmod(0o755)
    (work / "input.py").write_text(MC_INPUT % {"fake": fake})

    env = os.environ.copy()
    env.update({"PYTHONPATH": str(ROOT), "AUTOMAG_CODE_PATH": str(ROOT),
                "AUTOMAG_PATH": str(work), "AUTOMAG_WORKDIR": str(work),
                "AUTOMAG_EXTERNAL_CODES": str(work / "missing.json")})
    result = subprocess.run(
        [sys.executable, str(ROOT / "3_monte_carlo" / "5_run_mc.py")],
        cwd=str(work), env=env, capture_output=True, text=True, timeout=180)
    if result.returncode != 0:
        failures.append("5_run_mc.py failed with a registered fake backend: "
                        + (result.stderr or result.stdout)[-600:])
    output = work / "output"
    if not output.is_file() or "2.0 0.5" not in output.read_text():
        failures.append("the Monte Carlo output was not copied to ./output")
    for marker in ("pre.log", "post.log", "ran.log"):
        if not (work / "vampire" / marker).is_file():
            failures.append(f"{marker} missing (pre/post commands not run)")

    # an unusable code must stop the run with an actionable message
    codes = work / "codes.json"
    codes.write_text(json.dumps(
        {"vampire": {"available": False, "error": "build failed"}}))
    (work / "input.py").write_text(
        "mc_backend = 'vampire'\nmc_timeout = 60\n")
    output.unlink()
    env["AUTOMAG_EXTERNAL_CODES"] = str(codes)
    result = subprocess.run(
        [sys.executable, str(ROOT / "3_monte_carlo" / "5_run_mc.py")],
        cwd=str(work), env=env, capture_output=True, text=True, timeout=180)
    combined = result.stdout + result.stderr
    if result.returncode == 0:
        failures.append("an unavailable external code did not stop the run")
    elif "./build.sh --external" not in combined:
        failures.append("the error does not mention './build.sh --external'")
    return failures


def check(project, output):
    failures = []
    if "skipping instance 'coll_prophet'" not in output:
        failures.append("'reuse' did not skip the finished coll instance")
    if "skipping instance 'mc_vasp'" not in output:
        failures.append("missing 2_coll data did not skip the monte_carlo "
                        "instance")
    if "2_coll results exist for calculator(s): prophet_spin" not in output:
        failures.append("no advisory note about the available 2_coll data")

    mae_dir = project / "4_mae"
    produced = sorted(mae_dir.glob("setting*_processed.vasp"))
    if not produced:
        failures.append("mae_vasp did not execute 4_mae/1_submit.py")
    inputs = mae_dir / "input.py"
    if not inputs.is_file() or "coll_calculator = 'prophet_spin'" not in inputs.read_text():
        failures.append("generated input.py does not carry coll_calculator")
    if "Reading the collinear results produced by: prophet_spin" not in output:
        failures.append("4_mae/1_submit.py did not read the prophet bundle")
    return failures


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--keep", action="store_true",
                        help="keep the temporary project for inspection")
    args = parser.parse_args(argv)

    project = Path(tempfile.mkdtemp(prefix="automag_smoke_"))
    try:
        build_project(project)
        env = os.environ.copy()
        env["PYTHONPATH"] = str(ROOT)
        result = subprocess.run(
            [sys.executable, str(ROOT / "run_automag.py"),
             "--steps", "submit"],
            cwd=str(project), env=env, capture_output=True, text=True,
            timeout=600)
        output = result.stdout + result.stderr
        if result.returncode != 0:
            print(output[-4000:])
            print(f"FAILED: the pipeline exited with {result.returncode}")
            return 1
        failures = check(project, output)
        failures += check_external_codes(project)
        if failures:
            print(output[-4000:])
            print("FAILED:")
            for item in failures:
                print(" -", item)
            return 1
        print("smoke test OK: cross-calculator reads, reuse/skip policies, "
              "external-code wiring")
        return 0
    finally:
        if args.keep:
            print(f"project kept in {project}")
        else:
            shutil.rmtree(project, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
