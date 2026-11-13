#!/usr/bin/env python
"""
automag.tools.mae_chain_test
============================

Runs the **complete 4_mae chain** (``1_submit`` → ``2_plot_results`` →
``3_submit`` → ``4_plot_results``) without any DFT code, and checks the
parameter-sweep (convergence study) support.

The "calculation" is a shell one-liner that writes a syntactically valid
``OSZICAR``/``OUTCAR`` pair with an energy derived from the calculation
directory, so the post-processing scripts really parse the numbers and produce
real output files.  The whole pipeline is then driven twice with
``run_automag.py --wait``:

* once with a single parameter set (the historical behaviour, filenames without
  a suffix);
* once with ``"encut": [400, 500]``, which must produce one MAE calculation per
  set (``..._EN400_S01``/``..._EN500_S02`` as separate directories), run them in
  parallel and write the convergence summary.

    source /path/to/miniconda3/bin/activate automag
    python tools/mae_chain_test.py            # exit code 0 on success
    python tools/mae_chain_test.py --keep     # keep the temporary projects

.. codeauthor:: Automag contributors
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FORMULA = "Fe4O6"

#: Write a valid VASP-style OSZICAR/OUTCAR pair.  ``$PWD`` decides the energy,
#: so every angle/parameter set gets a different number and the MAE is non-zero.
FAKE_CALC = (
    "sh -c 'seed=$(printf \"%s\" \"$PWD\" | cksum | tr -dc 0-9 | tail -c 6); "
    "e=$(awk -v s=\"$seed\" \"BEGIN{printf \\\"%.6f\\\", -100 + (s%1000)/100.0}\"); "
    "printf \"       N       E                     dE             d eps       "
    "ncg     rms          rms(c)\\n\" > OSZICAR; "
    "printf \"DAV:   1     %s    %s   -0.456E+03   1234   0.123E+00\\n\" "
    "\"$e\" \"$e\" >> OSZICAR; "
    "printf \"   1 F= %s E0= %s  d E =0.0  mag= 4.0 4.0 -4.0 -4.0\\n\" "
    "\"$e\" \"$e\" >> OSZICAR; "
    "echo \"General timing and accounting informations for this job:\" > OUTCAR'"
)

CONFIG = '''"""Temporary MAE project written by tools/mae_chain_test.py."""

PARAMS = {
    "xc": "PBE", "prec": "Accurate", "ncore": 1,
    "encut": %(encut)s,
    "ediff": 1e-6, "ismear": 1, "sigma": 0.1, "nelm": 20, "kpts": %(kpts)s,
    "lcharg": False, "lwave": False, "isym": 0,
    "ldau": True, "ldautype": 2, "ldaul": [2, -1], "ldauu": [4.0, 0],
    "ldauj": [0, 0], "ldauprint": 2, "lreal": "Auto",
}

PIPELINE = {
    "workdir": ".",
    "structure": "structure.vasp",
    "state_file": "state.sqlite",
    "poll_interval": 1,
    "wait_timeout": 900,
    "run": [
        {"id": "mae_vasp", "stage": "mae", "calculator": "vasp",
         "coll_calculator": "prophet_spin", "configuration": "afm1",
         "Nph": 2, "Nth": 2, "N_MAE": 2,
         "standardize_cell": True, "use_primitive_cell": True,
         "spin_values": {"Fe": [4]}, "params": PARAMS,
         "max_parallel": 4},
    ],
    "submit": {
        "use_fireworks": False,
        "auto_submit": True,
        "jobheader": "#!/bin/bash",
        "calculator_command": %(command)s,
        "environment_activate": "",
        "environment_deactivate": "",
        "parallel_over_configurations": True,
    },
}
'''


def build_project(project, encut="300", kpts="2"):
    structure = ROOT / "examples" / "vasp" / "structure.vasp"
    (project / "geometries").mkdir(parents=True)
    (project / "CalcFold").mkdir()
    shutil.copy2(structure, project / "structure.vasp")
    shutil.copy2(structure, project / "geometries" / "structure.vasp")

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
        "encut": encut, "kpts": kpts, "command": repr(FAKE_CALC)})


def run_pipeline(project, *args, timeout=1800):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT)
    return subprocess.run(
        [sys.executable, str(ROOT / "run_automag.py"), *args],
        cwd=str(project), env=env, capture_output=True, text=True,
        timeout=timeout)


def check_single_set(project, failures):
    result = run_pipeline(project, "--wait")
    output = result.stdout + result.stderr
    if result.returncode != 0:
        failures.append("the single-set MAE run failed: " + output[-2000:])
        return
    state_dirs = glob.glob(str(project / "CalcFold" / FORMULA / "vasp" /
                             "afm1" / "mae_*"))
    if len(state_dirs) != 1:
        failures.append(f"expected exactly one MAE directory, got "
                        f"{[os.path.basename(d) for d in state_dirs]}")
    checks = {
        "MAE_theta_phi_*.txt": "MAE_theta_phi_{}_U4.0_J0.0_K2_EN300.txt",
        "MAE_curve_*.txt": "MAE_curve_{}_U4.0_J0.0_K2_EN300.txt",
    }
    for pattern, name in checks.items():
        expected = name.format("afm1")
        if not (project / "4_mae" / expected).is_file():
            found = [os.path.basename(p) for p in project.glob("4_mae/" + pattern)]
            failures.append(f"{expected} was not produced (found: {found})")
    mae_file = project / "4_mae" / "MAE_theta_phi_afm1_U4.0_J0.0_K2_EN300.txt"
    if mae_file.is_file():
        text = mae_file.read_text()
        for key in ("MAE_x:", "MAE_z:", "MAE_y:", "MAE:"):
            if key not in text:
                failures.append(f"{key} is missing from {mae_file.name}")
    curve = project / "4_mae" / "MAE_curve_afm1_U4.0_J0.0_K2_EN300.txt"
    if curve.is_file():
        text = curve.read_text()
        for key in ("MAE =", "K1", "M0"):
            if key not in text:
                failures.append(f"{key} is missing from {curve.name}")


def check_sweep(project, failures):
    result = run_pipeline(project, "--once")
    output = result.stdout + result.stderr
    if result.returncode != 0:
        failures.append("the sweep --once failed: " + output[-2000:])
        return
    if "MAE convergence study over 2 parameter sets" not in output:
        failures.append("the parameter lists were not recognised as a sweep: "
                        + output[-1500:])
    result = run_pipeline(project, "--wait")
    output = result.stdout + result.stderr
    if result.returncode != 0:
        failures.append("the sweep --wait failed: " + output[-2500:])
        return

    base = project / "CalcFold" / FORMULA / "vasp" / "afm1"
    dirs = sorted(os.path.basename(path) for path in
                  glob.glob(str(base / "mae_*")))
    if len(dirs) != 2:
        failures.append(f"expected 2 MAE directories for 2 ENCUT values, "
                        f"got {dirs}")
        return
    if not any("_S01" in name for name in dirs) or not any(
            "_S02" in name for name in dirs):
        failures.append(f"the sweep directories are not tagged: {dirs}")
    for dirname in dirs:
        for suffix in ("_S01", "_S02"):
            if suffix in dirname:
                break
    tables = sorted(os.path.basename(path) for path in
                    project.glob("4_mae/MAE_theta_phi_afm1_S*.txt"))
    if len(tables) != 2:
        failures.append(f"expected 2 MAE_theta_phi tables, got {tables}")
    curves = sorted(os.path.basename(path) for path in
                    project.glob("4_mae/MAE_curve_afm1_S*.txt"))
    if len(curves) != 2:
        failures.append(f"expected 2 MAE_curve tables, got {curves}")
    summary = project / "4_mae" / "MAE_convergence_afm1.txt"
    if not summary.is_file():
        failures.append("the convergence summary of the theta/phi stage is "
                        "missing")
    else:
        text = summary.read_text()
        if "encut" not in text or "relative change" not in text:
            failures.append("the convergence summary does not list the swept "
                            "parameter and the change between the sets")
    curve_summary = project / "4_mae" / "MAE_curve_convergence_afm1.txt"
    if not curve_summary.is_file():
        failures.append("the convergence summary of the MAE curve stage is "
                        "missing")
    # every set must be complete
    for dirname in dirs:
        phth = glob.glob(str(base / dirname / "PhTh_*"))
        if len(phth) != 9:
            failures.append(f"{dirname} has {len(phth)} angle calculations, "
                            "expected 9")
        kdirs = glob.glob(str(base / dirname / "K_*_RtMAE_*"))
        if len(kdirs) != 3:
            failures.append(f"{dirname} has {len(kdirs)} MAE curve "
                            "calculations, expected 3")
    # the pipeline state must know about both sets
    print("  sweep directories:", dirs)
    for path in sorted(project.glob("4_mae/*convergence*.txt")):
        print("  " + path.name + ":")
        print("    " + path.read_text().strip().replace("\n", "\n    "))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--keep", action="store_true")
    args = parser.parse_args(argv)

    base = Path(tempfile.mkdtemp(prefix="automag_mae_"))
    failures = []
    try:
        single = base / "single"
        build_project(single)
        check_single_set(single, failures)

        sweep = base / "sweep"
        build_project(sweep, encut="[400, 500]")
        check_sweep(sweep, failures)

        if failures:
            print("FAILED:")
            for item in failures:
                print(" -", item)
            return 1
        print("mae chain test OK: full 1-4 stage chain, parameter sweep, "
              "convergence summary")
        return 0
    finally:
        if args.keep:
            print(f"projects kept in {base}")
        else:
            shutil.rmtree(base, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
