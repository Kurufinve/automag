"""
automag.3_monte_carlo.5_run_mc
==============================

Run the Monte Carlo backend (VAMPIRE or ESpinS) from a command in the input
file, so that ``monte_carlo.postprocess`` does not require an external step.

Backends
--------
``mc_backend = "vampire"`` (default)
    The script expects ``vampire/vamp.ucf`` (written by
    ``2_write_vampire_ucf.py``), copies the VAMPIRE input and material files
    into ``vampire/`` (defaults are the bundled files in
    ``3_monte_carlo/vampire_input``), runs ``vampire_run_command`` and copies
    the resulting ``output`` to the stage directory.

``mc_backend = "espins"``
    If ``espins/run_mc.sh`` does not exist yet, ``4_write_espins_mcin.py`` is
    executed first (it also runs the two ESpinS initialisation steps).  Then
    ``run_mc.sh`` is executed.  The result is copied to ``output``;
    ``espins_output_file`` selects the file name and
    ``espins_output_converter`` can be used for a custom output format.

.. codeauthor:: Automag contributors
"""

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from common.external import require, run as run_external, shell_command
from common.input_loader import load_input
from common.paths import code_root

DEFAULTS = {
    "mc_backend": "vampire",
    "mc_timeout": None,
    # None -> taken from external/codes.json (common.external), which also
    # carries the optional pre/post commands (module load, conda activate, ...).
    "vampire_run_command": None,
    "vampire_input_file": None,
    "vampire_material_file": None,
    "vampire_output_file": "output",
    "espins_run_command": None,
    "espins_output_file": None,
    "espins_output_converter": None,
    "external_codes": None,
}

settings = dict(DEFAULTS)
settings.update(load_input(sys.argv[1] if len(sys.argv) > 1 else None))

backend = str(settings.get("mc_backend", "vampire")).lower()
#: Per-run override of the external-code registry (see common.external).
overrides = settings.get("external_codes") or None


def run_command(command, cwd, label):
    """Run an external command, raising a clear error on failure."""
    if not command:
        raise ValueError(f"No command configured for {label}.")
    print(f"Running {label}: {command}")
    print(f"  cwd: {cwd}")
    subprocess.run(
        command,
        shell=True,
        cwd=str(cwd),
        check=True,
        executable="/bin/bash" if os.path.exists("/bin/bash") else None,
        timeout=settings.get("mc_timeout"),
    )


def resolve_file(value, default_relative):
    """Resolve a user-supplied path, falling back to a bundled default."""
    if value:
        path = Path(value).expanduser()
        if path.is_absolute():
            candidates = [path]
        else:
            candidates = [
                Path.cwd() / path,
                code_root() / path,
                code_root() / Path(default_relative).parent / path.name,
            ]
        for candidate in candidates:
            if candidate.is_file():
                return candidate
        tried = ", ".join(str(item) for item in candidates)
        raise FileNotFoundError(f"File '{value}' not found. Tried: {tried}")

    default = code_root() / default_relative
    if not default.is_file():
        raise FileNotFoundError(f"Default file not found: {default}")
    return default


def copy_output(source, destination="output"):
    shutil.copyfile(source, destination)
    print(f"Monte Carlo output copied to {Path(destination).resolve()}")


def run_vampire():
    vampire_dir = Path("vampire")
    unit_cell = vampire_dir / "vamp.ucf"
    if not unit_cell.is_file():
        raise FileNotFoundError(
            f"{unit_cell} not found. Run the 'submit' phase of monte_carlo "
            "(2_write_vampire_ucf.py) first."
        )

    input_file = resolve_file(
        settings.get("vampire_input_file"),
        "3_monte_carlo/vampire_input/input",
    )
    material_file = resolve_file(
        settings.get("vampire_material_file"),
        "3_monte_carlo/vampire_input/vamp.mat",
    )

    input_text = input_file.read_text()
    material_name = material_file.name
    match = re.search(r"material:file\s*=\s*(\S+)", input_text)
    if match:
        material_name = match.group(1)

    unit_cell_name = "vamp.ucf"
    match = re.search(r"material:unit-cell-file\s*=\s*(\S+)", input_text)
    if match:
        unit_cell_name = match.group(1)

    shutil.copyfile(input_file, vampire_dir / "input")
    shutil.copyfile(material_file, vampire_dir / material_name)
    if unit_cell_name != "vamp.ucf":
        shutil.copyfile(unit_cell, vampire_dir / unit_cell_name)

    require("vampire", "the Monte Carlo simulation (3_monte_carlo)",
            overrides=overrides)
    try:
        run_external(
            "vampire", cwd=vampire_dir, timeout=settings.get("mc_timeout"),
            overrides=overrides,
            override=settings.get("vampire_run_command"),
        )
    except subprocess.TimeoutExpired as error:
        raise SystemExit(
            "VAMPIRE did not finish within 'mc_timeout' "
            f"({settings.get('mc_timeout')} s).") from error

    output = Path(settings.get("vampire_output_file") or "output")
    if not output.is_absolute():
        output = vampire_dir / output
    if output.is_file() and output.stat().st_size > 0:
        copy_output(output)
    elif Path("output").is_file():
        print("Using existing ./output produced by the VAMPIRE run.")
    else:
        raise FileNotFoundError(
            f"VAMPIRE finished but no output file was found (looked for "
            f"{output} and ./output). Set 'vampire_output_file' if VAMPIRE "
            "writes a different file name."
        )


def espins_seedname(espins_dir):
    candidates = sorted(espins_dir.glob("*.inp1.mcin"))
    if candidates:
        return candidates[0].name[:-len(".inp1.mcin")]
    return None


def run_espins():
    espins_dir = Path("espins")
    run_mc = espins_dir / "run_mc.sh"

    if not run_mc.is_file():
        require("espins", "the ESpinS Monte Carlo (3_monte_carlo)",
                overrides=overrides)
        effective = settings.get("espins_run_command") or shell_command(
            "espins", overrides=overrides)
        if not effective:
            raise ValueError(
                "ESpinS backend selected but no ESpinS command is configured "
                "('espins_run_command' or external/codes.json)."
            )
        script = code_root() / "3_monte_carlo" / "4_write_espins_mcin.py"
        print("Generating ESpinS inputs (4_write_espins_mcin.py) ...")
        subprocess.run(
            [sys.executable, str(script)], cwd=str(Path.cwd()), check=True
        )

    if not run_mc.is_file():
        raise FileNotFoundError(
            f"{run_mc} was not created; check 'espins_run_command'."
        )

    run_command("bash run_mc.sh", espins_dir, "ESpinS")

    seedname = espins_seedname(espins_dir)
    candidates = []
    configured = settings.get("espins_output_file")
    if configured:
        candidates.append(Path(configured))
    candidates.append(Path("output"))
    if seedname:
        candidates.extend([
            Path(f"{seedname}.out.mc"),
            Path(f"{seedname}.out"),
            Path(f"{seedname}.mc"),
            Path(f"{seedname}.out.dat"),
        ])

    for candidate in candidates:
        path = candidate if candidate.is_absolute() else espins_dir / candidate
        if path.is_file() and path.stat().st_size > 0:
            copy_output(path)
            return

    converter = settings.get("espins_output_converter")
    if converter:
        run_command(converter, Path.cwd(), "ESpinS output converter")
        if Path("output").is_file() and Path("output").stat().st_size > 0:
            print("ESpinS output converted to ./output.")
            return
        raise FileNotFoundError(
            "The ESpinS output converter did not create ./output."
        )

    looked_for = ", ".join(str(item) for item in candidates)
    raise FileNotFoundError(
        f"ESpinS finished but no output file was found (looked for: "
        f"{looked_for}). Set 'espins_output_file' or "
        "'espins_output_converter'."
    )


if backend == "vampire":
    run_vampire()
elif backend == "espins":
    run_espins()
else:
    raise ValueError(
        f"Unknown mc_backend '{backend}'. Use 'vampire' or 'espins'."
    )
