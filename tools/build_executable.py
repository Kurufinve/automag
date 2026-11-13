#!/usr/bin/env python
"""
automag.tools.build_executable
==============================

Build a single-file executable that runs the whole pipeline:

    dist/run_automag                                (Linux / macOS)
    dist\\run_automag.exe                            (Windows)

The executable bundles the automag sources (stage scripts, ``common/``,
``calculators/``, ``pipeline/`` and the reference geometries) together with
Python and its dependencies.  A copy of an ``examples/`` directory can then be
run on a cluster without a Python installation or an automag checkout:

    AUTOMAG_EXE=/shared/bin/run_automag ./run_pipeline.sh

How the bundled tree is addressed at run time
---------------------------------------------
``common.paths.code_root()`` returns ``sys._MEIPASS`` inside the frozen build,
and ``run_automag.py`` supports two internal switches:

``--run-script <path>``
    run a stage script (the runner uses it instead of ``python <script>``);
``--run-helper <path>``
    run one of the ``common/`` helper scripts (the job scripts and the
    prophet_spin interface use it).

Paths passed to those switches are relative to the bundle root, because a
one-file build unpacks to a new temporary directory on every start.

Usage
-----
    source /path/to/miniconda3/bin/activate automag
    python tools/build_executable.py            # -> dist/run_automag
    python tools/build_executable.py --clean    # remove build/ and dist/

``--onefile`` is the default; pass ``--onedir`` for a directory bundle that
starts faster (useful while debugging).

What goes into the bundle: the stage scripts (``*.py`` only - the example data
inside some stage directories is skipped), ``common/``, ``calculators/``,
``pipeline/``, ``geometries/`` and the data files of pymatgen/ase.  Heavy
optional stacks (torch/triton/nvidia/tensorflow/...) are excluded, which keeps
the result around 500 MB instead of several GB; the consequence is that the
``prophet_spin`` calculator is not available in the frozen build unless you
pass ``--prophet`` (a multi-GB bundle).  Everything else - VASP, Quantum
ESPRESSO, FHI-aims, generic ASE calculators - works.

.. codeauthor:: Automag contributors
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: Modules that are imported dynamically and that PyInstaller cannot see.
#: The stage scripts are bundled as *data* (``--add-data``), so their imports
#: are invisible to the analysis - every module they need has to be listed
#: here, exactly as the scripts import it.
HIDDEN_IMPORTS = (
    # pymatgen (imported by the stage scripts and by common/results.py)
    "pymatgen",
    "pymatgen.core",
    "pymatgen.core.periodic_table",
    "pymatgen.core.structure",
    "pymatgen.core.composition",
    "pymatgen.analysis",
    "pymatgen.analysis.magnetism",
    "pymatgen.analysis.magnetism.analyzer",
    "pymatgen.command_line",
    "pymatgen.command_line.enumlib_caller",
    "pymatgen.electronic_structure",
    "pymatgen.electronic_structure.core",
    "pymatgen.io",
    "pymatgen.io.ase",
    "pymatgen.io.vasp",
    "pymatgen.io.vasp.inputs",
    "pymatgen.io.vasp.outputs",
    "pymatgen.io.vasp.sets",
    "pymatgen.symmetry",
    "pymatgen.symmetry.analyzer",
    "pymatgen.transformations",
    "pymatgen.transformations.standard_transformations",
    "pymatgen.util",
    "pymatgen.util.coord",
    # ase
    "ase",
    "ase.atoms",
    "ase.io",
    "ase.io.vasp",
    "ase.calculators",
    "ase.calculators.calculator",
    "ase.calculators.vasp",
    "ase.calculators.espresso",
    "ase.calculators.aims",
    # everything else
    "spglib",
    "matplotlib",
    "matplotlib.pyplot",
    "numpy",
)

#: Packages whose *data files* have to be copied into the bundle.  PyInstaller
#: ships no hook for pymatgen or ase, so without this the frozen build dies on
#: ``FileNotFoundError: .../pymatgen/core/periodic_table.json``.
COLLECT_DATA = ("pymatgen", "ase")

#: Packages whose *modules* must all be collected: ase.io.json needs ase.io.db,
#: which is imported lazily and is therefore not found automatically.
COLLECT_SUBMODULES_EXTRA = ("ase",)

#: Packages whose subpackages are pulled in wholesale.  Keep this empty:
#: ``--collect-submodules pymatgen`` drags in every optional dependency the
#: environment happens to have (torch, tensorflow, ...) and turns a ~80 MB
#: bundle into several GB.  HIDDEN_IMPORTS above is the precise list.

#: Stage directories: only their scripts are bundled, never the example data
#: some of them contain (``3_monte_carlo/Fe12O18/espins`` alone is 600 MB of
#: tracked ESpinS output).  Keep ``1_submission``/``vampire_input`` subfolders,
#: which the scripts read at run time.
STAGE_DIRECTORIES = (
    "0_conv_tests",
    "1_lin_response",
    "2_coll",
    "2_relax",
    "3_monte_carlo",
    "4_mae",
)

#: Subdirectories of a stage directory that *are* part of the bundle.
STAGE_SUBDIRECTORIES = ("vampire_input",)

#: Directories copied in full (code, not data).
CODE_DIRECTORIES = ("common", "calculators", "pipeline", "geometries")

#: Source files that PyInstaller's frozen bytecode breaks, shipped as plain
#: source and loaded by the matching runtime hook (see
#: ``tools/pyi_rth_scipy_stats.py`` for the explanation).  The value is a list
#: of ``(search, replace)`` text pairs applied to the installed file; if a file
#: is missing, or the text to replace is not found (a different scipy version),
#: the patch is skipped and the hook stays inactive.
PATCHED_SOURCES = {
    "scipy_patches/_distn_infrastructure.py": (
        "scipy/stats/_distn_infrastructure.py",
        (("""for obj in [s for s in dir() if s.startswith('_doc_')]:
    exec('del ' + obj)
del obj""",
          """for _name in [s for s in list(globals()) if s.startswith('_doc_')]:
    globals().pop(_name, None)"""),),
    ),
}

#: Runtime hooks: code that runs inside the frozen application before the
#: program starts.  ``pyi_rth_scipy_stats.py`` works around pyinstaller#7992
#: (``scipy.stats`` cannot be imported in a frozen build otherwise).
RUNTIME_HOOKS = ("tools/pyi_rth_scipy_stats.py",)

#: Heavy optional stacks.  They are in the build environment only because some
#: other package depends on them - ``prophet_spin`` pulls in ~5 GB of
#: torch + triton + CUDA libraries, which has no place in a DFT executable.
#: ``--prophet`` bundles prophet_spin and keeps them (see the flag's help).
EXCLUDED_MODULES = (
    "torch", "torchvision", "torchaudio", "triton", "nvidia",
    "tensorflow", "tensorboard", "keras", "jax", "jaxlib", "cupy",
    "onnxruntime", "transformers", "datasets", "IPython", "jedi",
    "notebook", "pytest",
)


def data_entries(bundle_models=True):
    """``(source, destination directory)`` pairs for ``--add-data``."""
    entries = []
    # Prophet checkpoints: without them the frozen build cannot evaluate the
    # model (common.paths.models_dir() looks for <bundle>/models).
    # External-code registry (enumlib / VAMPIRE / ESpinS): common/external.py
    # reads <bundle>/external/codes.json, so bundle it when it exists.
    registry = ROOT / "external" / "codes.json"
    if registry.is_file():
        entries.append((registry, "external"))
    models = ROOT / "models"
    if not bundle_models:
        print("Models are NOT bundled: set PROPHET_MODEL and PROPHET_HEAD at "
              "run time (or 'model_path' in the stage params).")
    elif models.is_dir() and any(models.glob("*.pt")):
        entries.append((models, "models"))
    else:
        print("WARNING: no models/*.pt found - the frozen build will not be "
              "able to run prophet_spin without an explicit model_path.")
    for name in CODE_DIRECTORIES:
        source = ROOT / name
        if source.is_dir():
            entries.append((source, name))
    for name in STAGE_DIRECTORIES:
        source = ROOT / name
        if not source.is_dir():
            continue
        for path in sorted(source.rglob("*")):
            if path.is_dir() or "__pycache__" in path.parts:
                continue
            relative = path.relative_to(ROOT)
            keep = (path.suffix == ".py"
                    or any(part in STAGE_SUBDIRECTORIES
                           for part in relative.parts))
            if not keep:
                continue
            # --add-data copies a file into a *directory*
            entries.append((path, str(relative.parent)))
    return entries


def write_patched_sources(destination):
    """Write the sources of :data:`PATCHED_SOURCES` into ``destination``.

    Returns the ``(source, destination directory)`` pairs for ``--add-data``.
    """
    import importlib.util

    entries = []
    for target, (module_path, replacements) in PATCHED_SOURCES.items():
        spec = importlib.util.find_spec(module_path.split("/")[0])
        if spec is None or not spec.submodule_search_locations:
            continue
        package_root = Path(next(iter(spec.submodule_search_locations)))
        source_path = package_root / Path(*module_path.split("/")[1:])
        if not source_path.is_file():
            continue
        text = source_path.read_text()
        for old_text, new_text in replacements:
            if old_text not in text:
                print(f"WARNING: {module_path}: nothing to patch "
                      f"(scipy changed?) - the frozen build may not import "
                      f"scipy.stats on some Python versions")
                text = None
                break
            text = text.replace(old_text, new_text, 1)
        if text is None:
            continue
        output = Path(destination) / target
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(text)
        print(f"patched {module_path} -> {output}")
        entries.append((output, str(Path(target).parent)))
    return entries


def build_command(args):
    command = [
        sys.executable, "-m", "PyInstaller",
        "--noconfirm",
        "--clean",
        "--onefile" if args.onefile else "--onedir",
        "--name", args.name,
        "--distpath", str(ROOT / "dist"),
        "--workpath", str(ROOT / "build" / "pyinstaller"),
        "--specpath", str(ROOT / "build"),
        "--paths", str(ROOT),
    ]
    for module in HIDDEN_IMPORTS:
        command += ["--hidden-import", module]
    for package in COLLECT_DATA:
        command += ["--collect-data", package]
    for package in COLLECT_SUBMODULES_EXTRA:
        command += ["--collect-submodules", package]
    if not args.prophet:
        for module in EXCLUDED_MODULES + ("prophet_spin",):
            command += ["--exclude-module", module]
    for source, destination in write_patched_sources(ROOT / "build" / "patches"):
        command += ["--add-data", f"{source}{os.pathsep}{destination}"]
    for hook in RUNTIME_HOOKS:
        path = ROOT / hook
        if path.is_file():
            command += ["--runtime-hook", str(path)]
    for source, destination in data_entries(bundle_models=args.models):
        # os.pathsep is ':' on POSIX and ';' on Windows - PyInstaller expects
        # the platform separator here.
        command += ["--add-data", f"{source}{os.pathsep}{destination}"]

    if args.icon:
        command += ["--icon", args.icon]
    if args.strip:
        command += ["--strip"]
    # prophet_spin is only bundled on request (see EXCLUDED_MODULES).
    if args.prophet and _importable("prophet_spin"):
        command += ["--hidden-import", "prophet_spin"]
        command += ["--collect-submodules", "prophet_spin"]
        # torch ships Python sources plus a compiled library tree; PyInstaller
        # needs the data files as well, otherwise importing it fails in the
        # frozen build.
        if _importable("torch"):
            # ``--collect-all`` (sources + data + binaries) is the recipe that
            # works for torch: partial collection breaks its lazy imports
            # ("cannot import name 'detect_fake_mode'", "NameError: name
            # 'name'" in torch/_numpy).
            command += ["--collect-all", "torch"]
            # torch._numpy builds ufuncs with exec() of generated code and
            # fails in a frozen build; automag and prophet_spin never use it.
            command += ["--exclude-module", "torch._numpy"]

    command.append(str(ROOT / "run_automag.py"))
    return command


def _importable(module):
    try:
        __import__(module)
    except Exception:  # noqa: BLE001 - optional dependency
        return False
    return True


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--name", default="run_automag")
    parser.add_argument("--onedir", dest="onefile", action="store_false",
                        help="build a directory bundle instead of a single file")
    parser.add_argument("--strip", action="store_true",
                        help="strip the binary (Linux/macOS only)")
    parser.add_argument("--icon", default=None)
    parser.add_argument("--clean", action="store_true",
                        help="remove build/ and dist/ and exit")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the PyInstaller command and exit")
    parser.add_argument("--models", action=argparse.BooleanOptionalAction,
                        default=True,
                        help=("bundle the Prophet checkpoints from models/ "
                              "(default).  With --no-models the bundle expects "
                              "PROPHET_MODEL / PROPHET_HEAD (or an explicit "
                              "'model_path' in the stage params) at run time"))
    parser.add_argument("--prophet", action=argparse.BooleanOptionalAction,
                        default=True,
                        help=("bundle the prophet_spin calculator with its "
                              "torch/CUDA stack (several GB); without it the "
                              "frozen build cannot run prophet_spin"))
    args = parser.parse_args(argv)

    if args.clean:
        for directory in (ROOT / "build", ROOT / "dist"):
            shutil.rmtree(directory, ignore_errors=True)
        print("Removed build/ and dist/.")
        return 0

    command = build_command(args)
    print(" ".join(command))
    if args.dry_run:
        return 0

    result = subprocess.run(command, cwd=str(ROOT))
    if result.returncode != 0:
        return result.returncode

    suffix = ".exe" if os.name == "nt" else ""
    output = ROOT / "dist" / (args.name + suffix)
    # --onedir puts the executable inside a directory of the same name, so the
    # launcher hint must point at <dist>/<name>/<name>, not at the directory.
    executable = output if args.onefile else output / (args.name + suffix)
    print(f"\nBuilt {output}")
    print(f"Executable: {executable}")
    print("Copy it wherever you like and run the examples with:")
    print(f"    AUTOMAG_EXE={executable} ./run_pipeline_compiled.sh")
    print(f"    AUTOMAG_EXE={executable} ./run_pipeline_slurm_compiled.sh")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
