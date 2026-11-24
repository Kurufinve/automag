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
RUNTIME_HOOKS = ("tools/pyi_rth_scipy_stats.py", "tools/pyi_rth_sources.py")

#: Heavy optional stacks.  They are in the build environment only because some
#: other package depends on them - ``prophet_spin`` pulls in ~5 GB of
#: torch + triton + CUDA libraries, which has no place in a DFT executable.
#: ``--prophet`` bundles prophet_spin and keeps them (see the flag's help).
#: Modules and files dropped by the slim build.  The model only needs CPU
#: inference, so the whole GPU toolchain (CUDA runtime, libtorch_cuda, triton,
#: the JIT compiler and llvmlite) goes away, along with development helpers.
SLIM_EXCLUDED_MODULES = (
    "nvidia", "triton", "llvmlite", "numba", "jedi", "IPython",
)

#: Paths (relative to the onedir bundle) removed after PyInstaller finishes.
#: Nothing inside torch/ itself: e3nn and torch_geometric JIT-compile parts of
#: the model at import time and need torch's helpers (bin, distributed,
#: testing, _dynamo, _inductor) and, more importantly, the Python sources.
SLIM_PRUNED_PATHS = (
    "nvidia", "triton", "llvmlite", "numba",
    "torch/include",            # C++ headers, ~64 MB
)


EXCLUDED_MODULES = (
    "nvidia", "triton", "torch._inductor", "torch._dynamo",
    "torch.distributed", "torch.testing", "torch.utils.tensorboard",
    "llvmlite", "numba", "jedi", "IPython",
)

#: Paths (relative to the onedir bundle) removed after PyInstaller finishes.
SLIM_PRUNED_PATHS = (
    "nvidia", "triton", "llvmlite", "numba",
    "torch/include",            # C++ headers, ~64 MB
    "torch/bin",                # protoc and friends
    "torch/lib/libtorch_cuda.so",
    "torch/lib/libtorch_cuda_linalg.so",
    "torch/lib/libtorch_nvshmem.so",
    "torch/distributed", "torch/testing", "torch/test",
    "torch/utils/tensorboard", "torch/_inductor", "torch/_dynamo",
)


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
    if args.slim and not args.prophet:
        # Only prune when torch is not bundled: prophet's JIT path needs
        # numba/llvmlite and the sources, so a slim torch build is fragile.
        for module in SLIM_EXCLUDED_MODULES:
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
    # PyInstaller's --strip also strips the OpenBLAS libraries that scipy and
    # numpy ship, which corrupts them ("ELF load command address/offset not
    # properly aligned").  strip_binaries() below does the same job selectively.

    # prophet_spin is only bundled on request (see EXCLUDED_MODULES).
    if args.prophet and _importable("prophet_spin"):
        command += ["--hidden-import", "prophet_spin"]
        command += ["--collect-submodules", "prophet_spin"]
        # torch ships Python sources plus a compiled library tree; PyInstaller
        # needs the data files as well, otherwise importing it fails in the
        # frozen build.
        # e3nn / torch_geometric / prophet(_spin) JIT-compile model code at
        # import time, so their .py sources must be present in the bundle.
        for package in ("e3nn", "torch_geometric", "prophet_spin", "prophet"):
            if _importable(package):
                command += ["--collect-all", package]
                # TorchScript reads the source through inspect.getsource(),
                # which needs the .py files on disk, not only the .pyc.
                command += ["--add-data",
                            f"{_package_dir(package)}:{package}"]
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


def _package_dir(package):
    """Directory of an importable package (for --add-data)."""
    module = __import__(package)
    return Path(module.__file__).parent


def copy_package_sources(output, package):
    """Copy the ``.py`` files of a package next to the frozen copy.

    ``inspect.getsource`` (used by torch's config module and by TorchScript in
    e3nn/torch_geometric) reads the file from disk; PyInstaller only ships the
    compiled modules inside the PYZ archive.  The sources are ~30 MB, the
    binaries are already collected by PyInstaller.
    """
    import shutil

    try:
        module = __import__(package)
    except ImportError:
        return
    source = Path(module.__file__).parent
    target = output / "_internal" / package
    copied = 0
    for path in source.rglob("*.py"):
        relative = path.relative_to(source)
        destination = target / relative
        if destination.exists():
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)
        copied += 1
    if copied:
        print(f"Copied {copied} source file(s) of {package} for inspect.getsource()")


def prune_torch_bin(output):
    """Keep only ``torch_shm_manager`` in ``torch/bin``.

    The wheel ships protoc and a dozen test executables (``test_api`` alone is
    15 MB); torch only needs the shared-memory manager at runtime.
    """
    binary = output / "_internal" / "torch" / "bin"
    if not binary.is_dir():
        return
    removed = 0
    for item in binary.iterdir():
        if item.is_file() and item.name != "torch_shm_manager":
            removed += item.stat().st_size
            item.unlink()
    if removed:
        print(f"Removed {removed // 1024 // 1024} MiB of torch/bin tools")


def strip_binaries(output):
    """Run ``strip --strip-unneeded`` over the bundled shared libraries.

    ``libtorch_cpu.so`` alone loses ~90 MB of symbols this way; without it the
    frozen torch cannot get near the size of the ``nuitka``-built USPEX core.
    """
    import shutil
    import subprocess

    strip = shutil.which("strip")
    if strip is None:
        print("strip not found, skipping symbol removal")
        return
    total = 0
    # Only torch's own libraries: stripping the NumPy/SciPy stack (OpenBLAS,
    # libgfortran) breaks the ELF segment alignment those wheels require.
    library_dir = output / "_internal" / "torch" / "lib"
    for path in sorted(library_dir.glob("*.so*")):
        if not path.is_file() or path.stat().st_size < 1_000_000:
            continue
        before = path.stat().st_size
        result = subprocess.run([strip, "--strip-unneeded", str(path)],
                                capture_output=True)
        if result.returncode == 0:
            total += before - path.stat().st_size
    if total:
        print(f"Stripped {total // 1024 // 1024} MiB from shared libraries")


def ensure_torch_bin(output):
    """PyInstaller does not collect ``torch/bin`` (e.g. ``torch_shm_manager``).

    torch >= 2.7 aborts at import time when that executable is missing, so copy
    it from the build environment's torch into the bundle.
    """
    import shutil

    try:
        import torch
    except ImportError:
        return
    source = Path(torch.__file__).parent / "bin"
    target = output / "_internal" / "torch" / "bin"
    if not source.is_dir():
        return
    target.mkdir(parents=True, exist_ok=True)
    for item in source.iterdir():
        if item.is_file():
            shutil.copy2(item, target / item.name)


def prune_bundle(output):
    """Delete the GPU stack and development files from a onedir bundle."""
    import shutil

    freed = 0
    for relative in SLIM_PRUNED_PATHS:
        target = output / "_internal" / relative
        if not target.exists():
            continue
        size = (sum(path.stat().st_size for path in target.rglob("*")
                    if path.is_file()) if target.is_dir()
                else target.stat().st_size)
        freed += size
        if target.is_dir():
            shutil.rmtree(target, ignore_errors=True)
        else:
            target.unlink(missing_ok=True)
    print(f"Slim build: removed {freed / 2**20:.0f} MiB "
          "(GPU stack and development files)")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--name", default="run_automag")
    parser.add_argument("--onedir", dest="onefile", action="store_false",
                        help="build a directory bundle instead of a single file")
    parser.add_argument("--strip", dest="strip", action="store_true",
                        default=True,
                        help="strip the binary (Linux/macOS only, --no-strip "
                             "disables it)")
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
    parser.add_argument("--no-strip", dest="strip", action="store_false",
                        help="do not run strip over the bundled shared libraries")
    parser.add_argument("--slim", action=argparse.BooleanOptionalAction,
                        default=True,
                        help=("drop the GPU stack and development files "
                              "(default): nvidia/, triton/, libtorch_cuda, "
                              "llvmlite, torch/include, ...  "
                              "--no-slim keeps them, needed for device=cuda"))
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
    if args.slim and not args.onefile:
        # The pruned paths are safe for a CPU-only torch as well (there is no
        # nvidia/ or triton/ in such a wheel); torch/bin is restored below.
        prune_bundle(output)
    if args.prophet and not args.onefile:
        ensure_torch_bin(output)
        prune_torch_bin(output)
        copy_package_sources(output, "torch")
    if args.strip and not args.onefile:
        strip_binaries(output)
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
