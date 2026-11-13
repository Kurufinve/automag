#!/usr/bin/env python
"""
automag.tools.build_external
============================

Download and compile the optional external codes **enumlib**, **ESpinS** and
**VAMPIRE**, then record where they are in ``<prefix>/codes.json``.

Automag reads that file at run time (:mod:`common.external`), so the stage
scripts do not need any hard-coded paths.

Two modes
---------
*Build* (default): clone each missing code into ``<prefix>/src/<name>``,
compile it, and record the result.  A code that fails to compile is recorded
with ``"available": false`` and reported as a warning - automag keeps working,
the corresponding step just tells the user to compile it by hand.

*Register* (``--<code>-path``): use an installation that already exists.
Nothing is downloaded or compiled, the path (and optional pre/post commands)
are written to ``codes.json`` instead.

Examples
--------
    # build everything under external/
    python tools/build_external.py

    # use pre-installed codes, with the commands needed on this cluster
    python tools/build_external.py \\
        --vampire-path /opt/vampire/vampire-serial \\
        --vampire-pre-command "module load intel" \\
        --vampire-post-command "module purge" \\
        --espins-path /opt/espins/mc.x \\
        --enumlib-path /opt/enumlib \\
        --skip vampire

    # only refresh the registry from a file
    python tools/build_external.py --config cluster_codes.json

The same options are available through ``./build.sh --external ...``.

``--config`` accepts a JSON/YAML file with the same keys as the entry format
of ``codes.json`` plus ``prefix``::

    {"vampire": {"path": "/opt/vampire/vampire-serial",
                 "pre_command": "module load intel",
                 "post_command": "module purge"}}

.. codeauthor:: Automag contributors
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: Build recipes.  ``steps`` and ``step_dirs`` are parallel lists; ``probe``
#: maps an entry key to the path (relative to the clone) that must exist.
RECIPES = {
    "enumlib": {
        "repo": "https://github.com/msg-byu/enumlib",
        "recursive": True,
        "step_dirs": ["symlib/src", "src"],
        "steps": [
            ["make", "F90=gfortran"],
            ["make", "enum.x", "makestr.x", "F90=gfortran"],
        ],
        "probe": {"enum": "src/enum.x", "makestr": "aux_src/makeStr.py"},
        "extra": {"python": sys.executable},
        "requires": ["gfortran", "make", "git"],
    },
    "vampire": {
        "repo": "https://github.com/richard-evans/vampire",
        "recursive": False,
        "step_dirs": ["."],
        "steps": [["make", "serial"]],
        "probe": {"command": "vampire-serial"},
        "requires": ["g++", "make", "git"],
    },
    "espins": {
        "repo": "https://github.com/nafiserb/ESpinS",
        "recursive": False,
        "step_dirs": ["."],
        "steps": [["make"]],
        "probe": {"command": "mc.x"},
        # ESpinS ships an Intel/MKL make.sys; it is not buildable with gfortran.
        "requires": ["make", "git"],
        "hint": ("ESpinS expects Intel Fortran + MKL (edit make.sys). "
                 "Compile it manually and pass --espins-path."),
    },
}


def warn(message):
    print(f"WARNING: {message}", file=sys.stderr)


def have(program):
    return shutil.which(program) is not None


def run_step(command, cwd, dry_run=False, env=None):
    print(f"  $ {' '.join(command)}   (in {cwd})")
    if dry_run:
        return True
    try:
        result = subprocess.run(command, cwd=str(cwd), env=env)
    except OSError as error:
        warn(str(error))
        return False
    return result.returncode == 0


def system_env():
    """Environment without the conda toolchain (it often breaks linking)."""
    env = {key: value for key, value in os.environ.items()
           if key not in ("CC", "CXX", "FC", "F77", "F90", "CFLAGS",
                          "CXXFLAGS", "FFLAGS", "LDFLAGS", "CPPFLAGS",
                          "CONDA_PREFIX", "CONDA_DEFAULT_ENV")}
    # System compilers first: a conda g++ may not match the system glibc.
    env["PATH"] = "/usr/bin:/bin:/usr/local/bin:" + env.get("PATH", "")
    return env


def clone(name, recipe, destination, dry_run=False, env=None):
    if destination.is_dir():
        print(f"{name}: using the existing sources in {destination}")
        return True
    destination.parent.mkdir(parents=True, exist_ok=True)
    command = ["git", "clone", "--depth", "1"]
    if recipe.get("recursive"):
        command += ["--recursive", "--shallow-submodules"]
    command += [recipe["repo"], str(destination)]
    print(f"{name}: cloning {recipe['repo']}")
    return run_step(command, ROOT, dry_run, env=env)


def compile_code(name, recipe, source, dry_run=False):
    """Run the build steps; return True when every probe exists.

    A step that fails with the current environment is retried once with the
    system toolchain: the conda compilers frequently fail to link against the
    system glibc (undefined references to ``memcpy@GLIBC_2.14`` and friends).
    """
    for directory, step in zip(recipe["step_dirs"], recipe["steps"]):
        target = source / directory
        if run_step(step, target, dry_run):
            continue
        warn(f"{name}: build step failed in {target}; retrying with the "
             "system toolchain")
        if run_step(step, target, dry_run, env=system_env()):
            continue
        warn(f"{name}: build step failed in {target}")
        return False
    return True


def verify(name, recipe, source):
    """Return the entry recorded for ``source``, or None if a probe is missing."""
    entry = {"available": True}
    for key, relative in recipe["probe"].items():
        path = (source / relative).resolve()
        if not path.exists():
            warn(f"{name}: expected {path} after the build")
            return None
        entry[key] = str(path)
    entry.update(recipe.get("extra") or {})
    return entry


def build_one(name, recipe, prefix, dry_run=False):
    missing = [program for program in recipe["requires"] if not have(program)]
    if missing:
        warn(f"{name}: missing build tool(s): {', '.join(missing)}")
        return _unavailable(name, recipe)
    source = prefix / "src" / name
    if not clone(name, recipe, source, dry_run):
        return _unavailable(name, recipe)
    if not compile_code(name, recipe, source, dry_run):
        return _unavailable(name, recipe)
    if dry_run:
        return {"available": True}
    entry = verify(name, recipe, source)
    if entry is None:
        return _unavailable(name, recipe)
    print(f"{name}: OK -> {entry}")
    return entry


def _unavailable(name, recipe):
    hint = recipe.get("hint")
    message = (f"{name} could not be built, so the corresponding automag step "
               "will not work. Compile it manually and register the result "
               f"with '--{name}-path' (or edit <prefix>/codes.json).")
    if hint:
        message += f" Hint: {hint}"
    warn(message)
    return {"available": False, "error": "build failed"}


def register_from_path(name, recipe, value, args):
    """Record an existing installation instead of building one."""
    path = Path(value).expanduser()
    if not path.exists():
        warn(f"{name}: '{path}' does not exist; recording it anyway")
    entry = {"available": True}

    if name == "enumlib":
        directory = path if path.is_dir() else path.parent
        candidates = {
            "enum": [directory / "enum.x", directory / "src" / "enum.x",
                     directory / "bin" / "enum.x"],
            "makestr": [directory / "makeStr.py",
                        directory / "aux_src" / "makeStr.py"],
        }
        for key, options in candidates.items():
            found = next((item for item in options if item.exists()), None)
            if found is None:
                # A bare directory may also be the bin dir itself.
                found = options[-1]
            entry[key] = str(found.resolve())
        entry["python"] = args.enumlib_python or sys.executable
    else:
        entry["command"] = str(path.resolve()) if path.exists() else str(path)

    entry["pre_command"] = getattr(args, f"{name}_pre_command") or ""
    entry["post_command"] = getattr(args, f"{name}_post_command") or ""
    entry["source"] = "user-provided"
    print(f"{name}: registered {entry}")
    return entry


def load_config(path):
    from common.input_loader import load_input  # noqa: PLC0415 - optional path

    data = load_input(path)
    if "PIPELINE" in data or "pipeline" in data:
        data = data.get("PIPELINE") or data.get("pipeline")
    if not isinstance(data, dict):
        raise SystemExit(f"{path} must contain a mapping of code settings.")
    return data


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Build or register enumlib / ESpinS / VAMPIRE for automag.")
    parser.add_argument("--prefix", default=str(ROOT / "external"),
                        help="install prefix (default: <automag>/external)")
    parser.add_argument("--config", default=None,
                        help="JSON/YAML file with per-code settings")
    parser.add_argument("--skip", action="append", default=[],
                        choices=sorted(RECIPES),
                        help="leave this code untouched (repeatable)")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the commands without running them")
    parser.add_argument("--offline", action="store_true",
                        help="never clone or compile, only rewrite codes.json")
    for name in RECIPES:
        parser.add_argument(f"--{name}-path", default=None,
                            help=f"existing {name} installation (skips the build)")
        parser.add_argument(f"--{name}-pre-command", default=None,
                            help=f"command to run before {name} (module load, ...)")
        parser.add_argument(f"--{name}-post-command", default=None,
                            help=f"command to run after {name}")
    parser.add_argument("--enumlib-python", default=None,
                        help="interpreter used for makeStr.py")
    args = parser.parse_args(argv)

    config = load_config(args.config) if args.config else {}
    prefix = Path(config.get("prefix") or args.prefix).expanduser()
    if not prefix.is_absolute():
        prefix = (ROOT / prefix).resolve()
    print(f"External codes prefix: {prefix}")

    entries = {}
    for name, recipe in RECIPES.items():
        if name in args.skip:
            print(f"{name}: skipped")
            continue
        settings = config.get(name) or {}
        explicit = getattr(args, f"{name}_path") or settings.get("path")
        if explicit:
            if settings.get("pre_command") and not getattr(
                    args, f"{name}_pre_command"):
                setattr(args, f"{name}_pre_command", settings["pre_command"])
            if settings.get("post_command") and not getattr(
                    args, f"{name}_post_command"):
                setattr(args, f"{name}_post_command", settings["post_command"])
            entries[name] = register_from_path(name, recipe, explicit, args)
            continue
        if args.offline:
            print(f"{name}: offline, not built")
            continue
        entries[name] = build_one(name, recipe, prefix, args.dry_run)

    if args.dry_run:
        print("(dry run: nothing written)")
        return 0

    target = prefix / "codes.json"
    existing = {}
    if target.is_file():
        try:
            existing = json.loads(target.read_text())
        except ValueError:
            existing = {}
    existing.update(entries)
    existing["generated_by"] = "tools/build_external.py"
    existing["prefix"] = str(prefix)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(existing, indent=2, sort_keys=True) + "\n")
    print(f"\nWrote {target}")

    broken = [name for name, entry in entries.items()
              if not entry.get("available", True)]
    if broken:
        print("\nThe following codes are NOT usable:", ", ".join(broken))
        print("Compile them manually and register the paths, e.g.:")
        for name in broken:
            print(f"    ./build.sh --external --{name}-path /path/to/{name}")
        print("Until then the corresponding pipeline steps will report an "
              "error / be skipped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
