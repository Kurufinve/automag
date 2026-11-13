"""
automag.common.external
=======================

Registry of the optional external codes: **enumlib**, **VAMPIRE** and
**ESpinS**.

Automag only needs their executables at run time, so the paths and the shell
commands around them are kept in one small JSON file instead of being
hard-coded in the stage scripts.

Where the settings come from (later wins)
-----------------------------------------
1. built-in defaults (``enum.x`` / ``makeStr.py``, ``vampire``, ``mc.x`` on
   ``PATH``), so an existing setup keeps working without any file;
2. ``<code root>/external/codes.json``, written by
   ``tools/build_external.py`` (``./build.sh --external``);
3. a file pointed at by ``AUTOMAG_EXTERNAL_CODES`` - either the JSON file or
   the directory that contains it;
4. the optional ``external_codes`` mapping in the pipeline input file, e.g.::

       PIPELINE = {"external_codes": {"vampire": {"command": "/opt/vampire"}}}

Entry format::

    "vampire": {
        "available": true,
        "command": "/opt/vampire/vampire-serial",
        "pre_command": "module load intel",   # optional, run first
        "post_command": "module purge"        # optional, run last
    }

``enumlib`` has two executables instead of ``command``: ``enum`` and
``makestr`` (``python`` runs the latter).

``shell_command("vampire", args="input")`` returns one shell line with the
pre/post commands around the executable; :func:`run` executes it with
``/bin/bash`` when a shell is needed (``source``, ``module``).

.. codeauthor:: Automag contributors
"""

import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

from common.paths import code_root

#: Names of the external codes known to automag.
CODES = ("enumlib", "vampire", "espins")

#: Name of the generated settings file.
FILE_NAME = "codes.json"

#: Settings used when no ``codes.json`` exists (historical behaviour: the
#: executables are expected on ``PATH``).
DEFAULT_ENTRIES = {
    "enumlib": {
        "available": True,
        "enum": "enum.x",
        "makestr": "makeStr.py",
        "python": "python",
        "pre_command": "",
        "post_command": "",
    },
    "vampire": {
        "available": True,
        "command": "vampire",
        "pre_command": "",
        "post_command": "",
    },
    "espins": {
        "available": True,
        "command": "mc.x",
        "pre_command": "",
        "post_command": "",
    },
}


def codes_path():
    """Return the ``codes.json`` path that :func:`registry` reads.

    Order: ``AUTOMAG_EXTERNAL_CODES``; ``<code root>/external/codes.json``;
    finally ``<dir of the executable>/external/codes.json``, which is where a
    frozen single-file build looks for it (``_MEIPASS`` is temporary).
    """
    value = os.environ.get("AUTOMAG_EXTERNAL_CODES")
    if value:
        path = Path(value).expanduser()
        return path / FILE_NAME if path.is_dir() else path

    candidates = [code_root() / "external" / FILE_NAME]
    if getattr(sys, "frozen", False):
        candidates.append(Path(sys.executable).resolve().parent / "external"
                          / FILE_NAME)
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return candidates[0]


def _file_entries():
    path = codes_path()
    if not path.is_file():
        return {}
    try:
        data = json.loads(path.read_text())
    except (OSError, ValueError) as error:  # noqa: BLE001 - keep running
        print(f"WARNING: cannot read {path}: {error}")
        return {}
    return {key: value for key, value in data.items()
            if isinstance(value, dict)}


def registry(overrides=None):
    """Return the merged settings of every known code."""
    stored = _file_entries()
    overrides = overrides or {}
    merged = {}
    for name in CODES:
        entry = dict(DEFAULT_ENTRIES[name])
        entry.update(stored.get(name) or {})
        if isinstance(overrides.get(name), dict):
            entry.update(overrides[name])
        merged[name] = entry
    return merged


def tool(name, overrides=None):
    """Return the merged settings of one code."""
    if name not in CODES:
        raise KeyError(f"Unknown external code '{name}'. Known: {', '.join(CODES)}.")
    return registry(overrides)[name]


def available(name, overrides=None):
    """False when the build step recorded the code as unusable."""
    return bool(tool(name, overrides).get("available", True))


def command(name, key="command", overrides=None):
    """Return the executable of ``name`` (``key`` selects enumlib's program)."""
    value = tool(name, overrides).get(key)
    if not value:
        raise ValueError(
            f"No '{key}' configured for the external code '{name}'. "
            f"Set it in {codes_path()} or in 'external_codes'.")
    return str(value)


def shell_command(name, key="command", args="", overrides=None, override=None):
    """Return a single shell line for ``name`` with its pre/post commands."""
    entry = tool(name, overrides)
    parts = []
    pre = str(entry.get("pre_command") or "").strip()
    exe = str(override if override is not None else (entry.get(key) or "")).strip()
    post = str(entry.get("post_command") or "").strip()
    if pre:
        parts.append(pre)
    if exe:
        parts.append(f"{exe} {args}".strip())
    if post:
        parts.append(post)
    return "; ".join(parts)


def require(name, purpose="", overrides=None):
    """Return the entry of ``name`` or raise a helpful error."""
    entry = tool(name, overrides)
    if entry.get("available", True):
        return entry
    raise RuntimeError(
        f"'{name}' is not available, but it is needed for "
        f"{purpose or 'this step'}.\n"
        "Build it with './build.sh --external', or point automag at an "
        f"existing installation through {codes_path()}, the "
        "AUTOMAG_EXTERNAL_CODES environment variable, or "
        "'\"external_codes\"' in the pipeline input file.")


def _bash():
    return "/bin/bash" if os.path.exists("/bin/bash") else None


def _executable(entry, key, override):
    """Return the argv prefix for an entry, honouring ``python`` for .py scripts."""
    exe = str(override if override is not None else (entry.get(key) or "")).strip()
    if not exe:
        raise ValueError(f"No '{key}' configured for this external code.")
    interpreter = str(entry.get("python") or "").strip()
    if interpreter and exe.lower().endswith(".py"):
        return shlex.split(interpreter) + [exe]
    return shlex.split(exe)


def run(name, key="command", args="", cwd=None, timeout=None, check=True,
        overrides=None, override=None):
    """Run an external code in ``cwd`` (with its pre/post commands).

    Without pre/post commands the executable is started directly, so a
    ``timeout`` kills the process itself rather than a wrapping shell.  A
    ``python`` entry is used to run ``.py`` programs such as enumlib's
    ``makeStr.py`` (whose shebang may point at a missing interpreter).
    """
    entry = require(name, overrides=overrides)
    pre = str(entry.get("pre_command") or "").strip()
    post = str(entry.get("post_command") or "").strip()
    argv = _executable(entry, key, override) + shlex.split(args)
    exe = " ".join(argv[:len(argv) - len(shlex.split(args))])

    if not pre and not post:
        print(f"Running {name}: {' '.join(argv)}")
        return subprocess.run(argv, cwd=cwd, timeout=timeout, check=check)

    line = "; ".join(filter(None, [pre, f"{exe} {args}".strip(), post]))
    print(f"Running {name}: {line}")
    return subprocess.run(line, shell=True, cwd=cwd, timeout=timeout,
                          check=check, executable=_bash())


def summary(overrides=None):
    """Return a short human-readable report (used by the pipeline banner)."""
    lines = []
    for name, entry in registry(overrides).items():
        if name == "enumlib":
            target = entry.get("enum") or "?"
        else:
            target = entry.get("command") or "?"
        state = "ok" if entry.get("available", True) else "UNAVAILABLE"
        extra = ""
        if entry.get("pre_command"):
            extra = f" [pre: {entry['pre_command']}]"
        lines.append(f"{name:8s} {state:12s} {target}{extra}")
    return "\n".join(lines)
