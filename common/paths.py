"""
automag.common.paths
====================

Central place for resolving the two paths the pipeline needs:

* the **code root** - the directory of the ``automag`` repository, where
  ``common/``, ``geometries/`` and the stage folders live;
* the **data root** - the directory of a concrete project, where ``CalcFold/``,
  ``geometries/`` and all stage outputs are written.

Historically both roles were played by the ``AUTOMAG_PATH`` environment
variable.  This is what prevented running a pipeline from an arbitrary
directory.  The functions below keep ``AUTOMAG_PATH`` working as a fallback
for the data root, so the old per-stage ``input.py`` workflow is unchanged,
while the new unified pipeline sets

* ``AUTOMAG_PATH``     -> code root
* ``AUTOMAG_WORKDIR``  -> data root

.. codeauthor:: Automag contributors
"""

import os
import sys
from pathlib import Path


def frozen() -> bool:
    """True when running from the single-file (PyInstaller) build."""
    return bool(getattr(sys, "frozen", False))


def code_root() -> Path:
    """Return the directory that contains the ``automag`` sources.

    In a frozen build the sources are unpacked to ``sys._MEIPASS``.
    """
    bundle = getattr(sys, "_MEIPASS", None)
    if frozen() and bundle:
        return Path(bundle).resolve()
    for variable in ("AUTOMAG_CODE_PATH", "AUTOMAG_PATH"):
        value = os.environ.get(variable)
        if value:
            path = Path(value).expanduser()
            if path.is_dir():
                return path.resolve()
    return Path(__file__).resolve().parents[1]


def helper_command(relative_path, args="") -> str:
    """Return the shell command that runs one of automag's helper scripts.

    In a source checkout the script is executed with the current ``python``.
    In a frozen build the helpers live inside the executable, so the
    executable itself is asked to run the script (``--run-helper``).
    """
    relative_path = str(relative_path).replace(os.sep, "/")
    if frozen():
        command = f'"{sys.executable}" --run-helper {relative_path}'
    else:
        command = f'python "{code_root() / relative_path}"'
    if args:
        command = f"{command} {args}"
    return command


PROPHET_MODEL_NAME = "prophet-v1-mag.pt"
PROPHET_HEAD_NAME = "prophet-v1-mag-head.pt"


def models_dir() -> Path:
    """Directory holding the checkpoints bundled with automag.

    In a frozen build PyInstaller copies ``--add-data models:models`` next to
    the sources (``sys._MEIPASS/models``), so the same relative path works for
    a source checkout and for ``dist/run_automag``.
    """
    return code_root() / "models"


def model_file(name: str) -> Path:
    """Absolute path of one of the Prophet model files.

    Resolution order:

    1. the environment (``PROPHET_MODEL`` for the energy model,
       ``PROPHET_HEAD`` for the moment head) - this is what a build made with
       ``--no-models`` relies on;
    2. the copy bundled with automag (``models/`` in a checkout,
       ``<bundle>/_internal/models`` in the frozen build).
    """
    variable = {"prophet-v1-mag.pt": "PROPHET_MODEL",
                "prophet-v1-mag-head.pt": "PROPHET_HEAD"}.get(name)
    if variable:
        value = os.environ.get(variable)
        if value:
            path = Path(value).expanduser()
            if path.is_file():
                return path.resolve()
            raise FileNotFoundError(
                f"{variable} points at {path}, which does not exist.")
    path = models_dir() / name
    if not path.is_file():
        hint = f"set {variable}" if variable else "pass an explicit 'model_path'"
        raise FileNotFoundError(
            f"Model file '{name}' is missing from {models_dir()} - {hint}, or "
            "put the checkpoints of https://huggingface.co/kairosmaterials "
            "into that directory."
        )
    return path


def data_root() -> Path:
    """Return the directory where calculation results are stored.

    ``AUTOMAG_WORKDIR`` takes precedence.  If it is not set, fall back to
    ``AUTOMAG_PATH`` (the historical behaviour) and finally to the code root.
    """
    for variable in ("AUTOMAG_WORKDIR", "AUTOMAG_PATH"):
        value = os.environ.get(variable)
        if value:
            return Path(value).expanduser().resolve()
    return code_root()


def calcfold() -> Path:
    """Directory holding raw calculation inputs/outputs."""
    return data_root() / "CalcFold"


def geometries_dir() -> Path:
    """Directory holding the POSCAR files of the current project."""
    return data_root() / "geometries"
