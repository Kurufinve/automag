"""PyInstaller runtime hook: make ``scipy.stats`` importable in a frozen build.

``scipy/stats/_distn_infrastructure.py`` ends its docstring bookkeeping with::

    for obj in [s for s in dir() if s.startswith('_doc_')]:
        exec('del ' + obj)
    del obj

Inside a PyInstaller bundle module-level ``dir()`` returns an empty list
(pyinstaller#7992, #8186), so the loop never runs and ``del obj`` raises
``NameError`` - which aborts the import of ``scipy.stats``, and with it
``scipy.signal``, ``pymatgen`` and every stage script.

``tools/build_executable.py`` therefore ships a patched copy of that one file
(the loop rewritten with ``globals().pop(..., None)``); this hook installs an
import finder that loads it from the bundle instead of PyInstaller's bytecode.
A normal (non-frozen) interpreter is unaffected - the finder is only active
when the patched source is present.
"""

import os
import sys

TARGET = "scipy.stats._distn_infrastructure"
PATCH_DIR = "scipy_patches"
PATCH_FILE = "_distn_infrastructure.py"


def _patch_path():
    """Location of the patched source inside the bundle (or ``None``)."""
    roots = [getattr(sys, "_MEIPASS", None),
             os.path.dirname(os.path.abspath(sys.executable))]
    for root in roots:
        if not root:
            continue
        candidate = os.path.join(root, PATCH_DIR, PATCH_FILE)
        if os.path.isfile(candidate):
            return candidate
    return None


class _PatchedSourceFinder:
    """Load the patched scipy module from source instead of the frozen PYZ."""

    def __init__(self, path):
        self._path = path

    def find_spec(self, fullname, path=None, target=None):
        if fullname != TARGET:
            return None
        import importlib.util

        try:
            sys.meta_path.remove(self)
        except ValueError:
            pass
        spec = importlib.util.spec_from_file_location(fullname, self._path)
        if spec is None or spec.loader is None:
            return None
        return spec


def install():
    path = _patch_path()
    if path is None:
        return
    if any(isinstance(finder, _PatchedSourceFinder) for finder in sys.meta_path):
        return
    sys.meta_path.insert(0, _PatchedSourceFinder(path))


install()