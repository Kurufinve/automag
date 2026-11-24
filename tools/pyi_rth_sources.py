"""Make ``inspect.getsource()`` work inside a frozen build.

PyInstaller loads modules from the PYZ archive, so ``linecache`` has no entry for
them and ``inspect.getsource`` raises ``OSError: could not get source code``.
Several libraries need the source at import time:

* ``torch.utils._config_module`` reads its own module to find the
  ``compile_ignored`` comments (this aborts ``import torch``);
* ``torch/_sources.py`` (TorchScript) and e3nn/torch_geometric JIT-compile
  functions defined in Python.

The build ships the ``.py`` files as data next to ``sys._MEIPASS``, so this hook
falls back to reading them from disk.  It is installed as a PyInstaller runtime
hook, i.e. before the application starts.
"""

import inspect
import linecache
import os
import sys
import types

_original_getsource = inspect.getsource
_original_getsourcelines = inspect.getsourcelines


def _module_path(obj):
    """Best guess of the .py file that defines ``obj``."""
    module = obj if isinstance(obj, types.ModuleType) else inspect.getmodule(obj)
    path = getattr(module, "__file__", None)
    if path is None:
        code = getattr(obj, "__code__", None)
        path = getattr(code, "co_filename", None)
    if not path or not path.endswith(".py"):
        return None
    if not os.path.isfile(path):
        # The build time path (<env>/lib/pythonX.Y/site-packages/torch/...) does
        # not exist in the bundle: map the part after site-packages/ onto
        # <bundle>/_internal/ (and onto <bundle>/ for the stage scripts).
        bundle = getattr(sys, "_MEIPASS", None)
        marker = "%ssite-packages%s" % (os.sep, os.sep)
        if bundle and marker in path:
            relative = path.split(marker, 1)[1]
            for root in (os.path.join(bundle, "_internal"), bundle):
                candidate = os.path.join(root, relative)
                if os.path.isfile(candidate):
                    return candidate
    return path if os.path.isfile(path) else None


def _read(path):
    linecache.checkcache(path)
    lines = linecache.getlines(path)
    if not lines:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            lines = handle.readlines()
        linecache.cache[path] = (len("\n".join(lines)), None, lines, path)
    return lines


def getsource(obj, *args, **kwargs):
    try:
        return _original_getsource(obj, *args, **kwargs)
    except (OSError, TypeError):
        path = _module_path(obj)
        if path is None:
            raise
        return "".join(_read(path))


def getsourcelines(obj, *args, **kwargs):
    try:
        return _original_getsourcelines(obj, *args, **kwargs)
    except (OSError, TypeError):
        path = _module_path(obj)
        if path is None:
            raise
        return _read(path), 1


inspect.getsource = getsource
inspect.getsourcelines = getsourcelines
