"""
automag.common.input_loader
===========================

Load the settings of a stage (or of the whole pipeline) from a file.

Three formats are accepted and tried in this fixed order:

1. ``<base>.py``      - executed, its public module-level names are returned;
2. ``<base>.yaml``    - ``yaml.safe_load`` (``.yml`` is also accepted);
3. ``<base>.json``    - ``json.load``.

So ``input.py`` wins over ``input.yaml``, which wins over ``input.json``.
``<base>`` defaults to ``input`` for stage settings and to ``input_automag``
for the controlling file.

The historical stage scripts used ``from input import *``, which only works for
a file called exactly ``input.py`` in the current directory.  They now call
:func:`apply_stage_input`, which keeps the same "bare globals" behaviour but
supports every format and an explicit path.

YAML/JSON files hold a flat mapping of settings::

    calculator: vasp
    struct_suffix: ""
    configuration: lowest
    params:
      encut: 600
      kpts: 4

Values that the pipeline itself writes as ``range(...)`` become plain lists -
the stage code only iterates over them.

.. codeauthor:: Automag contributors
"""

import importlib.util
import json
import sys
from pathlib import Path

#: Supported input file formats, in priority order.
SUFFIXES = (".py", ".yaml", ".yml", ".json")


def resolve_input(base="input", path=None):
    """Return the input file to use, honouring the priority above.

    :param base: file stem to look for in the current directory (default
        ``"input"``).
    :param path: explicit file (or stem) taken from the command line; when it
        does not exist with any supported suffix, :class:`FileNotFoundError` is
        raised.
    """
    if path is None:
        for suffix in SUFFIXES:
            candidate = Path.cwd() / f"{base}{suffix}"
            if candidate.is_file():
                return candidate.resolve()
        tried = ", ".join(f"{base}{suffix}" for suffix in SUFFIXES)
        raise FileNotFoundError(
            f"No input settings found: tried {tried} in {Path.cwd()}.")

    given = Path(path).expanduser()
    if not given.is_absolute():
        given = Path.cwd() / given
    if given.is_file():
        return given.resolve()
    if not given.suffix:
        for suffix in SUFFIXES:
            candidate = given.with_suffix(suffix)
            if candidate.is_file():
                return candidate.resolve()
    tried = ", ".join(str(given.with_suffix(suffix)) for suffix in SUFFIXES)
    raise FileNotFoundError(f"Input file not found: {given} (tried {tried}).")


def load_input(path=None, base="input"):
    """Return the settings of an input file as a dictionary.

    ``path`` may be a full name (``my_input.yaml``), a stem (``my_input``) or
    ``None`` to search for ``<base>.py`` / ``.yaml`` / ``.yml`` / ``.json`` in
    the current directory.
    """
    resolved = resolve_input(base=base, path=path)
    suffix = resolved.suffix.lower()
    if suffix == ".py":
        return _load_python(resolved)
    if suffix in (".yaml", ".yml"):
        return _as_mapping(_load_yaml(resolved), resolved)
    if suffix == ".json":
        return _as_mapping(_load_json(resolved), resolved)
    raise ValueError(
        f"Unsupported input format '{suffix}' ({resolved}). "
        f"Use one of: {', '.join(SUFFIXES)}.")


def apply_stage_input(namespace, argv=None, base="input"):
    """Load the stage settings and inject them into ``namespace``.

    Replaces the historical ``from input import *`` prologue::

        from common.input_loader import apply_stage_input
        apply_stage_input(globals())

    :param namespace: usually ``globals()`` of the stage script.
    :param argv: command line to inspect for an explicit input file
        (defaults to ``sys.argv``).
    """
    argv = list(sys.argv if argv is None else argv)
    explicit = argv[1] if len(argv) > 1 else None
    resolved = resolve_input(base=base, path=explicit)
    print(f"Using the input file {resolved}")
    settings = load_input(resolved)
    namespace.update(settings)
    return settings


def _load_python(path):
    spec = importlib.util.spec_from_file_location("_automag_stage_input", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load input file: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return {key: value for key, value in vars(module).items()
            if not key.startswith("_")}


def _load_yaml(path):
    try:
        import yaml
    except ImportError as error:  # pragma: no cover - depends on the env
        raise ImportError(
            f"Reading {path} requires PyYAML (pip install pyyaml).") from error
    with open(path, "rt") as handle:
        return yaml.safe_load(handle) or {}


def _load_json(path):
    with open(path, "rt") as handle:
        return json.load(handle) or {}


def _as_mapping(data, path):
    if not isinstance(data, dict):
        raise ValueError(
            f"{path} must contain a mapping of settings at the top level, "
            f"got {type(data).__name__}.")
    return data
