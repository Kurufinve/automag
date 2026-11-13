"""
automag.common.sweep
====================

Parameter lists in the ``params`` dictionary of a stage.

A value that is a *list* instead of a single number means "run the stage for
every value and compare the results" - the convergence study of the MAE stage::

    "params": {
        "encut": [400, 500, 600],      # three sets
        "kpts":  [10, 20],             # x2 sets -> six in total
        "sigma": 0.1,                  # fixed
        "ldauu": [5.0, 0.0],           # NOT a sweep: per-species list
    }

Only the keys in :data:`SWEEPABLE` are expanded.  Every other list (``ldauu``,
``ldauj``, ``ldaul``, ``setups``, ``saxis``, ``magmom``, an explicit
``kpts`` mesh, ...) keeps its meaning.

.. codeauthor:: Automag contributors
"""

import itertools
import json
import os
from dataclasses import dataclass, field

#: Keys whose list value means "run a convergence study".
SWEEPABLE = (
    "encut", "kpts", "sigma", "ediff", "ediffg", "prec", "nelm", "ncore",
    "lreal", "ismear", "kgamma", "kspacing", "lasph", "amix", "bmix",
    "nsw", "ibrion", "isif",
)

#: Keys that are lists for a different reason (never expanded).
NEVER_SWEPT = (
    "ldauu", "ldauj", "ldaul", "setups", "saxis", "magmom", "magmoms",
    "positions", "cells", "species", "spin_values", "magnetic_atoms",
    "kpts_values", "encut_values", "sigma_values", "perturbations",
)


def kpts_label(kpts):
    """Compact label of a k-point specification (scalar or explicit mesh)."""
    if isinstance(kpts, (list, tuple)):
        return "x".join(str(item) for item in kpts)
    return str(kpts)


def is_sweep(key, value):
    """True when ``params[key] = value`` asks for a convergence study."""
    if key in NEVER_SWEPT or key not in SWEEPABLE:
        return False
    if not isinstance(value, (list, tuple)):
        return False
    # An explicit k-point mesh is a list of three integers, not a sweep.
    if key == "kpts" and len(value) == 3 and all(
            isinstance(item, (int, float)) and not isinstance(item, bool)
            for item in value):
        return False
    return len(value) > 1


@dataclass
class SweepSet:
    """One parameter combination of a convergence study."""

    index: int
    params: dict
    values: dict = field(default_factory=dict)

    @property
    def tag(self):
        """Suffix that distinguishes the calculation directories of the set."""
        return f"_S{self.index + 1:02d}" if self.values else ""

    @property
    def label(self):
        if not self.values:
            return "default"
        return ", ".join(f"{key}={value}" for key, value in self.values.items())

    @property
    def ename(self):
        """Short, file-system friendly name of the swept combination."""
        if not self.values:
            return "default"
        parts = []
        for key, value in self.values.items():
            text = kpts_label(value)
            text = text.replace(".", "p").replace("-", "m")
            parts.append(f"{key}{text}")
        return "_".join(parts)


def sweep_keys(params):
    """Return the keys of ``params`` that are convergence studies."""
    return [key for key, value in params.items() if is_sweep(key, value)]


def expand_params(params, extra=None):
    """Expand the parameter lists of ``params`` into :class:`SweepSet` objects.

    ``extra`` is merged over every set (handy for per-stage overrides).  The
    result always has at least one entry; when nothing is swept, that single
    entry has an empty ``tag`` and keeps the historical file names.
    """
    params = dict(params or {})
    keys = sweep_keys(params)
    if not keys:
        merged = dict(params)
        if extra:
            merged.update(extra)
        return [SweepSet(index=0, params=merged)]

    values = [list(params[key]) for key in keys]
    sets = []
    for index, combination in enumerate(itertools.product(*values)):
        merged = dict(params)
        chosen = dict(zip(keys, combination))
        merged.update(chosen)
        if extra:
            merged.update(extra)
        sets.append(SweepSet(index=index, params=merged, values=chosen))
    return sets


#: File in the stage directory listing the parameter sets of the current run.
MANIFEST_NAME = "sweep_sets.json"


def write_manifest(workdir, directories):
    """Record the calculation directories of the current parameter sets.

    The post-processing scripts read this file instead of scanning the
    directory tree, so leftovers of an earlier (different) parameter set are
    never analysed together with the current results.
    """
    path = os.path.join(str(workdir), MANIFEST_NAME)
    with open(path, "wt") as handle:
        json.dump({"directories": [str(item) for item in directories]},
                  handle, indent=2)
    return path


def read_manifest(workdir):
    """Directories recorded by :func:`write_manifest` (``None`` when absent)."""
    path = os.path.join(str(workdir), MANIFEST_NAME)
    if not os.path.isfile(path):
        return None
    try:
        with open(path, "rt") as handle:
            data = json.load(handle)
    except (OSError, ValueError):
        return None
    directories = [str(item) for item in (data.get("directories") or [])
                   if os.path.isdir(item)]
    return directories or None


def read_parameters(directory):
    """``(parameters, tag, label)`` stored next to one parameter set."""
    path = os.path.join(str(directory), "parameters.json")
    if not os.path.isfile(path):
        return None, "", "default"
    try:
        with open(path, "rt") as handle:
            meta = json.load(handle)
    except (OSError, ValueError):
        return None, "", "default"
    return (meta.get("parameters"), meta.get("tag", ""),
            meta.get("label", "default"))


def summary_table(rows, swept_labels=None):
    """Render a small text table of ``rows`` (lists of strings)."""
    if not rows:
        return ""
    widths = [max(len(str(row[column])) for row in rows)
              for column in range(len(rows[0]))]
    lines = []
    for number, row in enumerate(rows):
        line = "  ".join(str(value).ljust(width)
                         for value, width in zip(row, widths))
        lines.append(line.rstrip())
        if number == 0:
            lines.append("  ".join("-" * width for width in widths))
    return "\n".join(lines)


def convergence_report(entries, quantity="MAE"):
    """Compare one quantity across parameter sets.

    ``entries`` is a list of ``(label, value, extra)`` triples; the function
    returns a text table plus the relative change between neighbouring sets,
    which is what a convergence check looks at.
    """
    if not entries:
        return ""
    rows = [["set", "parameters", quantity]]
    for label, value, _ in entries:
        rows.append([label[0], label[1],
                     "n/a" if value is None else f"{value:.6f}"])
    text = [summary_table(rows)]
    numeric = [(label, value) for label, value, _ in entries
               if value is not None]
    if len(numeric) > 1:
        text.append("")
        text.append(f"relative change of {quantity} between successive sets:")
        for (previous, before), (current, after) in zip(numeric, numeric[1:]):
            if before == 0:
                continue
            text.append(
                f"  {previous[1]} -> {current[1]}: "
                f"{100.0 * abs(after - before) / abs(before):.2f} %")
    return "\n".join(text)
