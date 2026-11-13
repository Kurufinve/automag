"""
automag.common.results
======================

Helpers for reading the artefacts produced by the ``2_coll`` stage.

A finished ``2_coll`` run stores everything in a *bundle* directory::

    2_coll/<instance>/<formula><suffix>_<calculator>/
        trials_<formula><suffix>/settingNNN.vasp
        trials_<formula><suffix>/configurationsNNN.txt
        <formula><suffix>_<calculator>_statesNNN.txt
        <formula><suffix>_<calculator>_energiesNNN.txt
        <formula><suffix>_<calculator>_settingNNN.vasp
        <formula><suffix>_<calculator>_most_stable_configuration.txt

Downstream stages therefore have to say *which calculator's* collinear results
they want (:func:`coll_bundle`), because a pipeline may mix calculators, e.g.
``2_coll`` with a spin-lattice model followed by ``MAE`` with VASP.

Keeping this logic in one place avoids re-implementing the (slightly
inconsistent) folder layout of the collinear stage in every consumer.

.. codeauthor:: Automag contributors
"""

import re
from pathlib import Path

import numpy as np
from pymatgen.core.structure import Structure


def coll_root(data_root, formula, struct_suffix=""):
    """Return the directory where legacy ``2_coll`` results for ``formula`` live."""
    root = Path(data_root) / "2_coll"
    for candidate in (root / f"{formula}{struct_suffix}", root):
        if candidate.is_dir():
            return candidate
    return root


def _newest(paths):
    paths = list(paths)
    if not paths:
        return None
    return max(paths, key=lambda item: item.stat().st_mtime)


def bundle_name(formula, struct_suffix="", calculator="vasp"):
    """Return the basename of the ``2_coll`` bundle for a calculator."""
    return f"{formula}{struct_suffix}_{calculator}"


def coll_bundle(data_root, formula, struct_suffix="", calculator="vasp"):
    """Return the newest ``.../2_coll/**/<formula><suffix>_<calculator>`` dir.

    The bundle may live directly in ``2_coll/`` (single instance) or in a
    per-instance sub-directory such as ``2_coll/coll_prophet/``.
    """
    root = Path(data_root) / "2_coll"
    names = (bundle_name(formula, struct_suffix, calculator),
             f"{formula}_{calculator}")
    if root.is_dir():
        for name in names:
            match = _newest(path for path in root.rglob(name)
                            if path.is_dir())
            if match is not None:
                return match
    raise FileNotFoundError(
        f"No '{names[0]}' directory found under {root}. Run "
        "2_coll/2_plot_results.py (the 'postprocess' phase) first.")


def available_coll_calculators(data_root, formula, struct_suffix=""):
    """Return the set of calculators that produced a ``2_coll`` bundle."""
    root = Path(data_root) / "2_coll"
    prefix = f"{formula}{struct_suffix}_"
    found = set()
    if not root.is_dir():
        return found
    for path in root.rglob(f"{prefix}*"):
        if not path.is_dir():
            continue
        name = path.name
        # Ignore the convergence-style artefacts (files are filtered above)
        found.add(name[len(prefix):])
    return found


def _search_roots(data_root, formula, struct_suffix, calculator):
    roots = []
    try:
        roots.append(coll_bundle(data_root, formula, struct_suffix, calculator))
    except FileNotFoundError:
        pass
    roots.append(coll_root(data_root, formula, struct_suffix))
    return roots


def most_stable_state(data_root, formula, struct_suffix="", calculator="vasp"):
    """Return the name of the most stable configuration found by ``2_coll``."""
    names = (
        f"{formula}{struct_suffix}_{calculator}_most_stable_configuration.txt",
        f"{formula}_{calculator}_most_stable_configuration.txt",
    )
    for root in _search_roots(data_root, formula, struct_suffix, calculator):
        for name in names:
            match = _newest(root.rglob(name))
            if match is not None:
                return match.read_text().strip()
    raise FileNotFoundError(
        f"No '{names[0]}' found under {coll_root(data_root, formula, struct_suffix)}. "
        "Run 2_coll/2_plot_results.py first."
    )


def find_trials_dir(data_root, formula, struct_suffix="", calculator="vasp"):
    """Locate the ``trials_*`` folder that stores the generated configurations."""
    names = (f"trials_{formula}{struct_suffix}", f"trials_{formula}")
    for root in _search_roots(data_root, formula, struct_suffix, calculator):
        if not root.is_dir():
            continue
        for name in names:
            match = _newest(path for path in root.rglob(name) if path.is_dir())
            if match is not None:
                return match
        match = _newest(path for path in root.rglob("trials_*") if path.is_dir())
        if match is not None:
            return match
    raise FileNotFoundError(
        f"No trials_* folder found for '{formula}' "
        f"({struct_suffix or 'no suffix'}, {calculator}).")


def resolve_configuration(configuration, data_root, formula, struct_suffix="",
                          calculator="vasp"):
    """Resolve a configuration name for a downstream stage.

    ``"lowest"`` (case-insensitive), ``None`` and an empty string mean "the
    most stable configuration found by 2_coll".  Any other value is returned
    unchanged, so existing inputs keep working.

    Raises :class:`FileNotFoundError` if ``"lowest"`` is requested but the
    2_coll post-processing has not produced a most-stable file yet.
    """
    if configuration is None:
        return most_stable_state(data_root, formula, struct_suffix, calculator)
    name = str(configuration).strip()
    if not name or name.lower() in {"lowest", "most_stable", "most-stable"}:
        return most_stable_state(data_root, formula, struct_suffix, calculator)
    return configuration


def find_state(data_root, formula, state, struct_suffix="", calculator="vasp"):
    """Return geometry and spin pattern for a named configuration.

    :return: dict with ``setting``, ``spins``, ``structure`` (Path) and
        ``trials`` (Path).
    """
    trials = find_trials_dir(data_root, formula, struct_suffix, calculator)
    for configurations_file in sorted(trials.glob("configurations*.txt")):
        match = re.search(r"configurations(\d+)", configurations_file.name)
        if match is None:
            continue
        setting = int(match.group(1))
        for line in configurations_file.read_text().splitlines():
            values = line.split()
            if values and values[0] == state:
                return {
                    "setting": setting,
                    "spins": [int(item) for item in values[1:]],
                    "structure": trials / f"setting{setting:03d}.vasp",
                    "trials": trials,
                }
    raise KeyError(f"Configuration '{state}' not found under {trials}.")


def final_magmoms(data_root, formula, state, struct_suffix="", calculator="vasp"):
    """Return the converged magnetic moments of ``state`` (or ``None``)."""
    path = (
        Path(data_root)
        / "CalcFold"
        / f"{formula}{struct_suffix}_singlepoint_{calculator}.txt"
    )
    if not path.is_file():
        return None
    lines = path.read_text().splitlines()
    for index, line in enumerate(lines):
        values = line.split()
        if values and values[0] == state:
            for following in lines[index + 1:index + 3]:
                match = re.search(r"final_magmoms=\[([^\]]*)\]", following)
                if match:
                    return np.array(
                        [float(item) for item in match.group(1).split()]
                    )
    return None


def initial_magmoms(structure, spins, spin_values=None, nonmagnetic=0.0,
                    magnetic_atoms=None):
    """Build scalar magnetic moments from a spin pattern and ``spin_values``.

    ``2_coll`` and ``4_mae`` store the pattern for the *magnetic* atoms only,
    while the geometry contains all atoms.  If ``len(spins)`` is smaller than
    the number of sites, the pattern is expanded over the magnetic species
    (taken from ``magnetic_atoms``, or from the keys of ``spin_values``).
    """
    spin_values = spin_values or {}
    symbols = [site.specie.symbol for site in structure]
    spins = list(spins)

    if len(spins) != len(symbols):
        if magnetic_atoms:
            magnetic = set(magnetic_atoms)
        else:
            magnetic = set(spin_values)
        mask = [symbol in magnetic for symbol in symbols]
        if sum(mask) != len(spins):
            raise ValueError(
                f"Spin pattern has {len(spins)} entries but the structure has "
                f"{len(symbols)} atoms ({sum(mask)} of them magnetic)."
            )
        expanded, index = [], 0
        for is_magnetic in mask:
            if is_magnetic:
                expanded.append(spins[index])
                index += 1
            else:
                expanded.append(0)
        spins = expanded

    magmoms = []
    for symbol, spin in zip(symbols, spins):
        if spin == 0:
            magmoms.append(0.0)
            continue
        if symbol in spin_values:
            value = abs(float(spin_values[symbol][0]))
        else:
            value = abs(float(nonmagnetic))
        magmoms.append(spin * value)
    return np.array(magmoms)


def load_state(data_root, formula, state, struct_suffix="", calculator="vasp",
               spin_values=None):
    """Return ``(info, structure, magmoms)`` for a configuration.

    Converged magnetic moments from ``CalcFold`` are preferred; if they are
    not available (e.g. the collinear jobs are still running) the initial
    moments are reconstructed from ``spin_values``.
    """
    info = find_state(data_root, formula, state, struct_suffix,
                      calculator=calculator)
    structure = Structure.from_file(str(info["structure"]))
    magmoms = final_magmoms(
        data_root, formula, state, struct_suffix=struct_suffix,
        calculator=calculator
    )
    if magmoms is None or len(magmoms) != len(structure):
        magmoms = initial_magmoms(structure, info["spins"], spin_values)
    return info, structure, magmoms
