"""
automag.common.magnetic
=======================

Magnetic moments and magnetic-configuration enumeration.

Two jobs
--------
1. **Defaults.**  ``spin_values`` and ``magnetic_atoms`` are optional in the
   input file.  When they are missing, the moments come from
   :data:`pymatgen.analysis.magnetism.analyzer.DEFAULT_MAGMOMS` (falling back
   to the ``magmom`` site property of the structure, then to 0, i.e.
   non-magnetic).  This is what :func:`resolve_spin_values` and
   :func:`resolve_magnetic_elements` do.

2. **Enumeration backends.**  ``2_coll`` can generate the magnetic
   configurations either with the original hand-written enumlib driver or with
   :class:`pymatgen.analysis.magnetism.analyzer.MagneticStructureEnumerator`
   (``enumeration_mode = "pymatgen"``).  The latter understands the usual
   strategies (ferromagnetic, antiferromagnetic, antiferromagnetic_by_motif,
   ferrimagnetic_by_motif, ferrimagnetic_by_species, nonmagnetic) and returns
   :class:`~pymatgen.core.structure.Structure` objects whose spin is encoded in
   the species (``Fe,spin=np.float64(-5.0)``).

Both backends return the same three lists used by ``2_coll/1_submit.py``
(``lattices``, ``coordinates``, ``configurations``), so everything downstream
(``2_plot_results.py``, ``3_monte_carlo``, ``4_mae``) is unchanged.

.. codeauthor:: Automag contributors
"""

import os
import re
from pathlib import Path

import numpy as np
from pymatgen.analysis.magnetism.analyzer import DEFAULT_MAGMOMS
from pymatgen.core.structure import Structure

#: Strategies understood by ``MagneticStructureEnumerator``.
PYMATGEN_STRATEGIES = (
    "ferromagnetic",
    "antiferromagnetic",
    "antiferromagnetic_by_motif",
    "ferrimagnetic_by_motif",
    "ferrimagnetic_by_species",
    "nonmagnetic",
)

DEFAULT_STRATEGIES = ("ferromagnetic", "antiferromagnetic")

_SPIN_PATTERN = re.compile(r"spin\s*=\s*(?:np\.)?(?:float\d*)?\(?\s*([+-]?\d*\.?\d+)")


def default_moment(symbol, overrides=None):
    """Return the default magnetic moment of an element (0 if not magnetic)."""
    if overrides and symbol in overrides:
        value = overrides[symbol]
        return _numeric(value[0] if isinstance(value, (list, tuple)) else value)
    return _numeric(DEFAULT_MAGMOMS.get(symbol, 0.0))


def element_symbol(site):
    """Return the element symbol of a site (``specie`` may be Element or Species)."""
    specie = getattr(site, "specie", None)
    if specie is not None:
        element = getattr(specie, "element", None)
        symbol = getattr(element if element is not None else specie, "symbol", None)
        if symbol:
            return str(symbol)
    elements = getattr(getattr(site, "species", None), "elements", None)
    if elements:
        return str(elements[0].symbol)
    raise ValueError(f"Cannot determine the element of the site {site!r}.")


def resolve_magnetic_elements(structure, spin_values=None, magnetic_atoms=None):
    """Return the symbols treated as magnetic.

    Explicit ``magnetic_atoms`` wins; otherwise every element with a non-zero
    moment (from ``spin_values`` or ``DEFAULT_MAGMOMS``) is magnetic.
    """
    if magnetic_atoms:
        return [str(symbol) for symbol in magnetic_atoms]
    known = set(spin_values or {})
    result = []
    for symbol in _structure_symbols(structure):
        if symbol in known:
            moment = _first(spin_values[symbol])
        else:
            moment = default_moment(symbol)
        if moment:
            result.append(symbol)
    return result


def resolve_spin_values(structure, spin_values=None, magnetic_atoms=None,
                        overrides=None):
    """Return ``{symbol: [value, ...]}`` for every magnetic element.

    ``spin_values`` entries may hold one value (the magnitude) or two (the two
    magnitudes of a species that is split into an up and a down sublattice).
    Missing entries are filled from ``DEFAULT_MAGMOMS`` - via
    ``magnetic_atoms`` when given, otherwise for every element of the
    structure with a tabulated moment.  Elements that are not magnetic are
    omitted, which is exactly what the callers expect: a symbol absent from
    the dictionary is non-magnetic.
    """
    resolved = {}
    for symbol, value in (spin_values or {}).items():
        values = [_numeric(item) for item in _as_list(value)]
        if not values:
            continue
        if len(values) > 2:
            raise ValueError(
                f"spin_values['{symbol}'] = {value}: at most 2 values are "
                "supported (one magnitude, or two for a split sublattice).")
        resolved[str(symbol)] = values

    for symbol in resolve_magnetic_elements(structure, resolved, magnetic_atoms):
        if symbol in resolved:
            continue
        moment = default_moment(symbol, overrides)
        if moment:
            resolved[symbol] = [moment]
    return resolved


def structure_moments(structure, spin_values=None, magnetic_atoms=None,
                      overrides=None):
    """Return one scalar moment per site (0 for the non-magnetic atoms)."""
    resolved = resolve_spin_values(structure, spin_values, magnetic_atoms,
                                   overrides)
    magmoms = []
    for site in structure:
        symbol = element_symbol(site)
        magmoms.append(resolved[symbol][0] if symbol in resolved else 0.0)
    return magmoms


# ---------------------------------------------------------------------------
# Enumeration with pymatgen's MagneticStructureEnumerator
# ---------------------------------------------------------------------------
def ensure_enumlib_on_path():
    """Put enumlib where pymatgen expects it (it searches ``PATH``)."""
    from common.external import tool

    entry = tool("enumlib")
    directories = []
    for key in ("enum", "makestr"):
        value = entry.get(key)
        if value:
            candidate = Path(str(value)).expanduser()
            if not candidate.is_absolute():
                continue
            directories.append(str(candidate.parent))
    if not directories:
        return
    current = os.environ.get("PATH", "")
    known = current.split(os.pathsep) if current else []
    prefix = [item for item in dict.fromkeys(directories) if item not in known]
    if prefix:
        os.environ["PATH"] = os.pathsep.join(prefix + known)


def site_spin(site):
    """Return the spin stored in a pymatgen species (or ``None``)."""
    specie = getattr(site, "specie", None)
    value = getattr(specie, "spin", None) if specie is not None else None
    if value is not None:
        try:
            return float(value)
        except (TypeError, ValueError):
            pass
    match = _SPIN_PATTERN.search(str(site.species_string))
    return float(match.group(1)) if match else None


def pymatgen_configurations(structure, spin_values=None,
                            magnetic_atoms=None, overrides=None,
                            strategies=None, max_orderings=None,
                            truncate_by_symmetry=True, automatic=True):
    """Enumerate the magnetic configurations with pymatgen.

    :return: ``(lattices, coordinates, configurations)`` in the layout used by
        ``2_coll/1_submit.py``: one entry per distinct setting (lattice +
        fractional coordinates) and, per setting, the spin patterns as lists of
        ``+1`` / ``0`` / ``-1``.
    """
    from pymatgen.analysis.magnetism.analyzer import MagneticStructureEnumerator

    ensure_enumlib_on_path()
    resolved = resolve_spin_values(structure, spin_values, magnetic_atoms,
                                   overrides)
    moments = structure_moments(structure, resolved, magnetic_atoms, overrides)
    magnetic = resolve_magnetic_elements(structure, resolved, magnetic_atoms)
    if not magnetic:
        raise ValueError(
            "No magnetic element found. Set 'spin_values' (e.g. "
            "{'Fe': [5]}) or 'magnetic_atoms' in the 2_coll settings.")

    strategies = tuple(strategies or DEFAULT_STRATEGIES)
    unknown = [item for item in strategies if item not in PYMATGEN_STRATEGIES]
    if unknown:
        raise ValueError(
            f"Unknown enumeration strategy: {unknown}. "
            f"Use one of {', '.join(PYMATGEN_STRATEGIES)}.")

    ordered = structure.copy()
    ordered.add_site_property("magmom", moments)
    enumerator = MagneticStructureEnumerator(
        ordered,
        default_magmoms={symbol: resolved[symbol][0] for symbol in resolved},
        strategies=strategies,
        automatic=automatic,
        truncate_by_symmetry=truncate_by_symmetry,
    )
    candidates = enumerator.ordered_structures
    origins = enumerator.ordered_structure_origins
    if max_orderings is not None:
        candidates = candidates[:int(max_orderings)]
        origins = origins[:int(max_orderings)]
    if not candidates:
        raise ValueError(
            "MagneticStructureEnumerator produced no configuration with "
            f"strategies {list(strategies)}.")

    lattices, coordinates, configurations = [], [], []
    for candidate, origin in zip(candidates, origins):
        # Spins are encoded in the species; a missing spin means non-magnetic.
        pattern = []
        for new_site, old_site in zip(candidate, structure):
            symbol = element_symbol(old_site)
            if symbol not in resolved:
                pattern.append(0)
                continue
            spin = site_spin(new_site)
            if spin is None:
                spin = resolved[symbol][0]
            pattern.append(int(np.sign(spin)))

        # A "setting" is one cell (lattice + site order); the configurations of
        # that cell differ only in their spin pattern.  Grouping by the lattice
        # alone keeps the same layout as the enumlib driver, where
        # ``configurationsNNN.txt`` lists all patterns of ``settingNNN.vasp``.
        key = _round_array(candidate.lattice.matrix)
        index = None
        for position, existing in enumerate(lattices):
            if _round_array(existing.matrix) == key:
                index = position
                break
        if index is None:
            lattices.append(candidate.lattice)
            coordinates.append(candidate.frac_coords)
            configurations.append([])
            index = len(lattices) - 1
        if pattern not in configurations[index]:
            configurations[index].append(pattern)

    return lattices, coordinates, configurations, list(origins)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _as_list(value):
    if isinstance(value, (list, tuple, np.ndarray)):
        return list(value)
    return [value]


def _numeric(value):
    """Keep integral moments integral (``4.0`` -> ``4``), which the stage
    scripts format with ``{:2d}`` and pass on as spin multiplicities."""
    number = float(value)
    return int(round(number)) if number.is_integer() else number


def _first(value):
    values = _as_list(value)
    return _numeric(values[0]) if values else 0


def _structure_symbols(structure):
    seen = []
    for site in structure:
        symbol = element_symbol(site)
        if symbol not in seen:
            seen.append(symbol)
    return seen


def _round_array(values):
    return tuple(tuple(round(float(item), 5) for item in row) for row in values)
