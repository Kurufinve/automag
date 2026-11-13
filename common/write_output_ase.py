"""
automag.common.write_output_ase
===============================

Generic post-run writer for calculators that expose a
:class:`~automag.calculators.base.CalculatorInterface`.

It appends one line to ``CalcFold/<formula><suffix>_<mode>_<calculator>.txt``
in the same format as the VASP ``write_output.py``, so the existing
post-processing scripts (``2_coll/2_plot_results.py``,
``0_conv_tests/2_plot_results.py``, ``3_monte_carlo``, ...) work unchanged.

Usage (inside a job script, cwd is the calculation directory)::

    python /path/to/automag/common/write_output_ase.py <calculator> <mode> <suffix>

``mode`` is ``None`` (or ``coll``) for a collinear single-point calculation,
or ``encut`` / ``kgrid`` for convergence tests.

.. codeauthor:: Automag contributors
"""

import os
import re
import sys
from pathlib import Path

# The script is executed in place from the automag checkout, so its parent
# directory is the code root even when the job runs in another directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from calculators import get_calculator


def _data_root():
    """Project directory the results belong to.

    A job runs inside ``<project>/CalcFold/<...>/singlepoint``, so walking up
    to the ``CalcFold`` parent is the most reliable answer - it does not depend
    on the environment the scheduler handed to the job (a job started by the
    pipeline may not carry AUTOMAG_WORKDIR, and falling back to the code root
    would write the results into the automag checkout).
    """
    current = os.path.abspath(os.getcwd())
    while True:
        if os.path.basename(current) == "CalcFold":
            return os.path.dirname(current)
        parent = os.path.dirname(current)
        if parent == current:
            break
        current = parent
    for variable in ("AUTOMAG_WORKDIR", "AUTOMAG_PATH"):
        value = os.environ.get(variable)
        if value:
            return value
    return os.getcwd()


def _initial_magmoms(directory):
    candidates = [
        Path(directory) / "initial_magmoms.txt",
        Path(directory).parent / "singlepoint" / "initial_magmoms.txt",
    ]
    for candidate in candidates:
        if candidate.is_file():
            text = candidate.read_text().strip()
            if text:
                return [float(item) for item in text.split()]
    return []


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    calculator = argv[0] if argv else "vasp"
    calc_mode = argv[1] if len(argv) > 1 else "None"
    struct_suffix = argv[2] if len(argv) > 2 else ""

    if calc_mode.lower() in {"none", "coll", ""}:
        calc_mode = None

    directory = os.getcwd()
    interface = get_calculator(calculator)
    atoms = interface.read_structure(directory)
    energy = interface.read_energy(directory)

    try:
        converged = bool(interface.is_converged(directory))
    except Exception:  # noqa: BLE001 - assume converged if unknown
        converged = True

    state = os.path.basename(os.path.dirname(directory))
    match = re.search(r"\d+", state)
    calcid = int(match.group()) if match else 0

    formula = atoms.get_chemical_formula(mode="metal")
    try:
        structure = AseAtomsAdaptor.get_structure(atoms)
        spacegroup = SpacegroupAnalyzer(structure).get_space_group_symbol()
    except Exception:  # noqa: BLE001 - fall back to a neutral label
        spacegroup = "P1"

    convergence = ("convergedconverged" if converged
                   else "NONCONVERGEDNONCONVERGED")
    output_line = (
        f"{state}       {calcid}  singlepoint={convergence}  {formula}   "
        f"{spacegroup:10s}  energy={energy}\n"
    )

    if calc_mode is None:
        initial = _initial_magmoms(directory)
        if not initial:
            try:
                initial = [float(item)
                           for item in atoms.get_initial_magnetic_moments()]
            except Exception:  # noqa: BLE001
                initial = [0.0] * len(atoms)

        try:
            final = interface.read_magnetic_moments(directory)
        except Exception:  # noqa: BLE001
            final = None
        if final is None:
            final = initial

        output_line += "          initial_magmoms={}  ".format(
            np.array2string(np.array(initial), 10000)
        )
        output_line += "final_magmoms={}\n".format(
            np.array2string(np.array(final), 10000)
        )

    if calc_mode is None:
        filename = f"{formula}{struct_suffix}_singlepoint_{calculator}.txt"
    else:
        filename = f"{formula}{struct_suffix}_{calc_mode}_{calculator}.txt"

    with open(os.path.join(_data_root(), "CalcFold", filename), "a") as handle:
        handle.write(output_line)

    print(f"Wrote {filename}: {output_line.strip()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
