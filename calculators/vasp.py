"""
automag.calculators.vasp
========================

Reference implementation of :class:`~automag.calculators.base.CalculatorInterface`
for VASP, built on :class:`ase.calculators.vasp.Vasp`.

.. codeauthor:: Automag contributors
"""

import os

import numpy as np
from ase import Atoms
from ase.io import read
from ase.calculators.vasp import Vasp

from common.paths import code_root
from calculators.base import CalculatorInterface


def read_single(path) -> Atoms:
    """Read one structure, unwrapping the list ASE may return."""
    result = read(path)
    if isinstance(result, list):
        result = result[-1]
    return result


class VaspInterface(CalculatorInterface):
    name = "vasp"
    supports_ldau = True
    supports_noncollinear = True
    supports_relaxation = True
    supports_charges = True
    output_marker = "OUTCAR"

    def run_command(self):
        return self.command or "mpirun vasp_std"

    def write_input(self, atoms, directory, params=None):
        merged = dict(self.params)
        merged.update(params or {})
        merged.pop("directory", None)
        merged.pop("atoms", None)
        os.makedirs(directory, exist_ok=True)
        calc = Vasp(atoms=atoms, directory=directory, **merged)
        calc.write_input(atoms)
        return calc

    def read_energy(self, directory):
        atoms = read_single(os.path.join(directory, "OUTCAR"))
        return float(atoms.get_potential_energy())

    def read_magnetic_moments(self, directory):
        atoms = read_single(os.path.join(directory, "OUTCAR"))
        try:
            return np.asarray(atoms.get_magnetic_moments())
        except Exception:  # noqa: BLE001 - moments are optional
            return np.zeros(len(atoms))

    def read_structure(self, directory):
        for name in ("CONTCAR", "OUTCAR"):
            path = os.path.join(directory, name)
            if os.path.isfile(path) and os.path.getsize(path) > 0:
                return read_single(path)
        raise FileNotFoundError(f"No CONTCAR/OUTCAR found in {directory}")

    def is_converged(self, directory):
        try:
            return bool(Vasp(directory=directory).read_convergence())
        except Exception:  # noqa: BLE001 - missing/incomplete output
            return False

    def result_writer(self):
        script = code_root() / "common" / "write_output.py"
        return str(script), None

    def charges_writer(self):
        script = code_root() / "common" / "write_charges.py"
        return str(script), None

    def magmoms_writer(self):
        script = code_root() / "common" / "write_magmoms_input.py"
        return str(script), None
