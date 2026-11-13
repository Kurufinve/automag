"""
automag.calculators.espresso
============================

Quantum ESPRESSO interface built on :class:`ase.calculators.espresso.Espresso`.

ASE writes ``espresso.pwi`` and reads ``espresso.pwo`` in the calculation
directory.  The ``input_data`` / ``pseudopotentials`` / ``kpts`` dictionaries
from the input file are passed to ASE unchanged.

Note on magnetic configurations: Quantum ESPRESSO's ``starting_magnetization``
is defined per atomic species, so antiferromagnetic order within one species is
only fully representable if the structure distinguishes the sublattices.  The
pipeline still writes the generated moments to the ASE ``Atoms`` object; ASE
derives the per-species values from them.

.. codeauthor:: Automag contributors
"""

import os

import numpy as np
from ase.io import read

from common.paths import code_root
from calculators.base import CalculatorInterface


class EspressoInterface(CalculatorInterface):
    name = "espresso"
    # Structural relaxation currently writes VASP-specific tags (IBRION/ISIF/
    # NSW/EDIFFG) from the ``relax_params`` setting, so only VASP may run the
    # ``relax`` stage.  Flip this flag once relax_params are calculator-aware.
    supports_relaxation = False

    INPUT = "espresso.pwi"
    OUTPUT = "espresso.pwo"
    output_marker = OUTPUT

    def write_input(self, atoms, directory, params=None):
        from ase.calculators.espresso import Espresso, EspressoProfile

        merged = dict(self.params)
        merged.update(params or {})
        for key in ("directory", "atoms", "profile", "command"):
            merged.pop(key, None)

        pseudo_dir = merged.pop("pseudo_dir", None)
        input_data = merged.get("input_data", {})
        if pseudo_dir is None:
            pseudo_dir = input_data.get("control", {}).get("pseudo_dir", ".")

        # ASE 3.29+ requires an explicit profile; the real executable command
        # is supplied through ``calculator_command`` in the job script.
        profile = EspressoProfile(command="pw.x", pseudo_dir=pseudo_dir)
        os.makedirs(directory, exist_ok=True)
        calc = Espresso(profile=profile, directory=str(directory), **merged)
        if hasattr(calc, "write_inputfiles"):
            calc.write_inputfiles(atoms, ["energy"])
        else:  # older ASE releases
            calc.write_input(atoms)
        return calc

    def run_command(self):
        if self.command:
            return self.command
        return "mpirun pw.x -in espresso.pwi > espresso.pwo"

    def _output_path(self, directory):
        path = os.path.join(directory, self.OUTPUT)
        if os.path.isfile(path):
            return path
        alternative = os.path.join(directory, "output")
        if os.path.isfile(alternative):
            return alternative
        raise FileNotFoundError(
            f"No Quantum ESPRESSO output found in {directory} "
            f"(expected {self.OUTPUT})."
        )

    def _read_atoms(self, directory):
        return read(self._output_path(directory), format="espresso-out")

    def read_energy(self, directory):
        return float(self._read_atoms(directory).get_potential_energy())

    def read_structure(self, directory):
        return self._read_atoms(directory)

    def read_magnetic_moments(self, directory):
        try:
            return np.asarray(self._read_atoms(directory).get_magnetic_moments())
        except Exception:  # noqa: BLE001 - not all outputs carry moments
            return None

    def is_converged(self, directory):
        try:
            with open(self._output_path(directory), "r", errors="ignore") as handle:
                return "JOB DONE" in handle.read()
        except OSError:
            return False

    def result_writer(self):
        script = code_root() / "common" / "write_output_ase.py"
        return str(script), None
