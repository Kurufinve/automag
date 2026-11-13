"""
automag.calculators.generic
===========================

A generic interface that turns *any* ASE calculator into a pipeline backend.

Two usage modes:

1. ``calculator = "ase:ase.calculators.espresso.Espresso"`` - the class is
   imported by dotted path.  Input writing and the run command work out of the
   box; reading results requires ``output_format`` or ``reader``.
2. ``calculator = "myname"`` after registering a subclass of
   :class:`GenericASEInterface` from a plugin module.

Reading results generically is impossible for arbitrary codes, so it is
opt-in through ``output_format`` (any format understood by ``ase.io.read``)
or ``reader`` (a callable or dotted path to one).  The writer executed inside
the job script can likewise be supplied through ``result_writer_path``.

.. codeauthor:: Automag contributors
"""

import importlib
import os
from typing import Any, cast

from ase import Atoms

from calculators.base import CalculatorInterface


def import_object(dotted_path):
    """Import ``package.module:object`` or ``package.module.object``."""
    if callable(dotted_path):
        return dotted_path
    if ":" in dotted_path:
        module_name, attribute = dotted_path.split(":", 1)
    else:
        module_name, _, attribute = dotted_path.rpartition(".")
    module = importlib.import_module(module_name)
    return getattr(module, attribute)


class GenericASEInterface(CalculatorInterface):
    name = "ase"

    def __init__(self, calculator_path=None, params=None, command=None,
                 output_format=None, output_file=None, reader=None,
                 result_writer_path=None, **kwargs):
        super().__init__(params=params, command=command, **kwargs)
        calculator_path = calculator_path or self.options.get("calculator_path")
        if calculator_path is None:
            raise ValueError(
                "GenericASEInterface needs 'calculator_path', e.g. "
                "calculator = 'ase:ase.calculators.espresso.Espresso'."
            )
        self.calculator_path = calculator_path
        self.calculator_class: Any = import_object(calculator_path)
        self.output_format = output_format or self.options.get("output_format")
        self.output_file = output_file or self.options.get("output_file")
        self.reader = reader or self.options.get("reader")
        self.result_writer_path = result_writer_path or self.options.get(
            "result_writer_path"
        )

    def _build_calculator(self, directory):
        params = dict(self.params)
        params.pop("directory", None)
        params.pop("atoms", None)
        try:
            return self.calculator_class(directory=str(directory), **params)
        except TypeError:
            # Most ASE calculators use ``label`` instead of ``directory``.
            params.setdefault("label", os.path.join(str(directory), "calc"))
            return self.calculator_class(**params)

    def write_input(self, atoms, directory, params=None):
        os.makedirs(directory, exist_ok=True)
        merged = dict(self.params)
        merged.update(params or {})
        # Rebuild with the merged parameters so per-stage overrides are honoured.
        original = self.params
        try:
            self.params = merged
            calc = self._build_calculator(directory)
        finally:
            self.params = original
        atoms = atoms.copy()
        atoms.calc = calc
        if hasattr(calc, "write_inputfiles"):
            calc.write_inputfiles(atoms, ["energy"])
        elif hasattr(calc, "write_input"):
            calc.write_input(atoms)
        else:
            raise NotImplementedError(
                f"{self.calculator_path} has neither write_inputfiles nor "
                "write_input."
            )
        return calc

    def _read_atoms(self, directory) -> Atoms:
        if self.reader is not None:
            return cast(Atoms, import_object(self.reader)(directory))
        if self.output_format is None:
            raise NotImplementedError(
                "Generic ASE interface cannot read results without "
                "'output_format' or 'reader'. Set one in the input file."
            )
        filename = self.output_file or os.path.join(directory, "output")
        from ase.io import read

        result = read(filename, format=self.output_format)
        if isinstance(result, list):
            result = result[-1]
        return result

    def read_energy(self, directory):
        return float(self._read_atoms(directory).get_potential_energy())

    def read_structure(self, directory):
        return self._read_atoms(directory)

    def result_writer(self):
        if self.result_writer_path:
            return self.result_writer_path, None
        return None
