"""
automag.calculators.base
========================

The small interface every external code has to implement in order to be
usable by the pipeline.

The interface deliberately mirrors what the pipeline actually needs:

* write the input files of a calculation for a given ``Atoms`` object;
* return the command that runs the code;
* read back energy / magnetic moments / structure / convergence status.

Adding a new code means writing one class with these methods and registering
it (see :mod:`automag.calculators`).  No existing pipeline or stage code has
to be touched.

.. codeauthor:: Automag contributors
"""
from __future__ import annotations


import os


class CalculatorInterface:
    """Base class for a code interface.

    Subclasses only need to override the methods they support.  Unsupported
    operations raise :class:`NotImplementedError` with a clear message, so the
    pipeline can report *which* code lacks *which* feature.
    """

    #: Name used in ``calculator = "..."`` in the input file.
    name = "base"

    #: Capability flags used by the stages to fail early and explicitly.
    supports_ldau = False
    supports_noncollinear = False
    supports_relaxation = False
    supports_charges = False

    #: File whose presence (and non-zero size) means "this calculation has
    #: produced a result".  Used to decide whether an individual calculation
    #: of a batch has to be resubmitted.  Defaults to the ASE convention.
    output_marker = "output"

    def __init__(self, params=None, command=None, **kwargs):
        self.params = dict(params or {})
        self.command = command
        # Free-form extra settings from the input file (e.g. a dotted path
        # to an ASE calculator class for GenericASEInterface).
        self.options = dict(kwargs)

    # ------------------------------------------------------------------
    # Writing inputs
    # ------------------------------------------------------------------
    def write_input(self, atoms, directory, params=None):
        """Write the input files for a single calculation of ``atoms``.

        :param atoms: :class:`ase.Atoms` object.
        :param directory: directory where the files must be written.
        :param params: calculation parameters; merged over ``self.params``.
        """
        raise NotImplementedError(
            f"The '{self.name}' calculator interface does not implement "
            "write_input()."
        )

    def run_command(self):
        """Return the shell command that runs one calculation."""
        if self.command:
            return self.command
        raise NotImplementedError(
            f"The '{self.name}' calculator interface has no run command. "
            "Set 'calculator_command' in the input file."
        )

    # ------------------------------------------------------------------
    # Reading results
    # ------------------------------------------------------------------
    def read_energy(self, directory):
        raise NotImplementedError(
            f"The '{self.name}' calculator interface cannot read energies."
        )

    def read_magnetic_moments(self, directory):
        raise NotImplementedError(
            f"The '{self.name}' calculator interface cannot read magnetic "
            "moments."
        )

    def read_structure(self, directory):
        raise NotImplementedError(
            f"The '{self.name}' calculator interface cannot read structures."
        )

    def is_converged(self, directory):
        raise NotImplementedError(
            f"The '{self.name}' calculator interface cannot check convergence."
        )

    # ------------------------------------------------------------------
    # Optional helper scripts executed inside the job script
    # ------------------------------------------------------------------
    def result_markers(self, directory):
        """Marker files that prove one calculation in ``directory`` finished.

        The pipeline closes a calculation (and frees its ``max_parallel``
        slot) as soon as one of these files exists.  The default is the VASP
        ``OUTCAR``; calculators that do not write one must override this -
        otherwise their calculations never look finished.
        """
        return [{"kind": "outcar",
                 "path": os.path.join(str(directory), "singlepoint",
                                      "OUTCAR")}]

    def result_writer(self) -> tuple[str, str | None] | None:
        """Return ``(script_path, arguments)`` for the post-run writer or None."""
        return None

    def charges_writer(self) -> tuple[str, str | None] | None:
        """Return ``(script_path, arguments)`` for the LDA+U charges writer."""
        return None

    def magmoms_writer(self) -> tuple[str, str | None] | None:
        """Return ``(script_path, arguments)`` for the recalc MAGMOM writer."""
        return None

    def prepare_directory(self, directory):
        """Create ``directory`` (if needed) and return it as a string."""
        os.makedirs(directory, exist_ok=True)
        return str(directory)
