"""
automag.pipeline.stages
=======================

Declarative description of the pipeline stages.

A stage is a directory of the repository plus a mapping from named *phases*
(``submit``, ``postprocess``, ...) to the scripts that implement them.

For the ``--wait`` mode the same table also describes the dependency graph:

* ``full_steps``  - phases executed by a complete run (``--wait``);
* ``depends``     - phases that must finish before a phase can start;
* ``wait_after``  - phases that submit calculations and therefore have to be
  followed by a poll until the calculations finish.

Adding a stage or a phase is a matter of editing this table - the runner is
generic.

.. codeauthor:: Automag contributors
"""
from __future__ import annotations


from dataclasses import dataclass, field


@dataclass(frozen=True)
class Stage:
    name: str
    directory: str
    steps: dict
    default_steps: tuple
    work_dir: str = ""
    full_steps: tuple = ()
    depends: dict = field(default_factory=dict)
    wait_after: tuple = ()
    required: tuple = ()
    description: str = ""

    def phases(self, selection=None, full=False):
        """Return the ordered list of phase names to execute."""
        if selection is None:
            selection = self.full_steps if full else self.default_steps
        if isinstance(selection, str):
            selection = [item.strip() for item in selection.split(",") if item.strip()]
        phases = []
        for item in selection:
            if item in self.steps:
                phase = item
            else:
                # Also accept a script name and map it to its phase.
                matches = [name for name, scripts in self.steps.items()
                           if item in scripts]
                if not matches:
                    known = ", ".join(self.steps)
                    raise ValueError(
                        f"Unknown phase '{item}' for stage '{self.name}'. "
                        f"Known phases: {known}."
                    )
                phase = matches[0]
            if phase not in phases:
                phases.append(phase)
        return phases

    def resolve_steps(self, selection=None, full=False):
        """Return the ordered list of scripts to run."""
        scripts = []
        for phase in self.phases(selection, full=full):
            for script in self.steps[phase]:
                if script not in scripts:
                    scripts.append(script)
        return scripts


STAGE_ORDER = (
    "conv_tests",
    "lin_response",
    "coll",
    "relax",
    "monte_carlo",
    "mae",
)


STAGES = {
    "conv_tests": Stage(
        name="conv_tests",
        directory="0_conv_tests",
        work_dir="0_conv",
        steps={
            "submit": ("1_submit.py",),
            "postprocess": ("2_plot_results.py",),
        },
        default_steps=("submit",),
        full_steps=("submit", "postprocess"),
        depends={"postprocess": (("conv_tests", "submit"),)},
        wait_after=("submit",),
        required=("mode", "params"),
        description="ENCUT / k-mesh convergence tests.",
    ),
    "lin_response": Stage(
        name="lin_response",
        directory="1_lin_response",
        work_dir="1_lin_response",
        steps={
            "submit": ("1_submit.py",),
            "postprocess": ("2_plot_results.py",),
        },
        default_steps=("submit",),
        full_steps=("submit", "postprocess"),
        depends={"postprocess": (("lin_response", "submit"),)},
        wait_after=("submit",),
        required=("dummy_atom", "dummy_position", "perturbations", "params"),
        description="Hubbard U from linear response.",
    ),
    "coll": Stage(
        name="coll",
        directory="2_coll",
        work_dir="2_coll",
        steps={
            "submit": ("1_submit.py",),
            "postprocess": ("2_plot_results.py",),
        },
        default_steps=("submit",),
        full_steps=("submit", "postprocess"),
        depends={"postprocess": (("coll", "submit"),)},
        wait_after=("submit",),
        required=("supercell_size", "params"),
        description="Search for the ground collinear magnetic state.",
    ),
    "relax": Stage(
        name="relax",
        directory="2_relax",
        work_dir="2_relax",
        steps={
            "submit": ("1_submit.py",),
            "postprocess": ("2_collect.py",),
        },
        default_steps=("submit",),
        full_steps=("submit", "postprocess"),
        depends={
            "submit": (("coll", "postprocess"),),
            "postprocess": (("relax", "submit"),),
        },
        wait_after=("submit",),
        required=("params",),
        description="Re-relax the most stable magnetic configuration from 2_coll.",
    ),
    "monte_carlo": Stage(
        name="monte_carlo",
        directory="3_monte_carlo",
        work_dir="3_monte_carlo",
        steps={
            "submit": (
                "1_coupling_constants.py",
                "2_write_vampire_ucf.py",
            ),
            "run": ("5_run_mc.py",),
            "espins": ("4_write_espins_mcin.py",),
            "postprocess": ("3_plot_results.py",),
        },
        default_steps=("submit",),
        full_steps=("submit", "run", "postprocess"),
        depends={
            "submit": (("coll", "postprocess"),),
            "run": (("monte_carlo", "submit"),),
            "espins": (("monte_carlo", "submit"),),
            "postprocess": (("monte_carlo", "run"),),
        },
        wait_after=(),
        required=("configuration", "cutoff_radius", "control_group_size"),
        description="Heisenberg coupling constants and Monte Carlo inputs.",
    ),
    "mae": Stage(
        name="mae",
        directory="4_mae",
        work_dir="4_mae",
        steps={
            "submit": ("1_submit.py",),
            "postprocess": ("2_plot_results.py", "3_submit.py"),
            "axes": ("2_plot_results.py",),
            "curve": ("3_submit.py",),
            "curve_postprocess": ("4_plot_results.py",),
        },
        default_steps=("submit",),
        full_steps=("submit", "axes", "curve", "curve_postprocess"),
        depends={
            "submit": (("coll", "postprocess"),),
            "axes": (("mae", "submit"),),
            "curve": (("mae", "axes"),),
            "curve_postprocess": (("mae", "curve"),),
            "postprocess": (("mae", "submit"),),
        },
        wait_after=("submit", "curve"),
        required=("configuration", "params", "Nph", "Nth", "N_MAE"),
        description="Magnetocrystalline anisotropy energy.",
    ),
}
