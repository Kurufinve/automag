"""
Example controlling file for the unified automag pipeline.

**Copy this file** into your own project directory, adjust the placeholders
(structure, SLURM header, calculator command, environment) and run::

    cd <project dir>
    python /path/to/automag/run_automag.py --wait

This file can live in any directory.  The only mandatory setting is
``structure``; everything else has a sensible default.  Relative paths are
resolved against the **launch directory**, so the results appear next to this
file and the automag checkout is never modified.

Stages
------
``conv_tests``    ENCUT / k-mesh convergence tests        (0_conv_tests)
``lin_response``  Hubbard U from linear response          (1_lin_response)
``coll``          ground collinear magnetic state         (2_coll)
``relax``         re-relax the most stable configuration  (2_relax)
``monte_carlo``   Heisenberg couplings / MC inputs         (3_monte_carlo)
``mae``           magnetocrystalline anisotropy           (4_mae)

``run`` selects which stages are executed, in the given order.  Each stage is
configured by the dictionary with the same name.  Settings under ``submit``
(job header, calculator, environment) are shared by every stage.

A stage may also be repeated with different calculators: write ``run`` as a
list of instance dicts (``id``, ``stage``, ``calculator`` plus stage settings).
See ``examples/mixed_qe_vasp/`` for a full example, and ``examples/README.md``
for the launchers (one blocking command, manual passes, SLURM driver job).

Phases
------
Every stage is split into phases (``submit``, ``postprocess``, ...).  By
default only ``submit`` is executed, so one command launches all calculations.
After the cluster jobs have finished run::

    python run_automag.py --steps postprocess

Use ``--steps submit,postprocess`` to do both in one go, and
``python run_automag.py --list-stages`` to see all phases.

For a fully unattended run use ``--wait``: the runner submits a batch of
calculations, polls until they finish, runs the analysis, and continues with
the dependent stages::

    python run_automag.py --wait
"""

PIPELINE = {
    # Directory where all results are written.  Relative paths are resolved
    # against the launch directory (not this file).  Keep "." and run the
    # pipeline from your own project directory: the automag checkout must stay
    # code + examples only, so a workdir inside it is rejected unless it lives
    # under examples/.
    "workdir": ".",

    # Input structure.  May be absolute, relative to this file, or a name in
    # the automag/geometries folder.
    "structure": "geometries/Fe2O3-alpha_conventional.vasp",

    # Arbitrary suffix distinguishing structures with the same composition.
    "struct_suffix": "",

    # Stages to run, in order.  With --wait the runner respects the stage
    # dependencies automatically: relax/monte_carlo/mae wait for the 2_coll
    # post-processing instead of starting too early.
    #
    # To combine calculators or repeat a stage, use the instance form instead
    # of bare stage names (see examples/mixed_*):
    #
    #   "run": [
    #       {"id": "coll_qe", "stage": "coll", "calculator": "espresso",
    #        "params": QE_PARAMS, "supercell_size": 1},
    #       {"id": "mae_vasp", "stage": "mae", "calculator": "vasp",
    #        "coll_calculator": "espresso",       # read the QE 2_coll results
    #        "on_missing": "skip", "reuse": True,
    #        "configuration": "lowest", "params": VASP_PARAMS,
    #        "Nph": 8, "Nth": 4, "N_MAE": 8},
    #   ]
    "run": ["conv_tests", "coll", "relax", "monte_carlo", "mae"],

    # External calculator plugins (module paths registered with
    # calculators.register).  Any ASE calculator can also be used directly with
    # calculator = "ase:package.module.Calculator".
    "calculator_plugins": [],

    # Paths/commands of the optional external codes (enumlib for 2_coll,
    # VAMPIRE/ESpinS for monte_carlo).  Normally written for you by
    #   ./build.sh --external [--<code>-path P]   -> external/codes.json
    # Use this key only to override entries for one run, e.g.
    #   "external_codes": {"vampire": {"command": "/opt/vampire/vampire-serial",
    #                                "pre_command": "module load intel",
    #                                "post_command": "module purge"}},
    # "external_codes": None,

    # --wait settings (used only with `run_automag.py ... --wait`).
    # Maximum time to wait for one phase, polling interval, and how long a
    # result file may stop growing before the run continues with the
    # calculations that did finish.
    "wait_timeout": 30 * 24 * 3600,
    "wait_interval": 60,
    "wait_stable_seconds": 600,
    # VAMPIRE/ESpinS are external: set True to block until
    # <workdir>/3_monte_carlo/output appears; otherwise that phase is skipped.
    "wait_for_monte_carlo_output": False,

    # ------------------------------------------------------------------
    # Settings shared by all stages
    # ------------------------------------------------------------------
    "submit": {
        "calculator": "vasp",
        "use_fireworks": False,
        "parallel_over_configurations": True,
        "jobheader": """#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=1-00:00:00
#SBATCH --partition=lenovo
#SBATCH --output=vasp-%j.out
#SBATCH --error=vasp-%j.error
""",
        "calculator_command": (
            "export UCX_TLS=ud,sm,self; module load vasp/6.4.3; mpirun vasp_std"
        ),
        "environment_activate": (
            "source /home/dpoletaev/miniconda3/bin/activate automag"
        ),
        "environment_deactivate": "conda deactivate",
    },

    # ------------------------------------------------------------------
    # 0_conv_tests
    # ------------------------------------------------------------------
    "conv_tests": {
        "mode": "kgrid",  # 'kgrid' or 'encut'
        "params": {
            "xc": "PBE",
            "setups": "recommended",
            "prec": "Accurate",
            "ncore": 4,
            "encut": 830,
            "ediff": 1e-6,
            "ismear": 1,
            "lcharg": False,
            "lwave": False,
        },
        "encut_values": range(600, 900, 50),
        "sigma_values": [0.05, 0.10, 0.15, 0.20],
        "kpts_values": range(20, 60, 10),
        # "magnetic_atoms": ["Fe"],
        # "configuration": 6 * [4.0] + 6 * [-4.0] + 18 * [0.0],
    },

    # ------------------------------------------------------------------
    # 1_lin_response (not in 'run' above - add it when you need U)
    # ------------------------------------------------------------------
    "lin_response": {
        "dummy_atom": "Zn",
        "dummy_position": 0,
        "perturbations": [-0.08, -0.05, -0.02, 0.02, 0.05, 0.08],
        "params": {
            "xc": "PBE",
            "setups": "recommended",
            "prec": "Accurate",
            "ncore": 2,
            "encut": 830,
            "ediff": 1e-6,
            "ismear": 0,
            "sigma": 0.2,
            "kpts": 20,
            "lmaxmix": 4,
            "nelm": 200,
        },
        # "magnetic_atoms": ["Fe"],
    },

    # ------------------------------------------------------------------
    # 2_coll
    # ------------------------------------------------------------------
    "coll": {
        "supercell_size": 1,
        "spin_values": {"Fe": [5]},
        "params": {
            "xc": "PBE",
            "setups": "recommended",
            "prec": "Accurate",
            "ncore": 4,
            "encut": 830,
            "ediff": 1e-6,
            "ismear": 1,
            "sigma": 0.1,
            "nelm": 200,
            "kpts": 20,
            "lmaxmix": 4,
            "lcharg": False,
            "lwave": False,
            "isym": 0,
            "ldau": True,
            "ldautype": 2,
            "ldaul": [2, -1, -1],
            "ldauu": [5.17, 0, 0],
            "ldauj": [0, 0, 0],
            "ldauprint": 2,
            "lreal": "Auto",
        },
        # "lower_cutoff": 1.7,
    },

    # ------------------------------------------------------------------
    # 2_relax: re-relax the most stable magnetic configuration from 2_coll
    # ------------------------------------------------------------------
    "relax": {
        # "configuration": "afm1",  # default: most stable from 2_coll
        "spin_values": {"Fe": [5]},
        "params": {
            "xc": "PBE",
            "setups": "recommended",
            "prec": "Accurate",
            "ncore": 4,
            "encut": 830,
            "ediff": 1e-6,
            "ismear": 1,
            "sigma": 0.1,
            "nelm": 200,
            "kpts": 20,
            "lmaxmix": 4,
            "lcharg": False,
            "lwave": False,
            "isym": 0,
            "ldau": True,
            "ldautype": 2,
            "ldaul": [2, -1, -1],
            "ldauu": [5.17, 0, 0],
            "ldauj": [0, 0, 0],
            "ldauprint": 2,
            "lreal": "Auto",
        },
        "relax_params": {
            "ibrion": 2,
            "isif": 3,
            "nsw": 100,
            "ediffg": -0.01,
        },
        # copy the relaxed geometry to <workdir>/geometries
        "update_structure": False,
    },

    # ------------------------------------------------------------------
    # 3_monte_carlo
    # ------------------------------------------------------------------
    "monte_carlo": {
        # "lowest" = the most stable configuration found by 2_coll
        "configuration": "lowest",
        "cutoff_radius": 4.1,
        "control_group_size": 0.4,
        # required so that 2_write_vampire_ucf.py can read the constants
        "append_coupling_constants": True,
        "magnetic_atoms": ["Fe"],

        # Backend run by the 'run' phase (5_run_mc.py): 'vampire' or 'espins'.
        "mc_backend": "vampire",
        # "vampire_run_command": "vampire",   # default: external/codes.json (./build.sh --external)
        # Defaults are the bundled files in 3_monte_carlo/vampire_input.
        # "vampire_input_file": "3_monte_carlo/vampire_input/input",
        # "vampire_material_file": "3_monte_carlo/vampire_input/vamp.mat",
        # "vampire_output_file": "output",
        # "mc_timeout": 7 * 24 * 3600,

        # ESpinS backend:
        # "espins_run_command": "module load intel2021/mkl/latest; /path/to/mc.x",
        # "espins_output_file": "result.dat",
        # "espins_output_converter": "python convert_espins.py",
    },

    # ------------------------------------------------------------------
    # 4_mae
    # ------------------------------------------------------------------
    "mae": {
        # "lowest" = the most stable configuration found by 2_coll
        "configuration": "lowest",
        "Nph": 20,
        "Nth": 10,
        "N_MAE": 20,
        "standardize_cell": True,
        "use_primitive_cell": True,
        "spin_values": {"Fe": [5]},
        "params": {
            "xc": "PBE",
            "setups": "recommended",
            "prec": "Accurate",
            "ncore": 4,
            "encut": 830,
            "ediff": 1e-6,
            "ismear": 1,
            "sigma": 0.1,
            "nelm": 200,
            "kpts": 20,
            "lmaxmix": 4,
            "lcharg": False,
            "lwave": False,
            "isym": 0,
            "ldau": True,
            "ldautype": 2,
            "ldaul": [2, -1, -1],
            "ldauu": [5.17, 0, 0],
            "ldauj": [0, 0, 0],
            "ldauprint": 2,
            "lreal": "Auto",
        },
    },
}
