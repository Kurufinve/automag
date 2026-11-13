"""
Complex scenario: Prophet on this machine and VASP on a cluster, for the
original cell and for the relaxed one.

    prophet, local     coll -> monte_carlo -> mae        (structure.vasp)
    VASP, cluster      coll -> monte_carlo -> mae -> relax_vasp
    VASP, cluster      coll -> monte_carlo -> mae        (structure_vasp_relaxed.vasp)
    prophet, local     coll -> monte_carlo -> mae        (structure_vasp_relaxed.vasp)

So both codes are evaluated for the original and for the relaxed cell, and the
two halves are independent: their results live under different ``struct_suffix``
values, so nothing overwrites anything:

    prophet, original   structure.vasp                 CalcFold/Fe4O6/prophet_spin/...
    VASP, original      structure_vasp.vasp            CalcFold/Fe4O6_vasp/vasp/...
    VASP, relaxed       structure_vasp_relaxed.vasp    CalcFold/Fe4O6_vasp_relaxed/vasp/...
    prophet, relaxed    structure_vasp_relaxed.vasp    CalcFold/Fe4O6_relaxed/prophet_spin/...

The relaxation is done by VASP: ``Prophet-V1-Mag`` provides forces, but
``2_relax`` starts a DFT-style relaxation (``ibrion``/``isif``/``nsw``), which
the model interface does not implement.  Prophet therefore re-uses the VASP
relaxed cell instead of relaxing it itself.

``relax_vasp`` uses ``"update_structure": True``, so ``2_relax`` copies the
relaxed geometry to ``geometries/structure_vasp_relaxed.vasp``; the instances
below that read the relaxed cell name that file through ``"structure"`` and
depend on ``relax_vasp``.  Every ``monte_carlo``/``mae`` instance names its
``coll_instance``, because each calculator has two coll instances.

The VASP instances set ``"host": "cluster"``: their step scripts and
calculations are started over ssh, so the project directory must be visible on
the cluster under the *same absolute path* and ``automag_home`` must point at
the automag checkout there (or ``automag_exe`` at a compiled binary).

Edit the ``hosts`` block and the SLURM header before running.  The layout below
is the one of ``examples/mixed_fhi_aims_vasp``: parameter dictionaries at the
top, an explicit instance list in the order of ``examples/vasp``.

The launchers in this directory cover every way of running the pipeline:

    ./run_pipeline.sh                  # one command: run everything and wait
    ./run_pipeline_once.sh             # one pass, for manual restarts
    ./run_pipeline_once_loop.sh        # the same passes, in a loop
    ./run_pipeline_compiled.sh         # the same, with a compiled automag
    ./run_pipeline_slurm.sh            # the supervisor as a SLURM job
    ./run_pipeline_slurm_compiled.sh   # ... with the compiled binary
"""

import os

# --------------------------------------------------------------- calculators
PROPHET_PARAMS = {
    # The checkpoints ship with automag (<repo>/models, also inside the
    # frozen bundle); PROPHET_MODEL / PROPHET_HEAD still override them.
    "model_path": os.environ.get("PROPHET_MODEL"),
    "head_path": os.environ.get("PROPHET_HEAD"),
    "device": os.environ.get("PROPHET_DEVICE", "cpu"),
    "predict_moments": True,
}

PROPHET_SUBMIT = {
    "use_fireworks": False,
    "parallel_over_configurations": True,
    "jobheader": "#!/bin/bash",          # no batch directives -> local process
    "host": "local",
}

VASP_PARAMS = {
    "xc": "PBE",
    "setups": "recommended",
    "prec": "Accurate",
    "ncore": 4,
    "encut": 600,
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
    "ldaul": [2, -1],
    "ldauu": [4.0, 0],
    "ldauj": [0, 0],
    "ldauprint": 2,
    "lreal": "Auto",
}

VASP_SUBMIT = {
    "use_fireworks": False,
    "parallel_over_configurations": True,
    "host": "cluster",
    "jobheader": """#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=1-00:00:00
#SBATCH --partition=YOUR_PARTITION
#SBATCH --output=vasp-%j.out
#SBATCH --error=vasp-%j.error
""",
    "calculator_command": "module load vasp/6.4.3; mpirun vasp_std",
    "environment_activate": (
        "source /path/to/miniconda3/etc/profile.d/conda.sh && "
        "conda activate automag"
    ),
    "environment_deactivate": "conda deactivate",
}

SPIN_VALUES = {"Fe": [4]}
MAGNETIC_ATOMS = ["Fe"]

MC_SETTINGS = {
    "configuration": "lowest",
    "cutoff_radius": 4.1,
    "control_group_size": 0.4,
    "append_coupling_constants": True,
    "magnetic_atoms": MAGNETIC_ATOMS,
    "mc_backend": "vampire",
}

MAE_SETTINGS = {
    "configuration": "lowest",
    # Small angular grid: every point is one energy evaluation per half.
    "Nph": 4,
    "Nth": 2,
    "N_MAE": 4,
    "standardize_cell": True,
    "use_primitive_cell": True,
    "spin_values": SPIN_VALUES,
}


PIPELINE = {
    "workdir": ".",
    # Base structure of the prophet half; the VASP half uses a copy with a
    # different name so that its relaxed geometry gets its own file.
    "structure": "structure.vasp",
    "struct_suffix": "",

    "hosts": {
        "local": {"transport": "local"},
        "cluster": {
            "transport": "ssh",
            "host": "login.cluster.example",     # <-- EDIT
            "user": "myuser",                    # <-- EDIT
            # "identity_file": "~/.ssh/id_ed25519",
            # "port": 22,
            "workdir": "/shared/myproject",      # same path as on this machine
            "automag_home": "/shared/automag",   # automag checkout on the cluster
            # "automag_exe": "/shared/bin/run_automag",
            "python": "python3",
            "scheduler": "slurm",                # slurm | pbs | none
            # "environment": ["module load vasp/6.4.3"],
        },
    },

    "execution": {"mode": "process", "max_parallel": 4, "max_attempts": 3},

    "wait_timeout": 30 * 24 * 3600,
    "wait_interval": 60,
    "wait_stable_seconds": 600,

    "run": [
        # -------------------------------------------- prophet, local, original
        {
            "id": "coll",
            "stage": "coll",
            "calculator": "prophet_spin",
            "description": "magnetic configurations of the input structure",
            "params": PROPHET_PARAMS,
            "supercell_size": 1,
            "spin_values": SPIN_VALUES,
            "magnetic_atoms": MAGNETIC_ATOMS,
            "struct_suffix": "",
            **PROPHET_SUBMIT,
        },
        {
            "id": "mc",
            "stage": "monte_carlo",
            "calculator": "prophet_spin",
            "coll_instance": "coll",
            "struct_suffix": "",
            **MC_SETTINGS,
            **PROPHET_SUBMIT,
        },
        {
            "id": "mae",
            "stage": "mae",
            "calculator": "prophet_spin",
            "coll_instance": "coll",
            "struct_suffix": "",
            **MAE_SETTINGS,
            "params": PROPHET_PARAMS,
            **PROPHET_SUBMIT,
        },

        # ---------------------------------------------- VASP, cluster, original
        {
            "id": "coll_vasp",
            "stage": "coll",
            "calculator": "vasp",
            "description": "same magnetic search with VASP on the cluster",
            "structure": "structure_vasp.vasp",
            "struct_suffix": "_vasp",
            "params": VASP_PARAMS,
            "supercell_size": 1,
            "spin_values": SPIN_VALUES,
            "magnetic_atoms": MAGNETIC_ATOMS,
            **VASP_SUBMIT,
        },
        {
            "id": "mc_vasp",
            "stage": "monte_carlo",
            "calculator": "vasp",
            "coll_instance": "coll_vasp",
            "struct_suffix": "_vasp",
            **MC_SETTINGS,
            **VASP_SUBMIT,
        },
        {
            "id": "mae_vasp",
            "stage": "mae",
            "calculator": "vasp",
            "coll_instance": "coll_vasp",
            "struct_suffix": "_vasp",
            **MAE_SETTINGS,
            "params": VASP_PARAMS,
            **VASP_SUBMIT,
        },
        {
            "id": "relax_vasp",
            "stage": "relax",
            "calculator": "vasp",
            "coll_instance": "coll_vasp",
            "struct_suffix": "_vasp",
            "description": "relaxation of the VASP ground state",
            "spin_values": SPIN_VALUES,
            # Copy the relaxed geometry into geometries/ so that the two
            # relaxed passes below can read it.
            "update_structure": True,
            "params": VASP_PARAMS,
            "relax_params": {
                "ibrion": 2,
                "isif": 3,
                "nsw": 100,
                "ediffg": -0.01,
            },
            **VASP_SUBMIT,
        },

        # ------------------------------------------------ VASP, cluster, relaxed
        {
            "id": "coll_vasp_relaxed",
            "stage": "coll",
            "calculator": "vasp",
            "description": "magnetic configurations of the relaxed structure",
            "structure": "structure_vasp_relaxed.vasp",
            "depends_on": ["relax_vasp"],
            "struct_suffix": "_vasp_relaxed",
            "params": VASP_PARAMS,
            "supercell_size": 1,
            "spin_values": SPIN_VALUES,
            "magnetic_atoms": MAGNETIC_ATOMS,
            **VASP_SUBMIT,
        },
        {
            "id": "mc_vasp_relaxed",
            "stage": "monte_carlo",
            "calculator": "vasp",
            "coll_instance": "coll_vasp_relaxed",
            "struct_suffix": "_vasp_relaxed",
            **MC_SETTINGS,
            **VASP_SUBMIT,
        },
        {
            "id": "mae_vasp_relaxed",
            "stage": "mae",
            "calculator": "vasp",
            "coll_instance": "coll_vasp_relaxed",
            "struct_suffix": "_vasp_relaxed",
            **MAE_SETTINGS,
            "params": VASP_PARAMS,
            **VASP_SUBMIT,
        },

        # -------------------------------- prophet, local, on the relaxed cell
        {
            "id": "coll_relaxed",
            "stage": "coll",
            "calculator": "prophet_spin",
            "description": "prophet on the VASP-relaxed structure",
            "structure": "structure_vasp_relaxed.vasp",
            "depends_on": ["relax_vasp"],
            "struct_suffix": "_relaxed",
            "params": PROPHET_PARAMS,
            "supercell_size": 1,
            "spin_values": SPIN_VALUES,
            "magnetic_atoms": MAGNETIC_ATOMS,
            **PROPHET_SUBMIT,
        },
        {
            "id": "mc_relaxed",
            "stage": "monte_carlo",
            "calculator": "prophet_spin",
            "coll_instance": "coll_relaxed",
            "struct_suffix": "_relaxed",
            **MC_SETTINGS,
            **PROPHET_SUBMIT,
        },
        {
            "id": "mae_relaxed",
            "stage": "mae",
            "calculator": "prophet_spin",
            "coll_instance": "coll_relaxed",
            "struct_suffix": "_relaxed",
            **MAE_SETTINGS,
            "params": PROPHET_PARAMS,
            **PROPHET_SUBMIT,
        },
    ],
}