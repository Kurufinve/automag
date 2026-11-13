"""
VASP example: the full magnetic workflow on a cluster.

Copy this directory anywhere on the cluster, edit the SLURM header and
``calculator_command`` below (and ``AUTOMAG_HOME`` in ``run_pipeline.sh``),
then run

    ./run_pipeline.sh

The pipeline runs convergence tests, the search for the ground collinear state,
a re-relaxation of it, Heisenberg couplings with a VAMPIRE run, and the
magnetocrystalline anisotropy.  Settings are grouped in dictionaries above
``PIPELINE`` and referenced by the explicit instance list, so a stage can be
run twice with different parameters (see ``examples/prophet_then_vasp``).


The launchers in this directory cover every way of running the pipeline:

    ./run_pipeline.sh                  # one command: run everything and wait
    ./run_pipeline_once.sh             # one pass, for manual restarts
    ./run_pipeline_once_loop.sh        # the same passes, in a loop
    ./run_pipeline_compiled.sh         # the same, with a compiled automag
    ./run_pipeline_slurm.sh            # the supervisor as a SLURM job
    ./run_pipeline_slurm_compiled.sh   # ... with the compiled binary

"""

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

CONV_SETTINGS = {
    "mode": "kgrid",
    "params": {
        "xc": "PBE",
        "setups": "recommended",
        "prec": "Accurate",
        "ncore": 4,
        "encut": 600,
        "ediff": 1e-6,
        "ismear": 1,
        "lcharg": False,
        "lwave": False,
    },
    # Lists instead of range() so that the YAML/JSON variants of this file can
    # hold exactly the same settings.
    "encut_values": [500, 550, 600, 650],
    "sigma_values": [0.05, 0.10],
    "kpts_values": [10, 20, 30],
}

COLL_SETTINGS = {
    # 2 (or more) is what makes the monte_carlo stage possible: a 10-atom cell
    # with 4 magnetic atoms gives 5 configurations = 2 equations for ~6
    # exchange couplings, which cannot be fitted.
    "supercell_size": 2,
    "spin_values": SPIN_VALUES,
    "magnetic_atoms": MAGNETIC_ATOMS,
    "params": VASP_PARAMS,
}

RELAX_SETTINGS = {
    "spin_values": SPIN_VALUES,
    "params": VASP_PARAMS,
    "relax_params": {
        "ibrion": 2,
        "isif": 3,
        "nsw": 100,
        "ediffg": -0.01,
    },
    "update_structure": False,
}

MC_SETTINGS = {
    # "lowest" = the most stable configuration found by 2_coll
    "configuration": "lowest",
    "cutoff_radius": 4.1,
    "control_group_size": 0.4,
    "append_coupling_constants": True,
    "magnetic_atoms": MAGNETIC_ATOMS,
    "mc_backend": "vampire",
    # "vampire_run_command": "vampire",   # default: external/codes.json
}

MAE_SETTINGS = {
    "configuration": "lowest",
    "Nph": 8,
    "Nth": 4,
    "N_MAE": 8,
    "standardize_cell": True,
    "use_primitive_cell": True,
    "spin_values": SPIN_VALUES,
    "params": VASP_PARAMS,
}


PIPELINE = {
    "workdir": ".",
    "structure": "structure.vasp",
    "struct_suffix": "",

    "wait_timeout": 30 * 24 * 3600,
    "wait_interval": 60,
    "wait_stable_seconds": 600,
    "wait_for_monte_carlo_output": False,

    "run": [
        {
            "id": "conv_tests",
            "stage": "conv_tests",
            "calculator": "vasp",
            "description": "ENCUT / sigma / k-point convergence",
            **CONV_SETTINGS,
            **VASP_SUBMIT,
        },
        {
            "id": "coll",
            "stage": "coll",
            "calculator": "vasp",
            "description": "collinear magnetic configurations",
            **COLL_SETTINGS,
            **VASP_SUBMIT,
        },
        {
            "id": "relax",
            "stage": "relax",
            "calculator": "vasp",
            "coll_instance": "coll",
            "description": "relaxation of the ground state",
            **RELAX_SETTINGS,
            **VASP_SUBMIT,
        },
        {
            "id": "monte_carlo",
            "stage": "monte_carlo",
            "calculator": "vasp",
            "coll_instance": "coll",
            "description": "Heisenberg couplings and VAMPIRE input",
            **MC_SETTINGS,
            **VASP_SUBMIT,
        },
        {
            "id": "mae",
            "stage": "mae",
            "calculator": "vasp",
            "coll_instance": "coll",
            "description": "magnetocrystalline anisotropy",
            **MAE_SETTINGS,
            **VASP_SUBMIT,
        },
    ],
}
