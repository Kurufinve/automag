"""
Mixed-calculator pipeline: FHI-aims for the magnetic search, VASP for MAE.

    coll_aims     -> 2_coll with FHI-aims
    mc_aims       -> couplings from the FHI-aims energies
    mae_vasp      -> MAE with VASP, reading the FHI-aims 2_coll bundle

``mae_vasp`` sets ``"coll_calculator": "fhi-aims"`` so that it picks the
geometry, the spin pattern and the converged magnetic moments from the
FHI-aims bundle while the MAE calculations themselves are run with VASP and
written under ``CalcFold/<formula>/vasp/...``.

FHI-aims writes the initial moment of every atom to ``geometry.in``, so the
generated antiferromagnetic configurations are reproduced exactly.

The launchers in this directory cover every way of running the pipeline:

    ./run_pipeline.sh                  # one command: run everything and wait
    ./run_pipeline_once.sh             # one pass, for manual restarts
    ./run_pipeline_once_loop.sh        # the same passes, in a loop
    ./run_pipeline_compiled.sh         # the same, with a compiled automag
    ./run_pipeline_slurm.sh            # the supervisor as a SLURM job
    ./run_pipeline_slurm_compiled.sh   # ... with the compiled binary
"""

import os

AIMS_PARAMS = {
    "xc": "pbe",
    "kpts": [4, 4, 4],
    # FHI-aims species defaults, e.g. .../species_defaults/defaults_2020/tight
    "species_dir": "/path/to/aims/species_defaults/defaults_2020/tight",
    "sc_accuracy_etot": 1e-6,
    "sc_accuracy_eev": 1e-3,
    "sc_accuracy_rho": 1e-6,
}

AIMS_SUBMIT = {
    "use_fireworks": False,
    "parallel_over_configurations": True,
    "jobheader": """#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=1-00:00:00
#SBATCH --partition=YOUR_PARTITION
#SBATCH --output=aims-%j.out
#SBATCH --error=aims-%j.error
""",
    # ASE writes geometry.in / control.in; FHI-aims writes aims.out.
    "calculator_command": "module load fhi-aims; mpirun aims.x",
    "environment_activate": (
        "source /path/to/miniconda3/etc/profile.d/conda.sh && "
        "conda activate automag"
    ),
    "environment_deactivate": "conda deactivate",
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
    # "vampire_run_command": "vampire",   # default: external/codes.json (./build.sh --external)
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
            "id": "coll_aims",
            "stage": "coll",
            "calculator": "fhi-aims",
            "params": AIMS_PARAMS,
            "supercell_size": 1,
            "spin_values": SPIN_VALUES,
            "magnetic_atoms": MAGNETIC_ATOMS,
            **AIMS_SUBMIT,
        },
        {
            "id": "mc_aims",
            "stage": "monte_carlo",
            "calculator": "fhi-aims",
            "on_missing": "skip",
            **MC_SETTINGS,
            **AIMS_SUBMIT,
        },
        {
            "id": "mae_vasp",
            "stage": "mae",
            "calculator": "vasp",
            "coll_calculator": "fhi-aims",
            "on_missing": "skip",
            "description": "MAE with VASP from the FHI-aims ground state",
            "params": VASP_PARAMS,
            "configuration": "lowest",
            "Nph": 8,
            "Nth": 4,
            "N_MAE": 8,
            "standardize_cell": True,
            "use_primitive_cell": True,
            "spin_values": SPIN_VALUES,
            **VASP_SUBMIT,
        },
    ],
}
