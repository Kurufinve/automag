"""
Mixed-calculator pipeline: Quantum ESPRESSO for the magnetic search, VASP for
the Monte Carlo analysis and MAE.

    coll_qe       -> 2_coll with Quantum ESPRESSO
    mc_qe         -> couplings from the QE energies
    mc_vasp       -> the same analysis labelled as a VASP run, reading the QE
                     bundle through "coll_calculator": "espresso"
    mae_vasp      -> MAE with VASP from the QE ground state

Note: QE's ``starting_magnetization`` is defined per atomic species, so
antiferromagnetic order within one species is only fully representable if the
structure distinguishes the sublattices.

The launchers in this directory cover every way of running the pipeline:

    ./run_pipeline.sh                  # one command: run everything and wait
    ./run_pipeline_once.sh             # one pass, for manual restarts
    ./run_pipeline_once_loop.sh        # the same passes, in a loop
    ./run_pipeline_compiled.sh         # the same, with a compiled automag
    ./run_pipeline_slurm.sh            # the supervisor as a SLURM job
    ./run_pipeline_slurm_compiled.sh   # ... with the compiled binary
"""

import os

QE_INPUT = {
    "control": {
        "calculation": "scf",
        "prefix": "automag",
        "outdir": "./tmp",
        # Put the UPF files here or use an absolute directory.
        "pseudo_dir": "/path/to/qe/pseudos",
    },
    "system": {
        "ecutwfc": 60,
        "ecutrho": 480,
        "occupations": "smearing",
        "smearing": "gaussian",
        "degauss": 0.01,
        "nspin": 2,
    },
    "electrons": {
        "conv_thr": 1e-8,
        "mixing_beta": 0.3,
    },
}

QE_PARAMS = {
    "input_data": QE_INPUT,
    "pseudopotentials": {
        "Fe": "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF",
        "O": "O.pbe-n-kjpaw_psl.1.0.0.UPF",
    },
    "kpts": [4, 4, 4],
}

QE_SUBMIT = {
    "use_fireworks": False,
    "parallel_over_configurations": True,
    "jobheader": """#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=1-00:00:00
#SBATCH --partition=YOUR_PARTITION
#SBATCH --output=qe-%j.out
#SBATCH --error=qe-%j.error
""",
    # ASE writes espresso.pwi; QE writes espresso.pwo.
    "calculator_command": (
        "module load quantum-espresso/7.2; "
        "mpirun pw.x -in espresso.pwi > espresso.pwo"
    ),
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
            "id": "coll_qe",
            "stage": "coll",
            "calculator": "espresso",
            "params": QE_PARAMS,
            "supercell_size": 1,
            "spin_values": SPIN_VALUES,
            "magnetic_atoms": MAGNETIC_ATOMS,
            **QE_SUBMIT,
        },
        {
            "id": "mc_qe",
            "stage": "monte_carlo",
            "calculator": "espresso",
            "on_missing": "skip",
            **MC_SETTINGS,
            **QE_SUBMIT,
        },
        {
            "id": "mc_vasp",
            "stage": "monte_carlo",
            "calculator": "vasp",
            "coll_calculator": "espresso",
            "on_missing": "skip",
            "description": "same QE energies, labelled as a VASP run",
            **MC_SETTINGS,
            **VASP_SUBMIT,
        },
        {
            "id": "mae_vasp",
            "stage": "mae",
            "calculator": "vasp",
            "coll_calculator": "espresso",
            "on_missing": "skip",
            "description": "MAE with VASP from the QE ground state",
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
