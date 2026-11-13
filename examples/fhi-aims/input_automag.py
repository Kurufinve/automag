"""
FHI-aims example: the calculator-independent stages.

    ./run_pipeline.sh

``coll`` evaluates the magnetic configurations, ``monte_carlo`` fits the
Heisenberg couplings and writes the VAMPIRE input.  Edit ``species_dir``, the
aims command and the SLURM header below.

FHI-aims writes the initial moment of every atom to ``geometry.in``, so
antiferromagnetic configurations are reproduced exactly - unlike QE, whose
``starting_magnetization`` is defined per species.


The launchers in this directory cover every way of running the pipeline:

    ./run_pipeline.sh                  # one command: run everything and wait
    ./run_pipeline_once.sh             # one pass, for manual restarts
    ./run_pipeline_once_loop.sh        # the same passes, in a loop
    ./run_pipeline_compiled.sh         # the same, with a compiled automag
    ./run_pipeline_slurm.sh            # the supervisor as a SLURM job
    ./run_pipeline_slurm_compiled.sh   # ... with the compiled binary

"""

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

SPIN_VALUES = {"Fe": [4]}
MAGNETIC_ATOMS = ["Fe"]

COLL_SETTINGS = {
    "supercell_size": 2,
    "spin_values": SPIN_VALUES,
    "magnetic_atoms": MAGNETIC_ATOMS,
    "params": AIMS_PARAMS,
}

MC_SETTINGS = {
    "configuration": "lowest",
    "cutoff_radius": 4.1,
    "control_group_size": 0.4,
    "append_coupling_constants": True,
    "magnetic_atoms": MAGNETIC_ATOMS,
    "mc_backend": "vampire",
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
            "id": "coll",
            "stage": "coll",
            "calculator": "fhi-aims",
            "description": "collinear magnetic configurations",
            **COLL_SETTINGS,
            **AIMS_SUBMIT,
        },
        {
            "id": "monte_carlo",
            "stage": "monte_carlo",
            "calculator": "fhi-aims",
            "coll_instance": "coll",
            "description": "Heisenberg couplings and VAMPIRE input",
            **MC_SETTINGS,
            **AIMS_SUBMIT,
        },
    ],
}
