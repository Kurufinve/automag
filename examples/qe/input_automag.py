"""
Quantum ESPRESSO example: the calculator-independent stages.

    ./run_pipeline.sh

``coll`` evaluates the magnetic configurations, ``monte_carlo`` fits the
Heisenberg couplings and writes the VAMPIRE input.  Edit ``pseudo_dir``, the
UPF file names, the QE command and the SLURM header below.

Note: QE's ``starting_magnetization`` is defined per atomic species, so
antiferromagnetic order within one species is only fully representable if the
structure distinguishes the sublattices.  The FHI-aims example does not have
this restriction.


The launchers in this directory cover every way of running the pipeline:

    ./run_pipeline.sh                  # one command: run everything and wait
    ./run_pipeline_once.sh             # one pass, for manual restarts
    ./run_pipeline_once_loop.sh        # the same passes, in a loop
    ./run_pipeline_compiled.sh         # the same, with a compiled automag
    ./run_pipeline_slurm.sh            # the supervisor as a SLURM job
    ./run_pipeline_slurm_compiled.sh   # ... with the compiled binary

"""

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

SPIN_VALUES = {"Fe": [4]}
MAGNETIC_ATOMS = ["Fe"]

COLL_SETTINGS = {
    "supercell_size": 2,
    "spin_values": SPIN_VALUES,
    "magnetic_atoms": MAGNETIC_ATOMS,
    "params": QE_PARAMS,
}

MC_SETTINGS = {
    # The Heisenberg fit needs more independent configurations than couplings:
    # with supercell_size = 1 (5 configurations, 4 magnetic atoms) the system is
    # underdetermined.  "lowest" = the most stable configuration from 2_coll.
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
            "calculator": "espresso",
            "description": "collinear magnetic configurations",
            **COLL_SETTINGS,
            **QE_SUBMIT,
        },
        {
            "id": "monte_carlo",
            "stage": "monte_carlo",
            "calculator": "espresso",
            "coll_instance": "coll",
            "description": "Heisenberg couplings and VAMPIRE input",
            **MC_SETTINGS,
            **QE_SUBMIT,
        },
    ],
}
