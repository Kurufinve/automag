"""
Remote (SLURM cluster) execution of the ``coll`` stage with Quantum ESPRESSO.

The ``2_coll`` calculations of this example run on the machine described in
``hosts["cluster"]`` instead of the machine that starts the pipeline: the
instance sets ``"host": "cluster"``, so automag runs the enumeration and the
``pw.x`` single points there - through ``ssh`` and, because the job header
carries ``#SBATCH`` directives, through that machine's scheduler.

Requirements (see also ``examples/remote_vasp``):

* the project directory (``workdir``) must be visible on the remote side at the
  same path (shared filesystem) - results are read back from there;
* an automag checkout or compiled binary must be reachable there
  (``automag_home`` / ``automag_exe``);
* password-less ``ssh`` (``BatchMode=yes``) and real ``#SBATCH`` directives in
  the job header.

Note: QE's ``starting_magnetization`` is defined per atomic species, so
antiferromagnetic order within one species is only fully representable if the
structure distinguishes the sublattices.
"""

import os

QE_INPUT = {
    "control": {
        "calculation": "scf",
        "prefix": "automag",
        "pseudo_dir": os.environ.get("ESPRESSO_PSEUDO", "/path/to/upf"),
        "tstress": False,
        "tprnfor": True,
    },
    "system": {
        "ecutwfc": 60,
        "ecutrho": 480,
        "occupations": "smearing",
        "smearing": "mv",
        "degauss": 0.02,
        "nspin": 2,
    },
    "electrons": {"conv_thr": 1e-8, "mixing_beta": 0.3},
}

PIPELINE = {
    "workdir": ".",
    "structure": "structure.vasp",
    "struct_suffix": "",

    "run": ["coll"],

    "wait_timeout": 30 * 24 * 3600,
    "wait_interval": 60,
    "wait_stable_seconds": 600,

    # ----------------------------------------------------------------- hosts
    "hosts": {
        "local": {"transport": "local"},
        "cluster": {
            "transport": "ssh",
            "host": "login.cluster.example",      # <-- EDIT
            "user": "myuser",                     # <-- EDIT
            # "identity_file": "~/.ssh/id_ed25519",
            "workdir": "/shared/myproject",       # project dir on that machine
            "automag_home": "/shared/automag",    # automag checkout there
            # "automag_exe": "/shared/bin/run_automag",   # ... or the binary
            "python": "python3",
            "scheduler": "slurm",
            # shell lines executed before every command on that host:
            # "environment": ["module load qe/7.3"],
        },
        # A wrapper variant that needs no second machine (same code path):
        # "cluster": {"transport": "shell", "wrapper": "sh -c",
        #             "workdir": ".", "scheduler": "none"},
    },

    "execution": {"mode": "process", "max_parallel": 4, "max_attempts": 3},

    # ------------------------------------------------------------------ run
    "run": [
        {
            "id": "coll_qe_remote",
            "stage": "coll",
            "calculator": "espresso",
            "host": "cluster",              # <-- the stage runs on the cluster
            "description": "magnetic search with Quantum ESPRESSO on the cluster",
            "supercell_size": 1,
            "spin_values": {"Fe": [4]},
            "magnetic_atoms": ["Fe"],
            "params": {
                "pseudopotentials": {
                    "Fe": "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF",
                    "O": "O.pbe-n-kjpaw_psl.1.0.0.UPF",
                },
                "input_data": QE_INPUT,
                "kpts": (6, 6, 6),
            },
        },
    ],

    "submit": {
        "calculator": "espresso",
        "use_fireworks": False,
        "parallel_over_configurations": True,
        "jobheader": """#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=1-00:00:00
#SBATCH --partition=YOUR_PARTITION
#SBATCH --output=qe-%j.out
#SBATCH --error=qe-%j.error
""",
        "calculator_command": "module load qe/7.3; mpirun pw.x -in espresso.pwi",
        "environment_activate": (
            "source /path/to/miniconda3/etc/profile.d/conda.sh && "
            "conda activate automag"
        ),
        "environment_deactivate": "conda deactivate",
    },
}
