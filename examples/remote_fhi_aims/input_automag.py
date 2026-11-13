"""
Remote (SLURM cluster) execution of the ``coll`` stage with FHI-aims.

The ``2_coll`` calculations of this example run on the machine described in
``hosts["cluster"]`` instead of the machine that starts the pipeline: the
instance sets ``"host": "cluster"``, so automag runs the enumeration and the
FHI-aims single points there - through ``ssh`` and, because the job header
carries ``#SBATCH`` directives, through that machine's scheduler.

Requirements (see also ``examples/remote_vasp``):

* the project directory (``workdir``) must be visible on the remote side at the
  same path (shared filesystem) - results are read back from there;
* an automag checkout or compiled binary must be reachable there
  (``automag_home`` / ``automag_exe``);
* password-less ``ssh`` (``BatchMode=yes``) and real ``#SBATCH`` directives in
  the job header.

FHI-aims writes the initial moment of every atom to ``geometry.in``, so the
generated antiferromagnetic configurations are reproduced exactly.
"""

AIMS_PARAMS = {
    "species_dir": "/path/to/aims/species_defaults/defaults_2020/light",
    "xc": "pbe",
    "charge": 0,
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
            # "environment": ["module load aims/240507"],
        },
        # A wrapper variant that needs no second machine (same code path):
        # "cluster": {"transport": "shell", "wrapper": "sh -c",
        #             "workdir": ".", "scheduler": "none"},
    },

    "execution": {"mode": "process", "max_parallel": 4, "max_attempts": 3},

    # ------------------------------------------------------------------ run
    "run": [
        {
            "id": "coll_aims_remote",
            "stage": "coll",
            "calculator": "fhi-aims",
            "host": "cluster",              # <-- the stage runs on the cluster
            "description": "magnetic search with FHI-aims on the remote cluster",
            "supercell_size": 1,
            "spin_values": {"Fe": [4]},
            "magnetic_atoms": ["Fe"],
            "params": AIMS_PARAMS,
        },
    ],

    "submit": {
        "calculator": "fhi-aims",
        "use_fireworks": False,
        "parallel_over_configurations": True,
        "jobheader": """#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=1-00:00:00
#SBATCH --partition=YOUR_PARTITION
#SBATCH --output=aims-%j.out
#SBATCH --error=aims-%j.error
""",
        "calculator_command": "module load aims/240507; mpirun aims.x",
        "environment_activate": (
            "source /path/to/miniconda3/etc/profile.d/conda.sh && "
            "conda activate automag"
        ),
        "environment_deactivate": "conda deactivate",
    },
}
