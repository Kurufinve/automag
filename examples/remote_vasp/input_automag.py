"""
Remote (SLURM cluster) execution of the ``coll`` stage with VASP.

The ``2_coll`` calculations of this example do **not** run on the machine that
starts the pipeline: the instance sets ``"host": "cluster"``, and the pipeline
copies/streams everything to the host described in ``hosts`` below, runs the
enumeration and the DFT single points there (through its scheduler) and reads
the results back from the shared filesystem.

What has to be true for this to work:

* the project directory (``workdir``) is visible on the remote machine at the
  same path - a shared/cluster filesystem, or a path that exists on both sides;
* an automag checkout (or the compiled ``run_automag``) is reachable there;
  ``automag_home`` / ``automag_exe`` below tell the pipeline where;
* password-less ``ssh`` works (the transport runs ``ssh -o BatchMode=yes``);
  ``identity_file`` may point at a key instead of using the agent;
* the job header of the stage carries real ``#SBATCH`` directives, so the
  calculations are submitted with ``sbatch`` on that machine.

Start with the ``shell`` transport variant commented out at the bottom: it runs
the same code path through a wrapper on this machine, which is the quickest way
to check the configuration without a second host.
"""

VASP_PARAMS = {
    "xc": "PBE",
    "prec": "Accurate",
    "ncore": 4,
    "encut": 600,
    "ediff": 1e-6,
    "ismear": 1,
    "sigma": 0.1,
    "nelm": 200,
    "kpts": 20,
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

PIPELINE = {
    "workdir": ".",
    "structure": "structure.vasp",
    "struct_suffix": "",

    "wait_timeout": 30 * 24 * 3600,
    "wait_interval": 60,
    "wait_stable_seconds": 600,

    # ----------------------------------------------------------------- hosts
    "hosts": {
        # this machine - the pipeline itself and the analysis run here
        "local": {"transport": "local"},
        # the cluster that runs the 2_coll calculations
        "cluster": {
            "transport": "ssh",
            "host": "login.cluster.example",      # <-- EDIT
            "user": "myuser",                     # <-- EDIT
            # "identity_file": "~/.ssh/id_ed25519",
            # "port": 22,
            "workdir": "/shared/myproject",       # project dir on that machine
            "automag_home": "/shared/automag",    # automag checkout there
            # "automag_exe": "/shared/bin/run_automag",   # ... or the binary
            "python": "python3",
            "scheduler": "slurm",                 # slurm | pbs | none
            # shell lines executed before every command on that host, e.g.
            # "environment": ["module load vasp/6.4.3"],
        },
        # A second variant that needs no second machine: the same commands are
        # run on this host through a wrapper.  Handy to validate the setup:
        # "cluster": {"transport": "shell", "wrapper": "sh -c",
        #             "workdir": ".", "scheduler": "none"},
    },

    "execution": {"mode": "process", "max_parallel": 4, "max_attempts": 3},

    # ------------------------------------------------------------------ run
    # Only the magnetic search: everything else would need its own DFT run
    # elsewhere.  Add more instances when you need them.
    "run": [
        {
            "id": "coll_vasp_remote",
            "stage": "coll",
            "calculator": "vasp",
            "host": "cluster",              # <-- the stage runs on the cluster
            "description": "magnetic search with VASP on the remote cluster",
            "supercell_size": 1,
            "spin_values": {"Fe": [4]},
            "magnetic_atoms": ["Fe"],
            "params": VASP_PARAMS,
        },
    ],

    "submit": {
        "calculator": "vasp",
        "use_fireworks": False,
        "parallel_over_configurations": True,
        "jobheader": """#!/bin/bash
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
    },
}
