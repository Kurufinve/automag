"""
Magnetic configurations from pymatgen's MagneticStructureEnumerator.

Instead of the hand-written enumlib driver, ``2_coll`` can delegate the
enumeration to
:class:`pymatgen.analysis.magnetism.analyzer.MagneticStructureEnumerator`
(``enumeration_mode = "pymatgen"``).  It uses the same enumlib executable but
adds the pymatgen strategies:

``ferromagnetic``, ``antiferromagnetic``, ``antiferromagnetic_by_motif``,
``ferrimagnetic_by_motif``, ``ferrimagnetic_by_species``, ``nonmagnetic``.

Notice that ``spin_values`` and ``magnetic_atoms`` are **not** set below: when
they are missing, automag takes the moments from
``pymatgen.analysis.magnetism.analyzer.DEFAULT_MAGMOMS`` (Fe -> 5, O is not
tabulated and therefore non-magnetic).  Set them explicitly to override the
magnitudes, e.g. ``"spin_values": {"Fe": [4]}``.

The generated ``2_coll/trials_*`` folder has the same layout as usual
(``settingNNN.vasp`` + ``configurationsNNN.txt``), so the rest of the pipeline
is unchanged.

Prerequisite: enumlib (``./build.sh --external`` or ``--enumlib-path``), which
pymatgen looks up on ``PATH`` - automag adds the configured directories.

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

PIPELINE = {
    "workdir": ".",
    "structure": "structure.vasp",
    "struct_suffix": "",

    "run": ["coll", "relax", "monte_carlo"],

    "wait_timeout": 30 * 24 * 3600,
    "wait_interval": 60,
    "wait_stable_seconds": 600,
    "wait_for_monte_carlo_output": False,

    "submit": {
        "calculator": "vasp",
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
    },

    "coll": {
        # "enumlib" (default, original driver) or "pymatgen".
        "enumeration_mode": "pymatgen",
        # Strategies of MagneticStructureEnumerator.
        "enumeration_strategies": [
            "ferromagnetic",
            "antiferromagnetic",
            "antiferromagnetic_by_motif",
        ],
        # Optional cap on the number of generated configurations.
        "max_orderings": 20,
        "supercell_size": 1,
        # spin_values / magnetic_atoms are optional (DEFAULT_MAGMOMS).
        "params": VASP_PARAMS,
    },

    "relax": {
        "params": VASP_PARAMS,
        "relax_params": {"ibrion": 2, "isif": 3, "nsw": 100, "ediffg": -0.01},
        "update_structure": False,
    },

    "monte_carlo": {
        "configuration": "lowest",
        "cutoff_radius": 4.1,
        "control_group_size": 0.4,
        "append_coupling_constants": True,
        "mc_backend": "vampire",
    },
}
