"""
Prophet for the magnetocrystalline anisotropy (``4_mae``).

The chain is the usual one - ``2_coll`` finds the ground collinear state,
``4_mae`` rotates the moments away from it - but every energy comes from the
Prophet-V1-Mag model instead of a DFT code:

    coll  -> 2_coll with prophet_spin (energy of the collinear configurations)
    mae   -> 4_mae with prophet_spin (spin rotations, MAE(theta, phi), MAE curve)

``4_mae`` passes the rotated moments as three numbers per atom, which the
interface writes to ``magmoms.json`` unchanged.

Limitation worth knowing: the released Prophet-V1-Mag weights are isotropic in
the moment direction - the energy of ``M || z`` and ``M || x`` differs by less
than 1e-6 eV - so the resulting MAE is numerically zero and ``MAE_curve`` is a
flat line.  The example shows the whole workflow and is the template for a
direction-sensitive model or a higher-precision run; with VASP (see
``examples/vasp``) the same stages give the physical anisotropy.

The checkpoints are bundled with automag (``models/``), so no absolute paths
are needed; ``PROPHET_MODEL``, ``PROPHET_HEAD`` and ``PROPHET_DEVICE``
override them.

The launchers in this directory cover every way of running the pipeline:

    ./run_pipeline.sh                  # one command: run everything and wait
    ./run_pipeline_once.sh             # one pass, for manual restarts
    ./run_pipeline_once_loop.sh        # the same passes, in a loop
    ./run_pipeline_compiled.sh         # the same, with a compiled automag
    ./run_pipeline_slurm.sh            # the supervisor as a SLURM job
    ./run_pipeline_slurm_compiled.sh   # ... with the compiled binary
"""

import os

PROPHET_PARAMS = {
    # The checkpoints ship with automag (<repo>/models, also inside the
    # frozen bundle); PROPHET_MODEL / PROPHET_HEAD still override them.
    "model_path": os.environ.get("PROPHET_MODEL"),
    "head_path": os.environ.get("PROPHET_HEAD"),
    "device": os.environ.get("PROPHET_DEVICE", "cpu"),
    # The MAE stage sets the moment directions itself, so no moment prediction.
    "predict_moments": False,
}

PROPHET_SUBMIT = {
    "use_fireworks": False,
    "parallel_over_configurations": True,
    "jobheader": "#!/bin/bash",
}

SPIN_VALUES = {"Fe": [4]}
MAGNETIC_ATOMS = ["Fe"]


PIPELINE = {
    "workdir": ".",
    "structure": "structure.vasp",
    "struct_suffix": "",

    "wait_timeout": 30 * 24 * 3600,
    "wait_interval": 30,
    "wait_stable_seconds": 60,

    "run": [
        {
            "id": "coll",
            "stage": "coll",
            "calculator": "prophet_spin",
            "description": "collinear magnetic configurations with prophet",
            "params": PROPHET_PARAMS,
            "supercell_size": 1,
            "spin_values": SPIN_VALUES,
            "magnetic_atoms": MAGNETIC_ATOMS,
            **PROPHET_SUBMIT,
        },
        {
            "id": "mae",
            "stage": "mae",
            "calculator": "prophet_spin",
            "description": "MAE from the prophet energies",
            # Read the ground state found by 2_coll.
            "configuration": "lowest",
            # A small angular grid keeps the example cheap; the model is
            # evaluated for every point, so cost ~ Nph * Nth * N_MAE.
            "Nph": 4,
            "Nth": 2,
            "N_MAE": 4,
            "standardize_cell": True,
            "use_primitive_cell": True,
            "spin_values": SPIN_VALUES,
            "params": PROPHET_PARAMS,
            **PROPHET_SUBMIT,
        },
    ],
}