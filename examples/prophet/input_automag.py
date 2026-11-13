"""
Prophet-V1-Mag instead of a DFT code: the magnetic search without DFT.

    ./run_pipeline.sh

``prophet_spin`` is an ASE calculator, so ``2_coll`` works without
VASP/QE/FHI-aims.  It needs ``torch`` + ``prophet_spin`` in the environment; the
checkpoint ships with automag (``models/prophet-v1-mag.pt``), and
``PROPHET_MODEL`` / ``PROPHET_DEVICE`` (``cpu`` by default) override it.  Inference on one core is slow, so give the job several cores.

Only ``coll`` is enabled: ``monte_carlo`` would run, but it degenerates for this
cell - 5 configurations give 2 independent equations, so the fit shrinks
``cutoff_radius`` to the nearest shell and returns a single coupling.  Raise
``supercell_size`` below (2-3) before enabling it.  ``relax``/``mae``/
``lin_response`` need capabilities the model interface does not provide.


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
    # False: the energy model is evaluated for the initial moments written by
    # 2_coll.  True would predict the moments first (see
    # examples/prophet_moments).
    "predict_moments": False,
}

PROPHET_SUBMIT = {
    "use_fireworks": False,
    "parallel_over_configurations": True,
    "jobheader": "#!/bin/bash",
}

SPIN_VALUES = {"Fe": [4]}
MAGNETIC_ATOMS = ["Fe"]

COLL_SETTINGS = {
    "supercell_size": 1,
    "spin_values": SPIN_VALUES,
    "magnetic_atoms": MAGNETIC_ATOMS,
    "params": PROPHET_PARAMS,
}


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
            **COLL_SETTINGS,
            **PROPHET_SUBMIT,
        },
    ],
}
