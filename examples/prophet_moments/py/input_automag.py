"""
Prophet with magnetic-moment prediction (two steps per calculation).

Prophet-V1-Mag comes in two pieces:

* ``prophet-v1-mag-head.pt`` - the moment head.  Given the structure and an
  *initial* spin configuration it predicts the magnetic moment of every atom;
* ``prophet-v1-mag.pt`` - the energy model.  Given the structure and the
  moments it returns the energy.

So one calculation of this example is

    1. predict the moments from the initial configuration (``MagmomPredictor``);
    2. evaluate the energy for the predicted moments (``SpinCalculator``).

The stages themselves do not change: ``2_coll`` enumerates the spin
configurations and writes the initial moments, and ``prophet_spin`` evaluates
each of them with the two-step runner ``common/run_prophet.py``.  The predicted
moments are stored in ``results.json`` next to the energy, so the post-processing
of ``2_coll`` reads them as the converged moments of the configuration.

The checkpoints are bundled with automag (``models/prophet-v1-mag.pt`` and
``models/prophet-v1-mag-head.pt``), so no absolute paths are needed; the
environment variables ``PROPHET_MODEL``, ``PROPHET_HEAD`` and
``PROPHET_DEVICE`` override them.  Inference runs on one node: the CPU
device is the portable choice, ``"cuda"`` needs a GPU build of torch.

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
    # The moment head is only used when predict_moments is True.
    "head_path": os.environ.get("PROPHET_HEAD"),
    "device": os.environ.get("PROPHET_DEVICE", "cpu"),
    "predict_moments": True,
    # How the moments are seeded from the initial configuration: "lut" uses the
    # per-element table bundled with the head, "zero" starts from scratch.
    "seed_mode": "lut",
    "magnetic_threshold": 0.5,
}

PROPHET_SUBMIT = {
    "use_fireworks": False,
    "parallel_over_configurations": True,
    # No #SBATCH/PBS directives: the calculations are plain processes on this
    # machine, so the pipeline starts them directly (no external code needed).
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
            "description": "magnetic configurations, moments predicted by prophet",
            "params": PROPHET_PARAMS,
            "supercell_size": 1,
            "spin_values": SPIN_VALUES,
            "magnetic_atoms": MAGNETIC_ATOMS,
            **PROPHET_SUBMIT,
        },
    ],
}
