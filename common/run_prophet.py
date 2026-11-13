"""
automag.common.run_prophet
==========================

Run the ``prophet_spin`` ASE calculator in a calculation directory.

The pipeline's ``ProphetSpinInterface`` writes ``atoms.json``, ``magmoms.json``
and ``params.json``; this script loads them, evaluates the model and writes
``results.json`` with the energy and the final moments (the initial moments are
kept if the model does not provide them).

With ``predict_moments`` in the parameters the evaluation is two-step: the
magnetic-moment head (``head_path``) predicts the moments from the initial spin
configuration, then the energy model is evaluated for those moments.  The
predicted list is stored as ``magmoms_predicted``.

.. codeauthor:: Automag contributors
"""

import json
import os
import sys

import numpy as np


def main():
    from ase.io import read
    from prophet_spin import SpinCalculator

    if not os.path.isfile("atoms.json"):
        raise FileNotFoundError("atoms.json not found in the calculation directory.")

    atoms = read("atoms.json")
    with open("magmoms.json", "rt") as handle:
        magmoms = json.load(handle)
    if magmoms:
        # A list of numbers = collinear moments, a list of triples = vectors.
        atoms.set_initial_magnetic_moments(np.asarray(magmoms, dtype=float))

    with open("params.json", "rt") as handle:
        params = json.load(handle)

    model_path = params.get("model_path")
    if not model_path:
        # The checkpoints ship with automag (models/ next to the sources or
        # inside the frozen bundle), so the explicit path is optional.
        from common.paths import PROPHET_MODEL_NAME, model_file

        model_path = str(model_file(PROPHET_MODEL_NAME))

    predicted = None
    if params.get("predict_moments"):
        # Step 1: the magnetic-moment head predicts the moments of every atom
        # from the initial spin configuration written by the stage.
        from prophet_spin import MagmomPredictor

        head_path = params.get("head_path")
        if not head_path:
            from common.paths import PROPHET_HEAD_NAME, model_file

            head_path = str(model_file(PROPHET_HEAD_NAME))
        if not head_path:
            raise ValueError(
                "predict_moments needs 'head_path' in the stage params, e.g. "
                "params = {'model_path': ..., 'head_path': ..., "
                "'predict_moments': True}."
            )
        predictor = MagmomPredictor(
            model_path,
            head_path=head_path,
            device=params.get("device", "cpu"),
            seed_mode=params.get("seed_mode", "lut"),
            magnetic_threshold=params.get("magnetic_threshold", 0.5),
            use_kernel=params.get("use_kernel"),
        )
        atoms.calc = predictor
        predicted = [float(item) for item in atoms.get_magnetic_moments()]
        print("prophet_spin predicted moments: "
              + " ".join(f"{item:.4f}" for item in predicted))
        # Step 2 evaluates the energy for the predicted moments.
        atoms.set_initial_magnetic_moments(predicted)
        atoms.calc = None

    calculator = SpinCalculator(
        model_path,
        device=params.get("device", "cpu"),
        use_kernel=params.get("use_kernel"),
    )
    atoms.calc = calculator
    energy = float(atoms.get_potential_energy())

    try:
        final_magmoms = [float(item) for item in atoms.get_magnetic_moments()]
    except Exception:  # noqa: BLE001 - model may not expose moments
        final_magmoms = [float(item) for item in (predicted or magmoms)]

    payload = {"energy": energy, "magmoms": final_magmoms}
    if predicted is not None:
        payload["magmoms_predicted"] = predicted
        payload["magmoms_initial"] = [float(item) for item in magmoms]

    with open("results.json", "wt") as handle:
        json.dump(payload, handle, indent=2)

    print(f"prophet_spin energy: {energy} eV")
    return 0


if __name__ == "__main__":
    sys.exit(main())
