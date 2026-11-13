"""
automag.calculators.prophet
===========================

Interface for the ``prophet_spin`` model from
https://github.com/kairosmaterial/prophet.

``prophet_spin.SpinCalculator`` is an ASE calculator that predicts energies
(and forces) from a trained model.  Because it runs in-process rather than as a
separate executable, the interface writes the structure, the initial moments
and the model settings to JSON files, and the run command executes
``common/run_prophet.py``, which loads the model and writes ``results.json``.

The checkpoints live in ``<automag>/models`` and are used by default, so a
stage only needs ``params = {"device": "cpu"}``.  Anything else can be given
explicitly::

    params = {
        "model_path": "/path/to/model",
        "device": "cpu",          # or "cuda"
        "use_kernel": None,
    }

Set ``predict_moments = True`` (plus ``head_path``) to run the two-step
variant: the magnetic-moment head predicts the moments of every atom from the
initial spin configuration, and the energy model is then evaluated for those
predicted moments::

    params = {
        "model_path": "/path/to/prophet-v1-mag.pt",
        "head_path": "/path/to/prophet-v1-mag-head.pt",
        "device": "cpu",
        "predict_moments": True,
    }

.. codeauthor:: Automag contributors
"""

import json
import os

import numpy as np
from ase.io import read, write

from common.paths import (PROPHET_HEAD_NAME, PROPHET_MODEL_NAME,
                         code_root, helper_command, model_file)
from calculators.base import CalculatorInterface


class ProphetSpinInterface(CalculatorInterface):
    name = "prophet_spin"
    supports_relaxation = False
    # The interface carries non-collinear moments (three numbers per atom), so
    # 4_mae can run.  Note that the released Prophet-V1-Mag weights are
    # isotropic in the moment direction: the energies of two orientations
    # differ by less than 1e-6 eV, so the resulting MAE is numerically zero.
    supports_noncollinear = True
    output_marker = "results.json"

    def write_input(self, atoms, directory, params=None):
        merged = dict(self.params)
        merged.update(params or {})
        os.makedirs(directory, exist_ok=True)

        write(os.path.join(directory, "atoms.json"), atoms, format="json")
        # 4_mae passes non-collinear moments as three numbers per atom, while
        # the collinear stages pass one; keep both shapes intact.
        moments = np.asarray(atoms.get_initial_magnetic_moments(), dtype=float)
        payload_magmoms = (moments.tolist() if moments.ndim == 2
                           else [float(item) for item in moments])
        with open(os.path.join(directory, "magmoms.json"), "wt") as handle:
            json.dump(payload_magmoms, handle)
        model_path = merged.get("model_path") or str(
            model_file(PROPHET_MODEL_NAME))
        payload = {
            "model_path": model_path,
            "device": merged.get("device", "cpu"),
            "use_kernel": merged.get("use_kernel"),
            # Optional two-step mode: predict the magnetic moments from the
            # initial configuration first, then evaluate the energy for the
            # predicted moments (see common/run_prophet.py).
            "predict_moments": bool(merged.get("predict_moments", False)),
            "head_path": (merged.get("head_path")
                          or (str(model_file(PROPHET_HEAD_NAME))
                              if merged.get("predict_moments") else None)),
            "seed_mode": merged.get("seed_mode", "lut"),
            "magnetic_threshold": merged.get("magnetic_threshold", 0.5),
        }
        with open(os.path.join(directory, "params.json"), "wt") as handle:
            json.dump(payload, handle, indent=2)

    def run_command(self):
        if self.command:
            return self.command
        return helper_command("common/run_prophet.py")

    def _results(self, directory):
        path = os.path.join(directory, "results.json")
        if not os.path.isfile(path):
            raise FileNotFoundError(f"No prophet_spin results in {directory}.")
        with open(path, "rt") as handle:
            return json.load(handle)

    def read_energy(self, directory):
        return float(self._results(directory)["energy"])

    def read_magnetic_moments(self, directory):
        magmoms = self._results(directory).get("magmoms")
        return None if magmoms is None else np.asarray(magmoms)

    def read_structure(self, directory):
        return read(os.path.join(directory, "atoms.json"))

    def is_converged(self, directory):
        return os.path.isfile(os.path.join(directory, "results.json"))

    def result_markers(self, directory):
        # run_prophet.py writes results.json into the calculation directory
        # (<jobdir>/singlepoint) and there is no OUTCAR to wait for.
        return [{"kind": "file",
                 "path": os.path.join(str(directory), "singlepoint",
                                      "results.json")}]

    def result_writer(self):
        script = code_root() / "common" / "write_output_ase.py"
        return str(script), None
