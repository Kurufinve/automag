"""
automag.2_relax.1_submit
========================

Submit a structural relaxation of the most stable magnetic configuration found
by the ``2_coll`` stage.

The configuration can be chosen explicitly with ``configuration``; otherwise
the most stable one reported by ``2_coll/2_plot_results.py`` is used.  The
converged magnetic moments from the collinear calculation are reused when they
are available, so the relaxation starts from the correct spin state.

Relaxation tags are taken from ``params`` and can be overridden through the
``relax_params`` dictionary in the input file (``ibrion``, ``isif``, ``nsw``,
``ediffg`` by default).

.. codeauthor:: Automag contributors
"""

import json
import sys

from common.input_loader import load_input
from common.magnetic import resolve_magnetic_elements, resolve_spin_values
from common.paths import data_root, geometries_dir
from common.results import load_state, resolve_configuration
from common.SubmitManual import SubmitManual
from pymatgen.core.structure import Structure

DEFAULTS = {
    "use_fireworks": False,
    "calculator": "vasp",
    "jobheader": "#!/bin/bash",
    "calculator_command": "mpirun vasp_std",
    "environment_activate": "source .venv/bin/activate",
    "environment_deactivate": "source .venv/bin/deactivate",
    "struct_suffix": "",
    "parallel_over_configurations": True,
    "spin_values": {},
    "magnetic_atoms": [],
}

settings = dict(DEFAULTS)
settings.update(load_input(sys.argv[1] if len(sys.argv) > 1 else None))

calculator = settings["calculator"]
coll_calculator = settings.get("coll_calculator") or calculator
struct_suffix = settings.get("struct_suffix", "")
poscar_file = settings["poscar_file"]
params = settings["params"]

path_to_poscar = str(geometries_dir() / poscar_file)
input_structure = Structure.from_file(path_to_poscar)
formula = input_structure.formula.replace(' ', '')

# The collinear results may have been produced with a different calculator
# than the one used for the relaxation (mixed-calculator pipeline).
state = resolve_configuration(
    settings.get("configuration"), data_root(), formula, struct_suffix,
    coll_calculator)
print(f'Relaxing the magnetic configuration: {state}')

spin_values = settings.get("spin_values", {})
magnetic_atoms = settings.get("magnetic_atoms", [])
# Optional: missing moments come from pymatgen's DEFAULT_MAGMOMS.
spin_values = resolve_spin_values(input_structure, spin_values,
                                 magnetic_atoms)
magnetic_atoms = resolve_magnetic_elements(input_structure, spin_values,
                                          magnetic_atoms)
info, structure, magmoms = load_state(
    data_root(), formula, state, struct_suffix=struct_suffix,
    calculator=coll_calculator, spin_values=spin_values)

# Relaxation parameters: the base calculation parameters plus overrides.
relax_params = dict(params)
relax_params.update(settings.get("relax_params", {}))
relax_params.setdefault('ibrion', 2)
relax_params.setdefault('isif', 3)
relax_params.setdefault('nsw', 100)
relax_params.setdefault('ediffg', -0.01)

run = SubmitManual(
    str(info['structure']), mode='relax', fix_params=relax_params,
    magmoms=list(magmoms), struct_suffix=struct_suffix, name=f'relax_{state}',
    calculator=calculator, jobheader=settings["jobheader"],
    calculator_command=settings["calculator_command"],
    environment_activate=settings["environment_activate"],
    environment_deactivate=settings["environment_deactivate"],
    parallel_over_configurations=settings["parallel_over_configurations"])
run.submit()

# Remember what was relaxed so that 2_collect.py does not have to guess.
info_out = {
    'state': state,
    'setting': info['setting'],
    'structure_file': str(info['structure']),
    'formula': formula,
    'struct_suffix': struct_suffix,
    'calculator': calculator,
    'coll_calculator': coll_calculator,
    'magmoms': [float(item) for item in magmoms],
    'relax_params': {
        key: value for key, value in relax_params.items()
        if isinstance(value, (int, float, str, bool))
    },
}
with open(f'relax_{state}.json', 'wt') as f:
    json.dump(info_out, f, indent=2)

print(f'Relaxation inputs written for configuration {state}.')
