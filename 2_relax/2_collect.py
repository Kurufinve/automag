"""
automag.2_relax.2_collect
=========================

Collect the result of the relaxation started by ``1_submit.py``.

It reads ``CONTCAR``/``OUTCAR`` of the relaxation, writes ``relaxed_<state>.vasp``
and ``relaxed_<state>.json`` into the stage directory, and - if
``update_structure = True`` - copies the relaxed geometry into
``<workdir>/geometries/`` under a ``_relaxed`` name.

.. codeauthor:: Automag contributors
"""

import glob
import json
import os
import sys
from pathlib import Path

from ase import Atoms
from ase.io import read

from common.input_loader import load_input
from common.paths import data_root, geometries_dir
from common.results import most_stable_state
from pymatgen.core.structure import Structure

DEFAULTS = {
    "calculator": "vasp",
    "struct_suffix": "",
    "update_structure": False,
}

settings = dict(DEFAULTS)
settings.update(load_input(sys.argv[1] if len(sys.argv) > 1 else None))

calculator = settings["calculator"]
coll_calculator = settings.get("coll_calculator") or calculator
struct_suffix = settings.get("struct_suffix", "")
update_structure = settings.get("update_structure", False)
poscar_file = settings["poscar_file"]


def read_single(path) -> Atoms:
    """Read one structure, unwrapping the list ASE may return."""
    result = read(path)
    if isinstance(result, list):
        result = result[-1]
    return result


path_to_poscar = str(geometries_dir() / poscar_file)
formula = Structure.from_file(path_to_poscar).formula.replace(' ', '')

# Determine which configuration was relaxed.
state = settings.get("configuration")
if state and str(state).strip().lower() in {"lowest", "most_stable", "most-stable"}:
    state = most_stable_state(data_root(), formula, struct_suffix,
                              coll_calculator)
if not state:
    candidates = sorted(glob.glob('relax_*.json'),
                        key=os.path.getmtime, reverse=True)
    if candidates:
        with open(candidates[0], 'rt') as f:
            state = json.load(f).get('state')
if not state:
    state = most_stable_state(data_root(), formula, struct_suffix,
                              coll_calculator)

relax_dir = (data_root() / 'CalcFold' / f'{formula}{struct_suffix}'
             / calculator / f'relax_{state}' / 'singlepoint')
outcar = relax_dir / 'OUTCAR'
contcar = relax_dir / 'CONTCAR'

if not outcar.is_file() or not contcar.is_file() or contcar.stat().st_size == 0:
    raise OSError(
        f'Relaxation output for {state} not found in {relax_dir}. '
        'Run 1_submit.py and wait for the calculation to finish.'
    )

outcar_text = outcar.read_text(errors='ignore')
finished = (
    'General timing and accounting informations' in outcar_text
    or 'Voluntary context switches' in outcar_text
    or 'reached required accuracy' in outcar_text
)
if not finished:
    raise OSError(f'The relaxation in {relax_dir} has not finished yet.')

atoms = read_single(str(contcar))
relaxed_name = f'relaxed_{state}.vasp'
atoms.write(relaxed_name, format='vasp')

# Energy and magnetic moments live in OUTCAR, not in CONTCAR.
energy = None
magmoms = []
try:
    result_atoms = read_single(str(outcar))
    energy = float(result_atoms.get_potential_energy())
    magmoms = [float(item) for item in result_atoms.get_magnetic_moments()]
except Exception:  # noqa: BLE001 - OUTCAR parsing is best-effort
    print(f'WARNING: could not read energy/magnetic moments from {outcar}')

result = {
    'state': state,
    'formula': formula,
    'struct_suffix': struct_suffix,
    'calculator': calculator,
    'structure_file': relaxed_name,
    'energy': energy,
    'magmoms': magmoms,
}
with open(f'relaxed_{state}.json', 'wt') as f:
    json.dump(result, f, indent=2)

if update_structure:
    destination = geometries_dir() / f'{Path(poscar_file).stem}_relaxed.vasp'
    atoms.write(str(destination), format='vasp')
    print(f'Relaxed structure copied to {destination}')

print(f'Relaxed structure written to {relaxed_name}')
if energy is not None:
    print(f'Energy: {energy:.6f} eV')
