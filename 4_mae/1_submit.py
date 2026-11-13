"""
automag.4_mae.1_submit
=======================

Script which submits calculation for MAE (magnetocrystalline anisotropy energy) calculation.

.. codeauthor:: Daniil Poletaev <poletaev.dan@gmail.com>
"""

# default values for some input variables
use_fireworks = True
calculator = 'vasp'
jobheader = """#!/bin/bash""" 
calculator_command = "mpirun vasp_std"
environment_activate = ""      # no activation by default: the
# pipeline runs the jobs in the environment it was started in
environment_deactivate = ""
parallel_over_configurations = True
struct_suffix = ''
# number of grid points for phi angle rotation (0-360 degrees)
Nph = 20
# number of grid points for theta angle rotation (0-180 degrees)
Nth = 10
# number of grid points for MAE curve (0-360 degrees)
N_MAE = 20

import os,sys
import json

cwd = os.getcwd()

from common.input_loader import apply_stage_input

apply_stage_input(globals())


import subprocess
import numpy as np
import shutil

from ase.io import read as ase_read

from itertools import product
from pymatgen.io.vasp import Poscar
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from datetime import datetime
import spglib

# The MAE modes (mae_theta_phi / mae_curve) are only implemented in the manual
# submitter; the Fireworks path has no MAE support.
from common.SubmitManual import SubmitManual
from common.results import find_state, final_magmoms, initial_magmoms, resolve_configuration
from common.magnetic import resolve_spin_values
from common.sweep import expand_params, write_manifest
if use_fireworks:
    print('MAE calculations are not implemented with Fireworks: using SubmitManual.')


def fmt_kpts(kpts):
    """Compact label for a k-point specification (scalar or explicit mesh)."""
    if isinstance(kpts, (list, tuple)):
        return 'x'.join(str(k) for k in kpts)
    return str(kpts)


def mae_directory(state_dir, params, tag):
    """Directory of one MAE parameter set (mirrors SubmitManual.create_mae_dir)."""
    index = next(i for i, x in enumerate(params['ldaul']) if x > 0)
    U = np.round(float(params['ldauu'][index]), 1)
    J = np.round(float(params['ldauj'][index]), 1)
    return os.path.join(
        state_dir, f"mae_U{U:.1f}_J{J:.1f}_K{fmt_kpts(params['kpts'])}"
                   f"_EN{params['encut']}{tag}")


def magnetization(fileOutCar):
    Mag = 'Nan'
    with open(fileOutCar) as file:#open("./z/OUTCAR") as file:
        for line in file:
            if "General timing and accounting informations for this job:" in line:  
                with open("./z/OSZICAR") as file:
                    for line in file:
                        if 'mag=' in line:
                            head, sep, tail = line.partition('mag=')
                            Mag = [float(j) for j in re.findall(r"[-+]?\d*\.?\d+|\d+", tail)]
                M0 = np.sqrt((Mag[0]**2)+(Mag[1]**2)+(Mag[2]**2))
                #print('total magnetization:', np.sqrt((Mag[0]**2)+(Mag[1]**2)+(Mag[2]**2)))
    return M0




# full path to poscar file
rel_path_to_poscar = '/geometries/' + poscar_file
path_to_automag = os.environ.get('AUTOMAG_PATH')
path_to_poscar = path_to_automag + rel_path_to_poscar
input_structure = Structure.from_file(path_to_poscar)

formula = input_structure.formula.replace(' ','')
path_to_automag = os.environ.get('AUTOMAG_PATH')
# The collinear results may have been produced with a different calculator
# than the one used for the MAE calculations (mixed-calculator pipeline).
coll_calculator = globals().get('coll_calculator') or calculator
configuration = resolve_configuration(
    configuration, path_to_automag, formula, struct_suffix, coll_calculator)
print(f'Reading the collinear results produced by: {coll_calculator}')


# finding configuration and settings

try:
    collinear = find_state(path_to_automag, formula, configuration,
                           struct_suffix, calculator=coll_calculator)
except (KeyError, FileNotFoundError) as error:
    print(f'WARNING: no 2_coll results for {formula}: {error}')
    collinear = None

magmom = collinear['spins'] if collinear else None
trial_structure_file = str(collinear['structure']) if collinear else None
setting = int(collinear['setting']) if collinear else 1
use_collinear_results = (magmom is not None) and os.path.isfile(trial_structure_file)

if use_collinear_results:
    # create a pymatgen Structure object
    pmg_structure = Structure.from_file(trial_structure_file)

    # Converged magnetic moments from the collinear run (may be absent if the
    # 2_coll jobs are still running); in that case reconstruct them from the
    # spin pattern.
    full_magmom = final_magmoms(path_to_automag, formula, configuration,
                                struct_suffix, coll_calculator)
    if full_magmom is None or len(full_magmom) != len(pmg_structure):
        full_magmom = initial_magmoms(pmg_structure, magmom,
                                      globals().get('spin_values', {}))
    full_magmom = [float(item) for item in full_magmom]

    print(f'Magnetic moments of all atoms for configuration {configuration}:')
    print(full_magmom)

    if len(full_magmom) != len(pmg_structure):
        raise IOError(
            f'The number of converged magnetic moments ({len(full_magmom)}) does not match '
            f'the number of atoms in {trial_structure_file} ({len(pmg_structure)}).'
        )
else:
    # No 2_coll results for this formula. Build the initial moments from
    # spin_values (like MAE.py's NAGMOM_std) and use the geometry file
    # directly. Only unambiguous orderings (ferro-/non-magnetic) are supported
    # this way; antiferromagnetic ones need the 2_coll data.
    if not (configuration.startswith('fm') or configuration == 'nm'):
        raise IOError(
            f'configuration {configuration} not found and cannot be '
            'generated from spin_values. Run 2_coll for this structure first.')
    print(f'WARNING: no 2_coll results for {formula}; using the geometry file '
          f'{poscar_file} and spin_values for {configuration}.')
    pmg_structure = Structure.from_file(path_to_poscar)
    spin_values = globals().get('spin_values', {})
    # Optional: fall back on pymatgen's DEFAULT_MAGMOMS before giving up.
    spin_values = resolve_spin_values(
        pmg_structure, spin_values, globals().get('magnetic_atoms'))
    nonmagnetic_magmom = float(globals().get('nonmagnetic_magmom', 0.0))
    if configuration == 'nm':
        full_magmom = [0.0] * len(pmg_structure)
    else:
        full_magmom = []
        for site in pmg_structure:
            symbol = site.specie.symbol
            full_magmom.append(float(spin_values[symbol][0]) if symbol in spin_values
                               else nonmagnetic_magmom)
    print(f'Initial magnetic moments for configuration {configuration}:')
    print(full_magmom)

pmg_structure.add_site_property("magmom", full_magmom)

# ---------------------------------------------------------------------------
# Choose the cell that will actually be calculated.
#
# MAE.py runs on ./POSCAR exactly as given. Reducing/standardizing the cell
# changes the Cartesian frame in which SAXIS (and therefore the easy/hard axes)
# is expressed and the volume used for the MJ/m^3 conversion. The choice is
# therefore explicit here and is written to a sidecar JSON file which scripts
# 2-4 read, so that the calculation, the axes and the normalization all use the
# *same* cell.
# ---------------------------------------------------------------------------
standardize_cell = globals().get('standardize_cell', True)
use_primitive_cell = globals().get('use_primitive_cell', True)

if not standardize_cell:
    processed_structure = pmg_structure.copy()
    cell_type = 'input_cell'
elif use_primitive_cell:
    processed_structure = pmg_structure.get_primitive_structure(tolerance=0.2, use_site_props=True)
    cell_type = 'primitive_cell'
else:
    raise NotImplementedError(
        'standardize_cell=True with use_primitive_cell=False is not supported: '
        'the converged magnetic moments cannot be mapped unambiguously onto the '
        'conventional cell. Use standardize_cell=False to keep the input cell, '
        'or use_primitive_cell=True.'
    )

if 'magmom' not in processed_structure.site_properties:
    raise IOError('The magnetic moments were lost while processing the cell.')

processed_structure_name = f'setting{setting:03d}_{configuration}_processed.vasp'
processed_structure.to(filename=processed_structure_name, fmt='POSCAR')

processed_magmom = [float(m) for m in processed_structure.site_properties['magmom']]
standardized_magmom = [(0, 0, m) for m in processed_magmom]

print(f'Calculation cell: {cell_type}, {len(processed_structure)} atoms, '
      f'V = {processed_structure.volume:.4f} A^3')
print(f'Magnetic moments of all atoms in the calculated cell for configuration {configuration}:')
print(standardized_magmom)

# sidecar file used by 2_plot_results.py / 3_submit.py / 4_plot_results.py
cell_info = {
    'structure_file': processed_structure_name,
    'configuration': configuration,
    'setting': setting,
    'formula': processed_structure.composition.reduced_formula.replace(' ', ''),
    'cell_type': cell_type,
    'n_atoms': len(processed_structure),
    'volume': float(processed_structure.volume),
    'magmoms': [list(m) for m in standardized_magmom],
}
cell_info_name = f'setting{setting:03d}_{configuration}_processed.json'
with open(cell_info_name, 'wt') as f:
    json.dump(cell_info, f, indent=2)
print(f'Cell information written to {cell_info_name}')

# ---------------------------------------------------------------------------
# Submitting the calculations.
#
# A parameter that is given as a list instead of a single number turns the
# stage into a convergence study: every combination becomes an independent
# MAE calculation (own directory, own reference calculation) and the pipeline
# runs the sets in parallel like any other independent calculation.
# ---------------------------------------------------------------------------
processed_atoms = ase_read(processed_structure_name)
calc_formula = processed_atoms.get_chemical_formula(mode='metal')
state_dir = os.path.join(path_to_automag, 'CalcFold',
                         f'{calc_formula}{struct_suffix}', calculator,
                         configuration)

sweep_sets = expand_params(params)
if len(sweep_sets) > 1:
    print(f'MAE convergence study over {len(sweep_sets)} parameter sets:')
    for sweep_set in sweep_sets:
        print(f'  {sweep_set.tag}: {sweep_set.label}')
elif sweep_sets[0].values:
    print(f'MAE convergence study: {sweep_sets[0].label}')

set_directories = []
for sweep_set in sweep_sets:
    run = SubmitManual(processed_structure_name, mode='mae_theta_phi',
                       fix_params=sweep_set.params, magmoms=standardized_magmom,
                       struct_suffix=struct_suffix, Nth = Nth, Nph = Nph,
                       name=configuration, calculator = calculator, jobheader = jobheader, calculator_command = calculator_command,
                       environment_activate = environment_activate, environment_deactivate = environment_deactivate,
                       parallel_over_configurations = parallel_over_configurations,
                       mae_tag=sweep_set.tag)
    run.submit()

    # Remember the settings of this set next to its calculations, so that the
    # post-processing scripts do not have to guess them.
    set_dir = mae_directory(state_dir, sweep_set.params, sweep_set.tag)
    os.makedirs(set_dir, exist_ok=True)
    set_directories.append(set_dir)
    with open(os.path.join(set_dir, 'parameters.json'), 'wt') as handle:
        json.dump({
            'tag': sweep_set.tag,
            'label': sweep_set.label,
            'parameters': sweep_set.params,
            'values': sweep_set.values,
            'configuration': configuration,
            'mae_dir': set_dir,
        }, handle, indent=2, default=str)

# Record which calculation directories belong to this run, so that the
# analysis scripts ignore leftovers of an earlier parameter set.
write_manifest(os.getcwd(), set_directories)




# if use_fireworks:
#     print('Calculation of MAE is not implemented with Fireworks!')
#     # run = SubmitFirework(f'setting{i + 1:03d}.vasp', mode='singlepoint', fix_params=params, magmoms=conf,
#     #                      name=state)
#     # run.submit()

# else:
#     run = SubmitManual(f'setting{i + 1:03d}.vasp', mode='singlepoint', fix_params=params, magmoms=conf, ntheta = Nth, nphi = Nph,
#                          name=state, calculator = calculator, jobheader = jobheader, calculator_command = calculator_command,
#                          environment_activate = environment_activate, environment_deactivate = environment_deactivate)
#     run.submit()
