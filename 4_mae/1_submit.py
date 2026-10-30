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
environment_activate = "source .venv/bin/activate"
environment_deactivate = "source .venv/bin/deactivate"
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

try:
    input_file = sys.argv[1]
    print(f'Using the {input_file} file as input')
    try:
        exec(f"from {input_file.split('.')[0]} import *")
    except:
        from input import *
except IndexError:
    print(f'Using the input.py file from folder: {cwd}')
    from input import *


import subprocess
import numpy as np
import shutil

from itertools import product
from pymatgen.io.vasp import Poscar
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from datetime import datetime
import spglib

# The MAE modes (mae_theta_phi / mae_curve) are only implemented in the manual
# submitter; the Fireworks path has no MAE support.
from common.SubmitManual import SubmitManual
if use_fireworks:
    print('MAE calculations are not implemented with Fireworks: using SubmitManual.')


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
# path to the folder with results from collinear calculations
path_to_coll = path_to_automag + '/2_coll/' + f'{formula}{struct_suffix}/{formula}_{calculator}/'
if not os.path.isdir(path_to_coll):
    path_to_coll = path_to_automag + '/2_coll/' + f'{formula}_{calculator}/'
    if not os.path.isdir(path_to_coll):
        path_to_coll = path_to_automag + '/2_coll/'

print(f'Path to the results from collinear calculations: {path_to_coll}')

path_to_trials = path_to_coll + f'trials_{formula}/'


# finding configuration and settings

setting = 1
magmom = None
while os.path.isfile(f'{path_to_trials}/configurations{setting:03d}.txt'):
    with open(f'{path_to_trials}/configurations{setting:03d}.txt', 'rt') as f:
        for line in f:
            values = line.split()
            if values[0] == configuration:
                magmom = [int(item) for item in values[1:]]
                break
        if magmom is not None:
            break
    setting += 1

trial_structure_file = f'{path_to_trials}/setting{setting:03d}.vasp'
use_collinear_results = (magmom is not None) and os.path.isfile(trial_structure_file)

if use_collinear_results:
    # create a pymatgen Structure object
    pmg_structure = Structure.from_file(trial_structure_file)
    ase_structure = pmg_structure.to_ase_atoms()

    # assigning magnetic moments to the atoms

    calc_results = f"{ase_structure.get_chemical_formula(mode='metal')}{struct_suffix}_singlepoint_{calculator}.txt"

    path_to_calc_results = path_to_automag + '/CalcFold/' + calc_results
    full_magmom = None
    print(f'Reading full magmoms for configuration {configuration} from file {path_to_calc_results}')
    with open(path_to_calc_results, 'rt') as f:
        lines = f.readlines()
        i = 0
        for line in lines:
            values = line.split()
            if values[0] == configuration:
                print(f'Found configuration: {values[0]}')
                magmoms_str = lines[i+1]
                break
            i += 1
        full_magmom_str = magmoms_str.split('final_magmoms=')[1].strip().strip('[]')
        number_magmoms = full_magmom_str.split()
        full_magmom = [float(num) for num in number_magmoms]

    if full_magmom is None:
        raise IOError(f'Magmoms for configuration {configuration} not found.')

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
            f'configuration {configuration} not found under {path_to_trials} and cannot be '
            'generated from spin_values. Run 2_coll for this structure first.'
        )
    print(f'WARNING: no 2_coll results for {formula}; using the geometry file '
          f'{poscar_file} and spin_values for {configuration}.')
    pmg_structure = Structure.from_file(path_to_poscar)
    spin_values = globals().get('spin_values', {})
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

# Submitting the calculations
run = SubmitManual(processed_structure_name, mode='mae_theta_phi', fix_params=params, magmoms=standardized_magmom,  
                     struct_suffix=struct_suffix, Nth = Nth, Nph = Nph,
                     name=configuration, calculator = calculator, jobheader = jobheader, calculator_command = calculator_command,
                     environment_activate = environment_activate, environment_deactivate = environment_deactivate, 
                     parallel_over_configurations = parallel_over_configurations)
run.submit()




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
