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

def fmt_kpts(kpts):
    """Compact label for a k-point specification (scalar or explicit mesh)."""
    if isinstance(kpts, (list, tuple)):
        return 'x'.join(str(k) for k in kpts)
    return str(kpts)


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

if magmom is None:
    print(f'No 2_coll configuration data for {configuration}.')

# # get materials from magmom
# materials = [0 if item >= 0 else 1 for item in magmom]

# ---------------------------------------------------------------------------
# Use exactly the same cell as 1_submit.py: it is described by the sidecar
# JSON file written there. Re-deriving the cell from 2_coll (the old behaviour)
# is only used for runs started before the sidecar existed.
# ---------------------------------------------------------------------------
cell_info = None
cell_info_name = None
for setting_number in range(1, 10):
    candidate = f'setting{setting_number:03d}_{configuration}_processed.json'
    if os.path.isfile(candidate):
        with open(candidate, 'rt') as f:
            cell_info = json.load(f)
        cell_info_name = candidate
        break

if cell_info is not None and cell_info.get('setting'):
    setting = int(cell_info['setting'])

if magmom is None and cell_info is None:
    raise IOError(f'configuration {configuration} not found and no processed-cell sidecar file is present.')

if cell_info is not None and os.path.isfile(cell_info['structure_file']):
    standardized_structure_name = cell_info['structure_file']
    standardized_magmom = [tuple(m) for m in cell_info['magmoms']]
    print(f'Using the calculation cell from {cell_info_name}: '
          f'{standardized_structure_name} ({cell_info["n_atoms"]} atoms, '
          f'V = {cell_info["volume"]:.4f} A^3)')
else:
    print(f'WARNING: {cell_info_name} not found; re-deriving the cell from 2_coll.')
    # create a pymatgen Structure object
    pmg_structure = Structure.from_file(f'{path_to_trials}/setting{setting:03d}.vasp')
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
    pmg_structure.add_site_property("magmom", full_magmom)

    standardize_cell = globals().get('standardize_cell', True)
    use_primitive_cell = globals().get('use_primitive_cell', True)
    if not standardize_cell:
        standardized_structure = pmg_structure.copy()
    elif use_primitive_cell:
        standardized_structure = pmg_structure.get_primitive_structure(tolerance=0.2, use_site_props=True)
    else:
        raise NotImplementedError(
            'standardize_cell=True with use_primitive_cell=False is not supported.'
        )

    standardized_structure_name = f'setting{setting:03d}_{configuration}_processed.vasp'
    standardized_structure.to(filename=standardized_structure_name, fmt='POSCAR')
    standardized_magmom = [(0, 0, float(m)) for m in standardized_structure.site_properties['magmom']]

U = np.round(float(params['ldauu'][next(i for i, x in enumerate(params['ldaul']) if x > 0)]),1)
J = np.round(float(params['ldauj'][next(i for i, x in enumerate(params['ldaul']) if x > 0)]),1)
encut = params['encut']
kpts = params['kpts']
mae_file = f'MAE_theta_phi_{configuration}_U{U:.1f}_J{J:.1f}_K{fmt_kpts(kpts)}_EN{encut}.txt'
if not os.path.isfile(mae_file):
    raise IOError(f'{mae_file} not found: run 2_plot_results.py first.')

MAE_x = MAE_y = MAE_z = None
with open(mae_file, 'rt') as f:
    for line in f:
        if line.startswith('MAE_x:'):
            MAE_x = np.array(json.loads(line[line.find('['):line.find(']')+1]))
        elif line.startswith('MAE_y:'):
            MAE_y = np.array(json.loads(line[line.find('['):line.find(']')+1]))
        elif line.startswith('MAE_z:'):
            MAE_z = np.array(json.loads(line[line.find('['):line.find(']')+1]))

if MAE_x is None or MAE_y is None or MAE_z is None:
    raise IOError(f'Could not read MAE_x/MAE_y/MAE_z from {mae_file}.')

print(f'MAE_x (easy axis) = {MAE_x}')
print(f'MAE_y (hard axis) = {MAE_y}')
print(f'MAE_z             = {MAE_z}')

# Submitting the calculations
run = SubmitManual(standardized_structure_name, mode='mae_curve', fix_params=params, magmoms=standardized_magmom,  
                     struct_suffix=struct_suffix, N_MAE=N_MAE, MAE_x = MAE_x, MAE_y = MAE_y, MAE_z = MAE_z,
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
