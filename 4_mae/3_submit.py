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

if use_fireworks:
    from common.SubmitFirework import SubmitFirework
else:
    from common.SubmitManual import SubmitManual


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
    raise IOError(f'configuration {configuration} not found.')

# # get materials from magmom
# materials = [0 if item >= 0 else 1 for item in magmom]

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

standardized_structure = pmg_structure.get_primitive_structure(tolerance=0.2, use_site_props=True)


standardized_structure_name = f'setting{setting:03d}_{configuration}_standardized.vasp'

standardized_structure.to(filename=standardized_structure_name,fmt='POSCAR')
standardized_magmom = list(standardized_structure.site_properties['magmom'])
ncl_magmom = [(0,0,m) for m in standardized_magmom]
standardized_magmom = ncl_magmom

print(f'Magnetic moments of all atoms in standardized cell for configuration {configuration}:')
print(standardized_magmom)

U = np.round(float(params['ldauu'][next(i for i, x in enumerate(params['ldaul']) if x > 0)]),1)
J = np.round(float(params['ldauj'][next(i for i, x in enumerate(params['ldaul']) if x > 0)]),1)
encut = params['encut']
kpts = params['kpts']
# mae_dir = os.path.join(state_dir, f'mae_U{U:.1f}_J{J:.1f}_K{kpts}_EN{encut}')
outputs = open(f'MAE_theta_phi_{configuration}_U{U:.1f}_J{J:.1f}_K{kpts}_EN{encut}.txt','r')
lines = outputs.readlines()
for line in lines:
    if ('MAE_x' in line) or ('MAE_y' in line) or ('MAE_z' in line):
        start = line.find('[')
        end = line.find(']')
        array_str = line[start:end+1]
        if 'MAE_x' in line:
            MAE_x = np.fromstring(array_str[1:-1], sep=', ')
        if 'MAE_y' in line:
            MAE_y = np.fromstring(array_str[1:-1], sep=', ')
        if 'MAE_z' in line:
            MAE_z = np.fromstring(array_str[1:-1], sep=', ')

print(f'MAE_x = {MAE_x}')
print(f'MAE_y = {MAE_y}')
print(f'MAE_z = {MAE_z}')

outputs.close()

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
