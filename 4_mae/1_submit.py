"""
automag.4_mae.1_submit
=======================

Script which submits calculation for MAE (magnetocrystalline anisotropy energy) calculation.

.. codeauthor:: Daniil Poletaev <poletaev.dan@gmail.com>
"""

# default values for some input variables
# MAE grid/curve submission always goes through the manual submitter (matching
# MAE.py); the fireworks path has no MAE implementation, so default to False.
use_fireworks = False
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
from datetime import datetime
import spglib

if use_fireworks:
    from common.SubmitFirework import SubmitFirework
else:
    from common.SubmitManual import SubmitManual

# Shared helpers (same ones used by the refactored MAE scripts) so that the
# configuration lookup, magmom mapping and cell selection stay in one place.
from mae_utils import load_collinear_configuration, process_structure


# full path to poscar file
rel_path_to_poscar = '/geometries/' + poscar_file
path_to_automag = os.environ.get('AUTOMAG_PATH')
path_to_poscar = path_to_automag + rel_path_to_poscar
input_structure = Structure.from_file(path_to_poscar)

formula = input_structure.formula.replace(" ", "")
# Load the configuration and its final moments from the collinear results.
structure_path, full_magmom, setting = load_collinear_configuration(
    path_to_automag=path_to_automag,
    formula=formula,
    configuration=configuration,
    struct_suffix=struct_suffix,
    calculator=calculator,
)

print(f'Magnetic moments of all atoms for configuration {configuration}:')
print(full_magmom)

# Cell selection and magmom mapping follow the shared options
# (standardize_cell / use_primitive_cell), exactly as 1_submit_refactored.py.
standardize_cell = globals().get('standardize_cell', True)
use_primitive_cell = globals().get('use_primitive_cell', True)

standardized_structure_name, standardized_magmom, standardized_structure = process_structure(
    structure_path,
    full_magmom,
    f'setting{setting:03d}_{configuration}_processed.vasp',
    symmetrize=standardize_cell,
    to_primitive=use_primitive_cell,
)

print(f'Magnetic moments of all atoms in standardized cell for configuration {configuration}:')
print(standardized_magmom)

# Submitting the calculations
run = SubmitManual(standardized_structure_name, mode='mae_theta_phi', fix_params=params, magmoms=standardized_magmom,  
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
