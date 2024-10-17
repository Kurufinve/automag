"""
automag.0_conv_tests.1_submit
=============================

Script which submits convergence tests.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""
# default values for some input variables
use_fireworks = False
calculator = 'vasp'
jobheader = """#!/bin/bash""" 
calculator_command = "mpirun vasp_std"
environment_activate = "source .venv/bin/activate"
environment_deactivate = "source .venv/bin/deactivate"

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

from pymatgen.core.structure import Structure

if use_fireworks:
    from common.SubmitFirework import SubmitFirework
else:
    from common.SubmitManual import SubmitManual

# full path to poscar file
rel_path_to_poscar = '/geometries/' + poscar_file
path_to_automag = os.environ.get('AUTOMAG_PATH')
path_to_poscar = path_to_automag + rel_path_to_poscar

# magnetic configuration to use for the convergence test
if 'configuration' not in globals():
    configuration = []
    structure = Structure.from_file(path_to_poscar)
    for atom in structure.species:
        if 'magnetic_atoms' not in globals():
            if atom.is_transition_metal:
                configuration.append(4.0)
            else:
                configuration.append(0.0)
        else:
            if atom in magnetic_atoms:
                configuration.append(4.0)
            else:
                configuration.append(0.0)

# convergence test w.r.t. encut
if mode == 'encut':
    if 'encut_values' not in globals():
        encut_values = range(500, 1010, 10)

    if use_fireworks:
        convtest = SubmitFirework(path_to_poscar, mode='encut', fix_params=params, magmoms=configuration,
                                  encut_values=encut_values)
        convtest.submit()
    else:
        convtest = SubmitManual(path_to_poscar, mode='encut', fix_params=params, magmoms=configuration,
                                encut_values=encut_values,
                                calculator = calculator, jobheader = jobheader, calculator_command = calculator_command,
                                environment_activate = environment_activate, environment_deactivate = environment_deactivate)                                  
        convtest.submit()

# convergence test w.r.t. sigma and kpts
if mode == 'kgrid':
    if 'sigma_values' not in globals():
        sigma_values = [item / 100 for item in range(5, 25, 5)]

    if 'kpts_values' not in globals():
        kpts_values = range(20, 110, 10)

    # submit calculations
    if use_fireworks:
        convtest = SubmitFirework(path_to_poscar, mode='kgrid', fix_params=params, magmoms=configuration,
                                  sigma_values=sigma_values, kpts_values=kpts_values)
        convtest.submit()

    else:
        convtest = SubmitManual(path_to_poscar, mode='kgrid', fix_params=params, magmoms=configuration,
                                sigma_values=sigma_values, kpts_values=kpts_values,
                                calculator = calculator, jobheader = jobheader, calculator_command = calculator_command,
                                environment_activate = environment_activate, environment_deactivate = environment_deactivate)
        convtest.submit()

