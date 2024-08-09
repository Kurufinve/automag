"""
automag.1_lin_response.1_submit
===============================

Script which submits linear response U calculations.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

# default values for some input variables
use_fireworks = True

from input import *

from ase.io import read
from pymatgen.core.structure import Structure

from common.SubmitFirework import SubmitFirework

# full path to poscar file
path_to_poscar = '../geometries/' + poscar_file

# create ase.Atoms object
atoms = read(path_to_poscar)

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

# submit calculations
if use_fireworks:
    run = SubmitFirework(path_to_poscar, mode='perturbations', fix_params=params, pert_values=perturbations,
                         magmoms=configuration, dummy_atom=dummy_atom, dummy_position=dummy_position)
    run.submit()
else:
    run = SubmitManual(path_to_poscar, mode='perturbations', fix_params=params, pert_values=perturbations,
                       magmoms=configuration, dummy_atom=dummy_atom, dummy_position=dummy_position, calculator = calculator, jobheader = jobheader, calculator_command = calculator_command,
                       environment_activate = environment_activate, environment_deactivate = environment_deactivate)
    run.submit()    
