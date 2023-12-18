"""
automag.3_monte_carlo.4_run_ESpinS.py
======================================

Script which runs ESpinS for Monte Carlo simulation.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com> and Daniil Poletaev <d.poletaev@skoltech.ru>
"""

from input import *

import os
import numpy as np

from ase.io import read
from pymatgen.core.structure import Structure

# raise IOError if no trials folder
if not os.path.isdir('../2_coll/trials'):
    raise IOError('No trials folder found in ../2_coll.')

setting = 1
magmom = None
while os.path.isfile(f'../2_coll/trials/configurations{setting:03d}.txt'):
    with open(f'../2_coll/trials/configurations{setting:03d}.txt', 'rt') as f:
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

# assigning default values for parameters of biquadratic and Dzyalozhinskii-Moriya interactions
if 'Ham_bij' not in globals():
    Ham_bij = False
else:
    if "Shells_bij" not in globals:
        Shells_bij = len(coupling_constants)
if 'Ham_dij' not in globals():
    Ham_dij = False
else:
    if "Shells_dij" not in globals:
        Shells_dij = len(coupling_constants)

# get materials from magmom
materials = [0 if item >= 0 else 1 for item in magmom]

# create a pymatgen Structure object
structure = Structure.from_file(f'../2_coll/trials/setting{setting:03d}.vasp')

# assigning seedname for ESpinS
atoms = read(f'../2_coll/trials/setting{setting:03d}.vasp')
chemical_formula = atoms.get_chemical_formula(mode='metal', empirical=True)

seedname = f"{chemical_formula}_{configuration}"

# find out which atoms are magnetic
for element in structure.composition.elements:
    if 'magnetic_atoms' not in globals():
        element.is_magnetic = element.is_transition_metal
    else:
        if element.name in magnetic_atoms:
            element.is_magnetic = True
        else:
            element.is_magnetic = False

non_magnetic_atoms = [element.symbol for element in structure.composition.elements if not element.is_magnetic]
structure.remove_species(non_magnetic_atoms)
center_indices, point_indices, offset_vectors, distances = structure.get_neighbor_list(cutoff_radius)

# get thresholds for comparing atomic distances
thresholds = [(a + b) / 2 for a, b in zip(distances_between_neighbors[:-1], distances_between_neighbors[1:])]
thresholds = [0.0] + thresholds + [100.0]

try:
    os.system('mkdir espins')
except:
    pass

# get translation vectors for structure
vecs = structure.lattice.matrix

# write ESpinS inp1 file for 1st initialization step
with open(f'espins/{seedname}.inp1.mcin', 'w') as f:
    # print unit cell vectors for ESpinS inp1 file
    f.write('Begin unit_cell_cart\n')
    f.write(' Ang\n')
    f.write(f'{vecs[0][0]:12.8f}  {vecs[0][1]:12.8f}  {vecs[0][2]:12.8f}\n')
    f.write(f'{vecs[1][0]:12.8f}  {vecs[1][1]:12.8f}  {vecs[1][2]:12.8f}\n')
    f.write(f'{vecs[2][0]:12.8f}  {vecs[2][1]:12.8f}  {vecs[2][2]:12.8f}\n')
    f.write('End unit_cell_cart\n')

    # print fractional coordinates for ESpinS inp1 file
    f.write('Begin atoms_frac\n')
    for i, site in enumerate(zip(structure.sites,magmom)):
        label = structure.sites[i].label
        fc = structure.sites[i].frac_coords
        m = magmom[i]
        f.write(f'{label} {fc[0]: 12.8f}   {fc[1]: 12.8f}   {fc[2]: 12.8f} {m: 5.3f}\n')
    f.write('End atoms_frac\n')


    f.write(f'Shells_jij = {len(coupling_constants)}')

# run 1st initialization step before mc simulation in ESpinS
os.system(f"cd espins; {espins_run_command} -inp1 {seedname}")

# read ESpinS inp2 file for 2nd initialization step for editing 
with open(f'espins/{seedname}.inp2.mcin', 'r') as f:
    inp2_old = f.readlines()

inp2_new = []
# find ?????? in inp2.mcin file and replace them with Jij parameters from input.py
write_jij = False
i_jij = 0 # starting index of jij
for line in inp2_old:
    if "Begin Parameters_Jij" in line:
        write_jij = True
        inp2_new.append(line)
    elif "End Parameters_Jij" in line:
        write_jij = False
        inp2_new.append(line)
    else:
        if write_jij:
            newline = line.replace('??????',str(coupling_constants[i_jij]/1.60218e-19))
            inp2_new.append(newline)
            i_jij += 1
        else:
            inp2_new.append(line)

# write ESpinS inp2 file for 2nd initialization step
with open(f'espins/{seedname}.inp2.mcin', 'w') as f:
    f.writelines(inp2_new)

# run 2nd initialization step before mc simulation in ESpinS
os.system(f"cd espins; {espins_run_command} -inp2 {seedname}")

# write the command for running the Monte-Carlo simulation that should be executed separately
with open(f'espins/run_mc.sh', 'w') as f:
    f.write('#!/bin/bash\n')
    f.write(f'{espins_run_command} {seedname}')
