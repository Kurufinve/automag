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

# create Structure and SymmetrizedStructure objects
structure = Structure.from_file(path_to_poscar)
analyzer = SpacegroupAnalyzer(structure)
symmetrized_structure = analyzer.get_symmetrized_structure()

# find out which atoms are magnetic
for element in structure.composition.elements:
    if element.name in spin_values:
        element.is_magnetic = True
    else:
        element.is_magnetic = False



# geometrical settings and respective lists of magnetic configurations
lattices = [structure.lattice]
coordinates = [structure.frac_coords]
configurations = [[]]

# get the multiplicities of each Wyckoff position
multiplicities = [len(item) for item in symmetrized_structure.equivalent_indices]

wyckoff_magmoms = []
equivalent_multipliers = []
for multiplicity, wyckoff in zip(multiplicities, symmetrized_structure.equivalent_sites):
    if wyckoff[0].specie.is_magnetic:
        wyckoff_magmoms.append([1, 0, -1])

        if len(spin_values[wyckoff[0].specie.name]) == 2:
            val1 = spin_values[wyckoff[0].specie.name][0]
            val2 = spin_values[wyckoff[0].specie.name][1]

            if len(equivalent_multipliers) == 0:
                equivalent_multipliers.append(np.repeat(val1, multiplicity))
                equivalent_multipliers.append(np.repeat(val2, multiplicity))
            else:
                new_equivalent_multipliers = []
                for multiplier in equivalent_multipliers:
                    new_equivalent_multipliers.append(np.append(multiplier, np.repeat(val1, multiplicity)))
                    new_equivalent_multipliers.append(np.append(multiplier, np.repeat(val2, multiplicity)))

                equivalent_multipliers = new_equivalent_multipliers

        elif len(spin_values[wyckoff[0].specie.name]) == 1:
            val1 = spin_values[wyckoff[0].specie.name][0]

            if len(equivalent_multipliers) == 0:
                equivalent_multipliers.append(np.repeat(val1, multiplicity))
            else:
                new_equivalent_multipliers = []
                for multiplier in equivalent_multipliers:
                    new_equivalent_multipliers.append(np.append(multiplier, np.repeat(val1, multiplicity)))

                equivalent_multipliers = new_equivalent_multipliers

        else:
            raise ValueError('Max 2 spin values for each magnetic species.')

    else:
        wyckoff_magmoms.append([0])

        if len(equivalent_multipliers) == 0:
            equivalent_multipliers.append(np.repeat(0, multiplicity))
        else:
            new_equivalent_multipliers = []
            for multiplier in equivalent_multipliers:
                new_equivalent_multipliers.append(np.append(multiplier, np.repeat(0, multiplicity)))

            equivalent_multipliers = new_equivalent_multipliers

# get all possible configurations without any splitting
for conf in product(*wyckoff_magmoms):
    configuration = np.repeat(conf, multiplicities)

    for mult in equivalent_multipliers:
        candidate_conf = np.multiply(configuration, mult).tolist()

        # if the configuration is not NM, get the index of the first non-zero element
        if len(np.nonzero(candidate_conf)[0]) > 0:
            first_nonzero_index = np.nonzero(candidate_conf)[0][0]

            # if the first non-zero spin is negative flip all spins
            if candidate_conf[first_nonzero_index] < 0:
                candidate_conf = [-item for item in candidate_conf]

        # add to list of configurations for the current settings
        if candidate_conf not in configurations[0]:
            configurations[0].append(candidate_conf)

# split all possible combinations of Wyckoff positions
splits = []
possibilities = [[0, 1] if len(item) > 1 else [0] for item in wyckoff_magmoms]
for split in product(*possibilities):
    if sum(split) != 0:
        splits.append(split)

for i, split in enumerate(splits):
    launch_enumlib(i + 1, split)

# merge the first and second settings if they are equivalent
if len(lattices) > 1:
    transformation_matrix = np.dot(lattices[1].matrix, np.linalg.inv(lattices[0].matrix))
    inv_transformation_matrix = np.linalg.inv(transformation_matrix)
    determinant = int(round(np.linalg.det(transformation_matrix), 6))

    if determinant == 1:
        del_flag = True
        origin_shift = np.around(coordinates[1][0] - np.dot(coordinates[0][0], inv_transformation_matrix), decimals=6)
        for coord1, coord2 in zip(coordinates[0], coordinates[1]):
            if not np.allclose((np.dot(coord1, inv_transformation_matrix) + origin_shift) % 1, coord2):
                del_flag = False

        if del_flag:
            del lattices[0]
            del coordinates[0]
            configurations[0].extend(configurations[1])
            del configurations[1]

# write output and submit calculations
fm_count = 1
afm_count = 1
fim_count = 1
original_ch_symbols = [atom.name for atom in structure.species]
for i, (lattice, frac_coords, confs) in enumerate(zip(lattices, coordinates, configurations)):
    magnification = len(frac_coords) // len(structure.frac_coords)
    ch_symbols = np.repeat(original_ch_symbols, magnification)
    setting = Structure(lattice, ch_symbols, frac_coords)
    setting.to(fmt='poscar', filename=f'setting{i + 1:03d}.vasp')
    mask = [item.is_magnetic for item in setting.species]

    for conf in confs:
        conf_array = np.array(conf)
        with open(f'configurations{i + 1:03d}.txt', 'a') as f:
            if np.sum(np.abs(conf)) == 0:
                state = 'nm'
            elif min(conf) >= 0:
                state = 'fm' + str(fm_count)
                fm_count += 1
            elif np.sum(conf) == 0:
                state = 'afm' + str(afm_count)
                afm_count += 1
            else:
                state = 'fim' + str(fim_count)
                fim_count += 1

            f.write(f'{state:>6s}  ')
            f.write(' '.join(f'{e:2d}' for e in conf_array[mask]))
            f.write('\n')

        if use_fireworks:
            run = SubmitFirework(f'setting{i + 1:03d}.vasp', mode='singlepoint', fix_params=params, magmoms=conf,
                                 name=state)
            run.submit()

        else:
            run = SubmitManual(f'setting{i + 1:03d}.vasp', mode='singlepoint', fix_params=params, magmoms=conf,
                                 name=state, calculator = calculator, jobheader = jobheader, calculator_command = calculator_command,
                                 environment_activate = environment_activate, environment_deactivate = environment_deactivate)
            run.submit()
